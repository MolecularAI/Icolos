import os
import tempfile
from tracemalloc import start
from typing import List
from pydantic import BaseModel
import pandas as pd
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.step_utils.sdconvert_util import SDConvertUtil
from icolos.core.step_utils.structcat_util import StructcatUtil
from icolos.core.step_utils.structconvert import StructConvert
from icolos.core.workflow_steps.step import StepBase
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepActiveLearningEnum, StepBaseEnum
from icolos.utils.enums.write_out_enums import WriteOutEnum
from icolos.utils.general.convenience_functions import nested_get
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from sklearn.gaussian_process.kernels import DotProduct
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor
from icolos.core.workflow_steps.active_learning.al_utils import greedy_acquisition
from modAL.models.learners import BayesianOptimizer, ActiveLearner
from icolos.core.workflow_steps.active_learning.models.ffnn import FeedForwardNet
from skorch.regressor import NeuralNetRegressor
import torch
from torch import nn

_IE = StepInitializationEnum()
_SALE = StepActiveLearningEnum()
_WE = WorkflowEnum()
_WOE = WriteOutEnum()


class ActiveLearningBase(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def _initialize_oracle_step_from_dict(self, step_conf: dict) -> StepBase:
        # note this is a bit of a hack to get around a circular import, we can't use the main util
        _STE = StepBaseEnum
        step_type = nested_get(step_conf, _STE.STEP_TYPE, default=None)
        step_type = None if step_type is None else step_type.upper()
        if step_type in _IE.STEP_INIT_DICT.keys():
            return _IE.STEP_INIT_DICT[step_type](**step_conf)
        else:
            raise ValueError(
                f"Backend for step {nested_get(step_conf, _STE.STEPID, '')} unknown."
            )

    def construct_fingerprints(self, library: pd.DataFrame):
        # add morgan FPs
        library[_SALE.MORGAN_FP] = library.apply(
            lambda x: np.array(
                GetMorganFingerprintAsBitVect(x[_SALE.MOLECULE], 2, nBits=2048),
                dtype=np.float32,
            ),
            axis=1,
        )

        library[_SALE.IDX] = [i for i in range(len(library))]

        return library

    def _initialize_learner(self):
        """
        Initializes a range of surrogate models
        """
        running_mode = self.settings.additional[_SALE.MODEL]
        if running_mode == "gpr":
            learner = BayesianOptimizer(
                estimator=GaussianProcessRegressor(
                    kernel=DotProduct(), normalize_y=True
                ),
                query_strategy=greedy_acquisition,
            )
        elif running_mode == "random_forest":
            learner = ActiveLearner(
                estimator=RandomForestRegressor(n_estimators=100),
                query_strategy=greedy_acquisition,
            )
        elif running_mode == "ffnn":
            device = "cuda" if torch.cuda.is_available() else "cpu"
            regressor = NeuralNetRegressor(
                FeedForwardNet,
                criterion=nn.MSELoss,
                optimizer=torch.optim.Adam,
                train_split=None,
                verbose=0,
                device=device,
                max_epochs=100,
                batch_size=1024,
            )
            learner = ActiveLearner(
                estimator=regressor,
                query_strategy=greedy_acquisition,
            )
        else:
            raise KeyError(f"running mode: {running_mode} not supported")
        return learner

    def _initialize_oracle(
        self, compound_list: List[pd.Series] = None, work_dir: str = None
    ) -> WorkFlow:
        """
        Initialize a workflow object with the attached steps initialized
        """
        # list of step configs
        base_oracle_config = self.settings.additional["oracle_config"]
        wf_config = {
            # inherit header settings from the parent workflow
            _WE.HEADER: self.get_workflow_object().header,
            _WE.STEPS: [],
        }
        oracle_wf = WorkFlow(**wf_config)
        oracle_steps = []
        for step in base_oracle_config:
            step = self._initialize_oracle_step_from_dict(step)

            step.set_workflow_object(oracle_wf)
            oracle_steps.append(step)

        if compound_list is not None:

            # manually attach the compound objects to the oracle's lead step
            # subsequent steps should take their input from the the previous step, as ususal.
            for idx, compound in enumerate(compound_list):
                cmp = Compound(name=str(idx), compound_number=idx)
                cmp.add_enumeration(
                    Enumeration(
                        compound_object=cmp,
                        smile=compound[_SALE.SMILES],
                        original_smile=compound[_SALE.SMILES],
                        molecule=compound[_SALE.MOLECULE],
                    )
                )
                oracle_steps[0].data.compounds.append(cmp)
            self._logger.log(
                f"first step loaded with {len(oracle_steps[0].data.compounds)} compounds",
                _LE.DEBUG,
            )
        for step in oracle_steps:
            oracle_wf.add_step(step)
        return oracle_wf

    def _run_oracle_wf(
        self, oracle_wf: WorkFlow, skip_init_input: bool = False, work_dir: str = None
    ):
        for idx, step in enumerate(oracle_wf._initialized_steps):
            # only write initial input
            if skip_init_input and idx == 0:
                print("Skipping generating input")
            else:

                # input has been generated for the lead step from the virtual lib
                step.generate_input()
            self._logger.log(
                f"Starting execution of oracle step: {step.step_id}", _LE.INFO
            )
            if work_dir is not None:
                step.work_dir = work_dir
            step.execute()
            self._logger.log(
                f"Processing write-out blocks for {step.step_id}.", _LE.DEBUG
            )
            step.process_write_out()

        return oracle_wf

    def query_oracle(
        self,
        compound_list: List[pd.Series],
        oracle_type: str = "docking",
        fragment_lib: pd.DataFrame = None,
    ) -> List[Compound]:
        """
        Interface function with the oracle method
        """
        if oracle_type == "docking":
            oracle_wf = self._initialize_oracle(compound_list)
            # we have a fully initialized step with the compounds loaded.  Execute them
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf, skip_init_input=True)

            # retrieve compounds from the final step
            final_compounds = oracle_wf._initialized_steps[-1].data.compounds
            return final_compounds
        elif oracle_type == "FEP":
            self._logger.log("querying FEP oracle", _LE.DEBUG)
            # do not pass scores, generate a tmpdir, create the NCAA database, create the mutations file and run the FEP job on AWS
            # create tmpdir
            orig_dir = os.getcwd()
            tmp_dir = tempfile.mkdtemp()
            print(tmp_dir)
            os.chdir(tmp_dir)

            # extract the relevant amino acids using compound indices from the ncaa library
            compound_index = [row.IDX for row in compound_list]
            # retrieve the fragment using the index of the enumerated compound
            frags = [fragment_lib.iloc[idx] for idx in compound_index]
            idx = 0
            start_letter = "A"
            line_stub = self.get_additional_setting("mut_res")
            # create mutations file simultaneously
            with open("mutations.txt", "w") as f:
                for frag in frags:
                    frag_idx = str(idx).zfill(2)
                    elem_count = {}
                    with Chem.PDBWriter(f"{start_letter}{frag_idx}.pdb") as writer:
                        mol = frag.Molecule
                        mi = Chem.AtomPDBResidueInfo()
                        for a in mol.GetAtoms():
                            n = elem_count.get(a.GetAtomicNum(), 0)
                            n += 1
                            elem_count[a.GetAtomicNum()] = n
                            n = min(n, 99)
                            mi = Chem.AtomPDBResidueInfo()
                            mi.SetResidueName(f"{start_letter}{frag_idx}")
                            mi.SetResidueNumber(1)
                            atom_name = "{0:>2s}{1:<2d}".format(a.GetSymbol(), n)
                            mi.SetIsHeteroAtom(True)
                            mi.SetName(atom_name)
                            mi.SetOccupancy(0.0)
                            mi.SetTempFactor(0.0)
                            a.SetMonomerInfo(mi)
                        print(Chem.rdPartialCharges.ComputeGasteigerCharges(mol))
                        # mi.SetResidueName(f"{start_letter}{frag_idx}")
                        # [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
                        # rename the molecule
                        mol.SetProp(_WOE.RDKIT_NAME, frag_idx)

                        writer.write(frag.Molecule)

                    f.write(f"{line_stub}->{start_letter}{frag_idx}\n")
                    idx += 1
                    # not bullet proof, but gives 200 unique residue IDs, which should be enough to be getting on with
                    if idx > 99:
                        start_letter = "B"

            file_list = [f for f in os.listdir(tmp_dir) if f.endswith("pdb")]
            # need to add the header line to each pdb
            for pdb_file in file_list:
                with open(pdb_file, "r") as f:
                    lines = f.readlines()
                lines = [
                    f"HEADER    LIGAND                                  06-JAN-06   {pdb_file.split('.')[0]}\n"
                ] + lines
                with open(pdb_file, "w") as f:
                    for line in lines:
                        f.write(line)

            # # now convert sdf to mae
            converter = StructConvert(prefix_execution="module load schrodinger")
            for pdb in file_list:
                converter.pdb2mae(pdb, f"{pdb.split('.')[0]}.mae")

            # now concatenate
            concatenator = StructcatUtil(
                backend="structcat", prefix_execution="module load schrodinger"
            )
            concatenator.concatenate(
                [f"{pdb.split('.')[0]}.mae" for pdb in file_list], "database_nsr.mae"
            )

            # need to add the residue information to the resulting mae file

            # tmpdir is prepared, now initialize the FEP+ step, no need to prepare input, this is constant for FEP
            oracle_wf = self._initialize_oracle()
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf, work_dir=tmp_dir)

            # now parse the log file from the fep step

            os.chdir(orig_dir)

    def _extract_final_scores(
        self, compounds: List[Compound], criteria: str, highest_is_best: bool = False
    ) -> np.ndarray:
        """
        Takes a list of compound objects from the oracle and extracts the best score based on the provided criteria
        """
        top_scores = []
        for comp in compounds:
            scores = []
            for enum in comp.get_enumerations():
                for conf in enum.get_conformers():
                    scores.append(float(conf._conformer.GetProp(criteria)))

            # if docking generated no conformers
            if not scores:
                scores.append(0.0)

            best_score = max(scores) if highest_is_best else min(scores)
            top_scores.append(best_score)

        return np.absolute(top_scores, dtype=np.float32)

    def check_additional(self, key, val=True) -> bool:

        if (
            key in self.settings.additional.keys()
            and self.settings.additional[key] == val
        ):
            return True
        return False
