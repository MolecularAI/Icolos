import json
from rdkit import Chem
import os
import string
import tempfile
from typing import List
from pydantic import BaseModel
import pandas as pd
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.compound import Compound, Conformer, Enumeration
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
from icolos.utils.execute_external.execute import Executor
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
        """Initialize an oracle step, note this is similar to the main workflow initialization method, but needs to be a separate initializer for the AL workflow

        :param dict step_conf: config for the step to be initialized
        :raises ValueError: Step ID must be recognised as an internal step
        :return StepBase: Initialized step
        """
        _STE = StepBaseEnum
        step_type = nested_get(step_conf, _STE.STEP_TYPE, default=None)
        step_type = None if step_type is None else step_type.upper()
        if step_type in _IE.STEP_INIT_DICT.keys():
            return _IE.STEP_INIT_DICT[step_type](**step_conf)
        else:
            raise ValueError(
                f"Backend for step {nested_get(step_conf, _STE.STEPID, '')} unknown."
            )

    def construct_fingerprints(self, library: pd.DataFrame) -> pd.DataFrameame:
        """Add morgan FP column to dataframe containing rdkit mols

        :param pd.DataFrame library: lib containing the mols in column named Molecules
        :return pd.DataFrame: Modified dataframe with the extra fingerprints col
        """
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

    def _initialize_learner(self) -> ActiveLearner:
        """Initialize the surrogate model specified in the config

        :raises KeyError: Unsopported model specified
        :return ActiveLearner: Initialized active learner, initialized with a surrogate model architecture and acquisition function
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

    def _initialize_oracle_workflow(self, compound_list: List[pd.Series]) -> WorkFlow:
        """Initialize the oracle workflow

        :param List[pd.Series] compound_list: List of rows from the library df containing comounds to query
        :return WorkFlow: Initialized Icolos workflow
        """
        # path to the json config for the oracle config
        oracle_conf = self.settings.additional["oracle_config"]
        with open(oracle_conf, "r") as f:
            wf_config = json.load(f)

        # manually attach the compound objects to the oracle's lead step
        with Chem.SDWriter("compounds.sdf") as writer:
            for idx, comp in enumerate(compound_list):
                mol = comp[_SALE.MOLECULE]
                mol.SetProp(_WOE.RDKIT_NAME, f"{idx}:0")
                mol.SetProp(_WOE.COMPOUND_NAME, f"{idx}:0")
                writer.write(mol)

        compound_dict = {
            "source": "compounds.sdf",
            "source_type": "file",
            "format": "SDF",
        }
        wf_config["workflow"]["steps"][0]["input"]["compounds"] = [compound_dict]
        # inherit header from main workflow
        header = self.get_workflow_object().header
        oracle_wf = WorkFlow(**wf_config["workflow"])
        oracle_wf.header = header
        oracle_steps = []
        for step_conf in oracle_wf.steps:
            step_conf = oracle_wf._update_global_variables(conf=step_conf)
            st = self._initialize_oracle_step_from_dict(step_conf)

            st.set_workflow_object(oracle_wf)
            oracle_steps.append(st)

        for step in oracle_steps:
            oracle_wf.add_step(step)
        return oracle_wf

    def _run_oracle_wf(self, oracle_wf: WorkFlow, work_dir: str = None) -> WorkFlow:
        """Run the oracle workflow through for one set of compounds

        :param WorkFlow oracle_wf: initialized oracle workflow to be run
        :param str work_dir: _description_, defaults to None
        :return WorkFlow: Return completed workflow with data attached
        """
        for idx, step in enumerate(oracle_wf._initialized_steps):
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
    ) -> np.ndarray:
        """Main interface method with the oracle, controls initialization and executino

        :param List[pd.Series] compound_list: list of rows from library df to be queried by the oracle
        :param str oracle_type: type of oracle to initialize, defaults to "docking"
        :param pd.DataFrame fragment_lib: optional fragment library for fep and resi scanning oracles, defaults to None
        :raises NotImplementedError: Handles unknown oracle types
        :return np.ndarray: return array of scores for each compound from the oracle
        """
        if oracle_type == "pmx_rbfe":
            # TODO: with the pmx oracle, I think it only makes sense to use star maps, then we can take the ddG values from the hub compounds vs some reference
            oracle_wf = self._initialize_oracle_workflow(compound_list)
            # we have a fully initialized step with the compounds loaded.  Execute them
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf, skip_init_input=True)

            final_compounds = [
                n.conformer
                for n in oracle_wf._initialized_steps[-1]
                .get_perturbation_map()
                .get_nodes()
            ]

        elif oracle_type == "docking":

            oracle_wf = self._initialize_oracle_workflow(compound_list)
            # we have a fully initialized step with the compounds loaded.  Execute them
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf, skip_init_input=True)

            final_compounds = oracle_wf._initialized_steps[-1].data.compounds

        elif oracle_type in ("protein_FEP", "residue_scanning"):
            # do not pass scores, generate a tmpdir, create the NCAA database, create the mutations file and run the FEP job on AWS
            # create tmpdir
            orig_dir = os.getcwd()
            tmp_dir = tempfile.mkdtemp()
            os.chdir(tmp_dir)

            # extract the relevant amino acids using compound indices from the ncaa library
            compound_index = [row.IDX for row in compound_list]
            # retrieve the fragment using the index of the enumerated compound
            frags = [fragment_lib.iloc[idx] for idx in compound_index]
            letter_strings = string.ascii_uppercase
            all_letters = []
            for char1 in letter_strings:
                for char2 in letter_strings:
                    for char3 in letter_strings:
                        all_letters.append(f"{char1}{char2}{char3}")
            line_stub = self._get_additional_setting("mut_res")
            # create mutations file simultaneously
            with open("mutations.txt", "w") as f, Chem.SDWriter(
                "database.sdf"
            ) as writer:
                for idx, frag in enumerate(frags):
                    mol = frag.Molecule
                    # rename the molecule
                    mol.SetProp(_WOE.RDKIT_NAME, f"{all_letters[idx]}")
                    writer.write(frag.Molecule)
                    f.write(f"{line_stub}->{all_letters[idx]}\n")
                    idx += 1

            command = "$SCHRODINGER/run python3 $NSR_LIB_SCRIPT database.sdf"

            executor = Executor(prefix_execution="ml schrodinger")
            executor.execute(command, arguments=[], check=True, location=tmp_dir)
            # now ncaa_nca.maegz will be in the tmpdir

            # execution is the same for both cases, just different oracle
            oracle_wf = self._initialize_oracle_workflow(compound_list=compound_list)
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf, work_dir=tmp_dir)
            final_compounds = oracle_wf._initialized_steps[-1].data.compounds
            print(final_compounds)
            os.chdir(orig_dir)

        else:
            raise NotImplementedError(f"Oracle type {oracle_type} not implemented")
        self._logger.log("Extracting final scores", _LE.DEBUG)

        scores = self._extract_final_scores(
            final_compounds, self.settings.additional[_SALE.CRITERIA]
        )
        return scores

    def _extract_final_scores(
        self, compounds: List[Compound], criteria: str, highest_is_best: bool = False
    ) -> np.ndarray:
        """Extract final scores from compounds

        :param List[Compound] compounds: compounds extracted from final step in workflow
        :param str criteria: tag to extract score by
        :param bool highest_is_best: control whether highest vals correspond to best scores, negative for docking scores, affinities etc, defaults to False
        :return np.ndarray: array containing final scores per compound
        """
        top_scores = []

        if isinstance(compounds[0], Compound):
            for comp in compounds:
                scores = []
                for enum in comp.get_enumerations():
                    scores.append(float(enum.get_molecule().GetProp(criteria)))

                # if docking generated no conformers
                if not scores:
                    scores.append(0.0)

                best_score = max(scores) if highest_is_best else min(scores)
                top_scores.append(best_score)
        elif isinstance(compounds[0], Conformer):
            for conf in compounds:
                scores.append(float(conf._conformer.GetProp(criteria)))

        return np.array(top_scores, dtype=np.float32)
