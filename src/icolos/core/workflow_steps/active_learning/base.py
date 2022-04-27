import json
import os
import shutil
import string
import tempfile
from typing import Callable, List
from pydantic import BaseModel
import pandas as pd
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.compound import Compound, Conformer, Enumeration
from icolos.core.workflow_steps.step import StepBase
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import StepActiveLearningEnum, StepBaseEnum
from icolos.utils.enums.write_out_enums import WriteOutEnum
from icolos.utils.general.convenience_functions import nested_get
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from rdkit import Chem
from rdkit.Chem.AllChem import (
    GetMorganFingerprintAsBitVect,
)
from sklearn.gaussian_process.kernels import DotProduct
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor
from icolos.core.workflow_steps.active_learning.al_utils import (
    expected_improvement,
    greedy_acquisition,
)
from modAL.models.learners import BayesianOptimizer, ActiveLearner
from icolos.core.workflow_steps.active_learning.models.ffnn import FeedForwardNet
from skorch.regressor import NeuralNetRegressor
from skorch.callbacks import EarlyStopping
from icolos.utils.execute_external.execute import Executor
import torch
from torch import nn

_IE = StepInitializationEnum()
_SALE = StepActiveLearningEnum()
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

    def construct_fingerprints(self, library: pd.DataFrame) -> pd.DataFrame:
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

    def _get_acquisition_function(self) -> Callable:
        """Return the acquisition function specified in the config

        :return Callable: Functino that executes the acquisition function
        """
        acq_fn = self._get_additional_setting(
            _SALE.ACQUISITION_FUNCTION, default=_SALE.GREEDY
        )
        if acq_fn.lower() == _SALE.GREEDY:
            return greedy_acquisition
        elif acq_fn.lower() == _SALE.EI:
            return expected_improvement
        else:
            raise ValueError(f"Acquisition method {acq_fn} is not supported!")

    def _initialize_learner(self) -> ActiveLearner:
        """Initialize the surrogate model specified in the config

        :raises KeyError: Unsopported model specified
        :return ActiveLearner: Initialized active learner, initialized with a surrogate model architecture and acquisition function
        """
        running_mode = self.settings.additional[_SALE.MODEL]
        acquisition_function = self._get_acquisition_function()
        if running_mode == "gpr":
            learner = BayesianOptimizer(
                estimator=GaussianProcessRegressor(
                    kernel=DotProduct(), normalize_y=True
                ),
                query_strategy=acquisition_function,
            )
        elif running_mode == "random_forest":
            learner = ActiveLearner(
                estimator=RandomForestRegressor(n_estimators=32),
                query_strategy=acquisition_function,
            )
        elif running_mode == "ffnn":
            device = "cuda" if torch.cuda.is_available() else "cpu"
            regressor = NeuralNetRegressor(
                FeedForwardNet,
                criterion=nn.MSELoss,
                optimizer=torch.optim.Adam,
                callbacks=[EarlyStopping(patience=5)],
                warm_start=False,
                verbose=2,
                lr=1e-4,
                device=device,
                max_epochs=500,
                batch_size=1024,
            )
            learner = ActiveLearner(
                estimator=regressor,
                query_strategy=acquisition_function,
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
        with Chem.SDWriter(os.path.join(self.work_dir, "compounds.sdf")) as writer:
            for idx, comp in enumerate(compound_list):
                mol = comp[_SALE.MOLECULE]
                try:
                    name = comp[_SALE.ID]
                    mol.SetProp("original_name", name)
                except KeyError:
                    pass
                mol.SetProp(_WOE.RDKIT_NAME, f"{idx}:0")
                mol.SetProp(_WOE.COMPOUND_NAME, f"{idx}:0")
                writer.write(mol)
        compound_dict = {
            "source": os.path.join(self.work_dir, "compounds.sdf"),
            "source_type": "file",
            "format": "SDF",
        }
        try:
            # if other compound are already specified, add another entry to the list
            wf_config["workflow"]["steps"][0]["input"]["compounds"].append(
                compound_dict
            )
        except KeyError:
            # if no input block in the oracle template, add a blank one first,
            wf_config["workflow"]["steps"][0]["input"] = {}
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
        round: int = 0,
    ) -> np.ndarray:
        """Main interface method with the oracle, controls initialization and executino

        :param List[pd.Series] compound_list: list of rows from library df to be queried by the oracle
        :param str oracle_type: type of oracle to initialize, defaults to "docking"
        :param pd.DataFrame fragment_lib: optional fragment library for fep and resi scanning oracles, defaults to None
        :raises NotImplementedError: Handles unknown oracle types
        :return np.ndarray: return array of scores for each compound from the oracle
        """
        criteria = self._get_additional_setting(_SALE.CRITERIA)
        if oracle_type == "pmx_rbfe":
            # test code, append the right tags to the oracle compounds

            # output_compounds = []
            # for i in range(len(compound_list)):
            #     # simulate some dropped conformers
            #     if i not in (2, 4):
            #         parent_compound = Compound(name=f"{i}", compound_number=i)
            #         parent_enum = Enumeration(enumeration_id=0)
            #         parent_enum.set_compound_object(parent_compound)
            #         conf = Conformer(Chem.MolFromSmiles("CCCC"), conformer_id=0)
            #         conf.get_molecule().SetProp(criteria, str(np.random.uniform(-5, 5)))
            #         conf.set_enumeration_object(parent_enum)
            #         output_compounds.append(conf)

            oracle_wf = self._initialize_oracle_workflow(compound_list)
            # # we have a fully initialized step with the compounds loaded.  Execute them
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf)

            output_compounds = [
                n.conformer
                for n in oracle_wf._initialized_steps[-1]
                .get_perturbation_map()
                .get_nodes()
            ]
            original_compounds = [
                mol
                for mol in Chem.SDMolSupplier(
                    os.path.join(self.work_dir, "compounds.sdf")
                )
            ]
            final_compounds = []
            for comp in original_compounds:
                # look through the final compounds and find the one with the right comp/enumeration
                comp_name = comp.GetProp(_WOE.RDKIT_NAME)
                match_found = False
                for out_comp in output_compounds:
                    if (
                        comp_name
                        == out_comp.get_enumeration_object().get_index_string()
                    ):
                        final_compounds.append(out_comp)
                        self._logger.log(
                            f"matched output conformer for mol {comp_name}", _LE.DEBUG
                        )
                        match_found = True
                if not match_found:
                    # the enumeration was dropped in the workflow, append the start compound with no tags, will be zeroed
                    final_compounds.append(Conformer(conformer=comp))
                    self._logger.log(
                        f"No conformer belonging to enumeration {comp_name} found!",
                        _LE.WARNING,
                    )

            # TODO: really this should be coupled to the oracle type and should not change
            final_scores = []
            for conf in final_compounds:
                try:
                    final_scores.append(float(conf.get_molecule().GetProp(criteria)))
                except KeyError:
                    self._logger.log(
                        f"FEP score was not attached to conformer!", _LE.WARNING
                    )
                    final_scores.append(0.0)
            # move the old dir to backup
            try:
                shutil.move("output", f"output_{round}")
            except Exception as e:
                print("Could not back up output files!, error was ", e)
            return np.array(final_scores, dtype=np.float32)
        elif oracle_type == "docking":

            oracle_wf = self._initialize_oracle_workflow(compound_list)
            # we have a fully initialized step with the compounds loaded.  Execute them
            oracle_wf = self._run_oracle_wf(oracle_wf=oracle_wf)

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
            os.chdir(orig_dir)

        else:
            raise NotImplementedError(f"Oracle type {oracle_type} not implemented")
        self._logger.log("Extracting final scores", _LE.DEBUG)

        scores = self._extract_final_scores_from_compounds(
            final_compounds, self.settings.additional[_SALE.CRITERIA]
        )
        return scores

    def _extract_final_scores_from_compounds(
        self, compounds: List[Compound], criteria: str, highest_is_best: bool = False
    ) -> np.ndarray:
        """Extract final scores from compounds

        :param List[Compound] compounds: compounds extracted from final step in workflow
        :param str criteria: tag to extract score by
        :param bool highest_is_best: control whether highest vals correspond to best scores, negative for docking scores, affinities etc, defaults to False
        :return np.ndarray: array containing final scores per compound
        """
        top_scores = []
        for comp in compounds:
            # extract top score per compound
            scores = []
            for enum in comp.get_enumerations():
                # if conformers attached ,extract from there
                if enum.get_conformers():
                    for conf in enum.get_conformers():
                        try:
                            conf_score = conf.get_molecule().GetProp(criteria)
                            scores.append(float(conf_score))
                        except KeyError:
                            scores.append(0.0)
                else:
                    scores = [0.0]

            best_score = max(scores) if highest_is_best else min(scores)
            top_scores.append(best_score)

        return np.array(top_scores, dtype=np.float32)
