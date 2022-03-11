import os
import tempfile
from typing import List
from pydantic import BaseModel
import pandas as pd
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.step import StepBase
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepActiveLearningEnum, StepBaseEnum
from icolos.utils.general.convenience_functions import nested_get
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import Mol
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

    def _initialize_oracle(self, compound_list: List[pd.Series]) -> WorkFlow:
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

        # manually attach the compound objects to the oracle's lead step
        # subsequent steps should take their input from the the previous step, as ususal.
        for idx, compound in enumerate(compound_list):
            cmp = Compound(compound_number=idx)
            cmp.add_enumeration(
                Enumeration(
                    compound_object=cmp,
                    smile=compound[_SALE.SMILES],
                    original_smile=compound[_SALE.SMILES],
                    molecule=compound[_SALE.MOLECULE],
                )
            )
            oracle_steps[0].data.compounds.append(cmp)
        for step in oracle_steps:
            oracle_wf.add_step(step)

        return oracle_wf

    def query_oracle(self, compound_list: List[pd.Series]) -> List[Compound]:
        """
        Interface function with the oracle method, in the most likely case this is ligprep + docking

        Takes the requested compounds and runs them through the oracle workflow, returns the final compounds with annotations

        Notes:
        This could be an arbitrarily complex workflow, but the only thing that's going to change is the compounds.
        """
        oracle_wf = self._initialize_oracle(compound_list)
        # we have a fully initialized step with the compounds loaded.  Execute them
        for idx, step in enumerate(oracle_wf._initialized_steps):
            if idx != 0:
                # input has been generated for the lead step from the virtual lib
                step.generate_input()
            self._logger.log(
                f"Starting execution of oracle step: {step.step_id}", _LE.INFO
            )
            step.execute()
            self._logger.log(
                f"Processing write-out blocks for {step.step_id}.", _LE.DEBUG
            )
            step.process_write_out()
        self._logger.log(
            f"Execution of {len(self._initialized_steps)} steps completed.", _LE.INFO
        )

        # retrieve compounds from the final step
        final_compounds = oracle_wf.steps[-1].data.compounds
        return final_compounds

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

        return np.absolute(top_scores)

    def check_additional(self, key, val=True) -> bool:

        if (
            key in self.settings.additional.keys()
            and self.settings.additional[key] == val
        ):
            return True
        return False
