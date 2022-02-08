from typing import List, Tuple
import os
import pickle

from modAL.models.learners import BayesianOptimizer, ActiveLearner
from pydantic.main import BaseModel

from sklearn.gaussian_process.kernels import DotProduct
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor

try:
    import torch
    from torch import nn
except:
    print("Warning: PyTorch is not installed in the environment!")
try:
    from skorch.regressor import NeuralNetRegressor
except:
    print("Warning: skorch is not installed in the environment!")


from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.active_learning.models.ffnn import FeedForwardNet
from icolos.core.workflow_steps.step import StepBase
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import (
    StepBaseEnum,
    StepActiveLearningEnum,
)

from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import PandasTools, Mol
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np


_SALE = StepActiveLearningEnum()


class StepActiveLearning(StepBase, BaseModel):
    """
    Class to run an active learning framework
    Primarily designed for building QSAR models using a physics based method (embedding + docking) as an oracle

    Takes the step conf for the oracle as an additional argument.  The step with these settings is run with the queried compounds at each stage of the active learning loop
    """

    def __init__(self, **data):
        super().__init__(**data)

    def check_additional(self, key, val=True) -> bool:

        if (
            key in self.settings.additional.keys()
            and self.settings.additional[key] == val
        ):
            return True
        return False

    def _initialize_oracle(self, compound_list: List[pd.Series]) -> List[StepBase]:
        # list of step configs
        base_oracle_config = self.settings.additional["oracle_config"]
        oracle_steps = []
        for step in base_oracle_config:
            oracle_steps.append(self._initialize_oracle_step_from_dict(step))

        # manually attach the compound objects to the oracle's lead step
        # subsequent steps should take their input from the the first step.
        for idx, compound in enumerate(compound_list):
            cmp = Compound(compound_number=idx)
            cmp.add_enumeration(
                Enumeration(
                    compound_object=cmp,
                    smile=compound[_SALE.SMILES],
                    molecule=compound[_SALE.MOLECULE],
                )
            )
            oracle_steps[0].data.compounds.append(cmp)

        return oracle_steps

    def query_oracle(self, compound_list: List[Mol]) -> List:
        """
        Interface function with the oracle method, in the most likely case this is ligprep + docking

        Takes the requested compounds and runs them through the oracle workflow, returns the final compounds with annotations

        Notes:
        This could be an arbitrarily complex workflow, but the only thing that's going to change is the compounds.
        """
        oracle_steps = self._initialize_oracle(compound_list)
        # we have a fully initialized step with the compounds loaded.  Execute them
        for idx, step in enumerate(oracle_steps):
            # for subsequent steps we will need to read in from the previous one
            if idx != 0:
                step.generate_input()
            step.execute()
            step.process_write_out()

        # retrieve compounds from the final step
        final_compounds = oracle_steps[-1].data.compounds
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

    def _parse_library(
        self, extract_property: bool = True
    ) -> Tuple[DataFrame, np.ndarray]:
        """
        Loads the library file from disk
        This should be a .sdf file containing the compounds to be screened
        """
        lib_path = self.settings.additional[_SALE.VIRTUAL_LIB]
        assert lib_path.endswith(".sdf")

        # hold the lib in a pandas df
        library = PandasTools.LoadSDF(
            lib_path,
            smilesName=_SALE.SMILES,
            molColName=_SALE.MOLECULE,
            includeFingerprints=True,
            removeHs=False,
            embedProps=True,
        )
        # need the morgan fingerprints in the df
        library[_SALE.MORGAN_FP] = library.apply(
            lambda x: np.array(
                GetMorganFingerprintAsBitVect(x[_SALE.MOLECULE], 2, nBits=2048),
                dtype=np.float32,
            ),
            axis=1,
        )
        scores = (
            np.absolute(
                pd.to_numeric(
                    library[self.settings.additional[_SALE.CRITERIA]].fillna(0)
                )
            )
            if extract_property
            else []
        )

        return library, scores

    def _initialize_learner(self):
        running_mode = self.settings.additional["running_mode"]
        if running_mode == "gpr":
            learner = BayesianOptimizer(
                estimator=GaussianProcessRegressor(
                    kernel=DotProduct(), normalize_y=True
                ),
                query_strategy=self.greedy_acquisition,
            )
        elif running_mode == "random_forest":
            learner = ActiveLearner(
                estimator=RandomForestRegressor(n_estimators=100),
                query_strategy=self.greedy_acquisition,
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
                batch_size=self.settings.additional[_SALE.BATCH_SIZE],
            )
            learner = ActiveLearner(
                estimator=regressor,
                query_strategy=self.greedy_acquisition,
            )
        else:
            raise KeyError(f"running mode: {running_mode} not supported")
        return learner

    def _run_learning_loop(
        self, learner, lib, debug: bool = False, top_1_idx: list = None
    ):
        rounds = self.settings.additional[_SALE.N_ROUNDS]
        n_instances = self.settings.additional[_SALE.BATCH_SIZE]
        queried_compound_idx = []
        fraction_top1_hits = []

        for idx in range(rounds):
            query_idx, _ = learner.query(
                np.array(lib[_SALE.MORGAN_FP]),
                n_instances=n_instances,
                previous_idx=queried_compound_idx,
            )

            queried_compound_idx += list(query_idx)
            query_compounds = [lib.iloc[int(idx)] for idx in query_idx]

            if not debug:
                self._logger.log(
                    f"Querying oracle with {len(query_compounds)} compounds", _LE.INFO
                )

                compounds = self.query_oracle(query_compounds)
                scores = self._extract_final_scores(
                    compounds, self.settings.additional[_SALE.CRITERIA]
                )
            else:
                # retrieve precalculated scores instead of docking
                self._logger.log("Retrieving scores from precomputed data...", _LE.INFO)
                scores = np.absolute(
                    [
                        float(
                            lib.iloc[int(idx)][self.settings.additional[_SALE.CRITERIA]]
                        )
                        for idx in query_idx
                    ],
                    dtype=np.float32,
                )
                scores = scores.reshape(-1, 1)
            self._logger.log("Fitting with new data...", _LE.INFO)
            new_data = np.array(
                [compound[_SALE.MORGAN_FP] for compound in query_compounds],
                dtype=np.float32,
            )
            learner.teach(
                new_data,
                scores,
            )
            # calculate percentage of top-1% compounds queried
            if top_1_idx is not None:
                # what fraction of top1% compoudnds has the model requested to evaluate
                hits_queried = (
                    np.in1d(np.unique(np.array(queried_compound_idx)), top_1_idx).sum()
                    / len(top_1_idx)
                ) * 100
                self._logger.log(
                    f"Fraction of top 1% hits queries by round {idx}: {hits_queried}%",
                    _LE.INFO,
                )
                fraction_top1_hits.append(hits_queried)
        self._logger.log(f"Hits queried at each epoch:\n{fraction_top1_hits}", _LE.INFO)

    def execute(self):
        tmp_dir = self._make_tmpdir()
        lib, scores = self._parse_library(
            extract_property=self.settings.additional["evaluate"]
        )
        if scores is not None:
            top_1_percent = int(0.01 * len(scores))
            top_1_idx = np.argpartition(scores, -top_1_percent)[-top_1_percent:]

        learner = self._initialize_learner()

        # if we are using a pre-calculated dataset for evaluation
        if self.check_additional("debug"):
            # (val_lib, val_docking_scores) = self._parse_library(
            #     extract_property=True
            # )

            # top_1_val = int(0.01 * len(val_docking_scores))
            # # extract indices of top1% of compounds
            # top_1_val = np.argpartition(scores, -top_1_val)[-top_1_val:]
            self._logger.log("Starting debug run...", _LE.DEBUG)
            self._run_learning_loop(
                learner=learner, lib=lib, debug=True, top_1_idx=list(top_1_idx)
            )
        else:
            self._run_learning_loop(
                learner=learner,
                lib=lib,
                debug=False,
                top_1_idx=list(top_1_idx),
            )

        # pickle the final model
        with open(os.path.join(tmp_dir, "model.pkl"), "wb") as f:
            pickle.dump(learner, f)

        self._parse_output(tmp_dir)
