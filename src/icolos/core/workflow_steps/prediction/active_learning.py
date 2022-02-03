from typing import List, Tuple
import os
import random
import pickle

from modAL.acquisition import max_EI, max_UCB
from modAL.uncertainty import uncertainty_sampling
from modAL.models.learners import BayesianOptimizer, ActiveLearner
from pydantic.main import BaseModel

from sklearn.gaussian_process.kernels import DotProduct
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.decomposition import PCA
from sklearn.exceptions import NotFittedError
import matplotlib.pyplot as plt

from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.step import StepBase
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import (
    StepBaseEnum,
    StepGlideEnum,
    StepActiveLearningEnum,
)

from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import PandasTools, Mol
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np
from sklearn.metrics import mean_squared_error, confusion_matrix
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum

from icolos.utils.general.convenience_functions import nested_get


_SGE = StepGlideEnum()
_SALE = StepActiveLearningEnum()
_IE = StepInitializationEnum()


class StepActiveLearning(StepBase, BaseModel):
    """
    Class to run an active learning framework
    Primarily designed for building QSAR models using a physics based method (embedding + docking) as an oracle

    Takes the step conf for the oracle as an additional argument.  The step with these settings is run with the queried compounds at each stage of the active learning loop
    """

    _pca: PCA = PCA()

    def __init__(self, **data):
        super().__init__(**data)

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
        # initialize the basic oracle, load the query compounds for evaluation
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

    # def _reverse_sigmoid(self, x: float) -> float:
    #     """
    #     Scales compounds in range [-14,0] to be [1,0]
    #     """
    #     return 1.0 / (1 + np.e ** (0.45 * x + 4))

    def _extract_final_scores(
        self, compounds: List[Compound], criteria: str, highest_is_best: bool = False
    ) -> List[float]:
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

        return list(np.absolute(top_scores))

    def _generate_library(self) -> DataFrame:
        """
        Loads the library file from disk
        This should be a .sdf file with the pre-embedded compounds from a library enumeration or such
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
                GetMorganFingerprintAsBitVect(x[_SALE.MOLECULE], 2, nBits=2048)
            ),
            axis=1,
        )

        return library

    # def _prepare_initial_data(
    #     self, lib: pd.DataFrame
    # ) -> Tuple[np.ndarray, List[float]]:
    #     initial_compound_idx = random.sample(
    #         range(len(lib)), int(self.settings.additional[_SALE.INIT_SAMPLES])
    #     )
    #     data_rows = [lib.iloc[idx] for idx in initial_compound_idx]
    #     compounds = np.array([row[_SALE.MORGAN_FP] for row in data_rows])
    #     # return annotated compound list
    #     self._logger.log("Computing initial datapoints", _LE.INFO)
    #     annotated_compounds = self.query_oracle(data_rows)

    #     activities = self._extract_final_scores(
    #         annotated_compounds, criteria=_SGE.GLIDE_DOCKING_SCORE
    #     )
    #     self._logger.log(f"initial data points {activities}", _LE.DEBUG)

    #     return compounds, activities

    def _prepare_validation_data(self) -> Tuple[list[float], List[float], pd.DataFrame]:
        """
        parses sdf file with results to dataframe, extract fingerprints + results
        """
        val_lib = PandasTools.LoadSDF(
            self.settings.additional[_SALE.VALIDATION_LIB],
            smilesName=_SALE.SMILES,
            molColName=_SALE.MOLECULE,
            includeFingerprints=True,
            removeHs=False,
            embedProps=True,
        )
        val_lib[_SALE.MORGAN_FP] = val_lib.apply(
            lambda x: np.array(
                GetMorganFingerprintAsBitVect(x[_SALE.MOLECULE], 2, nBits=2048)
            ),
            axis=1,
        )
        scores = list(
            np.absolute(
                pd.to_numeric(
                    val_lib[self.settings.additional[_SALE.CRITERIA]].fillna(0)
                )
            )
        )
        scores = [float(x) for x in scores]

        return list(val_lib[_SALE.MORGAN_FP]), scores, val_lib

    def greedy_acquisition(
        self,
        estimator: RandomForestRegressor,
        X: np.ndarray,
        n_instances: int,
        highest_is_best: bool = False,
    ) -> np.ndarray:
        """
        Implement greedy acquisition strategy, return the n_samples best scores
        """
        # try:
        #     predictions = estimator.predict(X)
        # except:
        #     self._logger.log(
        #         "Estimator is not fitted, defaulting to random predictions", _LE.INFO
        #     )
        #     #     # if not initialized, generate random docking scores (absolute)
        predictions = np.random.uniform(12, 0, len(X))

        # zero those predictions we've seen before
        # for idx in previous_idx:
        #     predictions[idx] = 0
        # smaller before n_instances, largest after
        sorted_preds = np.argpartition(predictions, -n_instances, axis=0)[-n_instances:]
        return sorted_preds

    def execute(self):
        tmp_dir = self._make_tmpdir()
        lib = self._generate_library()

        (
            val_compounds,
            val_scores,
            val_lib,
        ) = self._prepare_validation_data()

        # fit tsne embedding
        self._pca.fit(val_compounds)

        running_mode = self.settings.additional["running_mode"]
        if running_mode == "bayes_opt":
            learner = BayesianOptimizer(
                estimator=GaussianProcessRegressor(
                    kernel=DotProduct(), normalize_y=True
                ),
                query_strategy=self.greedy_acquisition,
            )
        elif running_mode == "active_learning":
            learner = ActiveLearner(
                estimator=RandomForestRegressor(n_estimators=1000),
                query_strategy=self.greedy_acquisition,
            )
        else:
            raise KeyError(f"running mode: {running_mode} not supported")

            # generate baseline performance
        try:
            val_predictions = learner.predict(val_compounds)
            mse = mean_squared_error(val_scores, val_predictions)
            self._logger.log(f"Baseline val set rmsd: {np.sqrt(mse)}", _LE.INFO)
        except:
            pass

        if (
            "debug" in self.settings.additional.keys()
            and self.settings.additional["debug"] == True
        ):
            self._logger.log("Starting debug run...", _LE.DEBUG)
            self._run_learning_loop(
                learner=learner,
                lib=val_lib,
                val_compounds=val_compounds,
                val_scores=val_scores,
                debug=True,
                tmp_dir=tmp_dir,
            )
        else:
            self._run_learning_loop(
                learner=learner,
                lib=lib,
                val_compounds=val_compounds,
                val_scores=val_scores,
                debug=False,
                tmp_dir=tmp_dir,
            )

        # pickle the final model
        with open(os.path.join(tmp_dir, "model.pkl"), "wb") as f:
            pickle.dump(learner, f)

        self._parse_output(tmp_dir)

    def _run_learning_loop(
        self,
        learner,
        lib,
        val_compounds,
        val_scores,
        debug: bool = False,
        tmp_dir=None,
    ):
        rounds = self.settings.additional[_SALE.N_ROUNDS]
        n_instances = self.settings.additional[_SALE.BATCH_SIZE]
        fig, axs = plt.subplots(nrows=5, ncols=5, figsize=(40, 40), squeeze=True)
        axs = axs.ravel()

        for idx in range(rounds):
            query_idx, _ = learner.query(
                list(lib[_SALE.MORGAN_FP]),
                n_instances=n_instances,
            )
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
                # retrieve precalculated scores
                scores = [
                    float(lib.iloc[int(idx)][self.settings.additional[_SALE.CRITERIA]])
                    for idx in query_idx
                ]
                scores = list(np.absolute(scores))
                self._logger.log(f"Debug scores: {scores}", _LE.DEBUG)

            learner.teach(
                np.array([compound[_SALE.MORGAN_FP] for compound in query_compounds]),
                scores,
            )
            # get the predictions
            val_predictions = learner.predict(val_compounds)
            mse = mean_squared_error(val_scores, val_predictions)
            self._logger.log(
                f"Round {idx+1} Validation set rmsd: {np.sqrt(mse)}", _LE.INFO
            )
            self._logger.log(f"Predictions: \n{val_predictions[:5]}", _LE.INFO)
            self._logger.log(f"Actual: \n{val_scores[:5]}", _LE.INFO)

            # produce tsne embedding
            emb = self._pca.transform(learner.X_training)
            axs[idx].scatter(emb[0], emb[1])
        if tmp_dir is not None:
            fig.savefig(os.path.join(tmp_dir, "embeddings.png"))

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
