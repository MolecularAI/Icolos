from typing import List
import os
import random
import pickle

from modAL.acquisition import max_EI
from modAL.models.learners import BayesianOptimizer
from pydantic.main import BaseModel

from sklearn.gaussian_process.kernels import WhiteKernel, RBF
from sklearn.gaussian_process import GaussianProcessRegressor

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
from sklearn.metrics import mean_squared_error
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
            # we probably want to filter these before the model sees them
            if not scores:
                scores.append(0.0)

            best_score = max(scores) if highest_is_best else min(scores)
            top_scores.append(best_score)

        return top_scores

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

    def _prepare_initial_data(self, lib: pd.DataFrame):
        initial_compound_idx = random.sample(
            range(len(lib)), int(self.settings.additional[_SALE.INIT_SAMPLES])
        )
        data_rows = [lib.iloc[idx] for idx in initial_compound_idx]
        # return annotated compound list
        annotated_compounds = self.query_oracle(data_rows)

        # extract top score per compound
        init_scores: List[float] = self._extract_final_scores(
            annotated_compounds, criteria=_SGE.GLIDE_DOCKING_SCORE
        )
        init_compounds = np.array([row[_SALE.MORGAN_FP] for row in data_rows])

        return init_compounds, init_scores

    def _prepare_validation_data(self):
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
        # need the morgan fingerprints in the df
        val_lib[_SALE.MORGAN_FP] = val_lib.apply(
            lambda x: np.array(
                GetMorganFingerprintAsBitVect(x[_SALE.MOLECULE], 2, nBits=2048)
            ),
            axis=1,
        )
        scores = list(
            pd.to_numeric(val_lib[self.settings.additional[_SALE.CRITERIA]].fillna(0))
        )
        scores = [float(x) for x in scores]
        return list(val_lib[_SALE.MORGAN_FP]), scores

    def _filter_oracle_results(
        self, compound_rows: List[pd.Series], scores: List[float]
    ):
        final_compounds, final_scores = [], []
        for cmp, score in zip(compound_rows, scores):
            if score != 0.0:
                final_compounds.append(cmp)
                final_scores.append(score)

        return final_compounds, final_scores

    def execute(self):
        tmp_dir = self._make_tmpdir()

        # TODO: Implement comittee model

        # start with sdf of pre-calculatd ligand embeddings for each full peptide in the library
        lib = self._generate_library()
        init_compounds, init_scores = self._prepare_initial_data(lib)
        # load validation set for later
        validation_compounds, validation_scores = self._prepare_validation_data()

        kernel = RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(
            noise_level=1, noise_level_bounds=(1e-10, 1e2)
        )
        learner = BayesianOptimizer(
            # estimator=GaussianProcessRegressor(kernel=kernel),
            estimator=GaussianProcessRegressor(kernel),
            query_strategy=max_EI,
            X_training=init_compounds,
            y_training=init_scores,
        )

        for idx in range(int(self.settings.additional[_SALE.N_ROUNDS])):
            # generate the requested points from the learner
            query_idx, _ = learner.query(
                list(lib[_SALE.MORGAN_FP]),
                n_instances=int(self.settings.additional[_SALE.BATCH_SIZE]),
            )
            # generate oracle input
            query_compounds = [lib.iloc[int(idx)] for idx in query_idx]
            # query oracle

            compounds = self.query_oracle(query_compounds)
            scores = self._extract_final_scores(
                compounds, self.settings.additional[_SALE.CRITERIA]
            )
            # some of the scores will be zero if they didn't dock, do we want to filter these out, only hand back those compounds with a non-zero score?
            query_compounds, scores = self._filter_oracle_results(
                query_compounds, scores
            )

            learner.teach(
                np.array([compound[_SALE.MORGAN_FP] for compound in query_compounds]),
                scores,
            )
            # need a held-out test set with docking scores already computed
            performance = learner.score(validation_compounds, validation_scores)
            self._logger.log(
                f"Round {idx +1}; val set correlation: {performance}", _LE.INFO
            )
            # get the predictions
            predictions = learner.predict(validation_compounds)
            mse = mean_squared_error(validation_scores, predictions)
            self._logger.log(f"Round {idx+1}; rmse: {np.sqrt(mse)}", _LE.INFO)

        # pickle the final model
        with open(os.path.join(tmp_dir, "model.pkl"), "wb") as f:
            pickle.dump(learner, f)

        self._parse_output(tmp_dir)

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
