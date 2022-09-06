import pandas as pd
import numpy as np
import random
from scipy.stats import norm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from icolos.core.workflow_steps.active_learning.surrogate_model import SurrogateModel

from icolos.utils.enums.step_enums import StepActiveLearningEnum

_SALE = StepActiveLearningEnum()


class AcquisitionFunction:
    """class to execute all acquisition function heuristics for prospective REINVENT"""

    # TODO: assumption is that "lower scores are better" Add a variable to control this

    def __init__(
        self, surrogate: SurrogateModel, function: str, acquisition_batch_size: int
    ):
        # TODO: need to add a boolean for "highest is best"
        self.surrogate = surrogate
        self.function = function.lower()
        self.acquisition_batch_size = acquisition_batch_size

    def partition_compounds(
        self, original_smiles, fingerprints_pool, epsilon=0.2, beta=2, pi_epsilon=0.01
    ):
        """partitions REINVENT's batch of SMILES into acquired and non-acquired based on the acquisition function"""
        if self.function == _SALE.RANDOM:
            df = pd.DataFrame({"smiles": original_smiles})
            shuffled = df.sample(frac=1, axis=0)
            acquired_smiles = list(shuffled["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(
                shuffled["smiles"][self.acquisition_batch_size :]
            )

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.GREEDY:
            """greedy"""
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            df = self._get_sorted_df(
                smiles=original_smiles, predictions=predictions, ascending=True
            )
            acquired_smiles = list(df["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size :])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.GREEDY_EPSILON:
            """greedy-epsilon"""
            # TODO: have an "acquisition function specific parameters" block in the icolos JSON for "epsilon" parameter
            if random.random() <= epsilon:
                df = pd.DataFrame({"smiles": original_smiles})
                shuffled = df.sample(frac=1, axis=0)
                acquired_smiles = list(
                    shuffled["smiles"][: self.acquisition_batch_size]
                )
                non_acquired_smiles = list(
                    shuffled["smiles"][self.acquisition_batch_size :]
                )

                return acquired_smiles, non_acquired_smiles
            else:
                fps = self._get_morgan_fingerprints(smiles=original_smiles)
                predictions = list(self.surrogate.predict(fps))
                df = self._get_sorted_df(
                    smiles=original_smiles, predictions=predictions, ascending=True
                )
                acquired_smiles = list(df["smiles"][: self.acquisition_batch_size])
                non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size :])

                return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.PI:
            """probability of improvement"""
            # extract the current best
            try:
                with open("current_best.txt", "r") as f:
                    current_best = float(f.readlines()[0])
            except FileNotFoundError:
                current_best = 0.0

            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            standard_deviations = self.surrogate.get_std(fps)

            probabilities = np.zeros(len(fps))
            for idx, (pred, std) in enumerate(zip(predictions, standard_deviations)):
                # if there is no uncertainty, then PI = 1 or 0
                # depending on whether the prediction is better than the current best
                if int(std) == 0:
                    # TODO: assumes lower value is better
                    probabilities[idx] = 1 if pred < current_best else 0

                else:
                    probabilities[idx] = norm.cdf(
                        (pred - current_best + pi_epsilon) / std
                    )

            df = self._get_sorted_df(
                smiles=original_smiles, predictions=probabilities, ascending=False
            )
            acquired_smiles = list(df["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size :])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.EI:
            """expected improvement"""
            # extract the current best
            try:
                with open("current_best.txt", "r") as f:
                    current_best = float(f.readlines()[0])
            except FileNotFoundError:
                current_best = 0.0

            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            standard_deviations = self.surrogate.get_std(fps)

            expected_improvements = np.zeros(len(fps))
            for idx, (pred, std) in enumerate(zip(predictions, standard_deviations)):
                # if there is no uncertainty, then EI is the difference between the prediction and the current best
                if int(std) == 0:
                    # TODO: assumes lower value is better
                    expected_improvements[idx] = pred - current_best
                else:
                    numerator = pred - current_best + pi_epsilon
                    Z = numerator / std
                    expected_improvements[idx] = (numerator * norm.cdf(Z)) + (
                        std * norm.pdf(Z)
                    )

            df = self._get_sorted_df(
                smiles=original_smiles,
                predictions=expected_improvements,
                ascending=True,
            )
            acquired_smiles = list(df["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size :])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.TS:
            """thompson sampling"""
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            standard_deviations = self.surrogate.get_std(fps)
            gaussian_predictions = np.random.default_rng().normal(
                predictions, standard_deviations
            )

            df = self._get_sorted_df(
                smiles=original_smiles, predictions=gaussian_predictions, ascending=True
            )
            acquired_smiles = list(df["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size :])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.UCB:
            """upper-confidence bound"""
            # TODO: have an "acquisition function specific parameters" block in the icolos JSON for "beta" parameter
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            standard_deviations = self.surrogate.get_std(fps)
            ucb = predictions + beta * standard_deviations

            df = self._get_sorted_df(
                smiles=original_smiles, predictions=ucb, ascending=True
            )
            acquired_smiles = list(df["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size :])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.TANIMOTO:
            """Tanimoto. SMILES most dissimilar, i.e., lowest Tanimoto similarity to the pooled SMILES are acquired"""
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            similarities = np.zeros(len(fps))

            for idx, fp in enumerate(fps):
                similarities[idx] = np.mean(
                    DataStructs.BulkTanimotoSimilarity(fp, fingerprints_pool)
                )

            # sort by ascending Tanimoto similarities
            sorted_idx = np.argsort(similarities)
            acquired_idx = sorted_idx[: self.acquisition_batch_size]
            non_acquired_idx = sorted_idx[self.acquisition_batch_size :]

            acquired_smiles = [original_smiles[idx] for idx in acquired_idx]
            non_acquired_smiles = [original_smiles[idx] for idx in non_acquired_idx]

            return list(acquired_smiles), list(non_acquired_smiles)

        if self.function == _SALE.TANIMOTO_TRUE_BELIEVER:
            """
            Tanimoto. SMILES most dissimilar, i.e., lowest Tanimoto similarity to the pooled SMILES are acquired.
            To discourage acquiring clusters of SMILES, i.e., from the same chemical space, SMILES are acquired
            iteratively. After acquiring a SMILES, it is added to the "pool" and the same calculation is performed
            to acquire the next SMILES. This process is repeated until acquisition batch size SMILES is acquired.
            """
            from copy import deepcopy

            augmented_fps_pool = deepcopy(fingerprints_pool)

            def single_acquire(acquired_smiles, smiles_left, augmented_fps_pool):

                fps = self._get_morgan_fingerprints(smiles=smiles_left)
                similarities = np.zeros(len(fps))
                for idx, fp in enumerate(fps):
                    similarities[idx] = np.max(
                        DataStructs.BulkTanimotoSimilarity(fp, augmented_fps_pool)
                    )

                # sort by ascending Tanimoto similarities
                sorted_idx = np.argsort(similarities)
                acquire = smiles_left[sorted_idx[0]]
                acquired_smiles.append(acquire)
                smiles_left.remove(acquire)

                augmented_fps_pool.extend(self._get_morgan_fingerprints([acquire]))

                return acquired_smiles, smiles_left, augmented_fps_pool

            acquired_smiles, smiles_left, augmented_fps_pool = single_acquire(
                [], original_smiles, augmented_fps_pool
            )
            for iteration in range(self.acquisition_batch_size - 1):
                acquired_smiles, smiles_left, augmented_fps_pool = single_acquire(
                    acquired_smiles, smiles_left, augmented_fps_pool
                )

            non_acquired_smiles = smiles_left

            return list(acquired_smiles), list(non_acquired_smiles)

    def _get_morgan_fingerprints(self, smiles) -> list:
        """returns morgan fingerprints for a list of smiles"""
        molecules = [Chem.MolFromSmiles(s) for s in smiles]
        return [
            Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
            for mol in molecules
        ]

    def _get_sorted_df(self, smiles, predictions, ascending):
        """
        returns a DataFrame sorted by prediction values. Can be ascending or descending,
        as specified in the method call. The sorted DataFrame will be sliced into acquired
        and non-acquired SMILES
        """
        # TODO: need "highest_is_best" to be set in class instantiation
        df = pd.DataFrame({"smiles": smiles, "predictions": predictions})

        return df.sort_values(by="predictions", ascending=ascending)
