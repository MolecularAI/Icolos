import pandas as pd
import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from icolos.core.workflow_steps.active_learning.surrogate_model import SurrogateModel

from icolos.utils.enums.step_enums import StepActiveLearningEnum

_SALE = StepActiveLearningEnum()


class AcquisitionFunction:
    """class to execute all acquisition function heuristics for prospective REINVENT"""

    def __init__(
        self, surrogate: SurrogateModel, function: str, acquisition_batch_size: int
    ):
        # TODO: need to add a boolean for "highest is best"
        self.surrogate = surrogate
        self.function = function.lower()
        self.acquisition_batch_size = acquisition_batch_size

    def partition_compounds(self, original_smiles, fingerprints_pool, epsilon=0.2, beta=2):
        # TODO: missing acquisition functions requiring storage of the "current optimal"
        if self.function == _SALE.RANDOM:
            df = pd.DataFrame({"smiles": original_smiles})
            shuffled = df.sample(frac=1, axis=0)
            acquired_smiles = list(shuffled["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(shuffled["smiles"][self.acquisition_batch_size :])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.GREEDY:
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            df = pd.DataFrame({"smiles": original_smiles, "predictions": predictions})
            df = df.sort_values(by="predictions", ascending=True)
            acquired_smiles = list(df["smiles"][:self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size:])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.GREEDY_EPSILON:
            # TODO: have an "acquisition function specific parameters" block for epsilon and beta parameters
            if random.random() <= epsilon:
                self.function = _SALE.RANDOM
                self.partition_compounds(original_smiles=original_smiles, fingerprints_pool=fingerprints_pool)
            else:
                fps = self._get_morgan_fingerprints(smiles=original_smiles)
                predictions = list(self.surrogate.predict(fps))
                df = pd.DataFrame({"smiles": original_smiles, "predictions": predictions})
                df = df.sort_values(by="predictions", ascending=True)
                acquired_smiles = list(df["smiles"][:self.acquisition_batch_size])
                non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size:])

                return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.PI:
            pass

        if self.function == _SALE.EI:
            pass

        if self.function == _SALE.TS:
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            standard_deviations = self.surrogate.get_std(fps)
            gaussian_predictions = np.random.default_rng().normal(predictions, standard_deviations)

            df = self._get_sorted_df(smiles=original_smiles, predictions=gaussian_predictions)
            acquired_smiles = list(df["smiles"][:self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size:])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.UCB:
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            predictions = self.surrogate.predict(fps)
            standard_deviations = self.surrogate.get_std(fps)
            ucb = predictions + beta * standard_deviations

            df = self._get_sorted_df(smiles=original_smiles, predictions=ucb)
            acquired_smiles = list(df["smiles"][:self.acquisition_batch_size])
            non_acquired_smiles = list(df["smiles"][self.acquisition_batch_size:])

            return acquired_smiles, non_acquired_smiles

        if self.function == _SALE.TANIMOTO:
            fps = self._get_morgan_fingerprints(smiles=original_smiles)
            similarities = np.zeros(len(fps))

            for idx, fp in enumerate(fps):
                similarities[idx] = np.mean(DataStructs.BulkTanimotoSimilarity(fp, fingerprints_pool))

            # sort by ascending Tanimoto similarities
            sorted_idx = np.argsort(similarities)
            acquired_idx = sorted_idx[:self.acquisition_batch_size]
            non_acquired_idx = sorted_idx[self.acquisition_batch_size:]

            acquired_smiles = [original_smiles[idx] for idx in acquired_idx]
            non_acquired_smiles = [original_smiles[idx] for idx in non_acquired_idx]

            return list(acquired_smiles), list(non_acquired_smiles)

        else:
            raise ValueError("Error: Invalid/Unsupported Acquisition Function.")

    def _get_morgan_fingerprints(self, smiles) -> list:
        """returns morgan fingerprints for a list of smiles"""
        molecules = [Chem.MolFromSmiles(s) for s in smiles]
        return [
            Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
            for mol in molecules
        ]

    def _get_sorted_df(self, smiles, predictions, highest_is_best=False):
        # TODO: need "highest_is_best" to be set in class instantiation
        df = pd.DataFrame({"smiles": smiles, "predictions": predictions})
        if highest_is_best:
            return df.sort_values(by="predictions", ascending=False)
        else:
            return df.sort_values(by="predictions", ascending=True)

"""
        if self.function == _SALE.PI:
            # probability of improvement
            if model == None or current_max == None:
                raise ValueError(f"Incorrect or missed parameters in acquire for AF='{AF}'")

            sd = model.get_std(X)
            # calculate probabilities
            means = y
            probs = norm.cdf((means - current_max + xi) / sd)
        
        if self.function == _SALE.EI:
            # expected improvement

            if model == None or current_max == None:
                raise ValueError(f"Incorrect or missed parameters in acquire for AF='{AF}'")

            means = y
            Y = (means - current_max + xi)
            sd = model.get_std(X)
            Z = Y / sd
            expected = Y * norm.cdf(Z) + sd * norm.pdf(Z)
"""
