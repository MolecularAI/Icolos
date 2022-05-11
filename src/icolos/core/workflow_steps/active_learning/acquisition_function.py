import pandas as pd
import numpy as np
from icolos.core.workflow_steps.step import _LE

from icolos.core.workflow_steps.active_learning.surrogate_model import SurrogateModel

from icolos.utils.enums.step_enums import StepActiveLearningEnum

_SALE = StepActiveLearningEnum()


class AcquisitionFunction:
    """class to execute all acquisition function heuristics for prospective REINVENT"""

    def __init__(
        self, surrogate: SurrogateModel, function: str, acquisition_batch_size: int
    ):
        self.surrogate = surrogate
        self.function = function.lower()
        self.acquisition_batch_size = acquisition_batch_size

    def partition_compounds(self, original_smiles):
        # TODO: only the random acquisition function is tested. The rest of the acquisition functions' code
        #  is from Christian and Clara
        #  Will integrate/modify that (commented-out) code soon
        if self.function == _SALE.RANDOM:
            df = pd.DataFrame({"smiles": original_smiles})
            shuffled = df.sample(frac=1, axis=0)
            acquired_smiles = list(shuffled["smiles"][: self.acquisition_batch_size])
            non_acquired_smiles = list(
                shuffled["smiles"][self.acquisition_batch_size :]
            )

            return acquired_smiles, non_acquired_smiles


# TODO: need to somehow keep track of "current best" as we are not storing it in memory, at least for now
"""
        if self.function == _SALE.GREEDY:
            data = pd.DataFrame({"Smiles": smiles, "Fingerprint": X, "Score": y, "GT": ground_truth})
            out = data.sort_values(by="Score", ascending=False)
            out_tuple = (out["Fingerprint"][:n].to_numpy(),
                         out["GT"][:n].to_numpy(),
                         out["Smiles"][:n].to_numpy(),
                         out["Fingerprint"][n:].to_numpy(),
                         out["Score"][n:].to_numpy(),
                         out["GT"][n:].to_numpy(),
                         out["Smiles"][n:].to_numpy(),
                         )

            return out_tuple

        if self.function == _SALE.GREEDY_EPSILON:
            # TODO: get "special parameters"

            if np.random.default_rng().random() < epsilon:
                return acquire(smiles, X, y, ground_truth, n=n, AF='random')
            else:
                data = pd.DataFrame({"Smiles": smiles, "Fingerprint": X, "Score": y, "GT": ground_truth})
                out = data.sort_values(by="Score", ascending=False)

                out_tuple = (out["Fingerprint"][:n].to_numpy(),
                             out["GT"][:n].to_numpy(),
                             out["Smiles"][:n].to_numpy(),
                             out["Fingerprint"][n:].to_numpy(),
                             out["Score"][n:].to_numpy(),
                             out["GT"][n:].to_numpy(),
                             out["Smiles"][n:].to_numpy(),
                             )

                return out_tuple

        if self.function == _SALE.PI:
            # probability of improvement
            if model == None or current_max == None:
                raise ValueError(f"Incorrect or missed parameters in acquire for AF='{AF}'")

            sd = model.get_std(X)
            # calculate probabilities
            means = y
            probs = norm.cdf((means - current_max + xi) / sd)

            out = pd.DataFrame({"Smiles": smiles, "Fingerprint": X, "Probability": probs, "GT": ground_truth,
                                "Score": means}).sort_values(by="Probability", ascending=False)
            out_tuple = (out["Fingerprint"][:n].to_numpy(),
                         out["GT"][:n].to_numpy(),
                         out["Smiles"][:n].to_numpy(),
                         out["Fingerprint"][n:].to_numpy(),
                         out["Score"][n:].to_numpy(),
                         out["GT"][n:].to_numpy(),
                         out["Smiles"][n:].to_numpy(),
                         )

            return out_tuple
        
        if self.function == _SALE.EI:
            # expected improvement

            if model == None or current_max == None:
                raise ValueError(f"Incorrect or missed parameters in acquire for AF='{AF}'")

            means = y
            Y = (means - current_max + xi)
            sd = model.get_std(X)
            Z = Y / sd
            expected = Y * norm.cdf(Z) + sd * norm.pdf(Z)
            out = pd.DataFrame(
                {"Smiles": smiles, "Fingerprint": X, "Expected improvement": expected, "GT": ground_truth,
                 "Score": means}).sort_values(by="Expected improvement", ascending=False)
            out_tuple = (out["Fingerprint"][:n].to_numpy(),
                         out["GT"][:n].to_numpy(),
                         out["Smiles"][:n].to_numpy(),
                         out["Fingerprint"][n:].to_numpy(),
                         out["Score"][n:].to_numpy(),
                         out["GT"][n:].to_numpy(),
                         out["Smiles"][n:].to_numpy(),
                         )

            return out_tuple
        
        if self.function == _SALE.TS:
            # thompson sampling
            means = y
            sd = model.get_std(X)

            scores = np.random.default_rng().normal(means, sd)
            out = pd.DataFrame({"Smiles": smiles, "Fingerprint": X, "Score": scores, "GT": ground_truth}).sort_values(
                by="Score", ascending=False)
            out_tuple = (out["Fingerprint"][:n].to_numpy(),
                         out["GT"][:n].to_numpy(),
                         out["Smiles"][:n].to_numpy(),
                         out["Fingerprint"][n:].to_numpy(),
                         out["Score"][n:].to_numpy(),
                         out["GT"][n:].to_numpy(),
                         out["Smiles"][n:].to_numpy(),
                         )

            return out_tuple
        
        if self.function == _SALE.UCB:
            if beta == None:
                raise ValueError(f"Incorrect or missed parameters in acquire for AF='{AF}'")

            # Calculate the mean and standard deviation
            means = y
            std = model.get_std(X)

            UCB = means + beta * std

            out = pd.DataFrame(
                {"Smiles": smiles, "Fingerprint": X, "UCB": UCB, "GT": ground_truth, "Score": y}).sort_values(by="UCB",
                                                                                                              ascending=False)
            out_tuple = (out["Fingerprint"][:n].to_numpy(),
                         out["GT"][:n].to_numpy(),
                         out["Smiles"][:n].to_numpy(),
                         out["Fingerprint"][n:].to_numpy(),
                         out["Score"][n:].to_numpy(),
                         out["GT"][n:].to_numpy(),
                         out["Smiles"][n:].to_numpy(),
                         )

            return out_tuple

        if self.function == _SALE.TANIMOTO:
            # tanimoto similarity
            if tanimoto_comparison == None:
                raise ValueError(f"Incorrect or missed parameters in acquire for AF='{AF}'")

            scores = np.zeros((len(X)))

            for i, target in enumerate(X):
                scores[i] = np.mean(DataStructs.BulkTanimotoSimilarity(target, tanimoto_comparison))

            # sort list ascending
            sorted_idx = np.argsort(scores)
            acq_idx = sorted_idx[:n]
            nacq_idx = sorted_idx[n:]

            acq_fps = [X[i] for i in acq_idx]
            nacq_fps = [X[i] for i in nacq_idx]

            acq_smiles = [smiles[i] for i in acq_idx]
            nacq_smiles = [smiles[i] for i in nacq_idx]

            gt = ground_truth.to_numpy()
            acq_scores = gt[acq_idx]
            nacq_scores = gt[nacq_idx]
            nacq_preds = y[nacq_idx]

            return acq_fps, acq_scores, acq_smiles, nacq_fps, nacq_preds, nacq_scores, nacq_smiles
"""
