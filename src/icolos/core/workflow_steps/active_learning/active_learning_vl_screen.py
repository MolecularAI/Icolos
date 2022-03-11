from typing import List, Tuple
import os
import pickle

from modAL.models.learners import ActiveLearner
from pydantic.main import BaseModel


from icolos.core.workflow_steps.active_learning.base import ActiveLearningBase


from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import (
    StepActiveLearningEnum,
)
from rdkit.Chem import PandasTools
from rdkit import Chem
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np
import matplotlib.pyplot as plt


_SALE = StepActiveLearningEnum()


class StepActiveLearning(ActiveLearningBase, BaseModel):
    """
    Class to run an active learning framework
    Primarily designed for building QSAR models using a physics based method (embedding + docking) as an oracle

    Takes the step conf for the oracle as an additional argument.  The step with these settings is run with the queried compounds at each stage of the active learning loop
    """

    def __init__(self, **data):
        super().__init__(**data)

    def _parse_library(self, criteria: str = None) -> Tuple[DataFrame, np.ndarray]:
        """
        Loads the library file from disk
        This should be a .sdf or smi file containing the compounds to be screened
        """
        lib_path = self.settings.additional[_SALE.VIRTUAL_LIB]
        if lib_path.endswith("sdf"):
            # hold the lib in a pandas df
            library = PandasTools.LoadSDF(
                lib_path,
                smilesName=_SALE.SMILES,
                molColName=_SALE.MOLECULE,
                includeFingerprints=True,
                removeHs=False,
                embedProps=True,
            )
        elif lib_path.endswith("smi"):
            mols, smiles = [], []
            suppl = Chem.rdmolfiles.SmilesMolSupplier(lib_path)
            for mol in suppl:
                mols.append(mol)
                smiles.append(Chem.MolToSmiles(mol))
                # pd.concat([library, entry])
            library = pd.DataFrame({_SALE.SMILES: smiles, _SALE.MOLECULE: mols})
            print(library.head())
            library = self.construct_fingerprints(library)
        scores = (
            np.absolute(pd.to_numeric(library[criteria].fillna(0)))
            if criteria is not None
            else []
        )

        return library, scores

    def _run_learning_loop_virtual_lib(
        self, learner: ActiveLearner, lib: pd.DataFrame, top_1_idx: list = None
    ) -> tuple[pd.DataFrame, List]:
        rounds = int(self.settings.additional[_SALE.N_ROUNDS])
        n_instances = int(self.settings.additional[_SALE.BATCH_SIZE])
        queried_compound_idx = []
        queried_compounds_by_batch = []
        fraction_top1_hits = []
        key = _SALE.MORGAN_FP
        X = np.array(list(lib[key]))
        for rnd in range(rounds):
            query_idx, _ = learner.query(
                X,
                n_instances=n_instances,
                previous_idx=queried_compound_idx,
            )
            queried_compound_idx += list(query_idx)
            if (rnd + 1) % 5 == 0:
                queried_compounds_by_batch.append([int(i) for i in query_idx])
            query_compounds = [lib.iloc[int(idx)] for idx in query_idx]

            if self.get_additional_setting(_SALE.EVALUATE, default=False):
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
            else:
                # get scores from oracle
                self._logger.log(
                    f"Querying oracle with {len(query_compounds)} compounds", _LE.INFO
                )

                compounds = self.query_oracle(query_compounds)
                scores = self._extract_final_scores(
                    compounds, self.settings.additional[_SALE.CRITERIA]
                )

            if self.settings.additional["model"] == "ffnn":
                scores = scores.reshape(-1, 1)

            self._logger.log("Fitting with new data...", _LE.INFO)
            new_data = np.array([compound[key] for compound in query_compounds])
            learner.teach(new_data, scores, only_new=False)
            # calculate percentage of top-1% compounds queried
            if top_1_idx is not None:
                # what fraction of top1% compoudnds has the model requested to evaluate
                hits_queried = (
                    np.in1d(np.unique(np.array(queried_compound_idx)), top_1_idx).sum()
                    / len(top_1_idx)
                ) * 100
                self._logger.log(
                    f"Fraction of top 1% hits queries by round {rnd}: {hits_queried}%",
                    _LE.INFO,
                )
                fraction_top1_hits.append(hits_queried)
        self._logger.log(f"Hits queried at each epoch:\n{fraction_top1_hits}", _LE.INFO)

        # return the enriched df
        return lib.iloc[queried_compound_idx], queried_compounds_by_batch

    def execute(self):
        tmp_dir = self._make_tmpdir()
        print(tmp_dir)
        criteria = (
            self.get_additional_setting(_SALE.CRITERIA)
            if self.get_additional_setting(_SALE.EVALUATE, default=False)
            else None
        )
        lib, scores = self._parse_library(criteria=criteria)
        if scores is not None:
            top_1_percent = int(0.01 * len(scores))
            top_1_idx = np.argpartition(scores, -top_1_percent)[-top_1_percent:]

        learner = self._initialize_learner()

        self._logger.log("Evaluating virtual library...", _LE.DEBUG)
        (
            enriched_lib,
            queried_compound_idx_by_batch,
        ) = self._run_learning_loop_virtual_lib(
            learner=learner, lib=lib, top_1_idx=list(top_1_idx)
        )

        # compare distributions
        print(lib[criteria].head())
        print(lib[criteria].astype(float).describe())
        print(enriched_lib[criteria].astype(float).describe())
        enriched_lib.to_csv(os.path.join(tmp_dir, "enriched_lib.csv"))
        bins = np.linspace(-12, -9, 30)
        fig, axs = plt.subplots(
            len(queried_compound_idx_by_batch) + 1,
            1,
            figsize=(10, len(queried_compound_idx_by_batch) * 4),
        )
        axs[0].set_title("Full library distribution + final enriched set")
        axs[0].hist(
            lib[criteria].astype(np.float32),
            label="full lib",
            bins=bins,
            density=True,
            alpha=0.5,
        )
        axs[0].hist(
            enriched_lib[criteria].astype(np.float32),
            label=f"enriched lib",
            bins=bins,
            density=True,
            alpha=0.5,
        )
        axs[0].axvline(
            lib[criteria].astype(np.float32).mean(),
            color="red",
            linestyle="dashed",
            linewidth=2,
            label="full lib mean",
        )
        axs[0].axvline(
            enriched_lib[criteria].astype(np.float32).mean(),
            color="green",
            linestyle="dashed",
            linewidth=2,
            label="enriched lib mean",
        )
        axs[0].legend()
        for axis in ["top", "bottom", "left", "right"]:
            axs[0].spines[axis].set_linewidth(4)
        for idx, batch_idx in enumerate(queried_compound_idx_by_batch):
            axs[idx + 1].hist(
                lib.iloc[batch_idx][criteria].astype(np.float32),
                label=f"batch {(idx+1) * 5}",
                bins=bins,
                density=True,
                alpha=0.5,
            )
            axs[idx + 1].axvline(
                lib.iloc[batch_idx][criteria].astype(np.float32).mean(),
                color="green",
                linestyle="dashed",
                linewidth=1,
                label="mean",
            )

            axs[idx + 1].legend()
            # axs[idx + 1].set_title(f"enriched lib batch {(idx + 1) * 5}")
        fig.savefig(os.path.join(tmp_dir, "dist.png"), dpi=300)

        # pickle the final model
        # with open(os.path.join(tmp_dir, "model.pkl"), "wb") as f:
        #     pickle.dump(learner, f)

        self._parse_output(tmp_dir)
