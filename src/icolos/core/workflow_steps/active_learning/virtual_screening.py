from typing import List, Tuple
import os
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
        else:
            raise ValueError(f"File {lib_path} must of of type smi or sdf")
        library = self.construct_fingerprints(library)
        scores = (
            np.absolute(pd.to_numeric(library[criteria].fillna(0)))
            if criteria is not None
            else []
        )

        return library, scores

    def _run_learning_loop_virtual_lib(
        self,
        learner: ActiveLearner,
        lib: pd.DataFrame,
        tmp_dir: str,
        top_1_idx: list = None,
        replica=0,
        fragment_lib: pd.DataFrame = None,
    ) -> None:

        n_instances = int(self.settings.additional[_SALE.FRACTION_PER_EPOCH] * len(lib))
        rounds = int(
            (
                self.get_additional_setting(_SALE.MAX_SAMPLED_FRACTION, default=0.25)
                * len(lib)
            )
            / n_instances
        )
        self._logger.log(
            f"Running {rounds} rounds of {n_instances} componds, replica {replica}",
            _LE.DEBUG,
        )
        queried_compound_idx = []
        queried_compounds_per_epoch = []
        fraction_top1_hits_per_epoch = []
        key = _SALE.MORGAN_FP
        X = np.array(list(lib[key]))
        oracle_type = self.get_additional_setting(_SALE.ORACLE_TYPE, default="docking")

        for rnd in range(rounds):
            warmup = rnd < self.get_additional_setting(_SALE.WARMUP, default=1)
            init_batch_size = int(
                n_instances
                * self.get_additional_setting(_SALE.INIT_SAMPLE_FACTOR, default=1)
            )
            self._logger.log(f"Generating {init_batch_size} samples", _LE.DEBUG)
            query_idx, _ = learner.query(
                X,
                n_instances=init_batch_size if rnd == 0 else n_instances,
                previous_idx=queried_compound_idx,
                warmup=warmup,
                epsilon=self.get_additional_setting(_SALE.EPSILON, default=0.0),
            )
            queried_compound_idx += list(query_idx)
            queried_compounds_per_epoch.append([int(i) for i in query_idx])
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

                compounds = self.query_oracle(
                    query_compounds, fragment_lib=fragment_lib, oracle_type=oracle_type
                )
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
                    np.isin(np.unique(queried_compound_idx), top_1_idx).sum()
                    / len(top_1_idx)
                ) * 100
                self._logger.log(
                    f"Fraction of top 1% hits queries by round {rnd}: {hits_queried}%",
                    _LE.INFO,
                )
                fraction_top1_hits_per_epoch.append(hits_queried)
                if hits_queried > 99:
                    break
            df = lib.iloc[queried_compound_idx]
            df.to_pickle(
                os.path.join(tmp_dir, f"enriched_lib_rep_{replica}round_{rnd+1}.pkl")
            )

            if (
                self.get_additional_setting(_SALE.DYNAMIC_STOP, default=False) is True
                and rnd > 4
            ):
                # average of last three top-1 fractions
                rolling_avg = np.mean(fraction_top1_hits_per_epoch[-5:-1])
                if fraction_top1_hits_per_epoch[-1] - rolling_avg < 0.01:
                    self._logger.log(
                        f"Fraction top1 hits converged to 0.01, stopping...", _LE.INFO
                    )
                    break

        # create a results dataframe to store per-epoch properties
        self._logger.log("Generating results dataframe...", _LE.DEBUG)

        # we have a list of compound IDs per epoch, we want the transpose
        queried_compounds_per_position = np.array(queried_compounds_per_epoch).T
        col_dict = {
            "epoch": [i for i in range(len(fraction_top1_hits_per_epoch))],
            "top_1%_per_epoch": fraction_top1_hits_per_epoch,
        }
        # for idx, row in enumerate(queried_compounds_per_position):
        #     col_dict[f"seq_idx_pos_{idx}"] = list(row)
        resuts_df = pd.DataFrame(col_dict)
        resuts_df.to_csv(os.path.join(tmp_dir, f"results_rep_{replica}.csv"))

        # return the enriched df
        enriched_lib = lib.iloc[queried_compound_idx]
        enriched_lib.to_pickle(os.path.join(tmp_dir, f"enriched_lib_rep_{replica}.pkl"))

    def execute(self):
        tmp_dir = self._make_tmpdir()
        criteria = (
            self.get_additional_setting(_SALE.CRITERIA)
            if self.get_additional_setting(_SALE.EVALUATE, default=False)
            else None
        )
        lib, scores = self._parse_library(criteria=criteria)
        lib.to_pickle(os.path.join(tmp_dir, "starting_lib.pkl"))
        if scores is not None:
            top_1_percent = int(0.01 * len(scores))
            top_1_idx = np.argpartition(scores, -top_1_percent)[-top_1_percent:]

        # load fragment lib if provided
        print("loading fragment library")
        fragments_libary = (
            PandasTools.LoadSDF(
                self.get_additional_setting(_SALE.FRAGMENTS),
                smilesName=_SALE.SMILES,
                molColName=_SALE.MOLECULE,
                includeFingerprints=True,
                removeHs=False,
                embedProps=True,
            )
            if self.get_additional_setting(_SALE.FRAGMENTS) is not None
            else None
        )

        print(fragments_libary.head())

        replicas = self.get_additional_setting(_SALE.REPLICAS, default=1)
        for replica in range(replicas):
            learner = self._initialize_learner()

            self._logger.log("Evaluating virtual library...", _LE.DEBUG)
            self._run_learning_loop_virtual_lib(
                learner=learner,
                lib=lib,
                fragment_lib=fragments_libary,
                tmp_dir=tmp_dir,
                top_1_idx=list(top_1_idx),
                replica=replica,
            )
