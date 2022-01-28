import pandas as pd
from typing import List, Tuple

from pydantic import BaseModel

from icolos.core.containers.compound import Conformer

from icolos.utils.enums.step_enums import StepClusteringEnum
from icolos.core.workflow_steps.step import _LE
from icolos.core.workflow_steps.calculation.base import StepCalculationBase

from sklearn.cluster import KMeans

_SC = StepClusteringEnum()


class StepClustering(StepCalculationBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # extend parameters
        if _SC.N_CLUSTERS not in self.settings.arguments.parameters.keys():
            self.settings.arguments.parameters[_SC.N_CLUSTERS] = 3
        if _SC.MAX_ITER not in self.settings.arguments.parameters.keys():
            self.settings.arguments.parameters[_SC.MAX_ITER] = 300
        if _SC.TOP_N_PER_SOLVENT not in self.settings.additional.keys():
            self.settings.additional[_SC.TOP_N_PER_SOLVENT] = 3

    def _get_nclusters_and_top_n(self, len_conformers: int) -> Tuple[int, int]:
        n_clusters = self.settings.arguments.parameters[_SC.N_CLUSTERS]
        if n_clusters > len_conformers:
            n_clusters = len_conformers
            self._logger.log(
                f"Set number of clusters to {n_clusters} because not enough observations were provided.",
                _LE.DEBUG,
            )
        top_n_per_solvent = self.settings.additional[_SC.TOP_N_PER_SOLVENT]
        if top_n_per_solvent > len_conformers:
            top_n_per_solvent = len_conformers
            self._logger.log(
                f'Set number of "top_N_per_solvent" to {top_n_per_solvent} because not enough observations were provided.',
                _LE.DEBUG,
            )
        return n_clusters, top_n_per_solvent

    def _generate_feature_dataframe(self, conformers: List[Conformer]) -> pd.DataFrame:
        features = self.settings.additional[_SC.FEATURES]
        df_features = pd.DataFrame(columns=features)
        for conf in conformers:
            new_row = {}
            for feature in features:
                new_row[feature] = float(conf.get_molecule().GetProp(feature))
            df_features = df_features.append(new_row, ignore_index=True)
        return df_features

    def _get_representative_conformers(
        self, cluster_set: List[Tuple[int, Conformer]]
    ) -> List[int]:
        # for each selection (e.g. solvent), obtain the N top conformers (note, that the input is already clustered)
        # also get rid of duplicates in the indices
        rep_indices = []
        for solvent_key in self.settings.additional[_SC.FREE_ENERGY_SOLVENT_TAGS]:
            conf_indices = [tuple_conf[0] for tuple_conf in cluster_set]
            solvent_dGs = [
                float(tuple_conf[1].get_molecule().GetProp(solvent_key))
                for tuple_conf in cluster_set
            ]

            # sort list of global indices for this cluster according to their free energy for this solvent
            # note: from lowest (most negative) -> highest
            conf_indices_sorted = [
                idx for _, idx in sorted(zip(solvent_dGs, conf_indices))
            ]
            rep_indices = (
                rep_indices
                + conf_indices_sorted[
                    0 : min(
                        len(conf_indices),
                        self.settings.additional[_SC.TOP_N_PER_SOLVENT],
                    )
                ]
            )
        return list(set(rep_indices))

    def _cluster_conformers(self, conformers: List[Conformer]) -> List[Conformer]:
        # make sure the number of clusters specified and "N top per solvent" are not higher than the compound number
        n_clusters, top_n_per_solvent = self._get_nclusters_and_top_n(
            len_conformers=len(conformers)
        )

        # initialize K-means instance
        kmeans = KMeans(
            n_clusters=n_clusters,
            max_iter=self.settings.arguments.parameters[_SC.MAX_ITER],
            init="k-means++",
            n_init=10,
            tol=1e-04,
            random_state=0,
        )

        # generate dataframe with selected properties
        df_features = self._generate_feature_dataframe(conformers=conformers)

        # predict cluster and assign to conformer
        cluster_labels = kmeans.fit_predict(df_features)
        keep_indices = []
        for cluster_label in range(n_clusters):
            # keep the "global" index to select the appropriate conformers later
            cluster_set = [
                (i, conformers[i])
                for i in range(len(conformers))
                if cluster_labels[i] == cluster_label
            ]
            keep_indices = keep_indices + self._get_representative_conformers(
                cluster_set=cluster_set
            )
        return [conformers[i] for i in range(len(conformers)) if i in keep_indices]

    def execute(self):
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if len(enumeration.get_conformers()) == 0:
                    continue

                number_conformers_before = len(enumeration)

                # cluster conformers on the enumeration level
                clustered_conformers = self._cluster_conformers(
                    conformers=enumeration.get_conformers()
                )

                # add clustered conformers to enumeration
                enumeration.clear_conformers()
                for conf in clustered_conformers:
                    enumeration.add_conformer(conformer=conf, auto_update=True)
                number_conformers_after = len(enumeration)
                self._logger.log(
                    f"Clustered {number_conformers_before} into {number_conformers_after} conformers for enumeration {enumeration.get_index_string()}.",
                    _LE.INFO,
                )
