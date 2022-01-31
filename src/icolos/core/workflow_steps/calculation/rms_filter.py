import pandas as pd
from typing import List
from pydantic import BaseModel

from icolos.core.containers.compound import Conformer

from icolos.utils.enums.step_enums import StepRMSFilterEnum
from icolos.core.workflow_steps.step import _LE
from icolos.core.workflow_steps.calculation.base import StepCalculationBase

_SRF = StepRMSFilterEnum()


class StepRMSFilter(StepCalculationBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # extend parameters
        if _SRF.THRESHOLD not in self.settings.additional.keys():
            self.settings.additional[_SRF.THRESHOLD] = 1
        if _SRF.METHOD not in self.settings.additional.keys():
            self.settings.additional[_SRF.METHOD] = _SRF.METHOD_ALIGNMOL
        if _SRF.ORDER_BY not in self.settings.additional.keys():
            self.settings.additional[_SRF.ORDER_BY] = None
        else:
            if _SRF.ORDER_ASCENDING not in self.settings.additional.keys():
                self._logger.log(
                    'Setting order ascending not specified, setting to "True" (default).',
                    _LE.WARNING,
                )
                self.settings.additional[_SRF.ORDER_ASCENDING] = False

    def _get_representative_indices(
        self, df_rms: pd.DataFrame, prop_values: List[float]
    ) -> List[int]:
        keep_indices = []
        prop_idx = list(zip(prop_values, list(range(len(prop_values)))))
        threshold = self.settings.additional[_SRF.THRESHOLD]
        while len(prop_idx) > 0:
            # get the best (according to the property) element's index, add it to the list and remove it from
            # the remaining ones
            if self.settings.additional[_SRF.ORDER_BY] is not None:
                prop_idx = [
                    (prop, idx)
                    for prop, idx in sorted(
                        prop_idx, reverse=self.settings.additional[_SRF.ORDER_ASCENDING]
                    )
                ]
            cur_best_idx = prop_idx[0][1]
            keep_indices.append(cur_best_idx)
            del prop_idx[0]

            # remove all, that are fulfilling the RMS threshold
            for i in reversed(range(len(prop_idx))):
                comp_idx = prop_idx[i][1]
                cur_rms = df_rms.iloc[cur_best_idx, comp_idx]
                if cur_rms <= threshold:
                    del prop_idx[i]
        return keep_indices

    def _filter_conformers(self, conformers: List[Conformer]) -> List[Conformer]:
        # to select the "best" conformers, here the property to use for ordering / ranking is specified
        order_by = self.settings.additional[_SRF.ORDER_BY]
        if order_by is not None:
            prop_values = self._get_property_values(conformers, order_by)
        else:
            prop_values = [None for _ in range(len(conformers))]

        # generate RMS matrix (NxN, where N is the number of conformers)
        df_rms = self._calculate_rms_matrix(conformers, self._get_rms_method())

        keep_indices = self._get_representative_indices(df_rms, prop_values)

        return [conformers[i] for i in range(len(conformers)) if i in keep_indices]

    def execute(self):
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if self._input_object_empty(enumeration):
                    continue

                number_conformers_before = len(enumeration)

                # filter conformers on the enumeration level
                filtered_conformers = self._filter_conformers(
                    conformers=enumeration.get_conformers()
                )

                # add filtered conformers to enumeration
                enumeration.clear_conformers()
                for conf in filtered_conformers:
                    enumeration.add_conformer(conformer=conf, auto_update=True)
                number_conformers_after = len(enumeration)
                self._logger.log(
                    f"Filtered {number_conformers_before} conformers down to {number_conformers_after} for enumeration {enumeration.get_index_string()}.",
                    _LE.INFO,
                )
