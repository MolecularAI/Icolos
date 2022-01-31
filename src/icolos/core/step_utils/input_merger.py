from copy import deepcopy
from typing import List, Dict
from pydantic import BaseModel

from icolos.core.containers.compound import Enumeration, Compound
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class StepMerge(BaseModel):
    compounds: bool = True
    enumerations: bool = False
    merge_compounds_by: str = _SBE.INPUT_MERGE_BY_NAME
    merge_enumerations_by: str = _SBE.INPUT_MERGE_BY_ID


class InputMerger:
    def __init__(self, config: StepMerge):
        self.config = config

    def _group_enumerations(
        self, enumerations: List[Enumeration], by
    ) -> Dict[str, List[Enumeration]]:
        if by == _SBE.INPUT_MERGE_BY_SMILE:
            grouped = {enumeration.get_smile(): [] for enumeration in enumerations}
            for enum in enumerations:
                grouped[enum.get_smile()].append(enum)
        elif by == _SBE.INPUT_MERGE_BY_ID:
            grouped = {
                str(enumeration.get_enumeration_id()): []
                for enumeration in enumerations
            }
            for enum in enumerations:
                grouped[str(enum.get_enumeration_id())].append(enum)
        else:
            raise NotImplementedError
        return grouped

    def _group_compounds(
        self, compounds: List[Compound], by
    ) -> Dict[str, List[Compound]]:
        if by == _SBE.INPUT_MERGE_BY_NAME:
            names = {compound.get_name(): [] for compound in compounds}
            for compound in compounds:
                names[compound.get_name()].append(compound)
        elif by == _SBE.INPUT_MERGE_BY_ID:
            names = {str(compound.get_compound_number()): [] for compound in compounds}
            for compound in compounds:
                names[str(compound.get_compound_number())].append(compound)
        else:
            raise NotImplementedError
        return names

    def _merge_enumerations(
        self, enumerations: List[Enumeration], by
    ) -> List[Enumeration]:
        list_result = []

        # note that if it has been grouped by ID, the first (arbitrary) smile is used
        for _, enum_list in self._group_enumerations(enumerations, by).items():
            enum_combined = deepcopy(enum_list[0])
            enum_combined.clear_conformers()
            for enum in enum_list:
                enum_combined.add_conformers(
                    deepcopy(enum.get_conformers()), auto_update=False
                )
            list_result.append(enum_combined)
        return list_result

    def unroll_compounds(self, compounds: list) -> List[Compound]:
        list_buffer = []
        for ele in compounds:
            if isinstance(ele, list):
                list_buffer = list_buffer + self.unroll_compounds(ele)
            elif isinstance(ele, Compound):
                list_buffer.append(ele)
        return list_buffer

    def merge(self, compounds: List[Compound]) -> List[Compound]:
        list_result = []

        # if selected, combined compounds into one depending on the strategy
        if self.config.compounds:
            dict_grouped = self._group_compounds(
                compounds, self.config.merge_compounds_by
            )
            number = 0
            for name, compound_list in dict_grouped.items():
                # add the enumerations of all compounds together but do NOT auto-update yet (because enumerations might
                # also be merged later on)
                comp_combined = Compound(name=name, compound_number=number)
                for comp in compound_list:
                    comp_combined.add_enumerations(
                        deepcopy(comp.get_enumerations()), auto_update=False
                    )

                # as merging of enumerations only makes sense when there was a compound merge, keep
                # it on that indentation level
                if self.config.enumerations:
                    enumerations = self._merge_enumerations(
                        deepcopy(comp_combined.get_enumerations()),
                        self.config.merge_enumerations_by,
                    )
                    comp_combined.clear_enumerations()
                    comp_combined.add_enumerations(enumerations, auto_update=False)

                # now, rename the enumerations and conformers
                comp_combined.reset_all_ids()
                comp_combined.update_all_relations()

                list_result.append(comp_combined)
                number += 1
        return list_result
