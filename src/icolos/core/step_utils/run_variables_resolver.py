from copy import deepcopy

from icolos.core.containers.compound import Enumeration, Compound, Conformer
from icolos.utils.enums.write_out_enums import RunVariablesEnum

_RVE = RunVariablesEnum()


class RunVariablesResolver:
    def __init__(self):
        pass

    def _replace(self, input_str: str, pattern: str, replacement) -> str:
        if replacement is not None:
            pattern = _RVE.PREFIX + pattern + _RVE.POSTFIX
            input_str = input_str.replace(pattern, str(replacement))
        return input_str

    def resolve_compound_level(self, input_str: str, comp: Compound) -> str:
        comp = deepcopy(comp)
        resolved_str = self._replace(
            input_str, _RVE.COMPOUND_ID, comp.get_compound_number()
        )
        resolved_str = self._replace(resolved_str, _RVE.COMPOUND_NAME, comp.get_name())
        return resolved_str

    def resolve_enumeration_level(self, input_str: str, enum: Enumeration) -> str:
        enum = deepcopy(enum)
        resolved_str = self._replace(
            input_str, _RVE.ENUMERATION_ID, enum.get_enumeration_id()
        )
        resolved_str = self._replace(
            resolved_str, _RVE.ENUMERATION_STRING, enum.get_index_string()
        )
        return resolved_str

    def resolve_conformer_level(self, input_str: str, conf: Conformer) -> str:
        conf = deepcopy(conf)
        resolved_str = self._replace(
            input_str, _RVE.CONFORMER_ID, conf.get_conformer_id()
        )
        resolved_str = self._replace(
            resolved_str, _RVE.CONFORMER_STRING, conf.get_index_string()
        )
        return resolved_str

    def resolve(self, input_str: str, input_object) -> str:
        if not isinstance(input_str, str):
            return input_str

        if isinstance(input_object, Conformer):
            input_str = self.resolve_compound_level(
                input_str, input_object.get_enumeration_object().get_compound_object()
            )
            input_str = self.resolve_enumeration_level(
                input_str, input_object.get_enumeration_object()
            )
            return self.resolve_conformer_level(input_str, input_object)
        elif isinstance(input_object, Enumeration):
            input_str = self.resolve_compound_level(
                input_str, input_object.get_compound_object()
            )
            return self.resolve_enumeration_level(input_str, input_object)
        elif isinstance(input_object, Compound):
            return self.resolve_compound_level(input_str, input_object)
        else:
            raise ValueError(f'Object of type "{type(input_object)}" is not supported.')
