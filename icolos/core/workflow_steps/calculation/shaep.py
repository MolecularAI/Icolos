from icolos.utils.execute_external.execute import Executor
from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import StepShaepEnum
from icolos.utils.enums.program_parameters import PantherEnum, ShaepEnum
from icolos.core.containers.compound import Conformer
import tempfile
from pydantic import BaseModel
import os

_SSE = StepShaepEnum()
_SE = ShaepEnum()
_PE = PantherEnum()


class StepShaep(StepCalculationBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=Executor)

    def _prepare_tmp_input_dir(self):
        tmp_dir = tempfile.mkdtemp()
        return tmp_dir

    def _execute_backend(self, conf_path: str, tmp_dir: str, ni_path: str):
        arguments = [
            os.path.join(tmp_dir, ni_path),
            conf_path,
            os.path.join(tmp_dir, _SE.OUTPUT_SIMILARITY),
        ]
        self._backend_executor.execute(
            command=_SE.SHAEP_EXECUTABLE, arguments=arguments, check=True
        )

    def _parse_output(self, tmp_dir: str, conformer: Conformer):
        with open(os.path.join(tmp_dir, _SE.OUTPUT_SIMILARITY), "r") as f:
            # TODO: add support for multiple input structures; ignore the names (all will be in one line), but from
            #       position 8 (index 7 in python) onwards, the shape and esp similarities are reported in the same
            #       order as the input, i.e. <7 other values> mol1_shape mol1_esp mol2_shape ...
            parts = f.readlines()[1].split("\t")
            conformer.get_molecule().SetProp(_SE.TAG_SHAPE_SIMILARITY, str(parts[7]))
            conformer.get_molecule().SetProp(_SE.TAG_ESP_SIMILARITY, str(parts[8]))

    def execute(self):
        number_rescored = 0
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if len(enumeration.get_conformers()) == 0:
                    self._logger.log(
                        f"Found no conformers for enumeration {enumeration} for compound {compound}.",
                        _LE.WARNING,
                    )
                    # we can still execute shaep at the enumeration level, if the compounds are correcty annotated they should be written out ok.  Will be slower though
                    # easiest for now is to add the enumeration mol object as a single conformer and run that through shaep
                    mol = enumeration.get_molecule()
                    conf = Conformer(conformer=mol)
                    enumeration.add_conformer(conf)

                # TODO: ShaEP allow batch execution for any number of compounds (parsing gets more difficult though)
                #       Implement that to avoid overhead from file system issues
                # TODO: Refactor and add comments
                for conformer in enumeration.get_conformers():
                    tmp_dir = self._prepare_tmp_input_dir()
                    conf_path = os.path.join(tmp_dir, _SE.CONFORMER_PATH)
                    ni_file = self.data.generic.get_files_by_extension("mol2")[0]
                    ni_file.write(tmp_dir)
                    conformer.write(conf_path)
                    self._execute_backend(conf_path, tmp_dir, ni_file.get_file_name())
                    self._parse_output(tmp_dir, conformer)
                    self._logger.log(
                        f"Finished shaep execution for conformer {enumeration.get_index_string()}.",
                        _LE.DEBUG,
                    )
                    number_rescored += 1
                    self._remove_temporary(tmp_dir)
        self._logger.log(f"Executed ShaEP for {number_rescored} conformers.", _LE.INFO)
