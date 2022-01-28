import os
from icolos.core.step_utils.structconvert import StructConvert
from pydantic import BaseModel
from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.enums.step_enums import StepDesmondEnum
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum

_SDE = StepDesmondEnum()
_SEE = SchrodingerExecutablesEnum()


class StepDesmondSetup(StepSchrodingerBase, BaseModel):
    """
    Run preprocessing step to generate system for Desmond simulation
    """

    _struct_converter: StructConvert = None

    class Config:
        underscore_attrs_are_private = True
        arbitrary_types_allowed = True

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(SchrodingerExecutor)
        self._check_backend_availability()
        self._struct_converter = StructConvert(
            binary_location=self.execution.binary_location,
            prefix_execution=self.execution.prefix_execution,
        )

    def execute(self):
        tmp_dir = self._make_tmpdir()

        # need to take a structure file, possibly preprocess if pdb
        # get the structure file and extract the file name
        structure = self.data.generic.get_argument_by_extension(
            "pdb", rtn_file_object=True
        )
        pdb_id = structure.get_file_name()
        structure.write(tmp_dir)
        # convert the pdb file to mae

        self._struct_converter.pdb2mae(
            os.path.join(tmp_dir, pdb_id),
            os.path.join(tmp_dir, "desmond_md_job_1.mae"),
        )

        defaults = {
            "-HOST": "localhost",
            "-JOBNAME": "desmond_md_job_1",
            "-m": "config.msj desmond_md_job_1.mae",
            "-o": "setup.cms",
        }
        arguments = self._parse_arguments(defaults)
        # compile and write the msj to the tmpdir
        config_dict = (
            self.settings.additional[_SDE.MSJ_FIELDS]
            if _SDE.MSJ_FIELDS in self.settings.additional.keys()
            else {}
        )

        self._write_config(tmp_dir, dict_=config_dict, file_name=_SDE.PREPROCESS_MSJ)

        # execute
        self._backend_executor.execute(
            command=_SEE.MULTISIM_EXEC,
            arguments=arguments,
            check=True,
            location=tmp_dir,
        )
        self._parse_output(tmp_dir)

        self._remove_temporary(tmp_dir)
