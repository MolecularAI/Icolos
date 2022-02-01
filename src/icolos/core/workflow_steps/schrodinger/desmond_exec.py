import os
from icolos.core.step_utils.structconvert import StructConvert
from pydantic import BaseModel
from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.enums.step_enums import StepDesmondEnum
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum

_SDE = StepDesmondEnum()
_SEE = SchrodingerExecutablesEnum()


class StepDesmondExec(StepSchrodingerBase, BaseModel):
    """
    Executes a full Desmond multisim workflow
    """

    class Config:
        underscore_attrs_are_private = True
        arbitrary_types_allowed = True

    _struct_converter: StructConvert = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(SchrodingerExecutor)
        self._check_backend_availability()
        self._struct_converter = StructConvert(
            binary_location=self.execution.binary_location,
            prefix_execution=self.execution.prefix_execution,
        )

    def execute(self):
        # takes in the cms file from the preprocessor and runs the full multisim workflow on it
        tmp_dir = self._make_tmpdir()
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

        preprocess_defaults = {
            "-HOST": "localhost",
            "-JOBNAME": "desmond_md_job_1",
            "-m": "config.msj desmond_md_job_1.mae",
            "-o": "setup.cms",
        }
        arguments = self._parse_arguments(preprocess_defaults)
        # compile and write the msj to the tmpdir
        config_dict = (
            self.settings.additional[_SDE.SETUP_MSJ_FIELDS]
            if _SDE.SETUP_MSJ_FIELDS in self.settings.additional.keys()
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

        exec_defaults = {
            "-HOST": "localhost",
            "-JOBNAME": "desmond_production",
            "-maxjob": "1",
            "-cpu": "1",
            "-m": _SDE.PRODUCTION_MSJ,
            "-c": _SDE.PRODUCTION_CFG,
            "-description": '"Molecular Dynamics" setup.cms',
            "-mode": "umbrella",
            "-PROJ": tmp_dir,
            "-o": "out.cms",
            "-lic": _SDE.TOKEN_STR,
        }

        msj_config_dict = (
            self.settings.additional[_SDE.MSJ_FIELDS]
            if _SDE.MSJ_FIELDS in self.settings.additional.keys()
            else {}
        )

        cfg_config_dict = (
            self.settings.additional[_SDE.CFG_FIELDS]
            if _SDE.CFG_FIELDS in self.settings.additional.keys()
            else {}
        )

        # write the config files: msj for the full workflow, and a cfg for the production sim
        self._write_config(tmp_dir, msj_config_dict, _SDE.PRODUCTION_MSJ)

        self._write_config(tmp_dir, cfg_config_dict, _SDE.PRODUCTION_CFG)

        arguments = self._parse_arguments(exec_defaults)

        self._backend_executor.execute(
            command=_SEE.MULTISIM_EXEC,
            arguments=arguments,
            check=True,
            location=tmp_dir,
        )

        self._parse_output(tmp_dir)

        self._remove_temporary(tmp_dir)
