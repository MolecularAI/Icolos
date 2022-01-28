from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.workflow_steps.gromacs.editconf import StepGMXEditConf
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_Editconf(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_STRUCTURE_FILE), "r") as f:
            self.structure = f.read()

    def test_editconf_run(self):
        step_conf = {
            _SBE.STEPID: "test_editconf",
            _SBE.STEP_TYPE: "editconf",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-d": "1.0",
                        "-bt": "dodecahedron",
                    },
                }
            },
        }

        step_editconf = StepGMXEditConf(**step_conf)
        step_editconf.data.generic.add_file(
            GenericData(
                file_name="structure.gro", file_data=self.structure, argument=True
            )
        )
        step_editconf.execute()
        out_path = os.path.join(self._test_dir, "structure.gro")
        step_editconf.write_generic_by_name(self._test_dir, "structure.gro")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 22377)
