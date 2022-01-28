from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.gromacs.solvate import StepGMXSolvate

_SBE = StepBaseEnum
_SGE = StepGromacsEnum()


class Test_Solvate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP, "r") as f:
            self.topol = f.read()
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO, "r") as f:
            self.structure = f.read()

    def test_solvate(self):
        step_conf = {
            _SBE.STEPID: "test_solvate",
            _SBE.STEP_TYPE: "solvate",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                }
            },
        }

        step_solvate = StepGMXSolvate(**step_conf)
        step_solvate.data.generic.add_file(
            GenericData(
                file_name="structure.gro", file_data=self.structure, argument=True
            )
        )
        step_solvate.data.generic.add_file(
            GenericData(file_name="topol.top", file_data=self.topol, argument=True)
        )

        step_solvate.execute()

        out_path = os.path.join(self._test_dir, "structure.gro")
        step_solvate.write_generic_by_extension(
            self._test_dir, _SGE.FIELD_KEY_STRUCTURE
        )
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 650000)
