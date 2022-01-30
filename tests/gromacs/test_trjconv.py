from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.gromacs.trjconv import StepGMXTrjconv

SGE = StepGromacsEnum()
SBE = StepBaseEnum


class Test_Trjconv(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC, "rb") as f:
            self.xtc = f.read()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR, "rb") as f:
            self.tpr = f.read()

    def test_trjconv(self):
        step_conf = {
            SBE.STEPID: "test_trjconv",
            SBE.STEP_TYPE: "trjconv",
            SBE.EXEC: {
                SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            SBE.SETTINGS: {
                SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"],
                SBE.SETTINGS_ADDITIONAL: {SBE.PIPE_INPUT: "echo -ne 1 0"},
            },
        }

        step_trjconv = StepGMXTrjconv(**step_conf)
        step_trjconv.data.generic.add_file(
            GenericData(file_name="structure.xtc", file_data=self.xtc, argument=True)
        )
        step_trjconv.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.tpr, argument=True)
        )
        step_trjconv.execute()
        out_path = os.path.join(self._test_dir, "structure.xtc")
        step_trjconv.write_generic_by_extension(self._test_dir, "xtc")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 607000)
