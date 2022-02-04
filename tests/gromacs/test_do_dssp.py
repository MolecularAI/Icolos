from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.workflow_steps.gromacs.do_dssp import StepGMXDoDSSP
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
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR), "rb") as f:
            self.structure = f.read()
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC), "rb") as f:
            self.traj = f.read()

    def test_editconf_run(self):
        step_conf = {
            _SBE.STEPID: "test_dssp",
            _SBE.STEP_TYPE: "dssp",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                }
            },
        }

        step_do_dssp = StepGMXDoDSSP(**step_conf)
        step_do_dssp.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.structure)
        )
        step_do_dssp.data.generic.add_file(
            GenericData(file_name="traj.xtc", file_data=self.traj)
        )
        step_do_dssp.execute()
        out_path = os.path.join(self._test_dir, "info.dat")
        step_do_dssp.write_generic_by_name(self._test_dir, "info.dat")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 22377)
