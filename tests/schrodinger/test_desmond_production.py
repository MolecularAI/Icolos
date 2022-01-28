from icolos.core.containers.generic import GenericData
import unittest
from icolos.core.workflow_steps.schrodinger.desmond_exec import StepDesmondExec
from icolos.utils.general.files_paths import attach_root_path
import os
from tests.tests_paths import PATHS_EXAMPLEDATA


from icolos.utils.enums.step_enums import StepBaseEnum, StepDesmondEnum

_SBE = StepBaseEnum
_SDE = StepDesmondEnum()


class Test_Desmond_Exec(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/schrodinger")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.DESMOND_SETUP_PDB), "rb") as f:
            self.pdb = f.read()

    def test_desmond_production(self):
        step_conf = {
            _SBE.STEPID: "test_desmond_setup",
            _SBE.STEP_TYPE: "desmond_preprocess",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws && $SCHRODINGER/jsc local-server-start"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {},
                _SBE.SETTINGS_ADDITIONAL: {_SDE.CFG_FIELDS: {"time": "1"}},
            },
        }

        step_desmond_exec = StepDesmondExec(**step_conf)
        step_desmond_exec.data.generic.add_file(
            GenericData(file_name="structure.pdb", file_data=self.pdb, argument=True)
        )
        step_desmond_exec.execute()

        out_path = os.path.join(self._test_dir, "out.cms")
        step_desmond_exec.data.generic.write_out_all_files(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 23587000)
