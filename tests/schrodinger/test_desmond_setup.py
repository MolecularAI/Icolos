from icolos.core.containers.generic import GenericData
import unittest
from icolos.core.workflow_steps.schrodinger.desmond_preprocessor import StepDesmondSetup
from icolos.utils.general.files_paths import attach_root_path
import os
from tests.tests_paths import PATHS_EXAMPLEDATA


from icolos.utils.enums.step_enums import StepBaseEnum, StepDesmondEnum

_SBE = StepBaseEnum
_SDE = StepDesmondEnum()


class Test_Desmond_Setup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/schrodinger")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.DESMOND_SETUP_PDB), "r") as f:
            self.pdb = f.read()

    def test_desmond_preprocess(self):
        step_conf = {
            _SBE.STEPID: "test_desmond_setup",
            _SBE.STEP_TYPE: "desmond_preprocess",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {},
                _SBE.SETTINGS_ADDITIONAL: {_SDE.MSJ_FIELDS: {}},
            },
        }

        step_desmond_preprocess = StepDesmondSetup(**step_conf)
        step_desmond_preprocess.data.generic.add_file(
            GenericData(file_name="structure.pdb", file_data=self.pdb, argument=True)
        )
        step_desmond_preprocess.execute()

        out_path = os.path.join(self._test_dir, "setup.cms")
        step_desmond_preprocess.data.generic.write_out_all_files(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 22560500)
