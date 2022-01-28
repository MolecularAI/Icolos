from icolos.core.containers.generic import GenericData
import unittest
from icolos.core.workflow_steps.structure_prediction.pdb_fixer import StepPdbFixer
from icolos.utils.general.files_paths import attach_root_path
import os
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class TestPdbFixer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/structure_prediction")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.PANTHER_RECEPTOR_PDB), "r") as f:
            self.pdb = f.read()

    def test_pdb_fixer_default(self):
        step_conf = {
            _SBE.STEPID: "01_PDB_FIXER",
            _SBE.STEP_TYPE: _SBE.STEP_PDB_FIXER,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {_SBE.SETTINGS_ARGUMENTS: {}, _SBE.SETTINGS_ADDITIONAL: {}},
        }
        step_pdb_fixer = StepPdbFixer(**step_conf)
        test_pdb = GenericData(file_name="test.pdb", file_data=self.pdb)
        step_pdb_fixer.data.generic.add_file(test_pdb)
        step_pdb_fixer.execute()

        out_path = os.path.join(self._test_dir, "test.pdb")
        step_pdb_fixer.write_generic_by_extension(path=self._test_dir, ext="pdb")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 738000)

    def test_pdb_fixer(self):
        step_conf = {
            _SBE.STEPID: "01_PDB_FIXER",
            _SBE.STEP_TYPE: _SBE.STEP_PDB_FIXER,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "--keep-heterogens": "water",
                        "--ph": "4.0",
                    },
                }
            },
        }
        step_pdb_fixer = StepPdbFixer(**step_conf)
        test_pdb = GenericData(file_name="test_2.pdb", file_data=self.pdb)
        step_pdb_fixer.data.generic.add_file(test_pdb)
        step_pdb_fixer.execute()

        out_path = os.path.join(self._test_dir, "test_2.pdb")
        step_pdb_fixer.write_generic_by_extension(path=self._test_dir, ext="pdb")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 710000)
