from icolos.core.containers.generic import GenericData
import unittest
from icolos.core.workflow_steps.structure_prediction.dssp import StepDSSP
from icolos.utils.general.files_paths import attach_root_path
import os
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class TestDSSP(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/structure_prediction")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.DSSP_PDB_1), "r") as f:
            self.pdb1 = f.read()

        with open(attach_root_path(PATHS_EXAMPLEDATA.DSSP_PDB_2), "r") as f:
            self.pdb2 = f.read()

        with open(attach_root_path(PATHS_EXAMPLEDATA.DSSP_PDB_3), "r") as f:
            self.pdb3 = f.read()

    def test_dssp(self):
        step_conf = {
            _SBE.STEPID: "01_DSSP",
            _SBE.STEP_TYPE: _SBE.STEP_DSSP,
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load DSSP"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"--output-format": "dssp"}
                },
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }

        step_dssp = StepDSSP(**step_conf)
        pdb1 = GenericData(file_name="test_1.pdb", file_data=self.pdb1)
        pdb2 = GenericData(file_name="test_2.pdb", file_data=self.pdb2)
        pdb3 = GenericData(file_name="test_3.pdb", file_data=self.pdb3)
        step_dssp.data.generic.add_files([pdb1, pdb2, pdb3])
        step_dssp.execute()

        out_path = os.path.join(self._test_dir, "dssp_output_test_1.txt")
        step_dssp.write_generic_by_name(self._test_dir, "dssp_output_test_1.txt")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 1234)
