import unittest

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard
from icolos.core.containers.generic import GenericData
from tests.tests_paths import (
    MAIN_CONFIG,
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
)
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum, StepPrepwizEnum
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_SPE = StepPrepwizEnum()


class Test_ResiScanning(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/residue_scanning")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(PATHS_1UYD.PDB_PATH, "r") as f:
            apo_1uyd = f.read()
        self.apo_1uyd = GenericData(file_name="test_structure.pdb", file_data=apo_1uyd)
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE, "r") as f:
            self.holo_1uyd = f.read()

    def test_prepwizard(self):
        step_conf = {
            _SBE.STEPID: "01_ligprep",
            _SBE.STEP_TYPE: _SBE.STEP_PREPWIZARD,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: MAIN_CONFIG["SCHRODINGER_MODULE"],
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_PARAMETERS: {}}
            },
        }

        prepwiz_step = StepPrepwizard(**step_conf)
        prepwiz_step.data.generic.add_file(self.apo_1uyd)
        prepwiz_step.execute()

        out_file = prepwiz_step.data.generic.get_files_by_extension("pdb")[0].get_data()
        out_path = os.path.join(self._test_dir, "test_out.pdb")
        with open(out_path, "w") as f:
            f.write(out_file)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 155800)
