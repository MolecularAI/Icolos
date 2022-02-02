from icolos.utils.enums.program_parameters import PantherEnum
import os
import unittest
from tests.tests_paths import PATHS_1UYD, PATHS_EXAMPLEDATA, MAIN_CONFIG

from icolos.utils.enums.step_enums import StepBaseEnum, StepPantherEnum
from icolos.core.workflow_steps.calculation.panther import StepPanther
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SPE = StepPantherEnum()
_PE = PantherEnum()


class Test_Panther(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/panther")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def test_panther_run(self):
        step_conf = {
            _SBE.STEPID: "01_panther",
            _SBE.STEP_TYPE: _SBE.STEP_PANTHER,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SPE.PANTHER_LOCATION: MAIN_CONFIG["PANTHER_LOCATION"],
                    _SPE.PANTHER_CONFIG_FILE: attach_root_path(
                        PATHS_EXAMPLEDATA.PANTHER_CONFIG
                    ),
                    _SPE.FIELDS: {"1-Pdb file": PATHS_EXAMPLEDATA.PANTHER_HOLO_PDB},
                }
            },
        }
        panther_step = StepPanther(**step_conf)
        panther_step.execute()

        # check we get the negative image back
        out_path = os.path.join(self._test_dir, "neg_image.mol2")
        panther_step.write_generic_by_extension(self._test_dir, "mol2")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 12598)
