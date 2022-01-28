import unittest
import os
from icolos.core.workflow_steps.schrodinger.fep_plus_setup import StepFepPlusSetup
from icolos.utils.enums.step_enums import StepBaseEnum, StepGlideEnum, StepFepPlusEnum
from tests.tests_paths import PATHS_1UYD
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_docked_ligands_as_conformers,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path, empty_output_dir

_SBE = StepBaseEnum
_SGE = StepGlideEnum()
_SFE = StepFepPlusEnum()


class Test_FepPlusSetup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/fep_plus")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.FEP_PLUS_DOCKING_PV, "rb") as f:
            self.poseviewer = f.read()
        self.mol1 = get_docked_ligands_as_conformers(
            PATHS_1UYD.LIG4_POSES, poseviewer=self.poseviewer
        )
        self.mol2 = get_ligands_as_compounds_with_conformers(
            PATHS_1UYD.LIG_SDF, poseviewer=self.poseviewer
        )
        empty_output_dir(self._test_dir)

    def test_fep_setup_with_xray(self):
        step_conf = {
            _SBE.STEPID: "test_fep_setup_with_xray",
            _SBE.STEP_TYPE: _SBE.STEP_FEP_PLUS_SETUP,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SFE.XRAY_STRUCTURES: PATHS_1UYD.XRAY_STRUCTURES
                },
            },
        }
        step_fep_plus_setup = StepFepPlusSetup(**step_conf)
        step_fep_plus_setup.data.compounds = self.mol2
        step_fep_plus_setup.execute()

        # now confirm that the map has been generated properly
        out_path = os.path.join(self._test_dir, "xray_test_out.fmp")
        step_fep_plus_setup.write_generic_by_extension(
            path=os.path.join(self._test_dir, "xray_test_out.fmp"),
            ext="fmp",
            join=False,
        )
        stat_inf = os.stat(out_path)
        self.assertAlmostEqual(stat_inf.st_size, 821966, delta=500)

    def test_fep_setup(self):
        step_conf = {
            _SBE.STEPID: "test_fep_setup",
            _SBE.STEP_TYPE: _SBE.STEP_FEP_PLUS_SETUP,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                }
            },
        }

        step_fep_plus_setup = StepFepPlusSetup(**step_conf)
        step_fep_plus_setup.data.compounds = self.mol1
        step_fep_plus_setup.execute()

        # now confirm that the map has been generated properly
        out_path = os.path.join(self._test_dir, "test_out.fmp")
        step_fep_plus_setup.write_generic_by_extension(
            path=os.path.join(self._test_dir, "test_out.fmp"), ext="fmp", join=False
        )
        stat_inf = os.stat(out_path)
        self.assertAlmostEqual(stat_inf.st_size, 848697, delta=500)
