import unittest
import os

from icolos.core.workflow_steps.confgen.omega import StepOmega

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.enums.program_parameters import OMEGAEnum

from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
    get_mol_as_Compound,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_CE = OMEGAEnum()


class Test_OMEGA_confgen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls._test_dir = attach_root_path("tests/junk/OMEGA")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        self._paracetamol_molecule = get_mol_as_Compound(
            PATHS_EXAMPLEDATA.PARACETAMOL_PATH
        )
        self._aspirin_molecule = get_mol_as_Compound(PATHS_EXAMPLEDATA.ASPIRIN_PATH)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_coordinate_generation_neutral(self):
        step_conf = {
            _SBE.STEPID: "01_conf_gen_omega",
            _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.CLASSIC_MAXCONFS: 50,
                        _CE.CLASSIC_RMS: 0.05,
                    },
                }
            },
        }
        omega_step = StepOmega(**step_conf)
        omega_step.data.compounds = [self._paracetamol_molecule]
        omega_step.execute()

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertEqual(len(omega_step.get_compounds()[0][0]), 1)

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "OMEGA_conformers_paracetamol.sdf")
        omega_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 1878)

    def test_coordinate_generation_neutral_high_RMS(self):
        step_conf = {
            _SBE.STEPID: "01_conf_gen_omega",
            _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.CLASSIC_MAXCONFS: 10,
                        _CE.CLASSIC_RMS: 0.7,
                    },
                }
            },
        }
        omega_step = StepOmega(**step_conf)
        omega_step.data.compounds = [self._paracetamol_molecule]
        omega_step.execute()

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertEqual(len(omega_step.get_compounds()[0][0]), 1)

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(
            self._test_dir, "OMEGA_conformers_paracetamol_highRMS.sdf"
        )
        omega_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 1878)

    def test_coordinate_generation_charged(self):
        step_conf = {
            _SBE.STEPID: "01_conf_gen_omega",
            _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.CLASSIC_MAXCONFS: 10,
                        _CE.CLASSIC_RMS: 0.0,
                    },
                }
            },
        }
        omega_step = StepOmega(**step_conf)
        omega_step.data.compounds = [self._aspirin_molecule]
        omega_step.execute()

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertEqual(len(omega_step.get_compounds()[0][0]), 10)

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "OMEGA_conformers_aspirin.sdf")
        omega_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 19574)
