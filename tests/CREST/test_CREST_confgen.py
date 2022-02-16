import unittest
import os

from icolos.core.workflow_steps.confgen.crest import StepCREST

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.enums.program_parameters import CrestEnum

from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    MAIN_CONFIG,
    export_unit_test_env_vars,
    get_mol_as_Compound,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_CE = CrestEnum()


class Test_CREST_confgen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/CREST")
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
            _SBE.STEPID: "01_conf_gen_crest",
            _SBE.STEP_TYPE: _SBE.STEP_CREST,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: None,
                _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["CREST_BINARY_LOCATION"],
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-niceprint"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.CREST_OPT: "normal",
                        _CE.CREST_G: "h2o",
                        _CE.CREST_RTHR: 0.5,
                        _CE.CREST_ETHR: 0.25,
                        _CE.CREST_EWIN: 8.0,
                        _CE.CREST_PTHR: 0.4,
                        _CE.CREST_BTHR: 0.02,
                    },
                }
            },
        }
        crest_step = StepCREST(**step_conf)
        crest_step.data.compounds = [self._paracetamol_molecule]
        crest_step.execute()

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertGreaterEqual(len(crest_step.get_compounds()[0][0]), 14)

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "CREST_conformers_paracetamol.sdf")
        crest_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 28000)

    def test_coordinate_generation_charged(self):
        step_conf = {
            _SBE.STEPID: "01_conf_gen_crest",
            _SBE.STEP_TYPE: _SBE.STEP_CREST,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: None,
                _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["CREST_BINARY_LOCATION"],
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-niceprint"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.CREST_OPT: "normal",
                        _CE.CREST_G: "h2o",
                        _CE.CREST_RTHR: 0.5,
                        _CE.CREST_ETHR: 0.25,
                        _CE.CREST_EWIN: 8.0,
                        _CE.CREST_PTHR: 0.4,
                        _CE.CREST_BTHR: 0.02,
                    },
                }
            },
        }

        # check number of conformers returned
        crest_step = StepCREST(**step_conf)
        crest_step.data.compounds = [self._aspirin_molecule]
        crest_step.execute()

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertGreaterEqual(len(crest_step.get_compounds()[0][0]), 2)

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "CREST_conformers_aspirin.sdf")
        crest_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        print(stat_inf.st_size)
        self.assertGreater(stat_inf.st_size, 3200)
