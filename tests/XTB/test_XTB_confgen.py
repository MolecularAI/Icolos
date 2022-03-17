import unittest
import os

from icolos.core.workflow_steps.confgen.xtb import StepXTB

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.enums.program_parameters import XTBEnum

from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    MAIN_CONFIG,
    export_unit_test_env_vars,
    get_mol_as_Compound,
    get_ligands_as_compounds_with_conformers,
    get_mol_as_Conformer,
)
from icolos.utils.general.files_paths import attach_root_path
import time

_SBE = StepBaseEnum
_CE = XTBEnum()


class Test_XTB_confgen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/XTB")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        self._paracetamol_molecule = get_mol_as_Compound(
            PATHS_EXAMPLEDATA.PARACETAMOL_PATH
        )
        self._aspirin_molecule = get_mol_as_Compound(PATHS_EXAMPLEDATA.ASPIRIN_PATH)
        self._medium_molecules = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.SMALL_MOLECULES_SDF_PATH
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_coordinate_generation(self):
        step_conf = {
            _SBE.STEPID: "01_conf_genXTB",
            _SBE.STEP_TYPE: _SBE.STEP_XTB,
            _SBE.EXEC: {
                _SBE.EXEC_BINARYLOCATION: attach_root_path(
                    os.path.join(MAIN_CONFIG["XTBHOME"])
                ),
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.XTB_OPT: "vtight",
                        _CE.XTB_GBSA: "h2o",
                    },
                }
            },
        }
        xtb_step = StepXTB(**step_conf)
        xtb_step.data.compounds = [self._paracetamol_molecule]
        confs = get_mol_as_Conformer(
            attach_root_path(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        )
        xtb_step.data.compounds[0][0].add_conformers(confs, auto_update=True)
        self.assertListEqual(
            list(
                xtb_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        xtb_step.execute()
        self.assertListEqual(
            list(
                xtb_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.4367, 13.0798, 24.5152],
        )

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertEqual(len(xtb_step.get_compounds()[0][0]), 11)

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(
            self._test_dir, "XTB_conformers_from_OMEGA_paracetamol.sdf"
        )
        xtb_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 17768)

    def test_single_core_execution(self):
        step_conf = {
            _SBE.STEPID: "01_conf_genXTB",
            _SBE.STEP_TYPE: _SBE.STEP_XTB,
            _SBE.EXEC: {
                _SBE.EXEC_BINARYLOCATION: attach_root_path(
                    os.path.join(MAIN_CONFIG["XTBHOME"])
                ),
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.XTB_OPT: "vtight",
                        _CE.XTB_GBSA: "h2o",
                    },
                }
            },
        }
        xtb_step = StepXTB(**step_conf)
        xtb_step.data.compounds = self._medium_molecules
        self.assertListEqual(
            list(
                xtb_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.8851, -1.0363, -0.1124],
        )
        xtb_step.execute()
        self.assertListEqual(
            list(
                xtb_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.8526, -0.9638, -0.1394],
        )

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertEqual(len(xtb_step.get_compounds()[0][0]), 1)
        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(
            self._test_dir, "XTB_conformers_from_OMEGA_paracetamol.sdf"
        )
        xtb_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 3967)

    def test_parallel_execution(self):
        step_conf = {
            _SBE.STEPID: "01_conf_genXTB",
            _SBE.STEP_TYPE: _SBE.STEP_XTB,
            _SBE.EXEC: {
                _SBE.EXEC_BINARYLOCATION: attach_root_path(
                    os.path.join(MAIN_CONFIG["XTBHOME"])
                ),
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 8},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _CE.XTB_OPT: "vtight",
                        _CE.XTB_GBSA: "h2o",
                    },
                }
            },
        }
        xtb_step = StepXTB(**step_conf)
        xtb_step.data.compounds = self._medium_molecules
        self.assertListEqual(
            list(
                xtb_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.8851, -1.0363, -0.1124],
        )
        t1 = time.time()
        xtb_step.execute()
        t2 = time.time()
        self.assertLess(t2 - t1, 4)
        self.assertListEqual(
            list(
                xtb_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.8526, -0.9638, -0.1394],
        )

        # check number of conformers returned (only one Compound with only one Enumeration)
        self.assertEqual(len(xtb_step.get_compounds()[0][0]), 1)
        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(
            self._test_dir, "XTB_conformers_from_OMEGA_paracetamol.sdf"
        )
        xtb_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 3967)
