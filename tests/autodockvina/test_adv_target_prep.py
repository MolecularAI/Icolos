import unittest
import os

from icolos.core.workflow_steps.autodockvina.target_preparation import (
    StepAutoDockVinaTargetPreparation,
)
from icolos.utils.enums.step_enums import (
    StepBaseEnum,
    StepAutoDockVinaTargetPreparationEnum,
)
from icolos.utils.general.files_paths import attach_root_path
from tests.tests_paths import PATHS_1UYD

_SBE = StepBaseEnum
_SAE = StepAutoDockVinaTargetPreparationEnum()


class Test_ADV_target_preparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/ADV_target_prep")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)
        cls.receptor_output_path = os.path.join(cls._test_dir, "ADV_receptor.pdbqt")

    def setUp(self):
        self.receptor_input_path = PATHS_1UYD.PDB_PATH
        self.reference_ligand_sdf_path = PATHS_1UYD.NATIVE_LIGAND_SDF
        self.reference_ligand_pdb_path = PATHS_1UYD.NATIVE_LIGAND_PDB

    @classmethod
    def tearDownClass(cls):
        pass

    def test_extract_box(self):
        step_conf = {
            _SBE.STEPID: "01_ADV",
            _SBE.STEP_TYPE: _SBE.STEP_AUTODOCKVINA_TARGET_PREPARATION,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SAE.INPUT_RECEPTOR_PDB: self.receptor_input_path,
                    _SAE.OUTPUT_RECEPTOR_PDBQT: self.receptor_output_path,
                    _SAE.EXTRACT_BOX: {
                        _SAE.EXTRACT_BOX_REFERENCE_LIGAND_PATH: self.reference_ligand_sdf_path,
                        _SAE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT: _SAE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_SDF,
                    },
                },
            },
        }

        adv_tp_step = StepAutoDockVinaTargetPreparation(**step_conf)
        x_coords, y_coords, z_coords = adv_tp_step._extract_box()

        self.assertEqual(len(x_coords), 28)
        self.assertListEqual([4.403, 5.122, 5.091], x_coords[:3])
        self.assertListEqual([15.528, 15.084, 13.786], y_coords[:3])
        self.assertListEqual([26.579, 25.453, 24.846], z_coords[:3])

    def test_target_preparation(self):
        step_conf = {
            _SBE.STEPID: "01_ADV",
            _SBE.STEP_TYPE: _SBE.STEP_AUTODOCKVINA_TARGET_PREPARATION,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SAE.INPUT_RECEPTOR_PDB: self.receptor_input_path,
                    _SAE.OUTPUT_RECEPTOR_PDBQT: self.receptor_output_path,
                    _SAE.EXTRACT_BOX: {
                        _SAE.EXTRACT_BOX_REFERENCE_LIGAND_PATH: self.reference_ligand_pdb_path,
                        _SAE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT: _SAE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_PDB,
                    },
                },
            },
        }

        adv_tp_step = StepAutoDockVinaTargetPreparation(**step_conf)
        adv_tp_step.execute()

        # check SDF write-out
        stat_inf = os.stat(self.receptor_output_path)
        self.assertGreater(stat_inf.st_size, 290000)
