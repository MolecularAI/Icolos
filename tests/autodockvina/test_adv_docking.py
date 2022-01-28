import os
import unittest

from icolos.core.workflow_steps.autodockvina.docking import StepAutoDockVina

from icolos.utils.enums.step_enums import StepBaseEnum, StepAutoDockVinaEnum
from icolos.utils.enums.program_parameters import AutoDockVinaEnum

from tests.tests_paths import PATHS_1UYD, get_1UYD_ligands_as_Compounds
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SAE = StepAutoDockVinaEnum()
_EE = AutoDockVinaEnum()


class Test_ADV_docking(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/ADV")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._1UYD_compounds = get_1UYD_ligands_as_Compounds(
            abs_path=PATHS_1UYD.LIGANDS
        )
        self.receptor_path = PATHS_1UYD.PDBQT_PATH

    def test_ADV_run(self):
        step_conf = {
            _SBE.STEPID: "01_ADV",
            _SBE.STEP_TYPE: _SBE.STEP_AUTODOCKVINA_DOCKING,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load AutoDock_Vina",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 4},
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SAE.CONFIGURATION: {
                        _SAE.ADV_SEARCH_SPACE: {
                            _SAE.ADV_SEARCH_SPACE_CENTER_X: 3.3,
                            _SAE.ADV_SEARCH_SPACE_CENTER_Y: 11.5,
                            _SAE.ADV_SEARCH_SPACE_CENTER_Z: 24.8,
                            _SAE.ADV_SEARCH_SPACE_SIZE_Y: 10,
                            _SAE.ADV_SEARCH_SPACE_SIZE_Z: 10,
                        },
                        _SAE.NUMBER_POSES: 2,
                        _SAE.ADV_RECEPTOR_PATH: self.receptor_path,
                    }
                },
            },
        }

        adv_step = StepAutoDockVina(**step_conf)
        adv_step.data.compounds = self._1UYD_compounds

        adv_step.execute()
        self.assertEqual(len(adv_step.get_compounds()), 15)
        self.assertEqual(len(adv_step.get_compounds()[0][0].get_conformers()), 2)
        self.assertListEqual(
            list(
                adv_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.953, 10.113, 25.16],
        )
        self.assertListEqual(
            list(
                adv_step.get_compounds()[14][0][1]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [3.682, 15.421, 26.244],
        )
        self.assertEqual(
            adv_step.get_compounds()[0][0][0]
            .get_molecule()
            .GetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE),
            "-9.1",
        )

        # check SDF write-out
        out_path = os.path.join(self._test_dir, "adv_docked.sdf")
        adv_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 105000)
