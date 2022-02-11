import os
import unittest

from icolos.core.workflow_steps.ccdc.docking import StepGold

from icolos.utils.enums.step_enums import StepBaseEnum, StepGoldEnum

from tests.tests_paths import (
    get_1UYD_ligands_as_Compounds,
    PATHS_1UYD,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SGE = StepGoldEnum()


class Test_Gold_docking(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/Gold")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._1UYD_compounds = get_1UYD_ligands_as_Compounds(
            abs_path=PATHS_1UYD.LIGANDS
        )[:3]
        self.cavity_file = attach_root_path(os.path.join("../IcolosData", PATHS_1UYD.GOLD_CAVITY_MOL2))
        self.receptor_path = attach_root_path(os.path.join("../IcolosData", PATHS_1UYD.GOLD_MOL2_PROTEIN))

    def test_config_generation(self):
        step_conf = {
            _SBE.STEPID: "01_Gold",
            _SBE.STEP_TYPE: _SBE.STEP_GOLD_DOCKING,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load ccdc"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.CONFIGURATION: {
                        _SGE.AUTOMATIC_SETTINGS: {
                            _SGE.AUTOSCALE: 0.5
                        },
                        _SGE.FLOOD_FILL: {
                            _SGE.CAVITY_FILE: "whatever"
                        },
                        _SGE.PROTEIN_DATA: {
                            _SGE.PROTEIN_DATAFILE: "string"
                        }
                    }
                },
            },
        }

        gold_step = StepGold(**step_conf)

        # check the two mandatory settings and one "overwrite"
        self.assertEqual(gold_step.gold_additional.configuration.flood_fill.cavity_file, "whatever")
        self.assertEqual(gold_step.gold_additional.configuration.protein_data.protein_datafile, "string")
        self.assertEqual(gold_step.gold_additional.configuration.automatic_settings.autoscale, 0.5)

        # check config file size
        config_path = os.path.join(self._test_dir, "test_config.cfg")
        gold_step.generate_config_file(config_path,
                                       ["/a/path/to/an/SDF/with/a/compound.sdf",
                                        "/another/path.sdf"])
        stat_inf = os.stat(config_path)
        self.assertGreater(stat_inf.st_size, 1400)

    def test_docking(self):
        step_conf = {
            _SBE.STEPID: "01_Gold",
            _SBE.STEP_TYPE: _SBE.STEP_GOLD_DOCKING,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load ccdc"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.CONFIGURATION: {
                        _SGE.AUTOMATIC_SETTINGS: {
                            _SGE.AUTOSCALE: 0.5
                        },
                        _SGE.FLOOD_FILL: {
                            _SGE.CAVITY_FILE: self.cavity_file
                        },
                        _SGE.PROTEIN_DATA: {
                            _SGE.PROTEIN_DATAFILE: self.receptor_path
                        }
                    }
                },
            },
        }

        gold_step = StepGold(**step_conf)
        gold_step.data.compounds = self._1UYD_compounds
        gold_step.execute()

        self.assertEqual(len(gold_step.get_compounds()), 3)
        self.assertEqual(len(list(
                gold_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            )), 3)

        # check SDF write-out
        out_path = os.path.join(self._test_dir, "gold_docked.sdf")
        gold_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 50000)
