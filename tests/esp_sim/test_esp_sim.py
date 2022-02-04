import unittest
from icolos.core.workflow_steps.calculation.electrostatics.esp_sim import StepEspSim

from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import export_unit_test_env_vars

_SBE = StepBaseEnum


class Test_EspSim(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        export_unit_test_env_vars()

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_initialize_compound_SDF(self):
        step_conf = {
            _SBE.STEPID: "01_esp_sim",
            _SBE.STEP_TYPE: _SBE.STEP_ESP_SIM,
            _SBE.EXEC: {
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 3},
            },
            _SBE.SETTINGS: {_SBE.SETTINGS_ADDITIONAL: {"ref_smiles": "C(C(C(=O)O)O)O"}},
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: "C1=CC=C(C=C1)C(C(=O)O)O",
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STRING,
                    }
                ]
            },
        }
        step_esp_sim = StepEspSim(**step_conf)
        step_esp_sim.generate_input()
        step_esp_sim.execute()

        esp_sim_score = [0.811]

        shape_sim_score = [0.642]

        for i in range(len(esp_sim_score)):
            self.assertEqual(
                round(
                    float(
                        step_esp_sim.data.compounds[i]
                        .get_enumerations()[0]
                        .get_conformers()[0]
                        .get_molecule()
                        .GetProp("esp_sim")
                    ),
                    ndigits=3,
                ),
                esp_sim_score[i],
            )
            self.assertEqual(
                round(
                    float(
                        step_esp_sim.data.compounds[i]
                        .get_enumerations()[0]
                        .get_conformers()[0]
                        .get_molecule()
                        .GetProp("shape_sim")
                    ),
                    ndigits=3,
                ),
                shape_sim_score[i],
            )
