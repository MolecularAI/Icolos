import unittest
from icolos.core.workflow_steps.calculation.electrostatics.esp_sim import StepEspSim

from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import export_unit_test_env_vars

_SBE = StepBaseEnum()


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

        esp_sim_score = [
            0.8112564566774974,
            0.7940316946620978,
            0.8157010968264732,
            0.6927039160490105,
            0.6709748529493742,
            0.3780220716995563,
            0.7933792682013576,
            0.7672803082385128,
        ]

        shape_sim_score = [
            0.6419844502036283,
            0.9525606469002695,
            0.5686465433300876,
            0.5986955029179539,
            0.5460218408736349,
            0.5232662864004803,
            0.8305164319248827,
            0.7283643892339544,
        ]

        for i in range(len(esp_sim_score)):
            self.assertEqual(
                step_esp_sim.data.compounds[i]
                .get_enumerations()[0]
                .get_conformers()[0]
                .get_molecule()
                .GetProp("esp_sim"),
                str(esp_sim_score[i]),
            )
            self.assertEqual(
                step_esp_sim.data.compounds[i]
                .get_enumerations()[0]
                .get_conformers()[0]
                .get_molecule()
                .GetProp("shape_sim"),
                str(shape_sim_score[i]),
            )
