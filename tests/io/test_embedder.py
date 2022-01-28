import unittest

from icolos.core.workflow_steps.io.embedder import StepEmbedding
from icolos.utils.enums.step_enums import StepBaseEnum, StepEmbeddingEnum

from tests.tests_paths import PATHS_EXAMPLEDATA

_SBE = StepBaseEnum
_SEE = StepEmbeddingEnum()


class Test_Embedder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self._SMI_path = PATHS_EXAMPLEDATA.MEDIUM_MOLECULES_SMI_PATH

    @classmethod
    def tearDownClass(cls):
        pass

    def test_embed_with_RDkit_no_protonation(self):
        step_conf = {
            _SBE.STEPID: "01_embed_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_EMBEDDING,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _SEE.RDKIT_PROTONATE: False,
                        _SEE.METHOD: _SEE.METHOD_RDKIT,
                    }
                }
            },
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: self._SMI_path,
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                        _SBE.INPUT_FORMAT: _SBE.FORMAT_SMI,
                    }
                ]
            },
        }
        init_step = StepEmbedding(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertListEqual(
            list(
                init_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.9314762660385534, 0.06628711293694872, 4.923008037397455],
        )
        self.assertListEqual(
            list(
                init_step.get_compounds()[1][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.6176730474256593, 0.37859007619202606, 0.6065857814585477],
        )
        self.assertEqual(
            22, init_step.get_compounds()[0][0].get_molecule().GetNumAtoms()
        )

    def test_embed_with_RDkit_protonation(self):
        step_conf = {
            _SBE.STEPID: "01_embed_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_EMBEDDING,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _SEE.RDKIT_PROTONATE: True,
                        _SEE.METHOD: _SEE.METHOD_RDKIT,
                    }
                }
            },
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: self._SMI_path,
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                        _SBE.INPUT_FORMAT: _SBE.FORMAT_SMI,
                    }
                ]
            },
        }
        init_step = StepEmbedding(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertListEqual(
            list(
                init_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.9314762660385534, 0.06628711293694872, 4.923008037397455],
        )
        self.assertListEqual(
            list(
                init_step.get_compounds()[1][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.6176730474256593, 0.37859007619202606, 0.6065857814585477],
        )
        self.assertEqual(
            41, init_step.get_compounds()[0][0].get_molecule().GetNumAtoms()
        )
        self.assertListEqual(
            list(
                init_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[40]
            ),
            [-3.576148794943472, -0.8051119546399829, -0.9424118920903588],
        )
