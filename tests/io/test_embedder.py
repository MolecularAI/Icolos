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
        self._SMI_path = PATHS_EXAMPLEDATA.SMALL_MOLECULES_SMI_PATH

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
            [-3.676505807445281, -0.44491027005777944, 0.9478288681339868],
        )
        self.assertListEqual(
            list(
                init_step.get_compounds()[1][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.1635386441070907, 0.781887672700409, 2.383883775746168],
        )
        self.assertEqual(
            11, init_step.get_compounds()[0][0].get_molecule().GetNumAtoms()
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
            [-3.676505807445281, -0.44491027005777944, 0.9478288681339868],
        )
        self.assertListEqual(
            list(
                init_step.get_compounds()[1][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.1635386441070907, 0.781887672700409, 2.383883775746168],
        )
        self.assertEqual(
            20, init_step.get_compounds()[0][0].get_molecule().GetNumAtoms()
        )
        self.assertListEqual(
            list(
                init_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[18]
            ),
            [3.294854134599972, -0.7589589232493622, -0.4334701745138959],
        )
