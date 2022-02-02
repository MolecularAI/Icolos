import unittest

from icolos.core.workflow_steps.io.initialize_compound import StepInitializeCompound
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars

_SBE = StepBaseEnum


class Test_InitializeCompound(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        export_unit_test_env_vars()

    def setUp(self):
        self._paracetamol_path = PATHS_EXAMPLEDATA.PARACETAMOL_PATH
        self._SMI_path = PATHS_EXAMPLEDATA.SMALL_MOLECULES_SMI_PATH
        self._JSON_path = PATHS_EXAMPLEDATA.SMALL_MOLECULES_JSON_PATH
        self._CSV_path = PATHS_EXAMPLEDATA.SMALL_MOLECULES_CSV_PATH
        self._CSV_path_semicolon = (
            PATHS_EXAMPLEDATA.SMALL_MOLECULES_CSV_PATH_DELIMITER_SEMICOLON
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_initialize_compound_SDF(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: self._paracetamol_path,
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                        _SBE.INPUT_FORMAT: _SBE.FORMAT_SDF,
                    }
                ]
            },
        }
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 1)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 1)

        self.assertListEqual(
            list(
                init_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.8851, -1.0363, -0.1124],
        )

    def test_initialize_compound_SMI(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
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
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertEqual(
            init_step.get_compounds()[0][0].get_smile(),
            "CC(=O)Nc1ccc(O)cc1",
        )
        self.assertEqual(init_step.get_compounds()[1].get_name(), "Aspirin")

    def test_initialize_compound_JSON(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: self._JSON_path,
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                        _SBE.INPUT_FORMAT: _SBE.FORMAT_JSON,
                    }
                ]
            },
        }
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 3)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertEqual(
            init_step.get_compounds()[0][0].get_smile(),
            "O=C(C)Oc1ccccc1C(=O)O",
        )
        self.assertEqual(init_step.get_compounds()[1].get_name(), "paracetamol")

    def test_initialize_compound_smile(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: "abc:CN(C)CCn1cc(c2ccc(F)cc2)c(n1)n3cccc3;CN(C)CCn1cc(c2ccc(F)cc2)c(n1)n3cccc3",
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STRING,
                    }
                ]
            },
        }
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertEqual(
            init_step.get_compounds()[0][0].get_smile(),
            "CN(C)CCn1cc(c2ccc(F)cc2)c(n1)n3cccc3",
        )
        self.assertEqual(init_step.get_compounds()[0].get_name(), "abc")
        self.assertEqual(init_step.get_compounds()[1].get_name(), "1")
        self.assertEqual(init_step.get_compounds()[1].get_compound_number(), 1)

    def test_initialize_compound_smile_enforceIDs(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: "abc:CN(C)CCn1cc(c2ccc(F)cc2)c(n1)n3cccc3;CN(C)CCn1cc(c2ccc(F)cc2)c(n1)n3cccc3",
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STRING,
                        _SBE.INPUT_ENFORCE_IDS: {
                            _SBE.INPUT_ENFORCE_COMPOUND_IDS: ["3", 1],
                            _SBE.INPUT_ENFORCE_ENUMERATION_IDS: [10, 4],
                        },
                    }
                ],
            },
        }
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertEqual(
            init_step.get_compounds()[0][0].get_smile(),
            "CN(C)CCn1cc(c2ccc(F)cc2)c(n1)n3cccc3",
        )
        self.assertEqual(init_step.get_compounds()[0].get_name(), "abc")
        self.assertEqual(init_step.get_compounds()[1].get_name(), "1")
        self.assertEqual(init_step.get_compounds()[0].get_compound_number(), 3)
        self.assertEqual(init_step.get_compounds()[1][0].get_enumeration_id(), 4)

    def test_initialize_compound_CSV(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: self._CSV_path,
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                        _SBE.INPUT_CSV_COLUMNS: {
                            _SBE.INPUT_CSV_SMILES_COLUMN: "SMILES"
                        },
                        _SBE.INPUT_FORMAT: _SBE.FORMAT_CSV,
                    }
                ]
            },
        }
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)

        self.assertEqual(
            init_step.get_compounds()[0][0].get_smile(), "O=C(C)Oc1ccccc1C(=O)O"
        )
        self.assertEqual(init_step.get_compounds()[1].get_name(), "1")

    def test_initialize_compound_CSV_extended_options(self):
        step_conf = {
            _SBE.STEPID: "01_load_molecule",
            _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
            _SBE.INPUT: {
                _SBE.INPUT_COMPOUNDS: [
                    {
                        _SBE.INPUT_SOURCE: self._CSV_path_semicolon,
                        _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                        _SBE.INPUT_CSV_DELIMITER: ";",
                        _SBE.INPUT_CSV_COLUMNS: {
                            _SBE.INPUT_CSV_SMILES_COLUMN: "SMILES",
                            _SBE.INPUT_CSV_NAMES_COLUMN: "name",
                        },
                        _SBE.INPUT_FORMAT: _SBE.FORMAT_CSV,
                    }
                ]
            },
        }
        init_step = StepInitializeCompound(**step_conf)
        init_step.generate_input()
        init_step.execute()

        self.assertEqual(len(init_step.get_compounds()), 2)
        self.assertEqual(len(init_step.get_compounds()[0]), 1)
        self.assertEqual(len(init_step.get_compounds()[0][0]), 0)
        self.assertEqual(len(init_step.get_compounds()[1]), 1)

        self.assertEqual(
            init_step.get_compounds()[0][0].get_smile(), "O=C(C)Oc1ccccc1C(=O)O"
        )
        self.assertEqual(init_step.get_compounds()[1].get_name(), "Paracetamol")
