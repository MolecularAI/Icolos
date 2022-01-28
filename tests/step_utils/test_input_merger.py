import unittest

from icolos.core.step_utils.input_merger import InputMerger, StepMerge
from icolos.core.containers.compound import Conformer, Enumeration, Compound

from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class Test_InputMerger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        # comp1 has 2 enumerations, one with 2 and one with 3 conformers
        comp1 = Compound(name="test_molecule", compound_number=0)
        comp1_enum1 = Enumeration(smile="abc", molecule=None, enumeration_id=1)
        comp1_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum2 = Enumeration(smile="def", molecule=None, enumeration_id=2)
        comp1_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1.add_enumeration(comp1_enum1, auto_update=False)
        comp1.add_enumeration(comp1_enum2, auto_update=False)

        # comp2 has 3 enumerations, one with 1, one with 3 and one with 4 conformers
        comp2 = Compound(name="test_molecule_new", compound_number=0)
        comp2_enum1 = Enumeration(smile="kk", molecule=None, enumeration_id=0)
        comp2_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum2 = Enumeration(smile="abc", molecule=None, enumeration_id=1)
        comp2_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum3 = Enumeration(smile="xyz", molecule=None, enumeration_id=2)
        comp2_enum3.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum3.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum3.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2_enum3.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp2.add_enumeration(comp2_enum1, auto_update=False)
        comp2.add_enumeration(comp2_enum2, auto_update=False)
        comp2.add_enumeration(comp2_enum3, auto_update=False)

        # comp3 has 1 enumeration, with 2 conformers (and a different number and name)
        comp3 = Compound(name="test_molecule", compound_number=1)
        comp3_enum1 = Enumeration(smile="abc", molecule=None, enumeration_id=0)
        comp3_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp3_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp3_enum2 = Enumeration(smile="def", molecule=None, enumeration_id=1)
        comp3_enum2.add_conformer(Conformer(conformer_id=0), auto_update=False)
        comp3_enum2.add_conformer(Conformer(conformer_id=0), auto_update=False)
        comp3_enum2.add_conformer(Conformer(conformer_id=0), auto_update=False)
        comp3.add_enumeration(comp3_enum1, auto_update=False)
        comp3.add_enumeration(comp3_enum2, auto_update=False)
        self.list_compounds = [comp1, comp2, comp3]

    @classmethod
    def tearDownClass(cls):
        pass

    def test_merging_by_name_compound(self):
        conf = {
            _SBE.INPUT_MERGE_COMPOUNDS: True,
            _SBE.INPUT_MERGE_COMPOUNDS_BY: _SBE.INPUT_MERGE_BY_NAME,
            _SBE.INPUT_MERGE_ENUMERATIONS: False,
        }
        conf = StepMerge(**conf)
        merger = InputMerger(conf)
        list_compounds = merger.merge(self.list_compounds)

        self.assertEqual(len(list_compounds), 2)
        self.assertEqual(len(list_compounds[0].get_enumerations()), 4)
        self.assertEqual(len(list_compounds[1].get_enumerations()), 3)

        self.assertListEqual(
            [c.get_name() for c in list_compounds],
            ["test_molecule", "test_molecule_new"],
        )
        self.assertListEqual(
            [
                conf.get_index_string()
                for c in list_compounds
                for e in c.get_enumerations()
                for conf in e.get_conformers()
            ],
            [
                "0:0:0",
                "0:0:1",
                "0:1:0",
                "0:1:1",
                "0:1:2",
                "0:2:0",
                "0:2:1",
                "0:3:0",
                "0:3:1",
                "0:3:2",
                "1:0:0",
                "1:0:1",
                "1:1:0",
                "1:1:1",
                "1:1:2",
                "1:2:0",
                "1:2:1",
                "1:2:2",
                "1:2:3",
            ],
        )

    def test_merging_by_id_compound(self):
        conf = {
            _SBE.INPUT_MERGE_COMPOUNDS: True,
            _SBE.INPUT_MERGE_COMPOUNDS_BY: _SBE.INPUT_MERGE_BY_ID,
            _SBE.INPUT_MERGE_ENUMERATIONS: False,
        }
        conf = StepMerge(**conf)
        merger = InputMerger(conf)
        list_compounds = merger.merge(self.list_compounds)

        self.assertEqual(len(list_compounds), 2)
        self.assertEqual(len(list_compounds[0].get_enumerations()), 5)
        self.assertEqual(len(list_compounds[1].get_enumerations()), 2)

        self.assertListEqual([c.get_name() for c in list_compounds], ["0", "1"])

        self.assertListEqual(
            [
                conf.get_index_string()
                for c in list_compounds
                for e in c.get_enumerations()
                for conf in e.get_conformers()
            ],
            [
                "0:0:0",
                "0:0:1",
                "0:1:0",
                "0:1:1",
                "0:1:2",
                "0:2:0",
                "0:2:1",
                "0:3:0",
                "0:3:1",
                "0:3:2",
                "0:4:0",
                "0:4:1",
                "0:4:2",
                "0:4:3",
                "1:0:0",
                "1:0:1",
                "1:1:0",
                "1:1:1",
                "1:1:2",
            ],
        )

    def test_merging_by_name_compound_enumeration_smile(self):
        conf = {
            _SBE.INPUT_MERGE_COMPOUNDS: True,
            _SBE.INPUT_MERGE_COMPOUNDS_BY: _SBE.INPUT_MERGE_BY_NAME,
            _SBE.INPUT_MERGE_ENUMERATIONS: True,
            _SBE.INPUT_MERGE_ENUMERATIONS_BY: _SBE.INPUT_MERGE_BY_SMILE,
        }
        conf = StepMerge(**conf)
        merger = InputMerger(conf)
        list_compounds = merger.merge(self.list_compounds)

        self.assertEqual(len(list_compounds), 2)
        self.assertEqual(len(list_compounds[0].get_enumerations()), 2)
        self.assertEqual(len(list_compounds[1].get_enumerations()), 3)

        self.assertListEqual(
            [c.get_name() for c in list_compounds],
            ["test_molecule", "test_molecule_new"],
        )
        self.assertListEqual(
            [
                conf.get_index_string()
                for c in list_compounds
                for e in c.get_enumerations()
                for conf in e.get_conformers()
            ],
            [
                "0:0:0",
                "0:0:1",
                "0:0:2",
                "0:0:3",
                "0:1:0",
                "0:1:1",
                "0:1:2",
                "0:1:3",
                "0:1:4",
                "0:1:5",
                "1:0:0",
                "1:0:1",
                "1:1:0",
                "1:1:1",
                "1:1:2",
                "1:2:0",
                "1:2:1",
                "1:2:2",
                "1:2:3",
            ],
        )
        self.assertListEqual(
            [e.get_smile() for c in list_compounds for e in c.get_enumerations()],
            ["abc", "def", "kk", "abc", "xyz"],
        )

    def test_merging_by_name_compound_enumeration_id(self):
        conf = {
            _SBE.INPUT_MERGE_COMPOUNDS: True,
            _SBE.INPUT_MERGE_COMPOUNDS_BY: _SBE.INPUT_MERGE_BY_NAME,
            _SBE.INPUT_MERGE_ENUMERATIONS: True,
            _SBE.INPUT_MERGE_ENUMERATIONS_BY: _SBE.INPUT_MERGE_BY_ID,
        }
        conf = StepMerge(**conf)
        merger = InputMerger(conf)
        list_compounds = merger.merge(self.list_compounds)

        self.assertEqual(len(list_compounds), 2)
        self.assertEqual(len(list_compounds[0].get_enumerations()), 3)
        self.assertEqual(len(list_compounds[1].get_enumerations()), 3)

        self.assertListEqual(
            [c.get_name() for c in list_compounds],
            ["test_molecule", "test_molecule_new"],
        )
        self.assertListEqual(
            [
                conf.get_index_string()
                for c in list_compounds
                for e in c.get_enumerations()
                for conf in e.get_conformers()
            ],
            [
                "0:0:0",
                "0:0:1",
                "0:0:2",
                "0:0:3",
                "0:0:4",
                "0:1:0",
                "0:1:1",
                "0:1:2",
                "0:2:0",
                "0:2:1",
                "1:0:0",
                "1:0:1",
                "1:1:0",
                "1:1:1",
                "1:1:2",
                "1:2:0",
                "1:2:1",
                "1:2:2",
                "1:2:3",
            ],
        )
        self.assertListEqual(
            [e.get_smile() for c in list_compounds for e in c.get_enumerations()],
            ["abc", "def", "abc", "kk", "abc", "xyz"],
        )
