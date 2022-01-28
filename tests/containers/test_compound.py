import unittest
import os
from copy import deepcopy
from rdkit import Chem

from icolos.core.containers.compound import Conformer, Enumeration, Compound

from icolos.utils.enums.compound_enums import (
    CompoundContainerEnum,
    EnumerationContainerEnum,
)

from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.general.files_paths import attach_root_path


class Test_Compound(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._CC = CompoundContainerEnum()
        cls._EC = EnumerationContainerEnum()

        cls._test_dir = attach_root_path("tests/junk/Compound")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        comp = Compound(name="test_molecule", compound_number=0)
        enum1 = Enumeration(smile="", molecule=None)
        self.e1_conf1 = Conformer(conformer_id=0)
        self.e1_conf2 = Conformer(conformer_id=2)
        enum1.add_conformer(self.e1_conf1, auto_update=True)
        enum1.add_conformer(self.e1_conf2, auto_update=True)
        enum2 = Enumeration(smile="", molecule=None)
        self.e2_conf1 = Conformer(conformer_id=1)
        self.e2_conf2 = Conformer(conformer_id=3)
        self.e2_conf3 = Conformer(conformer_id=5)
        enum2.add_conformer(self.e2_conf1, auto_update=False)
        enum2.add_conformer(self.e2_conf2, auto_update=False)
        enum2.add_conformer(self.e2_conf3, auto_update=False)
        enum3 = Enumeration(smile="CCC", molecule=None, enumeration_id=4)
        self.e3_conf1 = Conformer(conformer_id=0)
        enum3.add_conformer(self.e3_conf1, auto_update=True)
        comp.add_enumeration(enumeration=enum1, auto_update=True)
        comp.add_enumeration(enumeration=enum2, auto_update=True)
        comp.add_enumeration(enumeration=enum3, auto_update=False)
        self.comp = comp
        self.enum1 = enum1
        self.enum2 = enum2
        self.enum3 = enum3

        _paracetamol_path = attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_PATH)
        mol_supplier = Chem.SDMolSupplier(_paracetamol_path, removeHs=False)
        for mol in mol_supplier:
            self._paracetamol_molecule = mol
        _aspirin_path = attach_root_path(PATHS_EXAMPLEDATA.ASPIRIN_PATH)
        mol_supplier = Chem.SDMolSupplier(_aspirin_path, removeHs=False)
        for mol in mol_supplier:
            self._aspirin_molecule = mol

    @classmethod
    def tearDownClass(cls):
        pass

    def test_general_handling(self):
        # Enumeration
        self.assertEqual(len(self.comp), 3)
        l_enums = self.comp.get_enumerations()
        self.assertEqual(l_enums[0].get_compound_object(), self.comp)
        self.assertEqual(l_enums[1].get_enumeration_id(), 1)
        self.assertIsNone(l_enums[2].get_compound_object())
        self.assertEqual(l_enums[2].get_enumeration_id(), 4)
        self.assertEqual(self.comp[2].get_enumeration_id(), 4)

        self.assertRaises(IndexError, self.comp.find_enumeration, 3)
        self.assertEqual(
            self.comp.find_enumeration(enumeration_id=4).get_smile(), "CCC"
        )

        self.assertListEqual([0, 1, 4], self.comp.get_enumeration_ids())
        self.comp.reset_enumeration_ids()
        self.assertListEqual([0, 1, 2], self.comp.get_enumeration_ids())

        # Conformer
        self.assertEqual(len(self.comp.find_enumeration(1)), 3)
        self.assertEqual(self.comp[1][1].get_conformer_id(), 3)
        self.assertListEqual([0, 1], self.comp[0].get_conformer_ids())
        self.assertListEqual([1, 3, 5], self.comp[1].get_conformer_ids())

        # Deletion
        self.comp[1].clear_conformers()
        self.assertEqual(len(self.comp[1]), 0)
        self.comp.clear_enumerations()
        self.assertEqual(len(self.comp), 0)

    def test_cloning_and_resetting(self):
        comp_clone = deepcopy(self.comp)
        comp_clone[0].set_enumeration_id(10)
        self.assertListEqual([0, 1, 4], self.comp.get_enumeration_ids())
        self.assertListEqual([10, 1, 4], comp_clone.get_enumeration_ids())

        all_conf_ids = []
        for enum in self.comp:
            for conf in enum:
                all_conf_ids.append(conf.get_conformer_id())
        self.assertListEqual([0, 1, 1, 3, 5, 0], all_conf_ids)

        comp_clone.reset_all_ids()
        all_conf_ids = []
        for enum in comp_clone:
            for conf in enum:
                all_conf_ids.append(conf.get_conformer_id())
        self.assertListEqual([0, 1, 0, 1, 2, 0], all_conf_ids)
