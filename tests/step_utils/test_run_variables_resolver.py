import unittest

from icolos.core.containers.compound import Conformer, Enumeration, Compound
from icolos.core.step_utils.run_variables_resolver import RunVariablesResolver
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class Test_RunVariablesResolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.resolver = RunVariablesResolver()

    def setUp(self):
        # comp1 has 2 enumerations, one with 2 and one with 3 conformers
        comp1 = Compound(name="test_molecule", compound_number=0)
        comp1_enum1 = Enumeration(
            smile="abc", molecule=None, enumeration_id=1, compound_object=comp1
        )
        comp1_enum1.add_conformer(
            Conformer(conformer_id=0, enumeration_object=comp1_enum1), auto_update=True
        )
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
        comp2_enum3.add_conformer(
            Conformer(conformer_id=0, enumeration_object=comp2_enum3), auto_update=True
        )
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
        comp3_enum2 = Enumeration(
            smile="def", molecule=None, enumeration_id=1, compound_object=comp3
        )
        comp3_enum2.add_conformer(Conformer(conformer_id=0), auto_update=False)
        comp3_enum2.add_conformer(Conformer(conformer_id=0), auto_update=False)
        comp3_enum2.add_conformer(Conformer(conformer_id=0), auto_update=False)
        comp3.add_enumeration(comp3_enum1, auto_update=False)
        comp3.add_enumeration(comp3_enum2, auto_update=False)
        self.list_compounds = [comp1, comp2, comp3]

    @classmethod
    def tearDownClass(cls):
        pass

    def test_compound_replacements(self):
        inp = "/a/path/to/nowhere/[compound_id]/[compound_id]/compound_id/whatever/[compound_name]"
        self.assertEqual(
            self.resolver.resolve_compound_level(inp, self.list_compounds[0]),
            "/a/path/to/nowhere/0/0/compound_id/whatever/test_molecule",
        )
        self.assertEqual(
            self.resolver.resolve_compound_level(inp, self.list_compounds[1]),
            "/a/path/to/nowhere/0/0/compound_id/whatever/test_molecule_new",
        )
        self.assertEqual(
            self.resolver.resolve_compound_level(inp, self.list_compounds[2]),
            "/a/path/to/nowhere/1/1/compound_id/whatever/test_molecule",
        )

        # test what happens, when no replacement is done
        inp = "/a/string/withouttreplacement"
        self.assertEqual(
            self.resolver.resolve_compound_level(inp, self.list_compounds[0]), inp
        )

    def test_enumeration_replacements(self):
        inp = "/a/path/to/nowhere/[compound_id]/[enumeration_id]/[enumeration_string]/whatever/[enumeration_id]"
        self.assertEqual(
            self.resolver.resolve_enumeration_level(inp, self.list_compounds[0][0]),
            "/a/path/to/nowhere/[compound_id]/1/0:1/whatever/1",
        )
        self.assertEqual(
            self.resolver.resolve_enumeration_level(inp, self.list_compounds[0][1]),
            "/a/path/to/nowhere/[compound_id]/2/:2/whatever/2",
        )
        self.assertEqual(
            self.resolver.resolve_enumeration_level(inp, self.list_compounds[2][1]),
            "/a/path/to/nowhere/[compound_id]/1/1:1/whatever/1",
        )

        # test what happens, when no replacement is done
        inp = "/a/string/withouttreplacement"
        self.assertEqual(
            self.resolver.resolve_enumeration_level(inp, self.list_compounds[0][0]), inp
        )

    def test_conformer_replacements(self):
        inp = "/a/path/[conformer_string]to/nowhere/[compound_id]/[conformer_id]/[enumeration_string]/whatever/[conformer_id]"
        self.assertEqual(
            self.resolver.resolve_conformer_level(inp, self.list_compounds[0][0][0]),
            "/a/path/0:1:0to/nowhere/[compound_id]/0/[enumeration_string]/whatever/0",
        )
        self.assertEqual(
            self.resolver.resolve_conformer_level(inp, self.list_compounds[0][0][1]),
            "/a/path/0:1:1to/nowhere/[compound_id]/1/[enumeration_string]/whatever/1",
        )
        self.assertEqual(
            self.resolver.resolve_conformer_level(inp, self.list_compounds[2][0][1]),
            "/a/path/:0:1to/nowhere/[compound_id]/1/[enumeration_string]/whatever/1",
        )
        self.assertEqual(
            self.resolver.resolve_conformer_level(inp, self.list_compounds[1][2][0]),
            "/a/path/:2:0to/nowhere/[compound_id]/0/[enumeration_string]/whatever/0",
        )

        # test what happens, when no replacement is done
        inp = "/a/string/withouttreplacement"
        self.assertEqual(
            self.resolver.resolve_conformer_level(inp, self.list_compounds[0][0][0]),
            inp,
        )

    def test_resolve(self):
        inp = "/a/path/[conformer_string]to/nowhere/[compound_id]/[conformer_id]/[enumeration_string]/whatever/[compound_name]"
        self.assertEqual(
            self.resolver.resolve(inp, self.list_compounds[0][0][0]),
            "/a/path/0:1:0to/nowhere/0/0/0:1/whatever/test_molecule",
        )
        self.assertEqual(
            self.resolver.resolve(inp, self.list_compounds[0][0]),
            "/a/path/[conformer_string]to/nowhere/0/[conformer_id]/0:1/whatever/test_molecule",
        )
        self.assertEqual(
            self.resolver.resolve(inp, self.list_compounds[0]),
            "/a/path/[conformer_string]to/nowhere/0/[conformer_id]/[enumeration_string]/whatever/test_molecule",
        )

        # fails for cases where the linking conformer -> enumeration -> compound is not established
        try:
            self.resolver.resolve(inp, self.list_compounds[2][0][1])
        except Exception as e:
            self.assertEqual(
                e.__str__(), "'NoneType' object has no attribute 'get_compound_number'"
            )

        # test what happens, when no replacement is done
        inp = "/a/string/withouttreplacement"
        self.assertEqual(self.resolver.resolve(inp, self.list_compounds[0][0][0]), inp)
