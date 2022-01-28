from icolos.core.containers.perturbation_map import PerturbationMap
import unittest
import os
from icolos.core.containers.generic import GenericData
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    construct_full_compound_object,
)
from icolos.utils.general.files_paths import attach_root_path


class Test_PerturbationMap(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:

        cls._test_dir = attach_root_path("tests/junk/perturbation_map")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        compounds = construct_full_compound_object(
            PATHS_EXAMPLEDATA.FEP_PLUS_LIGANDS,
        )
        with open(PATHS_EXAMPLEDATA.FEP_PLUS_PROTEIN, "r") as f:
            data = f.read()
        protein = GenericData(file_name="protein.pdb", file_data=data)
        p_map = PerturbationMap(compounds=compounds, protein=protein)

        p_map.parse_map_file(PATHS_EXAMPLEDATA.FEP_PLUS_MAP_LOG)
        self.p_map = p_map

    def test_perturbation_map(self):
        self.assertEqual(len(self.p_map.nodes), 38)
        self.assertEqual(len(self.p_map.edges), 62)
        self.assertEqual(
            self.p_map.nodes[5].get_conformer().get_enumeration_object().get_smile(),
            "[H]c1nc(N([H])c2c([H])c(C(=O)N([H])[H])c([H])c(N([H])S(=O)(=O)C([H])([H])[H])c2[H])nc(N([H])c2c(Cl)c([H])c([H])c3c2OC([H])([H])O3)c1[H]",
        )

    def test_vis_map(self):
        self.p_map.visualise_perturbation_map(self._test_dir)
        filepath = os.path.join(self._test_dir, "vmap.html")
        stat_inf = os.stat(filepath)
        self.assertGreater(stat_inf.st_size, 13300)
