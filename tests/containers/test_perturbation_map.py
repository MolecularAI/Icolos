from icolos.core.containers.perturbation_map import PerturbationMap
import unittest
import os
from icolos.core.containers.generic import GenericData
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    construct_full_compound_object,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path


class Test_PerturbationMap(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:

        cls._test_dir = attach_root_path("tests/junk/perturbation_map")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS,
        )
        with open(PATHS_EXAMPLEDATA.PMX_TNKS_PROTEIN, "r") as f:
            data = f.read()
        protein = GenericData(file_name="protein.pdb", file_data=data)
        p_map = PerturbationMap(compounds=compounds, protein=protein)

        p_map.parse_map_file(PATHS_EXAMPLEDATA.PMX_TNKS_MAP)
        self.p_map = p_map

    def test_perturbation_map(self):
        self.assertEqual(len(self.p_map.nodes), 2)
        self.assertEqual(len(self.p_map.edges), 1)
        self.assertEqual(
            self.p_map.nodes[1].get_conformer().get_enumeration_object().get_smile(),
            "[H]OC1([H])c2c(c(Oc3c([H])c(F)c([H])c(C#N)c3[H])c([H])c([H])c2S(=O)(=O)C([H])([H])C([H])([H])[H])C([H])([H])C1(F)F",
        )

    def test_vis_map(self):
        self.p_map.visualise_perturbation_map(self._test_dir)
        filepath = os.path.join(self._test_dir, "vmap.html")
        stat_inf = os.stat(filepath)
        self.assertGreater(stat_inf.st_size, 2600)
