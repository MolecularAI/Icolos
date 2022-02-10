import os
import unittest

from icolos.core.workflow_steps.ccdc.docking import StepGold

from icolos.utils.enums.step_enums import StepBaseEnum, StepGoldEnum

from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_1UYD_ligands_as_Compounds,
    PATHS_1UYD,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SGE = StepGoldEnum()


class Test_Gold_docking(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/ADV")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._1UYD_compounds = get_1UYD_ligands_as_Compounds(
            abs_path=PATHS_EXAMPLEDATA.PARACETAMOL_PATH
        )
        self.receptor_path = PATHS_1UYD.PDBQT_PATH

    def test_config_generation(self):
        step_conf = {
            _SBE.STEPID: "01_Gold",
            _SBE.STEP_TYPE: _SBE.STEP_GOLD_DOCKING,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load ccdc"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.CONFIGURATION: {
                        _SGE.DATA_FILES: {
                            _SGE.LIGAND_DATA_FILE: ["nope"]
                        },
                        _SGE.FLOOD_FILL: {
                            _SGE.CAVITY_FILE: "whatever"
                        },
                        _SGE.PROTEIN_DATA: {
                            _SGE.PROTEIN_DATA_FILE: "string"
                        }
                    }
                },
            },
        }

        gold_step = StepGold(**step_conf)
        gold_step.data.compounds = self._1UYD_compounds

        config_path = os.path.join(self._test_dir, "test_config.json")
        gold_step.generate_config_file(config_path)

    def test_docking(self):

        adv_step.execute()
        self.assertEqual(len(adv_step.get_compounds()), 1)
        self.assertEqual(len(adv_step.get_compounds()[0][0].get_conformers()), 2)
        self.assertListEqual(
            list(
                adv_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.305, 11.464, 24.663],
        )
        self.assertEqual(
            adv_step.get_compounds()[0][0][0]
            .get_molecule()
            .GetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE),
            "-6.0",
        )

        # check SDF write-out
        out_path = os.path.join(self._test_dir, "adv_docked.sdf")
        adv_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 3500)


"""import unittest
import os
import rdkit.Chem as Chem
from tests.tests_paths import MAIN_CONFIG
if "CSDHOME" in MAIN_CONFIG:
    os.environ["CSDHOME"] = MAIN_CONFIG["CSDHOME"]

from dockstream.core.Gold.Gold_docker import Gold, GoldParameters, GoldFitnessFunction, GoldResponseValue
from dockstream.utils.enums.Gold_enums import GoldDockingConfigurationEnum, GoldLigandPreparationEnum

from tests.tests_paths import PATHS_1UYD, PATH_GOLD_EXAMPLES
from dockstream.utils.files_paths import attach_root_path
from dockstream.core.ligand.ligand import Ligand
from dockstream.utils.smiles import to_smiles

from dockstream.utils.files_paths import lines_in_file


class Test_Gold_backend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._CE = GoldDockingConfigurationEnum()
        cls._LP = GoldLigandPreparationEnum()

        # specify absolute paths to the input file(s)
        cls.target_path = attach_root_path(PATH_GOLD_EXAMPLES.TARGETFILE)
        cls.ligands_sdf = attach_root_path(PATHS_1UYD.LIGANDS_WITH_ENUMERATION_SDF)

    def setUp(self):
        list_ligands = []
        for mol in Chem.SDMolSupplier(self.ligands_sdf, removeHs=False):
            name = mol.GetProp("_Name").split(':')
            lig_number = int(name[0])
            enumeration = int(name[1])
            list_ligands.append(Ligand(smile=to_smiles(mol),
                                       original_smile=to_smiles(mol),
                                       ligand_number=lig_number,
                                       enumeration=enumeration,
                                       molecule=mol,
                                       mol_type=self._LP.TYPE_RDKIT))
        self.ligands = list_ligands

    def test_Gold_docking_fitness(self):
        docker = Gold(
            input_pools=["Corina_pool"],
            parameters=GoldParameters(
                prefix_execution="module load ccdc/2020.3.0",
                receptor_paths=[self.target_path],
                fitness_function=GoldFitnessFunction.PLP,
                diverse_solutions=(False, None, None),
                early_termination=True,
                autoscale=100.0,
                ndocks=3,
            )
        )
        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(63, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertTrue(all(float(x) > 50 for x in list(result.iloc[::3, :]["score"])))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # test score retrieval
        self.assertEqual(17, len(docker.get_scores(best_only=True)))
        all_scores = docker.get_scores(best_only=False)
        self.assertEqual(63, len(all_scores))
        self.assertTrue(all_scores[0] > all_scores[1] > all_scores[2])
        self.assertTrue(all_scores[3] > all_scores[4] > all_scores[5])
        self.assertTrue(all_scores[6] > all_scores[7] > all_scores[8])

    def test_Gold_docking_value(self):
        docker = Gold(
            input_pools=["Corina_pool"],
            parameters=GoldParameters(
                prefix_execution="module load ccdc/2020.3.0",
                receptor_paths=[self.target_path],
                fitness_function=GoldFitnessFunction.PLP,
                diverse_solutions=(False, None, None),
                early_termination=True,
                autoscale=100.0,
                ndocks=3,
                response_value=GoldResponseValue.VALUE
            )
        )
        docker.add_molecules(molecules=self.ligands)
        docker.dock()
        result = docker.get_result()

        # test dataframe output
        self.assertEqual(63, result.shape[0])
        self.assertEqual(7, result.shape[1])
        self.assertTrue(all(float(x) < -40 for x in list(result.iloc[::3, :]["score"])))

        # test ligand retrieval
        list_docked_ligands = docker.get_docked_ligands()
        del list_docked_ligands[0]
        list_docked_ligands = docker.get_docked_ligands()
        self.assertEqual(21, len(list_docked_ligands))

        # write out poses and check length
        out_path = attach_root_path("tests/junk/Gold_backend_docked_all.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertGreater(lines_in_file(out_path), 13000)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_enumeration.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertGreater(lines_in_file(out_path), 4000)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_ligand.sdf")
        docker.write_docked_ligands(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertGreater(lines_in_file(out_path), 3000)

        # write out dataframe and check length
        out_path = attach_root_path("tests/junk/Gold_backend_docked_all.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_ALL)
        self.assertEqual(lines_in_file(out_path), 64)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_enumeration.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERENUMERATION)
        self.assertEqual(lines_in_file(out_path), 22)
        out_path = attach_root_path("tests/junk/Gold_backend_docked_best_per_ligand.csv")
        docker.write_result(path=out_path, mode=self._CE.OUTPUT_MODE_BESTPERLIGAND)
        self.assertEqual(lines_in_file(out_path), 18)"""
