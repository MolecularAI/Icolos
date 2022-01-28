import unittest
import os
from icolos.core.workflow_steps.pmx.ligandHybrid import StepPMXligandHybrid
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.utils.enums.program_parameters import PMXEnum, PMXLigandHybridEnum
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    create_test_dir,
    get_ligands_as_compounds_with_conformers,
    export_unit_test_env_vars,
)
from icolos.utils.general.files_paths import attach_root_path


_SBE = StepBaseEnum
_PE = PMXEnum()
_PHE = PMXLigandHybridEnum()


class Test_PMXligandHybrid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/test_ligandHybrid")
        # if not os.path.isdir(cls._test_dir):
        #     os.makedirs(cls._test_dir)
        create_test_dir(PATHS_EXAMPLEDATA.LIGAND_HYBRID_TEST_DIR, cls._test_dir)

    def setUp(self):
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS
        )
        p_map = PerturbationMap(compounds=self.compounds)
        p_map.parse_map_file(file_path=PATHS_EXAMPLEDATA.PMX_TNKS_MAP)
        self.p_map = p_map

        export_unit_test_env_vars()

    # def tearDown(self):
    #     shutil.rmtree(self._test_dir)

    def test_build_hybrid_topology_and_structure(self):
        merged_itp_path = os.path.join(
            self._test_dir, "0ec09ef_4afa8f9/hybridStrTop/merged.itp"
        )

        step_conf = {
            _SBE.STEPID: "ligand_hybrid",
            _SBE.STEP_TYPE: "pmx_ligandHybrid",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }

        step_ligand_hybrid = StepPMXligandHybrid(**step_conf)
        step_ligand_hybrid.work_dir = self._test_dir
        step_ligand_hybrid._workflow_object = WorkFlow()
        step_ligand_hybrid._workflow_object.workflow_data.perturbation_map = self.p_map
        step_ligand_hybrid.execute()

        stat_inf = os.stat(merged_itp_path)
        self.assertGreater(stat_inf.st_size, 39400)
