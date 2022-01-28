import unittest
import os
from icolos.core.workflow_steps.pmx.box_water_ions import StepPMXBoxWaterIons
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    create_test_dir,
    export_unit_test_env_vars,
    get_ligands_as_compounds_with_conformers,
    MAIN_CONFIG,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.containers.perturbation_map import PerturbationMap

_SBE = StepBaseEnum


class Test_PMXBoxWaterIons(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/test_box_water_ions")

        create_test_dir(PATHS_EXAMPLEDATA.BOX_WATER_IONS_TEST_DIR, cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        # initialise the map object for the two test ligands
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS
        )
        p_map = PerturbationMap(compounds=self.compounds)
        p_map.parse_map_file(file_path=PATHS_EXAMPLEDATA.PMX_TNKS_MAP)
        self.p_map = p_map

    # def tearDown(self):
    #     shutil.rmtree(self._test_dir)

    def test_box_water_ions(self):
        conf = {
            _SBE.STEPID: "01_PMX_BOX_WATER_IONS",
            _SBE.STEP_TYPE: _SBE.STEP_PMX_BOX_WATER_IONS,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }
        step = StepPMXBoxWaterIons(**conf)
        step.data.compounds = self.compounds
        step.work_dir = self._test_dir
        step._workflow_object = WorkFlow()
        step._workflow_object.workflow_data.perturbation_map = self.p_map
        step.execute()
        stat_inf = os.stat(
            os.path.join(self._test_dir, "0ec09ef_4afa8f9/water/tpr.tpr")
        )
        self.assertGreater(stat_inf.st_size, 167000)

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0ec09ef_4afa8f9/protein/tpr.tpr")
        )
        self.assertGreater(stat_inf.st_size, 1317800)
