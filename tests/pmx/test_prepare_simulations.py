import unittest
import os
from icolos.core.workflow_steps.pmx.prepare_simulations import StepPMXPrepareSimulations
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.composite_agents.workflow import WorkFlow
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    MAIN_CONFIG,
    export_unit_test_env_vars,
    create_test_dir,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.containers.perturbation_map import PerturbationMap

_SBE = StepBaseEnum


class Test_PMXPrepareSimulations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/test_prepare_simulations")

        create_test_dir(PATHS_EXAMPLEDATA.PREPARE_SIMULATIONS_TEST_DIR, cls._test_dir)

    def setUp(self):
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.FEP_PLUS_LIGANDS
        )
        p_map = PerturbationMap(compounds=self.compounds)
        p_map.parse_map_file(file_path=PATHS_EXAMPLEDATA.FEP_PLUS_MAP_LOG_SINGLE_EDGE)
        self.p_map = p_map

        export_unit_test_env_vars()

    def test_prepare_simulations(self):

        step_conf = {
            _SBE.STEPID: "prepare_simulations",
            _SBE.STEP_TYPE: "pmx_prepare_simulations",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["PMX"]["CLI_ENTRYPOINT"],
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {"sim_type": "em"},
            },
        }

        step_prepare_simulations = StepPMXPrepareSimulations(**step_conf)
        step_prepare_simulations.work_dir = self._test_dir
        step_prepare_simulations._workflow_object = WorkFlow()
        step_prepare_simulations._workflow_object.workflow_data.perturbation_map = (
            self.p_map
        )
        step_prepare_simulations.execute()

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0cd4b47_4f2ffa1/water/stateA/run1/em/tpr.tpr")
        )

        self.assertGreater(stat_inf.st_size, 213300)

        stat_inf = os.stat(
            os.path.join(
                self._test_dir, "0cd4b47_4f2ffa1/protein/stateB/run3/em/tpr.tpr"
            )
        )
        self.assertGreater(stat_inf.st_size, 3501000)
