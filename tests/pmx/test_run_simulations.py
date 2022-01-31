import unittest
import os
from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.pmx.run_simulations import StepPMXRunSimulations
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    create_test_dir,
    MAIN_CONFIG,
    export_unit_test_env_vars,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.core.composite_agents.workflow import WorkFlow

_SBE = StepBaseEnum


class Test_PMXRunSimulations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/run_simulations_test_dir")
        create_test_dir(PATHS_EXAMPLEDATA.RUN_SIMULATIONS_TEST_DIR, cls._test_dir)

    def setUp(self):
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.FEP_PLUS_LIGANDS
        )
        with open(PATHS_EXAMPLEDATA.FEP_PLUS_OTHER_PROTEIN, "r") as f:
            data = f.read()
        self.protein = GenericData(file_name="protein.pdb", file_data=data)
        p_map = PerturbationMap(compounds=self.compounds, protein=self.protein)
        p_map.parse_map_file(file_path=PATHS_EXAMPLEDATA.FEP_PLUS_MAP_LOG_MIN)
        self.p_map = p_map
        export_unit_test_env_vars()

    def test_run_simulations(self):
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

        step_run_simulations = StepPMXRunSimulations(**step_conf)
        step_run_simulations.work_dir = self._test_dir
        step_run_simulations._workflow_object = WorkFlow()
        step_run_simulations.get_workflow_object().workflow_data.perturbation_map = (
            self.p_map
        )
        step_run_simulations.execute()

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0cd4b47_4f2ffa1/water/stateB/run3/em/md.log")
        )

        self.assertEqual(stat_inf.st_size, 1347767)

        stat_inf = os.stat(
            os.path.join(
                self._test_dir, "0cd4b47_4f2ffa1/protein/stateB/run3/em/tpr.tpr"
            )
        )
        self.assertEqual(stat_inf.st_size, 3501084)

    def test_run_simulations_parallel(self):
        step_conf = {
            _SBE.STEPID: "prepare_simulations",
            _SBE.STEP_TYPE: "pmx_prepare_simulations",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["PMX"]["CLI_ENTRYPOINT"],
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 2},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {"sim_type": "em"},
            },
        }

        step_run_simulations = StepPMXRunSimulations(**step_conf)
        step_run_simulations.work_dir = self._test_dir
        step_run_simulations._workflow_object = WorkFlow()
        step_run_simulations.get_workflow_object().perturbation_map = self.p_map
        step_run_simulations.execute()

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0cd4b47_4f2ffa1/water/stateB/run3/em/md.log")
        )

        self.assertEqual(stat_inf.st_size, 1347767)

    stat_inf = os.stat(
        os.path.join(self._test_dir, "0cd4b47_4f2ffa1/protein/stateB/run3/em/tpr.tpr")
    )
    self.assertEqual(stat_inf.st_size, 3501084)
