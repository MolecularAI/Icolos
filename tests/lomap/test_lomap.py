import unittest
import os
from icolos.core.workflow_steps.calculation.lomap import StepLomap
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import (
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    get_docked_ligands_as_conformers,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.composite_agents.workflow import WorkFlow

_SBE = StepBaseEnum


class Test_Lomap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/lomap")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        # export_unit_test_env_vars()

    def setUp(self):
        # load in some compounds
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS
        )[:12]

    def test_generate_lomap_mcs_map(self):
        step_conf = {
            _SBE.STEPID: "test_lomap_mcs_map",
            _SBE.STEP_TYPE: _SBE.STEP_LOMAP,
            _SBE.EXEC: {
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 8}
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {"topology": "mcs"},
            },
        }
        w_flow = WorkFlow()
        step_lomap = StepLomap(**step_conf)
        step_lomap.set_workflow_object(w_flow)
        step_lomap.data.compounds = self.compounds
        step_lomap.execute()
        p_map = step_lomap.get_perturbation_map()
        self.assertEqual(len(p_map.edges), 16)
        self.assertEqual(len(p_map.nodes), 12)
