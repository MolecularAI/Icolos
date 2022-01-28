import unittest
import os
from icolos.core.workflow_steps.pmx.atomMapping import StepPMXatomMapping
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    create_test_dir,
    export_unit_test_env_vars,
    get_ligands_as_compounds_with_conformers,
)
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.general.files_paths import attach_root_path
from icolos.utils.enums.program_parameters import PMXEnum, PMXAtomMappingEnum

_SBE = StepBaseEnum
_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class Test_PMXatomMapping(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/test_atomMapping")
        create_test_dir(PATHS_EXAMPLEDATA.ATOM_MAPPING_TEST_DIR, cls._test_dir)

    def setUp(self):
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.FEP_PLUS_LIGANDS
        )
        p_map = PerturbationMap(compounds=self.compounds)
        p_map.parse_map_file(file_path=PATHS_EXAMPLEDATA.FEP_PLUS_MAP_LOG_SINGLE_EDGE)
        self.p_map = p_map

        export_unit_test_env_vars()

    def test_atomMapping(self):

        step_conf = {
            _SBE.STEPID: "atommapping",
            _SBE.STEP_TYPE: "pmx_atommapping",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }

        step_atom_mapping = StepPMXatomMapping(**step_conf)
        step_atom_mapping.work_dir = self._test_dir
        step_atom_mapping._workflow_object = WorkFlow()
        step_atom_mapping._workflow_object.workflow_data.perturbation_map = self.p_map
        step_atom_mapping.execute()

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0cd4b47_4f2ffa1/hybridStrTop/out_pdb1.pdb")
        )
        self.assertEqual(stat_inf.st_size, 4631)
