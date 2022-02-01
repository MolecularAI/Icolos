import unittest
import os
from icolos.core.workflow_steps.pmx.assemble_systems import StepPMXAssembleSystems
from icolos.core.containers.generic import GenericData
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    MAIN_CONFIG,
    create_test_dir,
    export_unit_test_env_vars,
    get_ligands_as_compounds_with_conformers,
)
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.utils.general.files_paths import attach_root_path
from icolos.utils.enums.program_parameters import PMXEnum, PMXAtomMappingEnum

_SBE = StepBaseEnum
_PE = PMXEnum()
_PAE = PMXAtomMappingEnum()


class Test_PMXAssembleSystems(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/test_assemble_systems")

        create_test_dir(PATHS_EXAMPLEDATA.ASSEMBLE_SYSTEMS_TEST_DIR, cls._test_dir)

    def setUp(self):
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS
        )
        p_map = PerturbationMap(compounds=self.compounds)
        p_map.parse_map_file(file_path=PATHS_EXAMPLEDATA.PMX_TNKS_MAP)
        self.p_map = p_map
        with open(PATHS_EXAMPLEDATA.PMX_TNKS_PROTEIN, "r") as f:
            data = f.read()
        self.protein = GenericData(file_name="protein.pdb", file_data=data)

        export_unit_test_env_vars()

    def test_assembleSystems(self):

        step_conf = {
            _SBE.STEPID: "assemble_systems",
            _SBE.STEP_TYPE: "pmx_assemble_systems",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["PMX"]["CLI_ENTRYPOINT"],
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

        step_assembleSystems = StepPMXAssembleSystems(**step_conf)
        step_assembleSystems.work_dir = self._test_dir
        step_assembleSystems._workflow_object = WorkFlow()
        step_assembleSystems._workflow_object.workflow_data.perturbation_map = (
            self.p_map
        )
        step_assembleSystems.data.generic.add_file(self.protein)
        step_assembleSystems.execute()

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0ec09ef_4afa8f9/complex/init.pdb")
        )
        self.assertGreater(stat_inf.st_size, 141300)

        stat_inf = os.stat(
            os.path.join(self._test_dir, "0ec09ef_4afa8f9/ligand/init.pdb")
        )
        self.assertGreater(stat_inf.st_size, 2800)
