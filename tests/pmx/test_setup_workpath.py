import unittest
import os
from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.pmx.setup_workpath import StepPMXSetup
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path
import shutil
from icolos.core.composite_agents.workflow import WorkFlow

_SBE = StepBaseEnum
_SGE = StepGromacsEnum()


class Test_PMX_setup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/test_setupWorkpath")
        if os.path.exists(cls._test_dir):
            shutil.rmtree(cls._test_dir)
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        self.compounds = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS
        )
        with open(PATHS_EXAMPLEDATA.PMX_TNKS_PROTEIN, "r") as f:
            data = f.read()
        self.protein = GenericData(file_name="protein.pdb", file_data=data)
        with open(PATHS_EXAMPLEDATA.PMX_TNKS_MAP, "r") as f:
            data = f.read()
        self.log_file = GenericData(
            file_name="map.log", file_data=data, extension="log"
        )

    def test_setup_workpath(self):
        step_conf = {
            _SBE.STEPID: "01_PMX_SETUP",
            _SBE.STEP_TYPE: _SBE.STEP_PMX_SETUP,
            _SBE.EXEC: {
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                }
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    # settings for protein parametrisation
                    "forcefield": "amber03",
                    "water": "tip3p",
                    _SGE.CHARGE_METHOD: "bcc",
                },
            },
        }

        step_setup = StepPMXSetup(**step_conf)
        step_setup.data.compounds = self.compounds
        step_setup.data.generic.add_file(self.protein)
        step_setup.data.generic.add_file(self.log_file)
        step_setup.data.generic.add_file(
            GenericData(
                file_name="mdp_files",
                extension="mdp",
                file_data=PATHS_EXAMPLEDATA.PMX_MDP_FILES,
            )
        )
        step_setup.work_dir = self._test_dir
        step_setup._workflow_object = WorkFlow()
        step_setup.execute()

        assert os.path.isdir(os.path.join(self._test_dir, "input"))
        assert os.path.isdir(
            os.path.join(self._test_dir, "0ec09ef_4afa8f9/ligand/stateA/run1/em")
        )
        # stat some of the ligand files and check they've been deposited in the right directory
