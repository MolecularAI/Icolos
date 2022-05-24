import shutil
import unittest
import os
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.pmx.mutate import StepPMXmutate
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path


_SBE = StepBaseEnum
_SGE = StepGromacsEnum()


class Test_PMXmutate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/mutate")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)
        if os.path.exists(cls._test_dir):
            shutil.rmtree(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.PMX_MUTATIONS_PROTEIN, "r") as f:
            data = f.read()
        self.system = GenericData(file_name="protein.pdb", file_data=data)
        with open(PATHS_EXAMPLEDATA.PMX_MUTATIONS_LIST, "r") as f:
            muts = f.read()
        self.muts = GenericData(file_name="mutations.mut", file_data=muts)

    def test_pmx_mutate(self):
        step_conf = {
            _SBE.STEPID: "01_PMX_SETUP",
            _SBE.STEP_TYPE: _SBE.STEP_PMX_MUTATE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    # settings for protein parametrisation
                    "forcefield": "amber14sbmut",
                    "water": "tip3p",
                },
            },
        }

        step_mutate = StepPMXmutate(**step_conf)
        step_mutate.data.generic.add_file(self.muts)
        step_mutate.data.generic.add_file(self.system)

        step_mutate.work_dir = self._test_dir
        step_mutate._workflow_object = WorkFlow()
        step_mutate.execute()
