import unittest
import os
from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.pmx.abfe import StepPMXabfe
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.composite_agents.workflow import WorkFlow
import shutil

_SBE = StepBaseEnum


class Test_PMXabfe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/abfe")
        if os.path.isdir(cls._test_dir):
            shutil.rmtree(cls._test_dir)
        os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.PMX_ABFE_INPUT_COMPLEX, "r") as f:
            data = f.read()
        self.protein = GenericData(file_name="complex.pdb", file_data=data)

    def test_pmx_abfe(self):
        step_conf = {
            _SBE.STEPID: "01_PMX_ABFE",
            _SBE.STEP_TYPE: _SBE.STEP_PMX_ABFE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_FLAGS: ["--build"]},
                _SBE.SETTINGS_ADDITIONAL: {
                    # settings for protein parametrisation
                    "forcefield": "amber03",
                    "water": "tip3p",
                },
            },
        }

        step_pmx_abfe = StepPMXabfe(**step_conf)
        step_pmx_abfe.data.generic.add_file(self.protein)

        step_pmx_abfe.work_dir = self._test_dir
        step_pmx_abfe._workflow_object = WorkFlow()
        step_pmx_abfe.execute()

        self.assertEqual(
            os.path.isfile(os.path.join(self._test_dir, "complex/genion.tpr")), True
        )

        stat_inf = os.stat(os.path.join(self._test_dir, "protein.gro"))
        self.assertGreater(stat_inf.st_size, 70000)
