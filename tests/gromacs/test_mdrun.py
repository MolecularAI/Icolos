from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.gromacs.mdrun import StepGMXMDrun


_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_MDrun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_TPR_FILE, "rb") as f:
            self.tpr = f.read()

    def test_mdrun(self):
        step_conf = {
            _SBE.STEPID: "test_mdrun",
            _SBE.STEP_TYPE: "mdrun",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
            },
        }

        step_mdrun = StepGMXMDrun(**step_conf)
        step_mdrun.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.tpr, argument=True)
        )
        step_mdrun.execute()

        out_path = os.path.join(self._test_dir, "structure.gro")
        step_mdrun.write_generic_by_extension(self._test_dir, "gro")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 874941)

    def test_mdrun_slurm(self):
        step_conf = {
            _SBE.STEPID: "test_mdrun",
            _SBE.STEP_TYPE: "mdrun",
            _SBE.EXEC: {
                _SBE.EXEC_RESOURCE: "slurm",
                _SBE.EXEC_JOB_CONTROL: {
                    _SBE.EXEC_JOB_CONTROL_PARTITION: "gpu",
                    _SBE.EXEC_JOB_CONTROL_GRES: "gpu:1",
                    _SBE.EXEC_JOB_CONTROL_MODULES: ["GROMACS/2020.3-fosscuda-2019a"],
                },
            },
        }

        step_mdrun = StepGMXMDrun(**step_conf)
        step_mdrun.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.tpr, argument=True)
        )
        step_mdrun.execute()

        out_path = os.path.join(self._test_dir, "structure.gro")
        step_mdrun.write_generic_by_extension(self._test_dir, "gro")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 874941)
