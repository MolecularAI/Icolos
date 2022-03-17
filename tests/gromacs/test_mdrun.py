from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.containers.gmx_state import GromacsState
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
        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR, "rb") as f:
            self.tpr = f.read()

    # def test_mdrun_external_tpr(self):
    #     step_conf = {
    #         _SBE.STEPID: "test_mdrun",
    #         _SBE.STEP_TYPE: "mdrun",
    #         _SBE.EXEC: {
    #             _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
    #         },
    #         _SBE.SETTINGS: {
    #             _SBE.SETTINGS_ARGUMENTS: {
    #                 _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-nsteps": 1000}
    #             }
    #         },
    #     }

    #     step_mdrun = StepGMXMDrun(**step_conf)
    #     topol = GromacsTopol()
    #     step_mdrun = StepGMXMDrun(**step_conf)
    #     wf = WorkFlow()
    #     wf.workflow_data.gmx_topol = topol
    #     step_mdrun.set_workflow_object(wf)
    #     step_mdrun.data.generic.add_file(
    #         GenericData(file_name=_SGE.STD_TPR, file_data=self.tpr, argument=True)
    #     )
    #     step_mdrun.execute()

    #     out_path = os.path.join(self._test_dir, "confout.gro")
    #     step_mdrun.write_generic_by_extension(self._test_dir, "gro")
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 3224400)

    def test_mdrun_internal_tpr(self):
        step_conf = {
            _SBE.STEPID: "test_mdrun",
            _SBE.STEP_TYPE: "mdrun",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-nsteps": 10000}
                }
            },
        }

        step_mdrun = StepGMXMDrun(**step_conf)
        topol = GromacsState()
        topol.tprs = [
            GenericData(file_name=_SGE.STD_TPR, file_data=self.tpr, argument=True)
        ]
        step_mdrun = StepGMXMDrun(**step_conf)
        wf = WorkFlow()
        wf.workflow_data.gmx_state = topol

        step_mdrun.set_workflow_object(wf)
        step_mdrun.execute()

        out_path = os.path.join(self._test_dir, "confout.gro")
        step_mdrun.write_generic_by_extension(self._test_dir, "gro")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 3224400)

    # def test_mdrun_slurm(self):
    #     step_conf = {
    #         _SBE.STEPID: "test_mdrun",
    #         _SBE.STEP_TYPE: "mdrun",
    #         _SBE.EXEC: {
    #             _SBE.EXEC_PLATFORM: "slurm",
    #             _SBE.EXEC_RESOURCES: {
    #                 _SBE.EXEC_RESOURCES_PARTITION: "gpu",
    #                 _SBE.EXEC_RESOURCES_GRES: "gpu:1",
    #                 _SBE.EXEC_RESOURCES_MODULES: [
    #                     "GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
    #                 ],
    #             },
    #         },
    #         _SBE.SETTINGS: {
    #             _SBE.SETTINGS_ARGUMENTS: {
    #                 _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-nsteps": 1000}
    #             }
    #         },
    #     }
    #     step_mdrun = StepGMXMDrun(**step_conf)
    #     topol = GromacsTopol()
    #     step_mdrun = StepGMXMDrun(**step_conf)
    #     wf = WorkFlow()
    #     wf.workflow_data.gmx_topol = topol
    #     step_mdrun.set_workflow_object(wf)
    #     step_mdrun.data.generic.add_file(
    #         GenericData(file_name=_SGE.STD_TPR, file_data=self.tpr, argument=True)
    #     )
    #     step_mdrun.execute()
    #     out_path = os.path.join(self._test_dir, "confout.gro")
    #     step_mdrun.write_generic_by_extension(self._test_dir, "gro")
    #     stat_inf = os.stat(out_path)
    #     self.assertEqual(stat_inf.st_size, 3224484)
