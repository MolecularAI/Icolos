from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.containers.gmx_state import GromacsState
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.gromacs.trjconv import StepGMXTrjconv

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_Trjconv(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC, "rb") as f:
            xtc = f.read()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR, "rb") as f:
            tpr = f.read()
        with open(
            attach_root_path(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO), "r"
        ) as f:
            struct = f.readlines()

        self.topol = GromacsState()
        self.topol.tprs = [GenericData(_SGE.STD_TPR, file_data=tpr)]
        self.topol.trajectories = [GenericData(_SGE.STD_XTC, file_data=xtc)]
        self.topol.structures = [GenericData(_SGE.STD_STRUCTURE, struct)]

    def test_trjconv(self):
        step_conf = {
            _SBE.STEPID: "test_trjconv",
            _SBE.STEP_TYPE: "trjconv",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"],
                _SBE.SETTINGS_ADDITIONAL: {_SBE.PIPE_INPUT: "echo -ne 1 0"},
            },
        }

        step_trjconv = StepGMXTrjconv(**step_conf)
        wf = WorkFlow()
        wf.workflow_data.gmx_state = self.topol
        step_trjconv.set_workflow_object(wf)
        step_trjconv.execute()
        out_path = os.path.join(self._test_dir, "traj.xtc")
        step_trjconv.get_topol().write_trajectory(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 607000)
