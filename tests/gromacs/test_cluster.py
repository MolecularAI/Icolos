from icolos.core.containers.generic import GenericData
from icolos.core.containers.gmx_state import GromacsState
from icolos.core.workflow_steps.gromacs.cluster import StepGMXCluster
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
SBE = StepBaseEnum


class Test_Cluster(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP), "r") as f:
            topol = f.readlines()
        with open(
            attach_root_path(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO), "r"
        ) as f:
            struct = f.readlines()

        self.topol = GromacsState()
        self.topol.structures = [GenericData(_SGE.STD_STRUCTURE, file_data=struct)]
        self.topol.top_lines = topol
        self.topol.set_tpr(path="", file=PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR)
        self.topol.set_trajectory(path="", file=PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC)

    def test_cluster(self):
        step_conf = {
            SBE.STEPID: "test_cluster",
            SBE.STEP_TYPE: "cluster",
            SBE.EXEC: {
                SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            SBE.SETTINGS: {
                SBE.SETTINGS_ARGUMENTS: {
                    SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-dt": "10",
                        "-n": "index.ndx",
                    },
                },
                SBE.SETTINGS_ADDITIONAL: {
                    SBE.PIPE_INPUT: "2 System",
                    _SGE.MAKE_NDX_COMMAND: "1 & a P",
                },
            },
        }

        step_cluster = StepGMXCluster(**step_conf)
        step_cluster.data.gmx_state = self.topol
        step_cluster.execute()
        out_path = os.path.join(self._test_dir, "clusters.pdb")
        step_cluster.write_generic_by_extension(self._test_dir, "pdb")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 7383600)
        out_path = os.path.join(self._test_dir, "cluster_id.xvg")
        step_cluster.write_generic_by_extension(self._test_dir, "xvg")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 700)
