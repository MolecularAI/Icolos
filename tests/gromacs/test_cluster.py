from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.gromacs.cluster import StepGMXCluster
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

SGE = StepGromacsEnum()
SBE = StepBaseEnum


class Test_Cluster(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC), "rb") as f:
            self.xtc = f.read()

        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR), "rb") as f:
            self.tpr = f.read()

        with open(
            attach_root_path(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO), "r"
        ) as f:
            self.structure = f.read()

    def test_cluster(self):
        step_conf = {
            SBE.STEPID: "test_cluster",
            SBE.STEP_TYPE: "cluster",
            SBE.EXEC: {
                SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
            },
            SBE.SETTINGS: {
                SBE.SETTINGS_ARGUMENTS: {
                    SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-dt": "1000",
                        "-n": "index.ndx",
                    },
                },
                SBE.SETTINGS_ADDITIONAL: {
                    SBE.PIPE_INPUT: "2 System",
                    SGE.MAKE_NDX_COMMAND: "1 & a P",
                },
            },
        }

        step_cluster = StepGMXCluster(**step_conf)
        step_cluster.data.generic.add_file(
            GenericData(file_name="tmp10249.xtc", file_data=self.xtc, argument=True)
        )
        step_cluster.data.generic.add_file(
            GenericData(file_name="tmp03942.tpr", file_data=self.tpr, argument=True)
        )
        step_cluster.data.generic.add_file(
            GenericData(
                file_name="structure.gro", file_data=self.structure, argument=True
            )
        )
        step_cluster.execute()
        out_path = os.path.join(self._test_dir, "clusters.pdb")
        step_cluster.write_generic_by_extension(self._test_dir, "pdb")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 2002553)
