from icolos.core.workflow_steps.gromacs.clusters_ts import StepClusterTS
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_ts_cluster(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_TS_CLUSTERS), "r") as f:
            self.data = f.read()

    def test_ts_cluster(self):
        step_conf = {
            _SBE.STEPID: "test_ts_cluster",
            _SBE.STEP_TYPE: "ts_cluster",
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load R"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "lengths": "10001",
                        "clustersNumber": "13",
                        "mdEngine": "GROMACS",
                    },
                }
            },
        }

        step_ts_cluster = StepClusterTS(**step_conf)
        step_ts_cluster.data.generic.add_file(
            GenericData(
                file_name="clusters_ts_example.xvg", file_data=self.data, argument=True
            )
        )

        step_ts_cluster.execute()

        out_path = os.path.join(self._test_dir, "clusters_ts.png")
        step_ts_cluster.write_generic_by_name(self._test_dir, "clusters_ts.png")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 36102)
