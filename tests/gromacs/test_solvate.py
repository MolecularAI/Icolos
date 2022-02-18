from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.containers.gromacs_topol import GromacsTopol
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.gromacs.solvate import StepGMXSolvate

_SBE = StepBaseEnum
_SGE = StepGromacsEnum()


class Test_Solvate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO, "r") as f:
            data = f.readlines()

        self.topol = GromacsTopol()
        self.topol.parse(PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP)
        self.topol.structures = [GenericData(_SGE.STD_STRUCTURE, file_data=data)]

    def test_solvate(self):
        step_conf = {
            _SBE.STEPID: "test_solvate",
            _SBE.STEP_TYPE: "solvate",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                }
            },
        }

        step_solvate = StepGMXSolvate(**step_conf)
        wf = WorkFlow()
        wf.workflow_data.gmx_topol = self.topol
        step_solvate.set_workflow_object(wf)
        step_solvate.execute()

        out_path = os.path.join(self._test_dir, "confout.gro")
        step_solvate.get_topol().write_structure(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 142800)
