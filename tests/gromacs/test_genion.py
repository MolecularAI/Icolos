from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
from icolos.core.containers.gromacs_topol import GromacsTopol
from icolos.core.workflow_steps.gromacs.genion import StepGMXGenion
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_Genion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP), "r") as f:
            topol = f.readlines()
        with open(attach_root_path(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR), "rb") as f:
            self.tpr = f.read()
        with open(
            attach_root_path(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO), "r"
        ) as f:
            struct = f.readlines()
        self.topol = GromacsTopol()
        self.topol.structure = struct
        self.topol.top_lines = topol

    def test_genion_run(self):
        step_conf = {
            _SBE.STEPID: "test_genion",
            _SBE.STEP_TYPE: "genion",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-neutral"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-pname": "NA",
                        "-nname": "CL",
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {_SBE.PIPE_INPUT: "3"},
            },
        }

        step_genion = StepGMXGenion(**step_conf)
        step_genion.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.tpr, argument=True)
        )
        wf = WorkFlow()
        wf.workflow_data.gmx_topol = self.topol
        step_genion.set_workflow_object(wf)
        step_genion.execute()

        out_path = os.path.join(self._test_dir, "confout.gro")
        step_genion.write_generic_by_name(self._test_dir, "confout.gro")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 2102900)
