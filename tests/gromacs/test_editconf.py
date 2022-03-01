from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.containers.gromacs_topol import GromacsTopol
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.workflow_steps.gromacs.editconf import StepGMXEditConf
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_Editconf(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(
            attach_root_path(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO), "r"
        ) as f:
            self.struct = f.readlines()
        self.topol = GromacsTopol()
        self.topol.structures = [GenericData(_SGE.STD_STRUCTURE, file_data=self.struct)]

    def test_editconf_wf_input(self):
        """
        Takes input from a workflow object
        """
        step_conf = {
            _SBE.STEPID: "test_editconf",
            _SBE.STEP_TYPE: "editconf",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-d": "1.0",
                        "-bt": "dodecahedron",
                    },
                }
            },
        }

        step_editconf = StepGMXEditConf(**step_conf)

        wf = WorkFlow()
        wf.workflow_data.gmx_topol = self.topol
        step_editconf.set_workflow_object(wf)

        step_editconf.execute()
        out_path = os.path.join(self._test_dir, "confout.gro")
        step_editconf.get_topol().write_structure(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 2102964)

    def test_editconf_external_input(self):
        step_conf = {
            _SBE.STEPID: "test_editconf",
            _SBE.STEP_TYPE: "editconf",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-d": "1.0",
                        "-bt": "dodecahedron",
                    },
                }
            },
        }

        step_editconf = StepGMXEditConf(**step_conf)

        wf = WorkFlow()
        wf.workflow_data.gmx_topol = self.topol
        step_editconf.set_workflow_object(wf)
        step_editconf.data.generic.add_file(
            GenericData(file_name=_SGE.STD_STRUCTURE, file_data=self.struct)
        )

        step_editconf.execute()
        out_path = os.path.join(self._test_dir, "confout.gro")
        step_editconf.get_topol().write_structure(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 2102964)
