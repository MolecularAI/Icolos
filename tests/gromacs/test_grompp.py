from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.core.workflow_steps.gromacs.grompp import StepGMXGrompp
from icolos.utils.general.files_paths import attach_root_path

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_Grompp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_GROMPP_INPUT_STRUCTURE, "r") as f:
            self.structure = f.read()
        with open(PATHS_EXAMPLEDATA.GROMACS_IONS_MDP, "r") as f:
            self.mdp = f.read()
        with open(PATHS_EXAMPLEDATA.GROMACS_GROMPP_TOPOL, "r") as f:
            self.topol = f.read()

    def test_grompp(self):
        step_conf = {
            _SBE.STEPID: "test_grompp",
            _SBE.STEP_TYPE: "grompp",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.FIELDS: {
                        "nsteps": 50,
                        "-nsteeps": 123,
                    },  # deliberate typo to check warning
                    _SGE.FORCEFIELD: "/projects/cc/mai/material/Icolos/forcefields/charmm36-feb2021.ff",
                    "-r": False,
                    _SGE.MAKE_NDX_COMMAND: "auto",
                },
            },
        }

        step_grompp = StepGMXGrompp(**step_conf)
        step_grompp.data.generic.add_file(
            GenericData(
                file_name="tmp029389.gro", file_data=self.structure, argument=True
            )
        )
        step_grompp.data.generic.add_file(
            GenericData(file_name="tmp03394.mdp", file_data=self.mdp, argument=True)
        )
        step_grompp.data.generic.add_file(
            GenericData(file_name="tmp91023.top", file_data=self.topol, argument=True)
        )

        step_grompp.execute()

        out_path = os.path.join(self._test_dir, "structure.tpr")
        step_grompp.write_generic_by_name(self._test_dir, "structure.tpr")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 596160)
