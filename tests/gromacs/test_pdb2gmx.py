import shutil
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum, StepOpenFFEnum
from tests.tests_paths import (
    MAIN_CONFIG,
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.gromacs.pdb2gmx import StepGMXPdb2gmx

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_SOFE = StepOpenFFEnum()


class Test_Pdb2gmx(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")

        if os.path.isdir(cls._test_dir):
            shutil.rmtree(cls._test_dir)
        os.makedirs(cls._test_dir, exist_ok=True)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_1UYD.PDB_PATH, "r") as f:
            self.structure = f.read()
        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_PDB, "r") as f:
            self.holo_structure = f.read()

    # def test_pdb2gmx_run(self):
    #     step_conf = {
    #         _SBE.STEPID: "test_pdb2gmx",
    #         _SBE.STEP_TYPE: "pdb2gmx",
    #         _SBE.EXEC: {
    #             _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
    #         },
    #         _SBE.SETTINGS: {
    #             _SBE.SETTINGS_ARGUMENTS: {
    #                 _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
    #                 _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                     "-water": "tip4p",
    #                     "-ff": "amber03",
    #                 },
    #             },
    #             _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
    #         },
    #     }

    #     step_pdb2gmx = StepGMXPdb2gmx(**step_conf)
    #     step_pdb2gmx.data.generic.add_file(
    #         GenericData(
    #             file_name="structure.pdb", file_data=self.structure, argument=True
    #         )
    #     )
    #     step_pdb2gmx.execute()
    #     out_path = os.path.join(self._test_dir, "confout.gro")
    #     step_pdb2gmx.write_generic_by_extension(
    #         self._test_dir, _SGE.FIELD_KEY_STRUCTURE
    #     )
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 22300)

    def test_lig_param_gaff(self):
        step_conf = {
            _SBE.STEPID: "test_pdb2gmx_lig_param",
            _SBE.STEP_TYPE: "pdb2gmx",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {},
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-water": "tip4p",
                        "-ff": "amber03",
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
            },
        }

        step_lig_param = StepGMXPdb2gmx(**step_conf)
        step_lig_param._workflow_object = WorkFlow()
        step_lig_param.data.generic.add_file(
            GenericData(
                file_name="tmp_whatever01923.pdb", file_data=self.holo_structure
            )
        )
        step_lig_param.execute()
        out_path = os.path.join(self._test_dir, "confout.gro")
        step_lig_param.write_generic_by_extension(
            self._test_dir, _SGE.FIELD_KEY_STRUCTURE
        )
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 73800)

        out_path = os.path.join(self._test_dir, "topol.top")
        step_lig_param.write_generic_by_extension(self._test_dir, _SGE.FIELD_KEY_TOPOL)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 73800)

    # def test_lig_param_openff(self):
    #     step_conf = {
    #         _SBE.STEPID: "test_pdb2gmx_lig_param",
    #         _SBE.STEP_TYPE: "pdb2gmx",
    #         _SBE.EXEC: {
    #             _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
    #         },
    #         _SBE.SETTINGS: {
    #             _SBE.SETTINGS_ADDITIONAL: {},
    #             _SBE.SETTINGS_ARGUMENTS: {
    #                 _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
    #                 _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                     "-water": "tip4p",
    #                     "-ff": "amber03",
    #                 },
    #             },
    #             _SBE.SETTINGS_ADDITIONAL: {
    #                 _SGE.CHARGE_METHOD: "gas",
    #                 _SGE.PARAM_METHOD: _SGE.OPENFF,
    #                 _SOFE.FORCEFIELD: os.path.join(
    #                     MAIN_CONFIG["OPENFF_FORCEFIELDS"], "openff-2.0.0.offxml"
    #                 ),
    #             },
    #         },
    #     }

    #     step_lig_param = StepGMXPdb2gmx(**step_conf)
    #     step_lig_param.data.generic.add_file(
    #         GenericData(
    #             file_name="tmp_whatever01923.pdb", file_data=self.holo_structure
    #         )
    #     )
    #     step_lig_param.execute()
    #     out_path = os.path.join(self._test_dir, "confout.gro")
    #     step_lig_param.write_generic_by_extension(
    #         self._test_dir, _SGE.FIELD_KEY_STRUCTURE
    #     )
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 73800)

    #     out_path = os.path.join(self._test_dir, "topol.top")
    #     step_lig_param.write_generic_by_extension(self._test_dir, _SGE.FIELD_KEY_TOPOL)
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 73800)
