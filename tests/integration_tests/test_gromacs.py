import unittest
import os
from icolos.loggers.steplogger import StepLogger
from icolos.utils.entry_point_functions.logging_helper_functions import (
    initialize_logging,
)
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from tests.tests_paths import (
    MAIN_CONFIG,
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum

_WE = WorkflowEnum()
_SBE = StepBaseEnum
_SGE = StepGromacsEnum()
_LE = LoggingConfigEnum()


class Test_GROMACS_MD(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/integration")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    # def test_MD_holo_fpocket(self):
    #     """
    #     run a full gromacs MD workflow from a pdb structure, then fpocket on the resulting trajectory
    #     MDPocket is run on the holo structure
    #     """

    #     conf = {
    #         _WE.HEADER: {
    #             _WE.ID: "gromacs_test",
    #             _WE.DESCRIPTION: "full md run with gromacs, with MDpocket run to extract descriptors for binding pocket",
    #             _WE.ENVIRONMENT: {
    #                 _WE.ENVIRONMENT_EXPORT: [
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_DD_COMMS",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_PME_PP_COMMS",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_FORCE_UPDATE_DEFAULT_GPU",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                 ]
    #             },
    #             _WE.GLOBAL_VARIABLES: {
    #                 "file_base": os.path.join(
    #                     MAIN_CONFIG["ICOLOS_TEST_DATA"], "gromacs/protein"
    #                 ),
    #                 "output_dir": attach_root_path("tests/junk/integration"),
    #             },
    #         },
    #         _WE.STEPS: [
    #             {
    #                 _SBE.STEPID: "01_pdb2gmx",
    #                 _SBE.STEP_TYPE: "pdb2gmx",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-water": "tip3p",
    #                             "-ff": "amber03",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: attach_root_path(
    #                                 PATHS_1UYD.PDB_PATH
    #                             ),
    #                             _SBE.INPUT_EXTENSION: "pdb",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "02_editconf",
    #                 _SBE.STEP_TYPE: "editconf",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-d": "1.0",
    #                             "-bt": "dodecahedron",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "03_solvate",
    #                 _SBE.STEP_TYPE: "solvate",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-cs": "spc216"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "04_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/ions.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "05_genion",
    #                 _SBE.STEP_TYPE: "genion",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-neutral"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-pname": "NA",
    #                             "-nname": "CL",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "pipe_input": "SOL",
    #                     },
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "06_grompp_eminim",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/minim.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "07_eminim_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "08_nvt_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": True,
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/nvt_equil.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "09_nvt_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "10_npt_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": True,
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/npt_equil.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "11_npt_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "12_prod_md_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                         "make_ndx_command": "auto",
    #                         "fields": {"nsteps": "5000"},
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/md.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "13_prod_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-nb": "gpu",
    #                             "-bonded": "gpu",
    #                             "-pme": "gpu",
    #                         },
    #                     }
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "12_prod_md_grompp",
    #                             _SBE.INPUT_EXTENSION: "tpr",
    #                         }
    #                     ]
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GMX: {
    #                             _SBE.WRITEOUT_GENERIC_KEY: "xtc, gro, tpr"
    #                         },
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}",
    #                         },
    #                     }
    #                 ],
    #             },
    #             {
    #                 _SBE.STEPID: "14_trjconv",
    #                 _SBE.STEP_TYPE: "trjconv",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"]
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {"pipe_input": "echo -ne 1 0"},
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/trjconv",
    #                         },
    #                     }
    #                 ],
    #             },
    #             {
    #                 _SBE.STEPID: "15_MDpocket",
    #                 _SBE.STEP_TYPE: "mdpocket",
    #                 _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load fpocket"},
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {}
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {"format": "gromacs"},
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "pdb"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}",
    #                             _SBE.WRITEOUT_DESTINATION_MODE: "dir",
    #                         },
    #                     }
    #                 ],
    #             },
    #         ],
    #     }
    #     export_unit_test_env_vars()
    #     wflow = WorkFlow(**conf)
    #     wflow.initialize()

    #     self.assertEqual(len(wflow.steps), 15)

    #     wflow.execute()

    #     out_path = os.path.join(self._test_dir, "traj.xtc")
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 62400)

    # def test_md_ligparam_external_tpr(self):
    #     """
    #     End to end gromacs workflow with ligand parametrisation and read-in of external tprs to "overwrite" existing files
    #     Checks both generic file handling and gromacs topol play nicely together
    #     """

    #     conf = {
    #         _WE.HEADER: {
    #             _WE.ID: "gromacs_test_ligparam",
    #             _WE.DESCRIPTION: "full md run with gromacs, with ligand parametrisation",
    #             _WE.ENVIRONMENT: {
    #                 _WE.ENVIRONMENT_EXPORT: [
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_DD_COMMS",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_PME_PP_COMMS",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_FORCE_UPDATE_DEFAULT_GPU",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                 ]
    #             },
    #             _WE.GLOBAL_VARIABLES: {
    #                 "output_dir": attach_root_path("tests/junk/integration"),
    #                 "file_base": PATHS_EXAMPLEDATA.GROMACS_PROTEIN_FILE_BASE,
    #             },
    #         },
    #         _WE.STEPS: [
    #             {
    #                 _SBE.STEPID: "01_pdb2gmx",
    #                 _SBE.STEP_TYPE: "pdb2gmx",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-water": "tip3p",
    #                             "-ff": "amber03",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: attach_root_path(
    #                                 PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE
    #                             ),
    #                             _SBE.INPUT_EXTENSION: "pdb",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "02_editconf",
    #                 _SBE.STEP_TYPE: "editconf",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-d": "1.2",
    #                             "-bt": "dodecahedron",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #                 _SBE.INPUT: {},
    #             },
    #             {
    #                 _SBE.STEPID: "03_solvate",
    #                 _SBE.STEP_TYPE: "solvate",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-cs": "spc216"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "04_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/ions.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "05_genion",
    #                 _SBE.STEP_TYPE: "genion",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-neutral"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-pname": "NA",
    #                             "-nname": "CL",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "pipe_input": "SOL",
    #                     },
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "06_grompp_eminim",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-maxwarn": 50},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/minim.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "07_eminim_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "08_nvt_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": True,
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/nvt_equil.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "09_nvt_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                     _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                 },
    #                 _SBE.SETTINGS_ADDITIONAL: {},
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "08_nvt_grompp",
    #                             _SBE.INPUT_EXTENSION: "tpr",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "10_npt_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": True,
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/npt_equil.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "11_npt_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                     _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "10_npt_grompp",
    #                             _SBE.INPUT_EXTENSION: "tpr",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "12_prod_md_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-n": "index.ndx",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                         "fields": {"nsteps": "10000"},
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/md.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "13_prod_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-nb": "gpu",
    #                             "-bonded": "gpu",
    #                             "-pme": "gpu",
    #                         },
    #                     }
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GMX: {
    #                             _SBE.WRITEOUT_GENERIC_KEY: "xtc, gro, tpr"
    #                         },
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}",
    #                         },
    #                     },
    #                     {
    #                         _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "log"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/md_0_1.log",
    #                         },
    #                     },
    #                 ],
    #             },
    #             {
    #                 _SBE.STEPID: "14_trjconv",
    #                 _SBE.STEP_TYPE: "trjconv",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-pbc": "mol",
    #                             "-n": "index.ndx",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "pipe_input": "Protein_Other System"
    #                     },
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "15_trjconv",
    #                 _SBE.STEP_TYPE: "trjconv",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-fit": "rot+trans",
    #                             "-n": "index.ndx",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "pipe_input": "Protein_Other System"
    #                     },
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/trjconv",
    #                         },
    #                     }
    #                 ],
    #             },
    #             {
    #                 _SBE.STEPID: "15_mmgbsa",
    #                 _SBE.STEP_TYPE: "gmx_mmpbsa",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA/1.4.3-foss-2021a-CUDA-11.3.1"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_FLAGS: []},
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "coupling_groups": "Protein Other",
    #                         "pipe_input": "Protein Other",
    #                         "forcefield": MAIN_CONFIG["FORCEFIELD"],
    #                     },
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "props"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}",
    #                         },
    #                     }
    #                 ],
    #             },
    #         ],
    #     }

    #     export_unit_test_env_vars()
    #     wflow = WorkFlow(**conf)
    #     wflow.initialize()

    #     self.assertEqual(len(wflow.steps), 16)

    #     wflow.execute()

    #     out_path = os.path.join(self._test_dir, "traj.xtc")
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 316516)

    #     out_path = os.path.join(self._test_dir, "ICOLOS_PROPS.dat")
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 200)

    # def test_md_dna_1_replica(self):
    #     """
    #     End to end gromacs workflow with ligand parametrisation and read-in of external tprs to "overwrite" existing files
    #     Checks both generic file handling and gromacs topol play nicely together
    #     """

    #     conf = {
    #         _WE.HEADER: {
    #             _WE.ID: "gromacs_test_ligparam",
    #             _WE.DESCRIPTION: "full md run with gromacs, with ligand parametrisation",
    #             _WE.ENVIRONMENT: {
    #                 _WE.ENVIRONMENT_EXPORT: [
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_DD_COMMS",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_PME_PP_COMMS",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                     {
    #                         _WE.ENVIRONMENT_EXPORT_KEY: "GMX_FORCE_UPDATE_DEFAULT_GPU",
    #                         _WE.ENVIRONMENT_EXPORT_VALUE: "True",
    #                     },
    #                 ]
    #             },
    #             _WE.GLOBAL_VARIABLES: {
    #                 "forcefield": "<path>/gmx_workflow/forcefields/amber14sb_OL15.ff",
    #                 "output_dir": attach_root_path("tests/junk/integration"),
    #                 "file_base": PATHS_EXAMPLEDATA.GROMACS_PROTEIN_FILE_BASE,
    #             },
    #         },
    #         _WE.STEPS: [
    #             {
    #                 _SBE.STEPID: "01_pdb2gmx",
    #                 _SBE.STEP_TYPE: "pdb2gmx",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-water": "tip3p",
    #                             "-ff": "amber03",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: attach_root_path(
    #                                 PATHS_EXAMPLEDATA.GROMACS_DNA_STRUCTURE
    #                             ),
    #                             _SBE.INPUT_EXTENSION: "pdb",
    #                         }
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "02_editconf",
    #                 _SBE.STEP_TYPE: "editconf",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-d": "1.2",
    #                             "-bt": "dodecahedron",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "03_solvate",
    #                 _SBE.STEP_TYPE: "solvate",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-cs": "spc216"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "04_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/ions.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "05_genion",
    #                 _SBE.STEP_TYPE: "genion",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-neutral"],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-pname": "NA",
    #                             "-nname": "CL",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "pipe_input": "SOL",
    #                     },
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "06_grompp_eminim",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/minim.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "07_eminim_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "08_nvt_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": True,
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/nvt_equil.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "09_nvt_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                     _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "10_npt_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": True,
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/npt_equil.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "11_npt_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                     _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
    #                     _SBE.SETTINGS_ADDITIONAL: {},
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "12_prod_md_grompp",
    #                 _SBE.STEP_TYPE: "grompp",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-n": "index.ndx",
    #                         },
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {
    #                         "-r": False,
    #                         "fields": {"nsteps": "5000"},
    #                         "make_ndx_command": "auto",
    #                     },
    #                 },
    #                 _SBE.INPUT: {
    #                     _SBE.INPUT_GENERIC: [
    #                         {
    #                             _SBE.INPUT_SOURCE: "{file_base}/md.mdp",
    #                             _SBE.INPUT_EXTENSION: "mdp",
    #                         },
    #                     ]
    #                 },
    #             },
    #             {
    #                 _SBE.STEPID: "13_prod_mdrun",
    #                 _SBE.STEP_TYPE: "mdrun",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                             "-nb": "gpu",
    #                             "-bonded": "gpu",
    #                             "-pme": "gpu",
    #                         },
    #                     }
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/md_0_1.xtc",
    #                             _SBE.STEP_TYPE: "file",
    #                             _SBE.WRITEOUT_DESTINATION_FORMAT: "TXT",
    #                         },
    #                     },
    #                     {
    #                         _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "log"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/md_0_1.log",
    #                             _SBE.STEP_TYPE: "file",
    #                             _SBE.WRITEOUT_DESTINATION_FORMAT: "TXT",
    #                         },
    #                     },
    #                     {
    #                         _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "gro"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/md_0_1.gro",
    #                             _SBE.STEP_TYPE: "file",
    #                             _SBE.WRITEOUT_DESTINATION_FORMAT: "TXT",
    #                         },
    #                     },
    #                     {
    #                         _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "tpr"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/md_0_1.tpr",
    #                             _SBE.STEP_TYPE: "file",
    #                             _SBE.WRITEOUT_DESTINATION_FORMAT: "TXT",
    #                         },
    #                     },
    #                 ],
    #             },
    #             {
    #                 _SBE.STEPID: "14_trjconv",
    #                 _SBE.STEP_TYPE: "trjconv",
    #                 _SBE.EXEC: {
    #                     _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
    #                 },
    #                 _SBE.SETTINGS: {
    #                     _SBE.SETTINGS_ARGUMENTS: {
    #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"]
    #                     },
    #                     _SBE.SETTINGS_ADDITIONAL: {"pipe_input": "echo -ne 1 0"},
    #                 },
    #                 _SBE.WRITEOUT: [
    #                     {
    #                         _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
    #                         _SBE.WRITEOUT_DESTINATION: {
    #                             _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/trjconv",
    #                         },
    #                     }
    #                 ],
    #             },
    #         ],
    #     }

    #     export_unit_test_env_vars()
    #     wflow = WorkFlow(**conf)
    #     wflow.initialize()

    #     self.assertEqual(len(wflow.steps), 14)

    #     wflow.execute()

    #     out_path = os.path.join(self._test_dir, "trjconv/traj.xtc")
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 44700)

    def test_md_complex_mmgbsa_hrex(self):
        """
        End to end gromacs workflow with ligand parametrisation and read-in of external tprs to "overwrite" existing files
        Checks both generic file handling and gromacs topol play nicely together
        """

        conf = {
            _WE.HEADER: {
                _WE.ID: "gromacs_test_ligparam",
                _WE.DESCRIPTION: "full md run with gromacs, with ligand parametrisation",
                _WE.ENVIRONMENT: {
                    _WE.ENVIRONMENT_EXPORT: [
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_DD_COMMS",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "True",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_PME_PP_COMMS",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "True",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMX_FORCE_UPDATE_DEFAULT_GPU",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "True",
                        },
                    ]
                },
                _WE.GLOBAL_VARIABLES: {
                    "forcefield": "<path>/gmx_workflow/forcefields/amber14sb_OL15.ff",
                    "output_dir": attach_root_path("tests/junk/integration"),
                    "file_base": PATHS_EXAMPLEDATA.GROMACS_PROTEIN_FILE_BASE,
                },
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "01_pdb2gmx",
                    _SBE.STEP_TYPE: "pdb2gmx",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-water": "tip3p",
                                "-ff": "amber03",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: attach_root_path(
                                    PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE
                                ),
                                _SBE.INPUT_EXTENSION: "pdb",
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "02_editconf",
                    _SBE.STEP_TYPE: "editconf",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-d": "1.2",
                                "-bt": "dodecahedron",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                    _SBE.INPUT: {},
                },
                {
                    _SBE.STEPID: "03_solvate",
                    _SBE.STEP_TYPE: "solvate",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-cs": "spc216"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
                {
                    _SBE.STEPID: "04_grompp",
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
                            "-r": False,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/ions.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "05_genion",
                    _SBE.STEP_TYPE: "genion",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-neutral"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-pname": "NA",
                                "-nname": "CL",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "pipe_input": "SOL",
                        },
                    },
                },
                {
                    _SBE.STEPID: "06_grompp_eminim",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-maxwarn": 50},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"-r": False, _SGE.REPLICAS: 4},
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/minim.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "07_eminim_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                        _SBE.EXEC_RESOURCES: {_SBE.EXEC_RESOURCES_TASKS: 4},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: True},
                    },
                },
                {
                    _SBE.STEPID: "08_nvt_grompp",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "-r": True,
                            "make_ndx_command": "auto",
                            _SGE.REPLICAS: 4,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/nvt_equil.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "09_nvt_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                        _SBE.EXEC_RESOURCES: {_SBE.EXEC_RESOURCES_TASKS: 4},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: True},
                        _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                    },
                },
                {
                    _SBE.STEPID: "10_npt_grompp",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "-r": True,
                            "make_ndx_command": "auto",
                            _SGE.REPLICAS: 4,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/npt_equil.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "11_npt_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                        _SBE.EXEC_RESOURCES: {_SBE.EXEC_RESOURCES_TASKS: 4},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                        _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: True},
                    },
                },
                {
                    _SBE.STEPID: "12_prod_md_grompp",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-n": "index.ndx",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "-r": False,
                            "fields": {"nsteps": "10000"},
                            "make_ndx_command": "auto",
                            _SGE.REPLICAS: 4,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/md.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "13_prod_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                        _SBE.EXEC_RESOURCES: {_SBE.EXEC_RESOURCES_TASKS: 4},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-nb": "gpu",
                                "-bonded": "gpu",
                                "-pme": "gpu",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: True},
                    },
                },
                {
                    _SBE.STEPID: "14_trjconv",
                    _SBE.STEP_TYPE: "trjconv",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"]
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "pipe_input": "Protein_Other System",
                            "-n": "index.ndx",
                        },
                    },
                },
                {
                    _SBE.STEPID: "15_trjconv",
                    _SBE.STEP_TYPE: "trjconv",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-fit": "rot+trans",
                                "-n": "index.ndx",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "pipe_input": "Protein_Other System"
                        },
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/processed",
                            },
                        }
                    ],
                },
                {
                    _SBE.STEPID: "15_mmgbsa",
                    _SBE.STEP_TYPE: "gmx_mmpbsa",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA/1.4.3-foss-2021a-CUDA-11.3.1"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_FLAGS: []},
                        _SBE.SETTINGS_ADDITIONAL: {
                            "coupling_groups": "Protein Other",
                            "pipe_input": "Protein Other",
                        },
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.WRITEOUT_GMX: {
                                _SBE.WRITEOUT_GENERIC_KEY: "props, xtc"
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/mmgbsa",
                            },
                        }
                    ],
                },
            ],
        }

        export_unit_test_env_vars()
        wflow = WorkFlow(**conf)
        log_conf = attach_root_path(_LE.PATH_CONFIG_DEBUG)
        _ = initialize_logging(log_conf_path=log_conf, workflow_conf=conf)
        wflow.initialize()

        self.assertEqual(len(wflow.steps), 15)

        wflow.execute()

        out_path = os.path.join(self._test_dir, "mmgbsa/ICOLOS_PROPS.dat")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 44700)

    def test_md_complex_mmgbsa_3_replicas(self):
        """
        End to end gromacs workflow with ligand parametrisation and read-in of external tprs to "overwrite" existing files
        Checks both generic file handling and gromacs topol play nicely together
        """

        conf = {
            _WE.HEADER: {
                _WE.ID: "gromacs_test_ligparam",
                _WE.DESCRIPTION: "full md run with gromacs, with ligand parametrisation",
                _WE.ENVIRONMENT: {
                    _WE.ENVIRONMENT_EXPORT: [
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_DD_COMMS",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "True",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMX_GPU_PME_PP_COMMS",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "True",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMX_FORCE_UPDATE_DEFAULT_GPU",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "True",
                        },
                    ]
                },
                _WE.GLOBAL_VARIABLES: {
                    "output_dir": attach_root_path("tests/junk/integration"),
                    "file_base": PATHS_EXAMPLEDATA.GROMACS_PROTEIN_FILE_BASE,
                },
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "01_pdb2gmx",
                    _SBE.STEP_TYPE: "pdb2gmx",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ignh"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-water": "tip3p",
                                "-ff": "amber03",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.CHARGE_METHOD: "gas"},
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: attach_root_path(
                                    PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE
                                ),
                                _SBE.INPUT_EXTENSION: "pdb",
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "02_editconf",
                    _SBE.STEP_TYPE: "editconf",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-d": "1.2",
                                "-bt": "dodecahedron",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                    _SBE.INPUT: {},
                },
                {
                    _SBE.STEPID: "03_solvate",
                    _SBE.STEP_TYPE: "solvate",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-cs": "spc216"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
                {
                    _SBE.STEPID: "04_grompp",
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
                            "-r": False,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/ions.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "05_genion",
                    _SBE.STEP_TYPE: "genion",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-neutral"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-pname": "NA",
                                "-nname": "CL",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "pipe_input": "SOL",
                        },
                    },
                },
                {
                    _SBE.STEPID: "06_grompp_eminim",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-maxwarn": 50},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"-r": False, _SGE.REPLICAS: 3},
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/minim.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "07_eminim_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: False},
                    },
                },
                {
                    _SBE.STEPID: "08_nvt_grompp",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "-r": True,
                            "make_ndx_command": "auto",
                            _SGE.REPLICAS: 3,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/nvt_equil.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "09_nvt_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                        _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                    },
                    _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: False},
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "08_nvt_grompp",
                                _SBE.INPUT_EXTENSION: "tpr",
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "10_npt_grompp",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-n": "index.ndx"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "-r": True,
                            "make_ndx_command": "auto",
                            _SGE.REPLICAS: 3,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/npt_equil.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "11_npt_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                        _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: False},
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "10_npt_grompp",
                                _SBE.INPUT_EXTENSION: "tpr",
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "12_prod_md_grompp",
                    _SBE.STEP_TYPE: "grompp",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-n": "index.ndx",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "-r": False,
                            "fields": {"nsteps": "5000"},
                            "make_ndx_command": "auto",
                            _SGE.REPLICAS: 3,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "{file_base}/md.mdp",
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "13_prod_mdrun",
                    _SBE.STEP_TYPE: "mdrun",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a",
                        _SBE.EXEC_RESOURCES: {_SBE.EXEC_RESOURCES_TASKS: 4},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-nb": "gpu",
                                "-bonded": "gpu",
                                "-pme": "gpu",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {_SGE.MULTIDIR: False},
                    },
                },
                {
                    _SBE.STEPID: "14_trjconv",
                    _SBE.STEP_TYPE: "trjconv",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-center"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-pbc": "mol",
                                "-n": "index.ndx",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "pipe_input": "Protein_Other System"
                        },
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}",
                            },
                        }
                    ],
                },
                {
                    _SBE.STEPID: "15_trjconv",
                    _SBE.STEP_TYPE: "trjconv",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-fit": "rot+trans",
                                "-n": "index.ndx",
                            }
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "pipe_input": "Protein_Other System"
                        },
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}",
                            },
                        }
                    ],
                },
                {
                    _SBE.STEPID: "15_mmgbsa",
                    _SBE.STEP_TYPE: "gmx_mmpbsa",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA/1.4.3-foss-2021a-CUDA-11.3.1"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_FLAGS: []},
                        _SBE.SETTINGS_ADDITIONAL: {
                            "coupling_groups": "Protein Other",
                            "pipe_input": "Protein Other",
                            "threads": 4,
                        },
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "dat"},
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "{output_dir}/scores.dat",
                            },
                        }
                    ],
                },
            ],
        }

        export_unit_test_env_vars()
        wflow = WorkFlow(**conf)
        log_conf = attach_root_path(_LE.PATH_CONFIG_DEBUG)
        _ = initialize_logging(log_conf_path=log_conf, workflow_conf=conf)
        wflow.initialize()

        self.assertEqual(len(wflow.steps), 16)

        wflow.execute()

        out_path = os.path.join(self._test_dir, "scores.dat")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 44700)
