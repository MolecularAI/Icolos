import unittest
import os
from icolos.core.composite_agents.workflow import WorkFlow
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
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum

_WE = WorkflowEnum()
_SBE = StepBaseEnum
_SGE = StepGromacsEnum()
_LE = LoggingConfigEnum()


class TestPMXrbfe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/integration/pmx")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def test_pmx_rbfe(self):
        conf = {
            _WE.HEADER: {
                _WE.ID: "pmx rbfe integration test",
                _WE.DESCRIPTION: "pmx rbfe run on single edge",
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
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "CONDA",
                            _WE.ENVIRONMENT_EXPORT_VALUE: MAIN_CONFIG["CONDA"],
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "PMX_PYTHON",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "${CONDA}/envs/pmx/bin/python",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "PMX",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "${CONDA}/envs/pmx/bin/pmx",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "GMXLIB",
                            _WE.ENVIRONMENT_EXPORT_VALUE: MAIN_CONFIG["PMX"]["GMXLIB"],
                        },
                    ]
                },
                _WE.GLOBAL_VARIABLES: {
                    "file_base": os.path.join(MAIN_CONFIG["ICOLOS_TEST_DATA"], "pmx")
                },
                _WE.GLOBAL_SETTINGS: {"single_directory": True},
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "01_pmx_setup",
                    _SBE.STEP_TYPE: "pmx_setup",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                        "parallelization": {"cores": 8},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            _SGE.CHARGE_METHOD: "gas",
                            "water": "tip3p",
                            "forcefield": "amber99sb-star-ildn-mut",
                            "replicas": 1,
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PMX_TNKS_PROTEIN,
                                _SBE.INPUT_EXTENSION: "pdb",
                            },
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PMX_TNKS_MAP,
                                _SBE.INPUT_EXTENSION: "log",
                            },
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PMX_MDP_FILES,
                                _SBE.INPUT_EXTENSION: "mdp",
                            },
                        ],
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PMX_TNKS_LIGANDS,
                                _SBE.INPUT_SOURCE_TYPE: "file",
                                _SBE.INPUT_FORMAT: "SDF",
                            }
                        ],
                        "work_dir": self._test_dir,
                    },
                },
                {
                    _SBE.STEPID: "02_pmx_atomMapping",
                    _SBE.STEP_TYPE: "pmx_atomMapping",
                    _SBE.EXEC: {
                        "parallelization": {"cores": 8},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
                {
                    _SBE.STEPID: "03_ligand_hybrid",
                    _SBE.STEP_TYPE: "pmx_ligandHybrid",
                    _SBE.EXEC: {
                        "parallelization": {"cores": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-cs": "spc216.gro"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
                {
                    _SBE.STEPID: "04_assemble_systems",
                    _SBE.STEP_TYPE: "pmx_assemble_systems",
                    _SBE.EXEC: {
                        "parallelization": {"cores": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                        _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["PMX"]["CLI_ENTRYPOINT"],
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
                {
                    _SBE.STEPID: "05-box_water_ions",
                    _SBE.STEP_TYPE: "pmx_box_water_ions",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"cores": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
                {
                    _SBE.STEPID: "06_prepare_simulations",
                    _SBE.STEP_TYPE: "pmx_prepare_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"cores": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "em"},
                    },
                },
                {
                    _SBE.STEPID: "06b_run_simulations",
                    _SBE.STEP_TYPE: "pmx_run_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"jobs": 2},
                        _SBE.EXEC_PLATFORM: "slurm",
                        _SBE.EXEC_RESOURCES: {
                            _SBE.EXEC_RESOURCES_MODULES: [
                                "GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
                            ],
                            _SBE.EXEC_RESOURCES_PARTITION: "gpu",
                            _SBE.EXEC_RESOURCES_GRES: "gpu:volta:1",
                            _SBE.EXEC_RESOURCES_CORES: "8",
                            _SBE.EXEC_RESOURCES_ADDITIONAL_LINES: [
                                'echo "hello, world!"'
                            ],
                        },
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "sim_type": "em",
                        },
                    },
                },
                {
                    _SBE.STEPID: "06c_prepare_simulations",
                    _SBE.STEP_TYPE: "pmx_prepare_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"jobs": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "nvt"},
                    },
                },
                {
                    _SBE.STEPID: "06d_run_simulations",
                    _SBE.STEP_TYPE: "pmx_run_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"jobs": 2},
                        _SBE.EXEC_PLATFORM: "slurm",
                        _SBE.EXEC_RESOURCES: {
                            _SBE.EXEC_RESOURCES_MODULES: [
                                "GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
                            ],
                            _SBE.EXEC_RESOURCES_PARTITION: "gpu",
                            _SBE.EXEC_RESOURCES_GRES: "gpu:volta:1",
                            _SBE.EXEC_RESOURCES_CORES: "8",
                        },
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "nvt"},
                    },
                },
                {
                    _SBE.STEPID: "08_prepare_simulations",
                    _SBE.STEP_TYPE: "pmx_prepare_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"cores": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "eq"},
                    },
                },
                {
                    _SBE.STEPID: "09_run_simulations",
                    _SBE.STEP_TYPE: "pmx_run_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"cores": 1},
                        _SBE.EXEC_PLATFORM: "slurm",
                        _SBE.EXEC_RESOURCES: {
                            _SBE.EXEC_RESOURCES_MODULES: [
                                "GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
                            ],
                            _SBE.EXEC_RESOURCES_PARTITION: "gpu",
                            _SBE.EXEC_RESOURCES_GRES: "gpu:volta:1",
                        },
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "eq"},
                    },
                },
                {
                    _SBE.STEPID: "10_prepare_simulations",
                    _SBE.STEP_TYPE: "pmx_prepare_transitions",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"cores": 8},
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2",
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "transitions"},
                    },
                },
                {
                    _SBE.STEPID: "11_run_simulations",
                    _SBE.STEP_TYPE: "pmx_run_simulations",
                    _SBE.EXEC: {
                        _SBE.EXEC_PARALLELIZATION: {"cores": 1},
                        _SBE.EXEC_PLATFORM: "slurm",
                        _SBE.EXEC_RESOURCES: {
                            _SBE.EXEC_RESOURCES_MODULES: [
                                "GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
                            ],
                            _SBE.EXEC_RESOURCES_PARTITION: "gpu",
                            _SBE.EXEC_RESOURCES_GRES: "gpu:volta:1",
                        },
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {"sim_type": "transitions"},
                    },
                },
                {
                    _SBE.STEPID: "12_prepare",
                    _SBE.STEP_TYPE: "pmx_run_analysis",
                    _SBE.EXEC: {_SBE.EXEC_PARALLELIZATION: {"cores": 8}},
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
                    },
                },
            ],
        }

        wflow = WorkFlow(**conf)
        wflow.initialize()

        log_conf = attach_root_path(_LE.PATH_CONFIG_DEBUG)
        _ = initialize_logging(log_conf_path=log_conf, workflow_conf=conf)
        self.assertEqual(len(wflow.steps), 14)
        wflow.execute()
        out_path = os.path.join(self._test_dir, "resultsAll.csv")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 223)
