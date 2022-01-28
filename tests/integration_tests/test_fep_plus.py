import unittest
import os
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGlideEnum, TokenGuardEnum

_WE = WorkflowEnum()
_SBE = StepBaseEnum
_SGE = StepGlideEnum()
_TE = TokenGuardEnum()


class TestFEPPlusWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls._test_dir = attach_root_path("tests/junk/integration")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_fep_plus_workflow(self):

        conf = {
            _WE.HEADER: {
                _WE.ID: "Docking/FEP+ combined workflow",
                _WE.DESCRIPTION: "test setup for FEP+ integration",
                _WE.ENVIRONMENT: {_WE.ENVIRONMENT_EXPORT: []},
                _WE.GLOBAL_VARIABLES: {
                    "smiles": "3,4-DIAMINOBENZOTRIFLUORIDE:Nc1ccc(cc1N)C(F)(F)F"
                },
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "initialization_smile",
                    _SBE.STEP_TYPE: "initialization",
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "{smiles}",
                                _SBE.INPUT_SOURCE_TYPE: "string",
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "Ligprep",
                    _SBE.STEP_TYPE: "ligprep",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws",
                        "parallelization": {"cores": 2, "max_length_sublists": 1},
                        "failure_policy": {"n_tries": 3},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-epik"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-ph": 7.0,
                                "-pht": 2.0,
                                "-s": 10,
                                "-bff": 14,
                                "-HOST": "localhost",
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "filter_file": {"Total_charge": "!= 0"}
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "initialization_smile",
                                _SBE.INPUT_SOURCE_TYPE: "step",
                            }
                        ]
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.INPUT_COMPOUNDS: {
                                _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_ENUMERATIONS
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "{entrypoint_dir}/ligprep_enums.sdf",
                                _SBE.STEP_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "SDF",
                            },
                        }
                    ],
                },
                {
                    _SBE.STEPID: "Glide",
                    _SBE.STEP_TYPE: "glide",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws",
                        _SBE.EXEC_PARALLELIZATION: {
                            _SBE.EXEC_PARALLELIZATION_CORES: 4,
                            _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                        },
                        _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 3},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-HOST": "localhost"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "configuration": {
                                "AMIDE_MODE": "trans",
                                "EXPANDED_SAMPLING": "True",
                                "GRIDFILE": [PATHS_EXAMPLEDATA.PRIME_COX2_GRID],
                                "NENHANCED_SAMPLING": "1",
                                "POSE_OUTTYPE": "poseviewer",
                                "POSES_PER_LIG": "1",
                                "POSTDOCK_NPOSE": "25",
                                "POSTDOCKSTRAIN": "True",
                                "PRECISION": "SP",
                                "REWARD_INTRA_HBONDS": "True",
                            }
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "Ligprep",
                                _SBE.INPUT_SOURCE_TYPE: "step",
                            }
                        ]
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.INPUT_COMPOUNDS: {"category": "conformers"},
                            "destination": {
                                "resource": "{entrypoint_dir}/tests/junk/docked_conformers_cox2_actives.sdf",
                                _SBE.STEP_TYPE: "file",
                                "format": "SDF",
                            },
                        }
                    ],
                },
                {
                    _SBE.STEPID: "FEP_plus_setup",
                    _SBE.STEP_TYPE: "fep_plus_setup",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws"
                    },
                    _SBE.SETTINGS: {},
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "Glide",
                                _SBE.INPUT_SOURCE_TYPE: "step",
                                "target_field": _SBE.INPUT_COMPOUNDS,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "FEP_plus_exec",
                    _SBE.STEP_TYPE: "fep_plus_exec",
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws"
                    },
                    _SBE.TOKEN_GUARD: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws",
                        _SBE.EXEC_BINARYLOCATION: "ssh 10.220.1.4 /opt/schrodinger/suite/installations/default",
                        _TE.TG_TOKEN_POOLS: {"FEP_GPGPU": 16},
                        "wait_interval_seconds": 30,
                        "wait_limit_seconds": 0,
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-JOBNAME": "test",
                                "-HOST": "fep-compute",
                            },
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "Glide",
                                _SBE.INPUT_SOURCE_TYPE: "step",
                                "target_field": _SBE.INPUT_COMPOUNDS,
                            }
                        ],
                        "generic": [
                            {_SBE.INPUT_SOURCE: "FEP_plus_setup", "extension": "fmp"}
                        ],
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.INPUT_COMPOUNDS: {
                                _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS,
                                _SBE.WRITEOUT_COMP_SELECTED_TAGS: [
                                    "dG",
                                    "docking_score",
                                ],
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION: os.path.join(
                                    self._test_dir, "fep_scored_conformers.csv"
                                ),
                                _SBE.STEP_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "CSV",
                            },
                        }
                    ],
                },
            ],
        }

        wflow = WorkFlow(**conf)
        wflow.initialize()
        wflow.execute()

        out_path = os.path.join(self._test_dir, "fep_scored_conformers.csv")
        stat_inf = os.stat(out_path)
        self.assertGreaterEqual(stat_inf.st_size, 4252)
