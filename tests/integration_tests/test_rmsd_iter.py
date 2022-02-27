import unittest
import os
from tests.tests_paths import PATHS_1UYD
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGlideEnum


_WE = WorkflowEnum()
_SBE = StepBaseEnum
_SGE = StepGlideEnum()


class TestDockingWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls._test_dir = attach_root_path("tests/junk/integration")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_iterator_workflow(self):
        """
        Runs the RMSD-corrected docking workflow using multiple xtb settings in parallel
        """

        conf = {
            _WE.HEADER: {
                _WE.ID: "RMSD_rescoring",
                _WE.DESCRIPTION: "Run RMSD rescoring on docking pose",
                _WE.ENVIRONMENT: {
                    _WE.ENVIRONMENT_EXPORT: [
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "OE_LICENSE",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "/opt/scp/software/oelicense/1.0/oe_license.seq1",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "XTBHOME",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "/projects/cc/mai/binaries/xtb-6.4.0",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "XTBPATH",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "${XTBHOME}/share/xtb",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "PATH",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "${PATH}:${XTBHOME}/bin",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "PKG_CONFIG_PATH",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "${PKG_CONFIG_PATH}:${XTBHOME}/lib64/pkgconfig",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "PARA_ARCH",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "MPI",
                        },
                        {
                            _WE.ENVIRONMENT_EXPORT_KEY: "PARNODES",
                            _WE.ENVIRONMENT_EXPORT_VALUE: "6",
                        },
                    ]
                },
                _WE.GLOBAL_VARIABLES: {
                    "smiles": "3,4-DIAMINOBENZOTRIFLUORIDE:Nc1ccc(cc1N)C(F)(F)F;aspirin:O=C(C)Oc1ccccc1C(=O)O"
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
                        "prefix_execution": "module load schrodinger/2021-2-js-aws",
                        "parallelization": {"cores": 4, "max_length_sublists": 1},
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
                                "-HOST": "cpu-only",
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
                },
                {
                    _SBE.STEPID: "Glide",
                    _SBE.STEP_TYPE: "glide",
                    _SBE.EXEC: {
                        "prefix_execution": "module load schrodinger/2021-2-js-aws",
                        "parallelization": {"cores": 4, "max_length_sublists": 1},
                        "failure_policy": {"n_tries": 3},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-HOST": "cpu-only"},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "configuration": {
                                "AMIDE_MODE": "trans",
                                "EXPANDED_SAMPLING": "True",
                                "GRIDFILE": [PATHS_1UYD.GRID_PATH],
                                "NENHANCED_SAMPLING": "1",
                                "POSE_OUTTYPE": "ligandlib_sd",
                                "POSES_PER_LIG": "3",
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
                            _SBE.INPUT_COMPOUNDS: {
                                _SBE.WRITEOUT_COMP_CATEGORY: "conformers"
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "tests/junk/integration/rmsd_rescoring_docked_conformers.sdf",
                                _SBE.STEP_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "SDF",
                            },
                        },
                        {
                            _SBE.INPUT_COMPOUNDS: {
                                _SBE.WRITEOUT_COMP_CATEGORY: "conformers",
                                _SBE.WRITEOUT_COMP_SELECTED_TAGS: [
                                    "docking_score",
                                    "grid_id",
                                ],
                                _SBE.WRITEOUT_COMP_AGGREGATION: {
                                    _SBE.WRITEOUT_COMP_AGGREGATION_MODE: "best_per_compound",
                                    _WE.ENVIRONMENT_EXPORT_KEY: "docking_score",
                                    _SBE.WRITEOUT_COMP_AGGREGATION_HIGHESTISBEST: False,
                                },
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: "tests/junk/integration/rmsd_rescoring_docked_conformers.csv",
                                _SBE.STEP_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "CSV",
                            },
                        },
                    ],
                },
                {
                    _SBE.STEPID: "compound_filter",
                    _SBE.STEP_TYPE: "data_manipulation",
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ADDITIONAL: {
                            "action": "filter",
                            "filter_level": _SBE.INPUT_COMPOUNDS,
                            "criteria": "docking_score",
                            "return_n": 1,
                            "highest_is_best": False,
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {_SBE.INPUT_SOURCE: "Glide", _SBE.INPUT_SOURCE_TYPE: "step"}
                        ]
                    },
                },
                {
                    _SBE.STEPID: "test_iterator",
                    _SBE.STEP_TYPE: "iterator",
                    "base_config": [
                        {
                            _SBE.STEPID: "xtb",
                            _SBE.STEP_TYPE: "xtb",
                            _SBE.EXEC: {
                                "binary_location": "/projects/cc/mai/binaries/xtb-6.4.0",
                                "parallelization": {"cores": 4},
                            },
                            _SBE.SETTINGS: {
                                _SBE.SETTINGS_ARGUMENTS: {
                                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                        "--gbsa": "h2o"
                                    },
                                }
                            },
                            _SBE.INPUT: {
                                _SBE.INPUT_COMPOUNDS: [
                                    {
                                        _SBE.INPUT_SOURCE: "compound_filter",
                                        _SBE.INPUT_SOURCE_TYPE: "step",
                                    }
                                ]
                            },
                            _SBE.WRITEOUT: [
                                {
                                    _SBE.INPUT_COMPOUNDS: {
                                        _SBE.WRITEOUT_COMP_CATEGORY: "conformers"
                                    },
                                    _SBE.WRITEOUT_DESTINATION: {
                                        _SBE.WRITEOUT_DESTINATION_RESOURCE: "tests/junk/rmsd_rescoring_xtb.sdf",
                                        _SBE.STEP_TYPE: "file",
                                        _SBE.WRITEOUT_DESTINATION_FORMAT: "SDF",
                                    },
                                }
                            ],
                        },
                        {
                            _SBE.STEPID: "data_manipulation",
                            _SBE.STEP_TYPE: "data_manipulation",
                            _SBE.SETTINGS: {
                                _SBE.SETTINGS_ADDITIONAL: {
                                    "action": "attach_conformers_as_extra",
                                    _SBE.INPUT_SOURCE: "xtb",
                                }
                            },
                            _SBE.INPUT: {
                                _SBE.INPUT_COMPOUNDS: [
                                    {
                                        _SBE.INPUT_SOURCE: "compound_filter",
                                        _SBE.INPUT_SOURCE_TYPE: "step",
                                    }
                                ]
                            },
                        },
                        {
                            _SBE.STEPID: "rmsd",
                            _SBE.STEP_TYPE: "rmsd",
                            _SBE.SETTINGS: {
                                _SBE.SETTINGS_ADDITIONAL: {"method": "alignmol"}
                            },
                            _SBE.INPUT: {
                                _SBE.INPUT_COMPOUNDS: [
                                    {
                                        _SBE.INPUT_SOURCE: "data_manipulation",
                                        _SBE.INPUT_SOURCE_TYPE: "step",
                                    }
                                ]
                            },
                            _SBE.WRITEOUT: [
                                {
                                    _SBE.INPUT_COMPOUNDS: {
                                        _SBE.WRITEOUT_COMP_CATEGORY: "conformers"
                                    },
                                    _SBE.WRITEOUT_DESTINATION: {
                                        _SBE.WRITEOUT_DESTINATION_RESOURCE: "tests/junk/integration/rmsd_rescoring.sdf",
                                        _SBE.STEP_TYPE: "file",
                                        _SBE.WRITEOUT_DESTINATION_FORMAT: "SDF",
                                    },
                                },
                                {
                                    _SBE.INPUT_COMPOUNDS: {
                                        _SBE.WRITEOUT_COMP_CATEGORY: "conformers",
                                        _SBE.WRITEOUT_COMP_SELECTED_TAGS: [
                                            "docking_score",
                                            "rmsd",
                                            "grid_id",
                                        ],
                                        _SBE.WRITEOUT_COMP_AGGREGATION: {
                                            _SBE.WRITEOUT_COMP_AGGREGATION_MODE: "best_per_compound",
                                            _WE.ENVIRONMENT_EXPORT_KEY: "docking_score",
                                        },
                                    },
                                    _SBE.WRITEOUT_DESTINATION: {
                                        _SBE.WRITEOUT_DESTINATION_RESOURCE: "tests/junk/integration/rmsd_rescoring.csv",
                                        _SBE.STEP_TYPE: "file",
                                        _SBE.WRITEOUT_DESTINATION_FORMAT: "CSV",
                                    },
                                },
                            ],
                        },
                    ],
                    "iter_settings": {
                        _SBE.SETTINGS: {
                            "xtb": {
                                _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                    "--opt": [
                                        "vtight",
                                        "vtight",
                                        "vtight",
                                        "vtight",
                                        "vtight",
                                        "vtight",
                                        "vtight",
                                        "tight",
                                    ]
                                }
                            }
                        },
                        "n_iters": 8,
                        "iter_mode": "single",
                        "parallelizer_settings": {
                            "parallelize": True,
                            "cores": 4,
                            "max_length_sublists": 3,
                        },
                    },
                },
            ],
        }

        wflow = WorkFlow(**conf)
        wflow.initialize()

        self.assertEqual(len(wflow.steps), 5)

        wflow.execute()

        out_path = os.path.join(self._test_dir, "run_0/rmsd_rescoring.csv")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 100)

        out_path = os.path.join(self._test_dir, "run_7/rmsd_rescoring.sdf")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 3893)
