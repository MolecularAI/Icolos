import unittest
import os
from tests.tests_paths import PATHS_EXAMPLEDATA, PATHS_1UYD
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

    def test_docking_workflow(self):

        conf = {
            _WE.HEADER: {
                _WE.ID: "NIBR",
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
                        _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-1-js-aws",
                        _SBE.EXEC_PARALLELIZATION: {
                            _SBE.EXEC_PARALLELIZATION_CORES: 4,
                            _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                        },
                        _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY: 3},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-epik"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                "-ph": 7.0,
                                "-pht": 1.0,
                                "-s": 1,
                                "-bff": 14,
                            },
                        },
                        _SBE.SETTINGS_ADDITIONAL: {},
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
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: f"{self._test_dir}/nibr_ligprep.sdf",
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
                            _SBE.EXEC_PARALLELIZATION_CORES: 8,
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
                                "GRIDFILE": [PATHS_1UYD.GRID_PATH],
                                "NENHANCED_SAMPLING": "1",
                                "POSE_OUTTYPE": "ligandlib_sd",
                                "POSES_PER_LIG": "15",
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
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: f"{self._test_dir}/tests/junk/nibr_glide.sdf",
                                _SBE.WRITEOUT_DESTINATION_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "SDF",
                            },
                        }
                    ],
                },
                {
                    _SBE.STEPID: "Shaep",
                    _SBE.STEP_TYPE: "shaep",
                    _SBE.EXEC: {_SBE.EXEC_BINARYLOCATION: "/projects/cc/mai/binaries"},
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PANTHER_NEGATIVE_IMAGE,
                                _SBE.INPUT_EXTENSION: "mol2",
                            }
                        ],
                        _SBE.INPUT_COMPOUNDS: [
                            {_SBE.INPUT_SOURCE: "Glide", _SBE.INPUT_SOURCE_TYPE: "step"}
                        ],
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.INPUT_COMPOUNDS: {
                                _SBE.WRITEOUT_COMP_CATEGORY: "conformers",
                                _SBE.WRITEOUT_COMP_SELECTED_TAGS: [
                                    "shape_similarity",
                                    "esp_similarity",
                                ],
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: f"{self._test_dir}/tests/junk/nibr_final_all.csv",
                                _SBE.STEP_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "CSV",
                            },
                        },
                        {
                            _SBE.INPUT_COMPOUNDS: {
                                _SBE.WRITEOUT_COMP_CATEGORY: "conformers",
                                _SBE.WRITEOUT_COMP_SELECTED_TAGS: [
                                    "shape_similarity",
                                    "esp_similarity",
                                ],
                                _SBE.WRITEOUT_COMP_AGGREGATION: {
                                    _SBE.WRITEOUT_COMP_AGGREGATION_MODE: "best_per_compound",
                                    _WE.ENVIRONMENT_EXPORT_KEY: "shape_similarity",
                                    _SBE.WRITEOUT_COMP_AGGREGATION_HIGHESTISBEST: True,
                                },
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: f"{self._test_dir}/nibr_final_bestpercompound.csv",
                                _SBE.STEP_TYPE: "file",
                                _SBE.WRITEOUT_DESTINATION_FORMAT: "CSV",
                            },
                        },
                    ],
                },
            ],
        }

        wflow = WorkFlow(**conf)
        wflow.initialize()

        self.assertEqual(len(wflow.steps), 4)

        wflow.execute()

        out_path = os.path.join(self._test_dir, "nibr_final_bestpercompound.csv")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 110)
