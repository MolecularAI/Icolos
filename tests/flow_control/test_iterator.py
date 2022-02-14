import unittest
from icolos.core.flow_control.iterator import StepIterator
from icolos.core.step_dispatch.dispatcher import StepDispatcher
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import (
    StepBaseEnum,
    StepGromacsEnum,
    StepTurbomoleEnum,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.utils.enums.program_parameters import TurbomoleEnum
import os
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    MAIN_CONFIG,
)

_SBE = StepBaseEnum
_TE = TurbomoleEnum()
_STE = StepTurbomoleEnum()
_SGE = StepGromacsEnum()


class TestIterator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/iterator")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self) -> None:
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO, "r") as f:
            self.structure = f.read()

        # def test_single_initialization(self):

        #     full_conf = {
        #         "base_config": [
        #             {
        #                 _SBE.STEPID: "01_turbomole",
        #                 _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
        #                 _SBE.EXEC: {
        #                     _SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73",
        #                     _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 1},
        #                 },
        #                 _SBE.SETTINGS: {
        #                     _SBE.SETTINGS_ARGUMENTS: {
        #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
        #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
        #                     },
        #                     _SBE.SETTINGS_ADDITIONAL: {
        #                         _TE.TM_CONFIG_DIR: MAIN_CONFIG["TURBOMOLE_CONFIG"],
        #                         _TE.TM_CONFIG_COSMO: os.path.join(
        #                             MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
        #                         ),
        #                         _STE.EXECUTION_MODE: _TE.TM_RIDFT,
        #                     },
        #                 },
        #             }
        #         ],
        #         "iter_settings": {
        #             "settings": {
        #                 "01_turbomole": {
        #                     _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
        #                     _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
        #                     _SBE.SETTINGS_ADDITIONAL: {
        #                         _TE.TM_CONFIG_BASENAME: [
        #                             "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
        #                             "blyp-ri-d3-def2-svp-int-coarse-charge2",
        #                             "some_other_spicy_functional",
        #                         ]
        #                     },
        #                 }
        #             },
        #             "iter_mode": "single",
        #             "n_iters": 3,  # for now this is manual, should match the number of settings to iterate over
        #         },
        #     }

        #     step_iterator = StepIterator(**full_conf)
        #     self.assertEqual(len(step_iterator.initialized_steps), 3)
        #     for i in step_iterator.initialized_steps:
        #         assert isinstance(i, StepBase)

        # def test_multi_iter_initialization(self):

        #     full_conf = {
        #         "base_config": [
        #             {
        #                 _SBE.STEPID: "01_turbomole",
        #                 _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
        #                 _SBE.EXEC: {
        #                     _SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73",
        #                     _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 1},
        #                 },
        #                 _SBE.SETTINGS: {
        #                     _SBE.SETTINGS_ARGUMENTS: {
        #                         _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
        #                         _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
        #                     },
        #                     _SBE.SETTINGS_ADDITIONAL: {
        #                         _TE.TM_CONFIG_DIR: MAIN_CONFIG["TURBOMOLE_CONFIG"],
        #                         _TE.TM_CONFIG_COSMO: os.path.join(
        #                             MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
        #                         ),
        #                         _TE.TM_CONFIG_BASENAME: "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
        #                         _STE.EXECUTION_MODE: _TE.TM_RIDFT,
        #                     },
        #                 },
        #             }
        #         ],
        #         "iter_settings": {
        #             # no changes in config, just run the same step through multiple times
        #             "iter_mode": "n_iters",
        #             "n_iters": 5,  # for now this is manual, should match the number of settings to iterate over
        #         },
        #     }

        # step_iterator = StepIterator(**full_conf)
        # self.assertEqual(len(step_iterator.initialized_steps), 5)
        # for i in step_iterator.initialized_steps:
        #     assert isinstance(i, StepBase)

    def test_n_iters_4_cores_gromacs(self):
        """
        Initialize 4 pdb2gmx jobs in separate workflows,
        """

        full_conf = {
            "base_config": [
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
                    _SBE.INPUT: {_SBE.INPUT_GENERIC: []},
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
                                _SBE.INPUT_SOURCE: os.path.join(
                                    MAIN_CONFIG["ICOLOS_TEST_DATA"],
                                    "gromacs/protein/ions.mdp",
                                ),
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
                    _SBE.INPUT: {
                        _SBE.INPUT_GENERIC: [
                            {
                                _SBE.INPUT_SOURCE: "04_grompp",
                                _SBE.INPUT_EXTENSION: "tpr",
                            },
                        ]
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.INPUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "top"},
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: f"{self._test_dir}/topol_out.top",
                            },
                        }
                    ],
                },
            ],
            "iter_settings": {
                "n_iters": 4,  # for now this is manual, should match the number of settings to iterate over
                "parallelizer_settings": {
                    "cores": 4,
                },
            },
        }

        step_gromacs_preprocess = StepIterator(**full_conf)

        # should return JobControl object
        assert isinstance(step_gromacs_preprocess.dispatcher, StepDispatcher)
        assert len(step_gromacs_preprocess.dispatcher.workflows) == 4

        # instantiate the Dispatcher object
        step_gromacs_preprocess.dispatcher.execute()
        out_path = os.path.join(self._test_dir, "run_3", "topol_out_0.top")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 1560)
