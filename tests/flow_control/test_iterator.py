import unittest
from icolos.core.flow_control.iterator import StepIterator
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import StepBaseEnum, StepTurbomoleEnum
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


class TestIterator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/iterator")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self) -> None:
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO, "r") as f:
            self.structure = f.read()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP, "r") as f:
            self.topol = f.read()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR, "rb") as f:
            self.tpr_file = f.read()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC, "rb") as f:
            self.xtc_file = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_POSRE, "r") as f:
            self.posre = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_LIG_ITP, "r") as f:
            self.lig_itp = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_LIG_POSRE, "r") as f:
            self.lig_posre = f.read()

    def test_single_initialization(self):

        full_conf = {
            "base_config": [
                {
                    _SBE.STEPID: "01_turbomole",
                    _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73",
                        _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 1},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            _TE.TM_CONFIG_DIR: MAIN_CONFIG["TURBOMOLE_CONFIG"],
                            _TE.TM_CONFIG_COSMO: os.path.join(
                                MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
                            ),
                            _STE.EXECUTION_MODE: _TE.TM_RIDFT,
                        },
                    },
                }
            ],
            "iter_settings": {
                "settings": {
                    "01_turbomole": {
                        _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                        _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        _SBE.SETTINGS_ADDITIONAL: {
                            _TE.TM_CONFIG_BASENAME: [
                                "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
                                "blyp-ri-d3-def2-svp-int-coarse-charge2",
                                "some_other_spicy_functional",
                            ]
                        },
                    }
                },
                "iter_mode": "single",
                "n_iters": 3,  # for now this is manual, should match the number of settings to iterate over
            },
        }

        step_iterator = StepIterator(**full_conf)
        self.assertEqual(len(step_iterator.initialized_steps), 3)
        for i in step_iterator.initialized_steps:
            assert isinstance(i, StepBase)

    def test_multi_iter_initialization(self):

        full_conf = {
            "base_config": [
                {
                    _SBE.STEPID: "01_turbomole",
                    _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73",
                        _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 1},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            _TE.TM_CONFIG_DIR: MAIN_CONFIG["TURBOMOLE_CONFIG"],
                            _TE.TM_CONFIG_COSMO: os.path.join(
                                MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
                            ),
                            _TE.TM_CONFIG_BASENAME: "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
                            _STE.EXECUTION_MODE: _TE.TM_RIDFT,
                        },
                    },
                }
            ],
            "iter_settings": {
                # no changes in config, just run the same step through multiple times
                "iter_mode": "n_iters",
                "n_iters": 5,  # for now this is manual, should match the number of settings to iterate over
            },
        }

        step_iterator = StepIterator(**full_conf)
        self.assertEqual(len(step_iterator.initialized_steps), 5)
        for i in step_iterator.initialized_steps:
            assert isinstance(i, StepBase)

    def test_single_initialization_parallel_execution(self):
        """
        Test running multiple steps in parallel
        """

        full_conf = {
            "base_config": [
                {
                    _SBE.STEPID: "test_mmgbsa",
                    _SBE.STEP_TYPE: _SBE.STEP_GMX_MMPBSA,
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2 && module load gmx_MMPBSA/1.3.3-fosscuda-2019a-Python-3.7.2"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            "make_ndx_command": "Protein Other",
                            "pipe_input": "Protein Other",
                        },
                    },
                }
            ],
            "iter_settings": {
                "n_iters": 4,  # for now this is manual, should match the number of settings to iterate over
                "parallelizer_settings": {
                    "parallelize": True,
                    "cores": 4,
                    "max_lenth_sublists": 1,
                },
            },
        }

        step_mmpbsa_job_control = StepIterator(**full_conf)
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="structure.gro", file_data=self.structure)
        # )
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="topol.top", file_data=self.topol)
        # )
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="structure.xtc", file_data=self.xtc_file)
        # )
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="structure.tpr", file_data=self.tpr_file)
        # )
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="posre.itp", file_data=self.posre)
        # )
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        # )
        # step_mmpbsa.data.generic.add_file(
        #     GenericData(file_name="posre_DMP:100.itp", file_data=self.lig_posre)
        # )

        # should return JobControl object
        assert isinstance(step_mmpbsa_job_control.initialized_workflows, StepBase)
        assert len(step_mmpbsa_job_control.initialized_workflows.workflows) == 4
        # TODO: there isn't really a good way to unit test this, it is a pain to load the data in to the individual steps
        # step_mmpbsa_job_control.initialized_steps.execute()
