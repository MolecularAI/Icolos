import unittest
from icolos.core.workflow_steps.active_learning.active_learning_reinvent import (
    StepActiveLearning,
)
from icolos.utils.enums.program_parameters import GlideEnum
import os
from icolos.utils.enums.step_enums import (
    StepActiveLearningEnum,
    StepBaseEnum,
    StepGlideEnum,
)
from icolos.utils.general.files_paths import attach_root_path
from tests.tests_paths import PATHS_1UYD, PATHS_EXAMPLEDATA

_SBE = StepBaseEnum
_EE = GlideEnum()
_SGE = StepGlideEnum()
_SALE = StepActiveLearningEnum()


class TestActiveLearning(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/active_learning")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self.ligands = PATHS_1UYD.LIGANDS
        self.receptor_path = PATHS_1UYD.GRID_PATH
        self.receptor_constraints_path = PATHS_1UYD.GRID_CONSTRAINTS_PATH
        self.receptor_path_COX2 = PATHS_EXAMPLEDATA.PRIME_COX2_GRID

    @classmethod
    def tearDownClass(cls):
        pass

    # def test_active_learning_docking(self):
    #     step_conf = {
    #         _SBE.STEPID: "01_active_learning",
    #         _SBE.STEP_TYPE: _SBE.STEP_ACTIVE_LEARNING,
    #         _SBE.EXEC: {},
    #         _SBE.SETTINGS: {
    #             _SBE.SETTINGS_ADDITIONAL: {
    #                 _SALE.VIRTUAL_LIB: self.ligands,
    #                 _SALE.N_ROUNDS: "2",
    #                 _SALE.INIT_SAMPLES: "2",
    #                 _SALE.BATCH_SIZE: "4",
    #                 _SALE.CRITERIA: _SGE.GLIDE_DOCKING_SCORE,
    #                 # config for embedding + docking
    #                 _SALE.ORACLE_CONFIG: [
    #                     {
    #                         _SBE.STEPID: "01_glide",
    #                         _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
    #                         _SBE.EXEC: {
    #                             _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
    #                             _SBE.EXEC_PARALLELIZATION: {
    #                                 _SBE.EXEC_PARALLELIZATION_CORES: 8,
    #                                 _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
    #                             },
    #                             _SBE.EXEC_FAILUREPOLICY: {
    #                                 _SBE.EXEC_FAILUREPOLICY_NTRIES: 1
    #                             },
    #                         },
    #                         _SBE.SETTINGS: {
    #                             _SBE.SETTINGS_ARGUMENTS: {
    #                                 _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
    #                                 _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
    #                                     _EE.GLIDE_HOST: "cpu-only"
    #                                 },
    #                             },
    #                             _SBE.SETTINGS_ADDITIONAL: {
    #                                 _SGE.CONFIGURATION: {
    #                                     _EE.GLIDE_AMIDE_MODE: "trans",
    #                                     _EE.GLIDE_EXPANDED_SAMPLING: "True",
    #                                     _EE.GLIDE_GRIDFILE: [self.receptor_path],
    #                                     _EE.GLIDE_NENHANCED_SAMPLING: "1",
    #                                     _EE.GLIDE_POSE_OUTTYPE: _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB,
    #                                     _EE.GLIDE_POSES_PER_LIG: "3",
    #                                     _EE.GLIDE_POSTDOCK_NPOSE: "25",
    #                                     _EE.GLIDE_POSTDOCKSTRAIN: "True",
    #                                     _EE.GLIDE_PRECISION: "SP",
    #                                     _EE.GLIDE_REWARD_INTRA_HBONDS: "True",
    #                                 }
    #                             },
    #                         },
    #                     },
    #                 ],
    #             },
    #         },
    #     }

    #     step_active_learning = StepActiveLearning(**step_conf)
    #     step_active_learning.execute()
    #     out_path = os.path.join(self._test_dir, "production_model.pkl")
    #     data = step_active_learning.data.generic.get_files_by_extension(ext="pkl")[
    #         0
    #     ].get_data()
    #     with open(out_path, "wb") as f:
    #         f.write(data)
    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 348000)

    def test_active_learning_soap_gpr(self):
        step_conf = {
            _SBE.STEPID: "01_active_learning",
            _SBE.STEP_TYPE: _SBE.STEP_ACTIVE_LEARNING,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SALE.VIRTUAL_LIB: self.ligands,
                    _SALE.N_ROUNDS: "4",
                    _SALE.INIT_SAMPLES: "2",
                    _SALE.BATCH_SIZE: "4",
                    _SALE.CRITERIA: _SGE.GLIDE_DOCKING_SCORE,
                    _SALE.RUNNING_MODE: "soap_gpr",
                    _SALE.EVALUATE: False,
                    # config for embedding + docking
                    _SALE.ORACLE_CONFIG: [
                        {
                            _SBE.STEPID: "01_glide",
                            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
                            _SBE.EXEC: {
                                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                                _SBE.EXEC_PARALLELIZATION: {
                                    _SBE.EXEC_PARALLELIZATION_CORES: 8,
                                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                                },
                                _SBE.EXEC_FAILUREPOLICY: {
                                    _SBE.EXEC_FAILUREPOLICY_NTRIES: 1
                                },
                            },
                            _SBE.SETTINGS: {
                                _SBE.SETTINGS_ARGUMENTS: {
                                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                        _EE.GLIDE_HOST: "localhost"
                                    },
                                },
                                _SBE.SETTINGS_ADDITIONAL: {
                                    _SGE.CONFIGURATION: {
                                        _EE.GLIDE_AMIDE_MODE: "trans",
                                        _EE.GLIDE_EXPANDED_SAMPLING: "True",
                                        _EE.GLIDE_GRIDFILE: [self.receptor_path],
                                        _EE.GLIDE_NENHANCED_SAMPLING: "1",
                                        _EE.GLIDE_POSE_OUTTYPE: _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB,
                                        _EE.GLIDE_POSES_PER_LIG: "3",
                                        _EE.GLIDE_POSTDOCK_NPOSE: "25",
                                        _EE.GLIDE_POSTDOCKSTRAIN: "True",
                                        _EE.GLIDE_PRECISION: "SP",
                                        _EE.GLIDE_REWARD_INTRA_HBONDS: "True",
                                    }
                                },
                            },
                        },
                    ],
                },
            },
        }

        step_active_learning = StepActiveLearning(**step_conf)
        step_active_learning.execute()
        out_path = os.path.join(self._test_dir, "production_model.pkl")
        data = step_active_learning.data.generic.get_files_by_extension(ext="pkl")[
            0
        ].get_data()
        with open(out_path, "wb") as f:
            f.write(data)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 348000)
