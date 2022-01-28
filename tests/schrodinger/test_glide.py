import os
import time
import unittest

from icolos.core.workflow_steps.schrodinger.glide import StepGlide

from icolos.utils.enums.step_enums import StepBaseEnum, TokenGuardEnum, StepGlideEnum
from icolos.utils.enums.program_parameters import GlideEnum

from tests.tests_paths import (
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    get_1UYD_ligands_as_Compounds,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SGE = StepGlideEnum()
_EE = GlideEnum()
_TE = TokenGuardEnum()


class Test_Glide(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/Glide")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._1UYD_compounds = get_1UYD_ligands_as_Compounds(
            abs_path=PATHS_1UYD.LIGANDS
        )
        self.receptor_path = PATHS_1UYD.GRID_PATH
        self.receptor_constraints_path = PATHS_1UYD.GRID_CONSTRAINTS_PATH
        self.receptor_path_COX2 = PATHS_EXAMPLEDATA.PRIME_COX2_GRID

    @classmethod
    def tearDownClass(cls):
        pass

    def test_Glide_run(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 4,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 2,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
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
        }

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = self._1UYD_compounds

        glide_step.execute()
        self.assertEqual(len(glide_step.get_compounds()), 15)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 3)
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-1.5198, 11.3439, 24.0245],
        )

        self.assertListEqual(
            list(
                glide_step.get_compounds()[14][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-2.1655, 12.4809, 24.137],
        )
        self.assertEqual(
            glide_step.get_compounds()[0][0][0]
            .get_molecule()
            .GetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE),
            "-8.4349",
        )
        self.assertEqual(
            glide_step.get_compounds()[0][0][1]
            .get_molecule()
            .GetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE),
            "-7.83118",
        )
        self.assertEqual(
            glide_step.get_compounds()[0][0][2]
            .get_molecule()
            .GetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE),
            "-6.0089",
        )

        # check SDF write-out
        out_path = os.path.join(self._test_dir, "glide_docked.sdf")
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 209000)

    def test_Glide_run_parallelization_1core_singleton(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 1,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
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
        }

        compounds = self._1UYD_compounds[:3]

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = compounds

        # execute on one core and put all in one list
        time_difference = time.time()
        glide_step.execute()
        time_difference = time.time() - time_difference
        self.assertGreater(time_difference, 100)
        self.assertGreater(325, time_difference)
        self.assertEqual(len(glide_step.get_compounds()), 3)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 3)
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-1.5198, 11.3439, 24.0245],
        )

        # check SDF write-out
        out_path = os.path.join(
            self._test_dir, "glide_docked_single_core_singleton_list.sdf"
        )
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 50500)

    def test_Glide_run_1_core(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 1},
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
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
        }

        compounds = self._1UYD_compounds[:3]

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = compounds

        # execute on one core and put all in one list
        time_difference = time.time()
        glide_step.execute()
        time_difference = time.time() - time_difference
        self.assertGreater(325, time_difference)
        self.assertEqual(len(glide_step.get_compounds()), 3)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 3)
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-1.5198, 11.3439, 24.0245],
        )

        # check SDF write-out
        out_path = os.path.join(
            self._test_dir, "glide_docked_merged_list_3compounds.sdf"
        )
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 50500)

    def test_Glide_run_parallelization_4cores(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 4},
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
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
        }

        compounds = self._1UYD_compounds[:3]

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = compounds

        # execute and put all in one list
        time_difference = time.time()
        glide_step.execute()
        time_difference = time.time() - time_difference
        self.assertGreater(150, time_difference)
        self.assertEqual(len(glide_step.get_compounds()), 3)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 3)
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-1.5198, 11.3439, 24.0245],
        )

        # check SDF write-out
        out_path = os.path.join(
            self._test_dir, "glide_docked_parallelized_3compounds.sdf"
        )
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 50500)

    def test_Glide_run_parallelization_4cores_in_file_usage(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 4},
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.CONFIGURATION: {
                        _EE.GLIDE_AMIDE_MODE: "trans",
                        _EE.GLIDE_EXPANDED_SAMPLING: "True",
                        _EE.GLIDE_GRIDFILE: [self.receptor_constraints_path],
                        _EE.GLIDE_NENHANCED_SAMPLING: "1",
                        _EE.GLIDE_POSE_OUTTYPE: _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB,
                        _EE.GLIDE_POSES_PER_LIG: "3",
                        _EE.GLIDE_POSTDOCK_NPOSE: "25",
                        _EE.GLIDE_POSTDOCKSTRAIN: "True",
                        _EE.GLIDE_PRECISION: "SP",
                        _EE.GLIDE_REWARD_INTRA_HBONDS: "True",
                    },
                    _SGE.MAESTRO_IN_FILE: {
                        _SGE.MAESTRO_IN_FILE_PATH: PATHS_EXAMPLEDATA.GLIDE_EXAMPLE_IN
                    },
                },
            },
        }

        compounds = self._1UYD_compounds[:3]

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = compounds
        glide_step.execute()

        # execute on one core and put all in one list
        self.assertEqual(len(glide_step.get_compounds()), 3)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 3)

        # would be [-2.5618, 10.8202, 25.2644] without constraints
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [3.1229, 4.5141, 24.8603],
        )

        # check SDF write-out
        out_path = os.path.join(
            self._test_dir, "glide_docked_parallelized_3compounds.sdf"
        )
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 50500)

    def test_Glide_run_parallelization_4cores_ensemble_docking(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 4},
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.CONFIGURATION: {
                        _EE.GLIDE_AMIDE_MODE: "trans",
                        _EE.GLIDE_EXPANDED_SAMPLING: "True",
                        _EE.GLIDE_GRIDFILE: [
                            self.receptor_path_COX2,
                            self.receptor_path,
                        ],
                        _EE.GLIDE_NENHANCED_SAMPLING: "1",
                        _EE.GLIDE_POSE_OUTTYPE: _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB,
                        _EE.GLIDE_POSES_PER_LIG: "3",
                        _EE.GLIDE_POSTDOCK_NPOSE: "25",
                        _EE.GLIDE_POSTDOCKSTRAIN: "True",
                        _EE.GLIDE_PRECISION: "SP",
                        _EE.GLIDE_REWARD_INTRA_HBONDS: "True",
                    },
                    _SBE.GRID_IDS: ["mygrid1", "mygrid2"],
                },
            },
        }

        compounds = self._1UYD_compounds[:3]

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = compounds

        # execute on one core and put all in one list
        glide_step.execute()
        self.assertEqual(len(glide_step.get_compounds()), 3)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 6)
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-1.5198, 11.3439, 24.0245],
        )
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][5]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [7.3776, 55.7005, 70.3807],
        )
        self.assertListEqual(
            ["mygrid2", "mygrid2", "mygrid1", "mygrid2", "mygrid1", "mygrid1"],
            [
                comp.get_molecule().GetProp(_SBE.ANNOTATION_GRID_ID)
                for comp in list(glide_step.get_compounds()[0][0])
            ],
        )

        # check SDF write-out
        out_path = os.path.join(
            self._test_dir, "glide_docked_parallelized_ensemble_docking.sdf"
        )
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 80000)

    def test_Glide_run_parallelization_poseviewer(self):
        step_conf = {
            _SBE.STEPID: "01_glide",
            _SBE.STEP_TYPE: _SBE.STEP_GLIDE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 4},
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_EE.GLIDE_HOST: "cpu-only"},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.CONFIGURATION: {
                        _EE.GLIDE_AMIDE_MODE: "trans",
                        _EE.GLIDE_EXPANDED_SAMPLING: "True",
                        _EE.GLIDE_GRIDFILE: [self.receptor_path],
                        _EE.GLIDE_NENHANCED_SAMPLING: "1",
                        _EE.GLIDE_POSE_OUTTYPE: _EE.GLIDE_POSE_OUTTYPE_POSEVIEWER,
                        _EE.GLIDE_POSES_PER_LIG: "3",
                        _EE.GLIDE_POSTDOCK_NPOSE: "25",
                        _EE.GLIDE_POSTDOCKSTRAIN: "True",
                        _EE.GLIDE_PRECISION: "SP",
                        _EE.GLIDE_REWARD_INTRA_HBONDS: "True",
                    }
                },
            },
        }

        compounds = self._1UYD_compounds[:3]

        glide_step = StepGlide(**step_conf)
        glide_step.data.compounds = compounds

        # execute on one core and put all in one list
        glide_step.execute()

        self.assertEqual(len(glide_step.get_compounds()), 3)
        self.assertEqual(len(glide_step.get_compounds()[0][0].get_conformers()), 3)
        self.assertListEqual(
            list(
                glide_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-1.5198, 11.3439, 24.0245],
        )

        # check SDF write-out
        out_path = os.path.join(
            self._test_dir, "glide_docked_parallelized_3compounds_pv.sdf"
        )
        glide_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 50500)
