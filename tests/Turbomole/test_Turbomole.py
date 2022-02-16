import unittest
import os

from icolos.core.workflow_steps.calculation.turbomole import StepTurbomole

from icolos.utils.enums.step_enums import StepBaseEnum, StepTurbomoleEnum
from icolos.utils.enums.program_parameters import TurbomoleEnum

from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
    get_mol_as_Compound,
    get_mol_as_Conformer,
    MAIN_CONFIG,
)
from icolos.utils.enums.compound_enums import ConformerContainerEnum
from icolos.utils.general.files_paths import attach_root_path
import time

_SBE = StepBaseEnum
_TE = TurbomoleEnum()
_COE = ConformerContainerEnum()
_STE = StepTurbomoleEnum()


class Test_Turbomole(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/Turbomole")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        # initialize a Compound with 1 Enumeration and 2 Conformers (done by OMEGA)
        _paracetamol_molecule = get_mol_as_Compound(PATHS_EXAMPLEDATA.PARACETAMOL_PATH)
        confs = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        _paracetamol_molecule[0].add_conformers(confs, auto_update=True)
        self._paracetamol_molecule = _paracetamol_molecule

    @classmethod
    def tearDownClass(cls):
        pass

    def test_Turbomole_run_ridft_single_core(self):
        step_conf = {
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
                    _TE.TM_CONFIG_BASENAME: "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
                    _TE.TM_CONFIG_COSMO: os.path.join(
                        MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
                    ),
                    _STE.EXECUTION_MODE: _TE.TM_RIDFT,
                },
            },
        }

        os.environ["PARA_ARCH"] = "MPI"
        os.environ["PARNODES"] = "4"
        tm_step = StepTurbomole(**step_conf)
        tm_step.data.compounds = [self._paracetamol_molecule]

        # conformer coordinates should not be touched by the execution
        self.assertListEqual(
            list(
                tm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [0.8785, 0.6004, -0.2173],
        )
        tm_step.execute()
        self.assertListEqual(
            list(
                tm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [0.8785, 0.6004, -0.2173],
        )
        cosmofile = tm_step.get_compounds()[0][0][0].get_extra_data()[
            _COE.EXTRA_DATA_COSMOFILE
        ]
        coordfile = tm_step.get_compounds()[0][0][0].get_extra_data()[
            _COE.EXTRA_DATA_COORDFILE
        ]
        self.assertTrue("basgrd points=   9806" in cosmofile[5])

        # check write-out
        out_path = os.path.join(self._test_dir, "paracetamol_conf1_CosmoFile")
        with open(out_path, "w") as f:
            f.writelines(cosmofile)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 132018)

        out_path = os.path.join(self._test_dir, "paracetamol_conf1_CoordFile")
        with open(out_path, "w") as f:
            f.writelines(coordfile)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 13544)

    def test_Turbomole_run_ridft_dual_core(self):
        step_conf = {
            _SBE.STEPID: "01_turbomole",
            _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 2},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _TE.TM_CONFIG_DIR: MAIN_CONFIG["TURBOMOLE_CONFIG"],
                    _TE.TM_CONFIG_BASENAME: "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
                    _TE.TM_CONFIG_COSMO: os.path.join(
                        MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
                    ),
                    _STE.EXECUTION_MODE: _TE.TM_RIDFT,
                },
            },
        }
        os.environ["PARA_ARCH"] = "MPI"
        os.environ["PARNODES"] = "4"

        tm_step = StepTurbomole(**step_conf)
        tm_step.data.compounds = [self._paracetamol_molecule]

        # conformer coordinates should not be touched by the execution
        self.assertListEqual(
            list(
                tm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        t1 = time.time()
        tm_step.execute()
        t2 = time.time()

        self.assertLess(t2 - t1, 50)
        self.assertListEqual(
            list(
                tm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [0.8785, 0.6004, -0.2173],
        )
        cosmofile = tm_step.get_compounds()[0][0][0].get_extra_data()[
            _COE.EXTRA_DATA_COSMOFILE
        ]
        coordfile = tm_step.get_compounds()[0][0][0].get_extra_data()[
            _COE.EXTRA_DATA_COORDFILE
        ]
        self.assertTrue("basgrd points=   9806" in cosmofile[5])

        # check write-out
        out_path = os.path.join(self._test_dir, "paracetamol_conf1_CosmoFile")
        with open(out_path, "w") as f:
            f.writelines(cosmofile)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 132018)

        out_path = os.path.join(self._test_dir, "paracetamole_conf1_CoordFile")
        with open(out_path, "w") as f:
            f.writelines(coordfile)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 13544)

    def test_Turbomole_run_jobex(self):
        step_conf = {
            _SBE.STEPID: "01_turbomole",
            _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73",
                _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 2},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-ri"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _TE.TM_JOBEX_C: 70,
                        _TE.TM_JOBEX_GCART: 3,
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _TE.TM_CONFIG_DIR: MAIN_CONFIG["TURBOMOLE_CONFIG"],
                    _TE.TM_CONFIG_BASENAME: "b97-3c-ri-d3-def2-mtzvp-int-charge",
                    _TE.TM_CONFIG_COSMO: os.path.join(
                        MAIN_CONFIG["TURBOMOLE_CONFIG"], "cosmoprep_eps80.tm"
                    ),
                    _STE.EXECUTION_MODE: _TE.TM_JOBEX,
                },
            },
        }

        os.environ["PARA_ARCH"] = "MPI"
        os.environ["PARNODES"] = "3"
        tm_step = StepTurbomole(**step_conf)
        tm_step.data.compounds = [self._paracetamol_molecule]

        # conformer coordinates should be touched by the execution (this is geo opt)
        self.assertListEqual(
            list(
                tm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        tm_step.execute()
        self.assertListEqual(
            list(
                tm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-0.7887, -0.0618, 0.1129],
        )
        cosmofile = tm_step.get_compounds()[0][0][0].get_extra_data()[
            _COE.EXTRA_DATA_COSMOFILE
        ]

        self.assertTrue("nspa=   92" in cosmofile[5])

        # check write-out
        out_path = os.path.join(self._test_dir, "paracetamol_conf1_CosmoFile_jobex")
        with open(out_path, "w") as f:
            f.writelines(cosmofile)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 115864)
