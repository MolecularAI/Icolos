import unittest
import os
from rdkit import Chem

from icolos.core.composite_agents.workflow import WorkFlow

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.program_parameters import OMEGAEnum
from icolos.utils.enums.program_parameters import XTBEnum
from icolos.utils.enums.program_parameters import CrestEnum
from icolos.utils.enums.program_parameters import TurbomoleEnum
from icolos.utils.enums.program_parameters import PantherEnum
from icolos.core.steps_utils import initialize_step_from_dict

from tests.tests_paths import PATHS_EXAMPLEDATA, MAIN_CONFIG
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_WE = WorkflowEnum()
_OE = OMEGAEnum()
_XE = XTBEnum()
_CE = CrestEnum()
_TE = TurbomoleEnum()
_PE = PantherEnum()


class Test_workflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/workflow")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        _paracetamol_path = PATHS_EXAMPLEDATA.PARACETAMOL_PATH
        mol_supplier = Chem.SDMolSupplier(_paracetamol_path, removeHs=False)
        for mol in mol_supplier:
            self._molecule = mol

        # TODO: move header variables to MAIN_CONFIG
        self._HEADER_EXPORT = {
            _WE.ENVIRONMENT_EXPORT: [
                {
                    _WE.ENVIRONMENT_EXPORT_KEY: "OE_LICENSE",
                    _WE.ENVIRONMENT_EXPORT_VALUE: "/opt/scp/software/oelicense/1.0/oe_license.seq1",
                },
                {
                    _WE.ENVIRONMENT_EXPORT_KEY: "XTBHOME",
                    _WE.ENVIRONMENT_EXPORT_VALUE: "/opt/scp/services/reinvent/Icolos/binaries/xtb-6.3.2",
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
            ]
        }

    @classmethod
    def tearDownClass(cls):
        pass

    def test_workflow_initialization(self):
        conf = {
            _WE.HEADER: {_WE.ID: "test_workflow", _WE.ENVIRONMENT: self._HEADER_EXPORT},
            _WE.STEPS: [
                {
                    _SBE.STEPID: "crest_confgen",
                    _SBE.STEP_TYPE: _SBE.STEP_CREST,
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: None,
                        _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["CREST_BINARY_LOCATION"],
                        _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
                    },
                },
                {
                    _SBE.STEPID: "omega_confgen",
                    _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
                    _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
                },
            ],
        }
        wflow = WorkFlow(**conf)
        wflow.initialize()
        self.assertEqual(len(wflow.steps), 2)
        wflow.add_step(
            initialize_step_from_dict(
                {
                    _SBE.STEPID: "omega_confgen2",
                    _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
                    _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
                }
            )
        )
        self.assertEqual(len(wflow.steps), 2)
        self.assertEqual(len(wflow.get_steps()), 3)

    def test_workflow_with_global_variables(self):
        out_path = os.path.join(self._test_dir, "global_variables_out.sdf")
        conf = {
            _WE.HEADER: {
                _WE.ID: "test_workflow",
                _WE.DESCRIPTION: "this is a test description",
                _WE.ENVIRONMENT: self._HEADER_EXPORT,
                _WE.GLOBAL_VARIABLES: {
                    "root_dir": attach_root_path(""),
                },
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "01_initialization",
                    _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PARACETAMOL_PATH,
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                                _SBE.INPUT_FORMAT: _SBE.FORMAT_SDF,
                            }
                        ]
                    },
                }
            ],
        }
        wflow = WorkFlow(**conf)
        wflow.initialize()
        wflow.execute()

    def test_workflow_execution(self):
        conf = {
            _WE.HEADER: {
                _WE.ID: "test_workflow",
                _WE.DESCRIPTION: "this is a test description",
                _WE.ENVIRONMENT: self._HEADER_EXPORT,
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "01a_initialization",
                    _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PARACETAMOL_PATH,
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                                _SBE.INPUT_FORMAT: _SBE.FORMAT_SDF,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "01b_initialization",
                    _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.ASPIRIN_PATH,
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                                _SBE.INPUT_FORMAT: _SBE.FORMAT_SDF,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "02_omega_confgen",
                    _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
                    _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                _OE.CLASSIC_MAXCONFS: 10,
                                _OE.CLASSIC_RMS: 0.0,
                            },
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "01a_initialization",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            },
                            {
                                _SBE.INPUT_SOURCE: "01b_initialization",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            },
                        ],
                        _SBE.INPUT_MERGE: {
                            _SBE.INPUT_MERGE_COMPOUNDS: True,
                            _SBE.INPUT_MERGE_COMPOUNDS_BY: "id",
                            _SBE.INPUT_MERGE_ENUMERATIONS: True,
                            _SBE.INPUT_MERGE_ENUMERATIONS_BY: "id",
                        },
                    },
                },
                {
                    _SBE.STEPID: "02_conf_gen_crest",
                    _SBE.STEP_TYPE: _SBE.STEP_CREST,
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: None,
                        _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["CREST_BINARY_LOCATION"],
                        _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: ["-niceprint"],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                _CE.CREST_OPT: "normal",
                                _CE.CREST_G: "h2o",
                                _CE.CREST_RTHR: 0.5,
                                _CE.CREST_ETHR: 0.25,
                                _CE.CREST_EWIN: 8.0,
                                _CE.CREST_PTHR: 0.4,
                                _CE.CREST_BTHR: 0.02,
                            },
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "01a_initialization",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "01_conf_genXTB",
                    _SBE.STEP_TYPE: _SBE.STEP_XTB,
                    _SBE.EXEC: {
                        _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["XTBHOME"],
                        _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                _XE.XTB_OPT: "vtight",
                                _XE.XTB_GBSA: "h2o",
                            },
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "02_omega_confgen",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            },
                            {
                                _SBE.INPUT_SOURCE: "02_conf_gen_crest",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            },
                        ],
                        _SBE.INPUT_MERGE: {
                            _SBE.INPUT_MERGE_COMPOUNDS: True,
                            _SBE.INPUT_MERGE_COMPOUNDS_BY: "id",
                            _SBE.INPUT_MERGE_ENUMERATIONS: True,
                            _SBE.INPUT_MERGE_ENUMERATIONS_BY: "id",
                        },
                    },
                },
            ],
        }
        wflow = WorkFlow(**conf)
        wflow.initialize()
        wflow.execute()

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "02a_omega_confgen.sdf")
        wflow[2].write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreaterEqual(stat_inf.st_size, 1800)
        out_path = os.path.join(self._test_dir, "02b_crest_confgen.sdf")
        wflow[3].write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreaterEqual(stat_inf.st_size, 47156)
        out_path = os.path.join(self._test_dir, "03_XTB_from_omega_and_crest.sdf")
        wflow[4].write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreaterEqual(stat_inf.st_size, 52807)

    def test_ePSA_workflow_execution(self):
        conf = {
            _WE.HEADER: {
                _WE.ID: "test_workflow",
                _WE.DESCRIPTION: "this is a test description",
                _WE.ENVIRONMENT: self._HEADER_EXPORT,
            },
            _WE.STEPS: [
                {
                    _SBE.STEPID: "01_initialization_paracetamol",
                    _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.PARACETAMOL_PATH,
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                                _SBE.INPUT_FORMAT: _SBE.FORMAT_SDF,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "01_initialization_aspirin",
                    _SBE.STEP_TYPE: _SBE.STEP_INITIALIZATION,
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: PATHS_EXAMPLEDATA.ASPIRIN_PATH,
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_FILE,
                                _SBE.INPUT_FORMAT: _SBE.FORMAT_SDF,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "02_omega_confgen",
                    _SBE.STEP_TYPE: _SBE.STEP_OMEGA,
                    _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load omega"},
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                _OE.CLASSIC_MAXCONFS: 200,
                                _OE.CLASSIC_RMS: 0.0,
                                _OE.CLASSIC_CANON_ORDER: "false",
                            },
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "01_initialization_paracetamol",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            },
                            {
                                _SBE.INPUT_SOURCE: "01_initialization_aspirin",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            },
                        ]
                    },
                },
                {
                    _SBE.STEPID: "03_conf_optXTB",
                    _SBE.STEP_TYPE: _SBE.STEP_XTB,
                    _SBE.EXEC: {
                        _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["XTBHOME"],
                        _SBE.EXEC_PARALLELIZATION: {_SBE.EXEC_PARALLELIZATION_CORES: 7},
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                                _XE.XTB_OPT: "vtight",
                                _XE.XTB_GBSA: "h2o",
                            },
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "02_omega_confgen",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "04_turbomole",
                    _SBE.STEP_TYPE: _SBE.STEP_TURBOMOLE,
                    _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load turbomole/73"},
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        },
                        _SBE.SETTINGS_ADDITIONAL: {
                            _TE.TM_CONFIG_DIR: "/projects/cc/mai/material/Icolos/turbomole_config",
                            _TE.TM_CONFIG_BASENAME: "b97-3c-ri-d3-def2-mtzvp-int-nosym-charge",
                            _TE.TM_CONFIG_COSMO: "/projects/cc/mai/material/Icolos/turbomole_config/cosmoprep_eps80.tm",
                        },
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "03_conf_optXTB",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            }
                        ]
                    },
                },
                {
                    _SBE.STEPID: "05_cosmo",
                    _SBE.STEP_TYPE: _SBE.STEP_COSMO,
                    _SBE.EXEC: {
                        _SBE.EXEC_PREFIXEXECUTION: "module load COSMOtherm/19.0.4"
                    },
                    _SBE.SETTINGS: {
                        _SBE.SETTINGS_ARGUMENTS: {
                            _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                            _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                        }
                    },
                    _SBE.INPUT: {
                        _SBE.INPUT_COMPOUNDS: [
                            {
                                _SBE.INPUT_SOURCE: "04_turbomole",
                                _SBE.INPUT_SOURCE_TYPE: _SBE.INPUT_SOURCE_TYPE_STEP,
                            }
                        ]
                    },
                    _SBE.WRITEOUT: [
                        {
                            _SBE.WRITEOUT_COMP: {
                                _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                            },
                            _SBE.WRITEOUT_DESTINATION: {
                                _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                                _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_SDF,
                                _SBE.WRITEOUT_DESTINATION_RESOURCE: os.path.join(
                                    self._test_dir, "05_cosmo_ePSA_workflow.sdf"
                                ),
                            },
                        }
                    ],
                },
            ],
        }
        wflow = WorkFlow(**conf)
        wflow.initialize()
        wflow.execute()

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "02_omega_confgen.sdf")
        wflow.find_step_by_step_id("02_omega_confgen").write_conformers(out_path)
        self.assertGreater(8200, os.stat(out_path).st_size)
        out_path = os.path.join(self._test_dir, "03_conf_optXTB.sdf")
        wflow.find_step_by_step_id("03_conf_optXTB").write_conformers(out_path)
        self.assertGreater(8200, os.stat(out_path).st_size)
        out_path = os.path.join(self._test_dir, "04_turbomole.sdf")
        wflow.find_step_by_step_id("04_turbomole").write_conformers(out_path)
        self.assertGreater(82008, os.stat(out_path).st_size)
        out_path = os.path.join(self._test_dir, "05_cosmo_ePSA_workflow.sdf")
        self.assertGreater(12500, os.stat(out_path).st_size)
