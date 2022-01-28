import unittest

from icolos.core.workflow_steps.schrodinger.ligprep import StepLigprep
from icolos.utils.enums.step_enums import StepBaseEnum, TokenGuardEnum, StepLigprepEnum
from icolos.utils.enums.program_parameters import LigprepEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_mol_as_Compound,
    get_test_Compounds_without_molecules,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_LSBE = StepLigprepEnum()
_LIE = LigprepEnum()
_TE = TokenGuardEnum()


class Test_Ligprep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        self._paracetamol_molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_PATH), compound_number=0
        )
        self._aspirin_molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.ASPIRIN_PATH), compound_number=1
        )
        self._Aspirin = get_test_Compounds_without_molecules(compound_numbers=[2])[
            "Aspirin"
        ]

    @classmethod
    def tearDownClass(cls):
        pass

    def test_LigPrep_run(self):
        step_conf = {
            _SBE.STEPID: "01_ligprep",
            _SBE.STEP_TYPE: _SBE.STEP_LIGPREP,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2020-4",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 2,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 2,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 2},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _LIE.LIGPREP_F: "/a/path/to/be/ignored/filter.txt"
                    }
                }
            },
        }

        ligprep_step = StepLigprep(**step_conf)
        ligprep_step.data.compounds = [
            self._paracetamol_molecule,
            self._aspirin_molecule,
            self._Aspirin,
        ]

        ligprep_step.execute()
        self.assertEqual(
            ["0:0", "1:0", "2:0"],
            [
                enum.get_index_string()
                for comp in ligprep_step.get_compounds()
                for enum in comp
            ],
        )
        self.assertEqual(
            [
                "[H]c1c([H])c(Cl)c2c(=O)nc(N([H])C(=O)C([H])([H])[H])sc2c1[H]",
                "[H]c1c([H])c(Cl)c2c(=O)nc(C(=O)[O-])sc2c1[H]",
                "O=C(C)Oc1ccccc1C(=O)O",
            ],
            [
                enum.get_original_smile()
                for comp in ligprep_step.get_compounds()
                for enum in comp
            ],
        )
        self.assertEqual(
            [
                "[H]c1c([H])c(Cl)c2c(=O)nc(N([H])C(=O)C([H])([H])[H])sc2c1[H]",
                "[H]OC(=O)c1nc(=O)c2c(Cl)c([H])c([H])c([H])c2s1",
                "[H]OC(=O)c1c([H])c([H])c([H])c([H])c1OC(=O)C([H])([H])[H]",
            ],
            [
                enum.get_smile()
                for comp in ligprep_step.get_compounds()
                for enum in comp
            ],
        )
        self.assertListEqual(
            list(
                ligprep_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-4.9037, 3.0725, 2.0034],
        )
        self.assertListEqual(
            list(
                ligprep_step.get_compounds()[1][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-4.8794, 3.0688, -2.0104],
        )
        self.assertListEqual(
            list(
                ligprep_step.get_compounds()[2][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [0.0243, 2.4719, -0.3164],
        )

    def test_LigPrep_run_EPIK_stereo_filtering(self):
        step_conf = {
            _SBE.STEPID: "01_ligprep",
            _SBE.STEP_TYPE: _SBE.STEP_LIGPREP,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2020-4",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 2,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 2,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 2},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [_LIE.LIGPREP_EPIK],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _LIE.LIGPREP_PH: 7.0,
                        _LIE.LIGPREP_PHT: 2.0,
                        _LIE.LIGPREP_S: 10,
                        _LIE.LIGPREP_BFF: 14,
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {_LSBE.FILTER_FILE: {"Total_charge": "!= 0"}},
            },
        }

        ligprep_step = StepLigprep(**step_conf)
        ligprep_step.data.compounds = [
            self._paracetamol_molecule,
            self._aspirin_molecule,
            self._Aspirin,
        ]

        ligprep_step.execute()
        self.assertEqual(
            ["0:0", "0:1", "1:0"],
            [
                enum.get_index_string()
                for comp in ligprep_step.get_compounds()
                for enum in comp
            ],
        )
        self.assertEqual(
            [
                "[H]c1c([H])c(Cl)c2c(=O)nc(N([H])C(=O)C([H])([H])[H])sc2c1[H]",
                "[H]c1c([H])c(Cl)c2c(=O)nc(N([H])C(=O)C([H])([H])[H])sc2c1[H]",
                "[H]c1c([H])c(Cl)c2c(=O)nc(C(=O)[O-])sc2c1[H]",
            ],
            [
                enum.get_original_smile()
                for comp in ligprep_step.get_compounds()
                for enum in comp
            ],
        )
        self.assertEqual(
            [
                "[H]c1c([H])c(Cl)c2c(=O)n([H])/c(=N\\C(=O)C([H])([H])[H])sc2c1[H]",
                "[H]c1c([H])c(Cl)c2c(=O)nc(N([H])C(=O)C([H])([H])[H])sc2c1[H]",
                "[H]c1c([H])c(Cl)c2c(=O)[n+]([H])c(C(=O)[O-])sc2c1[H]",
            ],
            [
                enum.get_smile()
                for comp in ligprep_step.get_compounds()
                for enum in comp
            ],
        )
        self.assertListEqual(
            list(
                ligprep_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-4.7828, 5.0389, -2.1622],
        )
        self.assertListEqual(
            list(
                ligprep_step.get_compounds()[0][1]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-4.9037, 3.0725, 2.0034],
        )
        self.assertListEqual(
            list(
                ligprep_step.get_compounds()[1][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-5.2155, 3.215, -1.1152],
        )
