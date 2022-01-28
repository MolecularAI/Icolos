import unittest
import os

from icolos.core.workflow_steps.schrodinger.macromodel import StepMacromodel

from icolos.utils.enums.step_enums import StepBaseEnum, TokenGuardEnum
from icolos.utils.enums.program_parameters import MacromodelEnum

from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Compound
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_CE = MacromodelEnum()
_TE = TokenGuardEnum()


class Test_Macromodel_confgen(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/MacroModel")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._paracetamol_molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_PATH)
        )
        self._aspirin_molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.ASPIRIN_PATH)
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_MacroModel_run(self):
        step_conf = {
            _SBE.STEPID: "01_macromodel",
            _SBE.STEP_TYPE: _SBE.STEP_MACROMODEL,
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2020-4"},
            _TE.TG: {
                _TE.TG_PREFIX_EXECUTION: "module load schrodinger/2020-4",
                _TE.TG_TOKEN_POOLS: {"MMOD_MACROMODEL": 2},
                _TE.TG_WAIT_INTERVAL_SECONDS: 30,
                _TE.TG_WAIT_LIMIT_SECONDS: 900,
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [_CE.MACROMODEL_WAIT],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {_CE.MACROMODEL_NJOBS: 2},
                }
            },
        }

        mm_step = StepMacromodel(**step_conf)
        mm_step.data.compounds = [self._paracetamol_molecule]

        # conformer coordinates should not be touched by the execution
        self.assertListEqual(
            list(
                mm_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-3.8276, -1.0625, 0.3279],
        )
        mm_step.execute()
        self.assertListEqual(
            list(
                mm_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-3.8276, -1.0625, 0.3279],
        )
        self.assertEqual(len(mm_step.get_compounds()[0][0].get_conformers()), 10)
        self.assertEqual(
            list(
                mm_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [-4.2269, -0.441, 0.2359],
        )

        # check write-out
        out_path = os.path.join(self._test_dir, "macromodel_output_file.sdf")
        mm_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 25637)
