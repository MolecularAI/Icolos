import os
import unittest

from icolos.core.workflow_steps.calculation.jazzy import StepJazzy
from icolos.utils.enums.program_parameters import JazzyEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_mol_as_Compound,
    get_mol_as_Conformer,
    MAIN_CONFIG,
)

from icolos.utils.enums.step_enums import StepBaseEnum, StepJazzyEnum
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SJE = StepJazzyEnum()
_JE = JazzyEnum()


class Test_Jazzy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/jazzy")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._paracetamol_molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_PATH)
        )

    def test_jazzy_run(self):
        step_conf = {
            _SBE.STEPID: "01_jazzy",
            _SBE.STEP_TYPE: _SBE.STEP_JAZZY,
            _SBE.EXEC: {_SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["JAZZY_LOCATION"]},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {},
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }
        jazzy_step = StepJazzy(**step_conf)
        jazzy_step.data.compounds = [self._paracetamol_molecule]
        confs = get_mol_as_Conformer(
            attach_root_path(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        )
        jazzy_step.data.compounds[0][0].add_conformers(confs, auto_update=True)
        self.assertListEqual(
            list(
                jazzy_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        jazzy_step.execute()
        self.assertListEqual(
            list(
                jazzy_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        self.assertEqual(
            jazzy_step.get_compounds()[0][0][0].get_molecule().GetProp(_JE.RESULT_SDX),
            "1.4151"
        )

        # check SDF write-out
        out_path = os.path.join(self._test_dir, "jazzy_paracetamol.sdf")
        jazzy_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 18225)
