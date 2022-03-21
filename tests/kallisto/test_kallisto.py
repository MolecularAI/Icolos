import os
import unittest

from icolos.core.workflow_steps.calculation.kallisto import StepKallisto
from icolos.utils.enums.program_parameters import KallistoEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_mol_as_Compound,
    get_mol_as_Conformer,
    MAIN_CONFIG,
)

from icolos.utils.enums.step_enums import StepBaseEnum, StepKallistoEnum
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SKE = StepKallistoEnum()
_KE = KallistoEnum()


class Test_Kallisto(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/kallisto")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._paracetamol_molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_PATH)
        )

    def test_kallisto_run(self):
        step_conf = {
            _SBE.STEPID: "01_kallisto",
            _SBE.STEP_TYPE: _SBE.STEP_KALLISTO,
            _SBE.EXEC: {_SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["KALLISTO_LOCATION"]},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {},
                _SBE.SETTINGS_ADDITIONAL: {_SKE.FEATURES: [_KE.ALP, _KE.CNS, _KE.VDW]},
            },
        }
        kallisto_step = StepKallisto(**step_conf)
        kallisto_step.data.compounds = [self._paracetamol_molecule]
        confs = get_mol_as_Conformer(
            attach_root_path(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        )
        kallisto_step.data.compounds[0][0].add_conformers(confs, auto_update=True)
        self.assertListEqual(
            list(
                kallisto_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        kallisto_step.execute()
        self.assertListEqual(
            list(
                kallisto_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [5.3347, 12.9328, 24.6745],
        )
        self.assertEqual(
            kallisto_step.get_compounds()[0][0][0].get_molecule().GetProp(_KE.VDW),
            "3.4378599199634587|3.4396538784100725|3.4382871194521347|3.4390089115175315|3.397830365486183|3.3935108635482236|3.360611990090132|3.2709073572038108|3.419724012319036|3.279653714976599|3.272947644270613",
        )
        self.assertEqual(
            kallisto_step.get_compounds()[0][0][1].get_molecule().GetProp(_KE.VDW),
            "3.4396309533823914|3.437787905571396|3.439055955216591|3.4382857970612846|3.3978687711274818|3.3934963890468453|3.360556672297951|3.2709909790501612|3.419779086892037|3.279703488664608|3.2729758490815994",
        )

        # check SDF write-out (including energy-as-tag annotation)
        out_path = os.path.join(self._test_dir, "kallisto_paracetamol.sdf")
        kallisto_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 23491)
