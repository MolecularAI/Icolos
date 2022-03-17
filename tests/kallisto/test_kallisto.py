import os
import unittest

from icolos.core.workflow_steps.calculation.kallisto import StepKallisto
from icolos.utils.enums.program_parameters import KallistoEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Compound, get_mol_as_Conformer, MAIN_CONFIG

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
            _SBE.EXEC: {
                _SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["KALLISTO_LOCATION"]
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {

                },
                _SBE.SETTINGS_ADDITIONAL: {
                }
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
