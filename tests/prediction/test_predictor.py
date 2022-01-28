import unittest
import os

from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.prediction.predictor import StepPredictor

from icolos.utils.enums.step_enums import StepBaseEnum, StepPredictorEnum

from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Conformer
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SPE = StepPredictorEnum()


class Test_Predictor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/Prediction")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._example_model_path = attach_root_path(PATHS_EXAMPLEDATA.EPSA_MODEL_PATH)
        self._example_mol_path = attach_root_path(
            PATHS_EXAMPLEDATA.EPSA_EXAMPLE_MOLECULE
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_predict_ePSA_with_descriptors(self):
        step_conf = {
            _SBE.STEPID: "01_predict_ePSA",
            _SBE.STEP_TYPE: _SBE.STEP_PREDICTION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_PARAMETERS: {}},
                _SBE.SETTINGS_ADDITIONAL: {
                    _SPE.MODEL_PATH: self._example_model_path,
                    _SPE.FEATURES: [
                        "bf_weighted_volume_boltzfactor_dmso",
                        "bf_weighted_area_boltzfactor_dmso",
                        "bf_weighted_HB_acc_boltzfactor_dmso",
                        "bf_weighted_HB_don_boltzfactor_dmso",
                        "bf_weighted_sigma2_boltzfactor_dmso",
                        "bf_weighted_Gsolv_meoh_boltzfactor_dmso",
                    ],
                    _SPE.NAME_PREDICTED: "pred_ePSA",
                },
            },
        }
        pred_step = StepPredictor(**step_conf)
        pred_step.get_compounds().append(Compound())
        pred_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformer = get_mol_as_Conformer(self._example_mol_path)
        pred_step.data.compounds[0][0].add_conformers(conformer, auto_update=True)
        pred_step.execute()

        self.assertEqual(len(pred_step.get_compounds()), 1)
        self.assertEqual(len(pred_step.get_compounds()[0]), 1)
        self.assertEqual(len(pred_step.get_compounds()[0][0]), 1)

        # check SDF write-out (including ePSA prediction as tag)
        out_path = os.path.join(self._test_dir, "ePSA_predicted_annotated.sdf")
        pred_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 4448)
