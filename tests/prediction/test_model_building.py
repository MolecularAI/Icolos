import json
import unittest
import os

from icolos.core.workflow_steps.prediction.model_building import StepModelBuilder
from icolos.utils.enums.program_parameters import ModelBuilderEnum

from icolos.utils.enums.step_enums import StepBaseEnum, StepModelBuilderEnum

from tests.tests_paths import PATHS_EXAMPLEDATA, load_SDF_docked, MAIN_CONFIG
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SME = ModelBuilderEnum()
_SMBE = StepModelBuilderEnum()


class Test_Model_Building(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/model_building")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        self._example_JSON = PATHS_EXAMPLEDATA.MODEL_BUILDER_EXAMPLE_JSON
        self._compounds = load_SDF_docked(
            PATHS_EXAMPLEDATA.MODEL_BUILDER_TEST_INPUT_SDF
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_build_model(self):
        step_conf = {
            _SBE.STEPID: "01_model_building",
            _SBE.STEP_TYPE: _SBE.STEP_PREDICTION,
            _SBE.EXEC: {
                _SBE.EXEC_BINARYLOCATION: " ".join(
                    [
                        MAIN_CONFIG["OPTUNA_AZ"]["ENVIRONMENT_PYTHON"],
                        MAIN_CONFIG["OPTUNA_AZ"]["ENTRY_POINT_LOCATION"],
                    ]
                )
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _SME.CONFIG: self._example_JSON,
                        _SME.BEST_BUILDCONFIG_OUTPATH: os.path.join(
                            self._test_dir, "buildconfig.json"
                        ),
                        _SME.BEST_MODEL_OUTPATH: os.path.join(
                            self._test_dir, "best_model_trial.pkl"
                        ),
                        _SME.MERGED_MODEL_OUTPATH: os.path.join(
                            self._test_dir, "production_model.pkl"
                        ),
                    }
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SMBE.DATA: {
                        _SMBE.DATA_INPUT_COLUMN: "original_smiles",
                        _SMBE.DATA_RESPONSE_COLUMN: _SBE.ANNOTATION_TAG_DOCKING_SCORE,
                    }
                },
            },
        }
        model_step = StepModelBuilder(**step_conf)
        model_step.data.compounds = self._compounds

        model_step.execute()

        # check, that the input data has been written as expected
        out_path = os.path.join(self._test_dir, "best_param.json")
        container = model_step.data.generic.get_files_by_extension(ext="json")[0]
        with open(out_path, "w") as f:
            json.dump(container.get_data(), f, indent=4)
        stat_inf = os.stat(out_path)
        self.assertEqual(_SMBE.TMP_OUTPUT_BEST_PARAMETERS, container.get_file_name())
        self.assertGreater(stat_inf.st_size, 800)

        # check, that a model has been produced
        # note, that the model's size strongly depends on the underlying algorithm / hyper-parameters chosen
        out_path = os.path.join(self._test_dir, "production_model.pkl")
        data = model_step.data.generic.get_files_by_extension(ext="pkl")[0].get_data()
        with open(out_path, "wb") as f:
            f.write(data)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 5000)
