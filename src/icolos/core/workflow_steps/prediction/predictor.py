import pickle
from copy import deepcopy

import numpy as np
from typing import List

from pydantic import BaseModel
from rdkit import Chem

from icolos.utils.general.icolos_exceptions import StepFailed, get_exception_message
from icolos.utils.enums.step_enums import StepPredictorEnum
from icolos.core.workflow_steps.io.base import StepIOBase
from icolos.core.workflow_steps.step import _LE

from icolos.utils.general.convenience_functions import *

_SPE = StepPredictorEnum()


class StepPredictor(StepIOBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    @classmethod
    def _load_scikit_model(cls, model_path: str):
        with open(model_path, "rb") as f:
            scikit_model = pickle.load(f)
            return scikit_model

    def _get_feature_values(
        self, conformer: Chem.Mol, feature_names: List[str]
    ) -> np.ndarray:
        list_values = []
        for feature in feature_names:
            try:
                list_values.append(float(conformer.GetProp(feature)))
            except KeyError as e:
                self._logger.log(
                    f"Could not find feature / property, error message: {get_exception_message(e)}",
                    _LE.ERROR,
                )
                raise e

        # cast list to 2D array
        return np.array([list_values])

    def execute(self):
        # get parameters
        parameters = deepcopy(self.settings.additional)
        model_path = nested_get(parameters, _SPE.MODEL_PATH, default=None)
        feature_names = nested_get(parameters, _SPE.FEATURES, default=None)
        name_predicted = nested_get(
            parameters, _SPE.NAME_PREDICTED, default=_SPE.NAME_PREDICTED_DEFAULT
        )

        # check parameters; model_path and features are mandatory
        if model_path is None or feature_names is None:
            message = f"Parameters {_SPE.MODEL_PATH} (path to model) and {_SPE.FEATURES} (list with features) have to be set - abort."
            self._logger.log(message, _LE.ERROR)
            raise StepFailed(message)
        if name_predicted == _SPE.NAME_PREDICTED_DEFAULT:
            self._logger.log(
                f"Name of predicted property not specified, using default value {_SPE.NAME_PREDICTED_DEFAULT} instead (not recommended).",
                _LE.WARNING,
            )

        # load model from file and predict endpoint
        model = self._load_scikit_model(model_path=model_path)
        predicted = 0
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                for conformer in enumeration.get_conformers():
                    if not self._input_object_valid(conformer):
                        continue

                    f_values = self._get_feature_values(
                        conformer=conformer.get_molecule(), feature_names=feature_names
                    )
                    conformer.get_molecule().SetProp(
                        name_predicted, str(model.predict(X=f_values)[0])
                    )
                    predicted += 1
        self._logger.log(
            f"Predicted {name_predicted} for {predicted} conformers in {len(self.get_compounds())} compounds.",
            _LE.INFO,
        )
