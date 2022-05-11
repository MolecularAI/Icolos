import pandas as pd
import numpy as np
from icolos.core.workflow_steps.step import _LE

from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR

from icolos.utils.enums.step_enums import (
    StepActiveLearningEnum,
)

_SALE = StepActiveLearningEnum()


class SurrogateModel:
    """class for the surrogate model used in prospective REINVENT"""

    # TODO: only random forests have been tested so far. Test SVR and maybe deep learning models
    def __init__(self, model_type: str):
        def initialize_RF():
            model = RandomForestRegressor(
                n_estimators=100, max_depth=8, n_jobs=-1, criterion="mse"
            )  # the rest of the parameters are left default

            return model

        def initialize_SVR():
            model = SVR(kernel="rbf", gamma="scale")
            return model

        self.type = model_type.lower()

        # random forest regressor
        if self.type == _SALE.RANDOM_FOREST_REGRESSOR:
            self.model = initialize_RF()
        # support vector regressor
        elif self.type == _SALE.SUPPORT_VECTOR_REGRESSOR:
            self.model = initialize_SVR()
        else:
            raise ValueError("Model not supported.")

    def fit(self, X_train, y_train):
        if self.type == _SALE.RANDOM_FOREST_REGRESSOR:
            # TODO: setting a seed would allow reproducibility, add if needed --> np.random.seed(0)
            self.model.fit(X_train, y_train)
        if self.type == _SALE.SUPPORT_VECTOR_REGRESSOR:
            self.model.fit(np.array(X_train), y_train)

    def predict(self, data):
        if (
            self.type == _SALE.RANDOM_FOREST_REGRESSOR
            or self.type == _SALE.SUPPORT_VECTOR_REGRESSOR
        ):
            return self.model.predict(data)

    def get_std(self, X):
        if self.type == _SALE.RANDOM_FOREST_REGRESSOR:
            individual_trees = self.model.estimators_
            subEstimates = np.array([tree.predict(X) for tree in individual_trees])
            return np.std(subEstimates, axis=0)

        elif self.type == _SALE.SUPPORT_VECTOR_REGRESSOR:
            raise NotImplementedError(
                "Standard deviation calculation not supported for SVR"
            )
