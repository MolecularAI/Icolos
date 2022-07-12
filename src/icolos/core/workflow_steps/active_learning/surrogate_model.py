import numpy as np

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
                # enforcing max-depth may cause trees to be unable to predict extreme values
                # due to under-representation in the training data
                n_estimators=100, n_jobs=-1, criterion="mse"
            )

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
        """train surrogate model from scratch. No online learning for random forest and support vectors"""
        if self.type == _SALE.RANDOM_FOREST_REGRESSOR:
            self.model.fit(X_train, y_train)
        if self.type == _SALE.SUPPORT_VECTOR_REGRESSOR:
            self.model.fit(np.array(X_train), y_train)

    def predict(self, data):
        """make predictions given query fingerprints"""
        if (
            self.type == _SALE.RANDOM_FOREST_REGRESSOR
            or self.type == _SALE.SUPPORT_VECTOR_REGRESSOR
        ):
            return self.model.predict(data)

    def get_std(self, X):
        """
        get standard deviation. This method is only implemented for random forests and is calculated
        using the prediction differences of each individual tree in the ensemble.
        """
        if self.type == _SALE.RANDOM_FOREST_REGRESSOR:
            individual_trees = self.model.estimators_
            subEstimates = np.array([tree.predict(X) for tree in individual_trees])
            return np.std(subEstimates, axis=0)

        elif self.type == _SALE.SUPPORT_VECTOR_REGRESSOR:
            raise NotImplementedError(
                "Standard deviation calculation not supported for SVR"
            )
