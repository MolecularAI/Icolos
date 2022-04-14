from typing import Callable, List
from pydantic import BaseModel
import numpy as np
from icolos.core.workflow_steps.step import _LE

# TODO: estimator can be either a pytorch model or a sklearn estimator, using skorch wrapping for now, but the skorch API seems fragile in some cases


class ActiveLearner(BaseModel):
    """
    Implements the base active learner
    :param fit_to_new: bool: control whether to train the model just on the new data, or on the entire dataset at each refitting step
    """

    def __init__(
        self,
        estimator,
        acquisition_function: Callable,
        fit_to_new: bool = True,
        X: np.ndarray = None,
        y: np.ndarray = None,
    ) -> None:
        self.estimator = estimator
        self.acquisition_function = acquisition_function
        self.fit_to_new = fit_to_new
        self.X = X
        self.y = y

    def teach(self, X: np.ndarray, y: np.ndarray):
        """
        Take new batch of data from the oracle, train the model
        """
        if self.fit_to_new:
            # just fit the model on the new data, works for tf models etc
            self.estimator.fit(X, y)
        else:
            self.append_new_data(X, y)
            self.estimator.fit(self.X, self.y)

    def append_new_data(self, X: np.ndarray, y: np.ndarray):
        if self.X is None:
            self.X = X
        else:
            # TODO: make sure dims match here
            self.X = np.concatenate(self.X, self)
        if self.y is None:
            self.y = y
        else:
            # TODO: make sure dims match here
            self.y = np.concatenate(self.y, self)

    def query(self, previous_idx: List[int], n_samples: int):
        """
        Query the model using the acquisition function, return the n_samples
        :param previous_idx: list of indices of previously requested samples
        """
        return self.acquisition_function(
            self.estimator, self.X, previous_idx, n_samples
        )
