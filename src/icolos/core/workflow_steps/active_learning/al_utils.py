from typing import List
import numpy as np
from icolos.loggers.steplogger import StepLogger
from icolos.core.workflow_steps.step import _LE

_logger = StepLogger()


def greedy_acquisition(
    estimator,
    X: np.ndarray,
    previous_idx: List[int],
    n_instances: int,
    warmup: bool = False,
    epsilon: float = 0.0,
) -> np.ndarray:
    """Implements greedy acquisition by querying model for top predicted points

    :param  estimator: SKLearn-type estimator to be queried
    :param np.ndarray X: array of fingerprints for each compound to be predicted by the model
    :param List[int] previous_idx: List of the previously queried indices to avoid repitition
    :param int n_instances: batch size to query
    :param bool warmup: Control warmup period in which random samples are generated, defaults to False
    :param float epsilon: enable greedy epsilon acquisition, defaults to 0.0
    :return np.ndarray: array of top predictions from the estimator
    """
    if warmup:
        _logger.log("Warmup epoch, using random sampling...", _LE.DEBUG)
        return np.array([np.random.randint(0, X.shape[0]) for _ in range(n_instances)])
    else:
        try:
            predictions = estimator.predict(X)

        except:
            _logger.log("estimator not fitted, using random sampling...", _LE.DEBUG)
            predictions = np.random.uniform(0, 12, X.shape[0])
            return np.array(
                [np.random.randint(0, X.shape[0]) for _ in range(n_instances)]
            )

    # zero those predictions we've seen before
    for idx in previous_idx:
        predictions[idx] = 0
    # smaller before n_instances, largest after
    sorted_preds = np.argpartition(predictions, -n_instances, axis=0)[-n_instances:]

    if epsilon > 0.0:
        # replace that fraction of the predictions with random samples
        n_replacements = int(n_instances * epsilon)

        # select random indices in sorted preds to select
        indices_to_replace = np.random.choice(
            [i for i in range(len(sorted_preds))], n_replacements, replace=False
        )
        # select random compound incides from the dataset to substitute into the returned predictions
        replacements = [np.random.randint(0, X.shape[0]) for _ in range(n_replacements)]
        _logger.log(
            f"Replacing {len(indices_to_replace)} predictions with random values",
            _LE.DEBUG,
        )
        for idx, val in zip(indices_to_replace, replacements):
            predictions[idx] = val
    return sorted_preds
