from typing import List
import numpy as np
from icolos.loggers.steplogger import StepLogger
from icolos.core.workflow_steps.step import _LE
from scipy.stats import norm
from modAL import ActiveLearner

_logger = StepLogger()


def greedy_acquisition(
    learner: ActiveLearner,
    X: np.ndarray,
    previous_idx: List[int],
    n_instances: int,
    warmup: bool = False,
    epsilon: float = 0.0,
) -> np.ndarray:
    """Implements greedy acquisition by querying model for top predicted points, note that typically we are dealing with affinity prediction so highest is not best.

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
        indices = range(0, X.shape[0])
        # sample indices without replacement
        return np.random.choice(indices, replace=False, size=n_instances)
    else:
        try:
            predictions = learner.predict(X)

        except:
            _logger.log(
                "estimator not fitted, selecting random compounds ...", _LE.DEBUG
            )

            return np.random.choice(indices, replace=False, size=n_instances)

    # zero those predictions we've seen before
    for idx in previous_idx:
        predictions[idx] = 0
    # smaller before n_instances, largest after, take the most negative
    sorted_preds = np.argpartition(predictions, n_instances, axis=0).reshape(-1)[
        :n_instances
    ]
    print("sorted preds", sorted_preds.shape, predictions[sorted_preds])

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


def expected_improvement(
    learner: ActiveLearner,
    X: np.ndarray,
    previous_idx: List[int],
    n_instances: int,
    warmup: bool = False,
    highest_is_best: bool = False,
    **kwargs
) -> np.ndarray:
    """Select new points to acquire based on expected improvement

    :param  estimator: SKLearn-type estimator to be queried
    :param np.ndarray X: array of fingerprints for each compound to be predicted by the model
    :param List[int] previous_idx: List of the previously queried indices to avoid repitition
    :param int n_instances: batch size to query
    :param bool warmup: Control warmup period in which random samples are generated, defaults to False
    :param float epsilon: enable greedy epsilon acquisition, defaults to 0.0
    :return np.ndarray: array of top predictions from the estimator
    """
    # if no stdev implemented (only rf + GP can do this as std), generate the standard deviation from repeated sampling
    # TODO: at the moment this will only work with RF
    if warmup:
        _logger.log("Warmup epoch, using random sampling...", _LE.DEBUG)
        return np.array([np.random.randint(0, X.shape[0]) for _ in range(n_instances)])

    estimator = learner.estimator
    subestimates = []
    estimators = estimator.estimators_

    for i in range(len(estimators)):
        estimator.estimators_ = estimators[0 : i + 1]
        subestimates.append(estimator.predict(X))
    subestimates = np.array(subestimates)
    stdev = np.std(subestimates, axis=0).reshape((-1))
    y_hat = np.array(learner.predict(X))
    mu = np.mean(y_hat)
    y_best = np.max(y_hat) if highest_is_best else np.min(y_hat)
    gamma = (mu - y_best) / stdev
    ei = stdev * gamma * norm.cdf(gamma) + stdev * norm.pdf(gamma)

    for idx in previous_idx:
        ei[idx] = 0
    # smaller before n_instances, largest after
    return np.argpartition(ei, n_instances, axis=0).reshape(-1)[:n_instances]
