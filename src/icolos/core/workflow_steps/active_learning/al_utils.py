from typing import List
import numpy as np
from icolos.loggers.steplogger import StepLogger
from icolos.core.workflow_steps.step import _LE
from scipy.stats import norm
from modAL import ActiveLearner
from skorch.regressor import NeuralNetRegressor
from numpy.linalg import norm
import torch
from torch import nn
from skorch.exceptions import NotInitializedError

_logger = StepLogger()


class FeedForwardNet(nn.Module):
    def __init__(self) -> None:
        super(FeedForwardNet, self).__init__()

        self.lin1 = nn.Linear(2048, 128)
        self.relu1 = nn.LeakyReLU()
        self.lin2 = nn.Linear(128, 128)
        self.relu2 = nn.LeakyReLU()
        self.lin3 = nn.Linear(128, 128)
        self.relu3 = nn.LeakyReLU()
        self.dropout = nn.Dropout()
        self.lin4 = nn.Linear(128, 1)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        x = self.relu1(self.lin1(inputs))
        x = self.relu2(self.lin2(x))
        x = self.relu3(self.lin3(x))
        x = self.lin4(x)
        return x

    def get_latent_embedding(self, inputs: torch.Tensor) -> torch.Tensor:
        """Return the latent feature representatino from the last hidden layer

        :param torch.Tensor inputs: input features
        :return _type_: _description_
        """
        with torch.no_grad():
            x = self.relu1(self.lin1(inputs))
            x = self.relu2(self.lin2(x))
            x = self.lin3(x)
        return x


def latent_distance(
    learner: ActiveLearner, X: np.ndarray, previous_idx: List[int], n_instances: int
):
    """Selects new points to acquire based on the UQ approach detailed here https://doi.org/10.1039/C9SC02298H

    :param ActiveLearner learner: trained active learning model to query
    :param np.ndarray X: array containing the fingerprints to predict values from
    """
    try:
        predictions = learner.predict(X)

    except NotInitializedError:
        _logger.log("estimator not fitted, selecting random compounds ...", _LE.DEBUG)

        indices = range(0, X.shape[0])
        return np.random.choice(indices, replace=False, size=n_instances)

    # compute pairwise euclidian distances of each predicted point to the known points, return those indices with the highest distance to known
    known_X = X[previous_idx]
    # unknown_X = np.delete(X, previous_idx, axis=0)
    # print("unknown X", unknown_X.shape)

    # need the embedding after the final layer
    # access the skorch object
    estimator: NeuralNetRegressor = learner.estimator
    known_latent_embedding = estimator.module_.get_latent_embedding(
        inputs=torch.from_numpy(known_X).to(estimator.device)
    )
    print("known_latent_embedding", known_latent_embedding.shape)
    # n_samples, emb_dim
    # include the full dataset here to retain indexing, some norms will be 0, but that is fine
    all_embedding = estimator.module_.get_latent_embedding(
        inputs=torch.from_numpy(torch.from_numpy(X)).to(estimator.device)
    )
    known_latent_embedding, all_embedding = (
        known_latent_embedding.cpu(),
        all_embedding.cpu(),
    )
    print("unknown_latent_embedding", all_embedding.shape)

    # now compute distance to each known embedding vector for each unknown point
    distances = []
    for idx, embedding in enumerate(all_embedding):
        # min(abs(x), axis=0)
        # return the distance to the closest point in the training set
        # Note, these distances are uncalibrated
        dist = np.linalg.norm(embedding - known_latent_embedding, ord=-1)
        if dist < 1:
            print(idx, dist)
        distances.append(dist)
    print(np.array(distances).shape)
    print(distances)
    # get indices of the largest
    sorted_preds = np.argpartition(np.array(distances), n_instances, axis=0).reshape(
        -1
    )[:n_instances]

    return sorted_preds


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

            indices = range(0, X.shape[0])
            return np.random.choice(indices, replace=False, size=n_instances)

    # zero those predictions we've seen before
    # TODO: remove this, just predict on the unseen variables
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
    **kwargs,
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
