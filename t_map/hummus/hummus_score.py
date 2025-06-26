from typing import Union, List, Tuple, Callable, Iterable, Set, Optional
from t_map.gene.gene import Gene
from enum import Enum, auto
import itertools

ScoreFn = Callable[[Union[List[Gene], List[Tuple[Gene, float]]]], None]


class ScoreTypes(Enum):
    TOP_K = auto()


class HummusScore:
    """Hummus Score is a module used for scoring a particular
    prioritization algorithm. Its main purpose is to be passed
    into an instance of Hummus to the `with_scoring` parameter.
    If provided during testing and training time, the Hummus
    module with additionally return a score function, which will
    accept a list of prediction with or without a score.

    Note: Hummus Score function returns None. It is the responsibility
    of Hummus Score to later provide results.
    """

    def __init__(self, score_type: ScoreTypes = ScoreTypes.TOP_K, k: int = 20):
        self._score_type = score_type
        self._k = k
        self._training = False
        self._scores: List[dict] = []
        self._training_scores: List[List[dict]] = []
        self._testing_scores: List[List[dict]] = []
        self._ground_truth: Union[None, Set[Gene]] = None

    def train(self) -> None:
        if len(self._scores) > 0:
            if self._training:
                self._training_scores.append(self._scores)
            else:
                self._testing_scores.append(self._scores)
            self._scores = []
        self._training = True

    def test(self) -> None:
        if len(self._scores) > 0:
            if self._training:
                self._training_scores.append(self._scores)
            else:
                self._testing_scores.append(self._scores)
            self._scores = []
        self._training = False

    def reset_ground_truth(self) -> None:
        self._ground_truth = None

    def score(self, heldout_genes: Iterable[Gene],
              predictions: Union[Iterable[Gene],
                                 Iterable[Tuple[Gene, float]]]) -> None:
        """

        Raises:
            Exception: If the requested score type does not exist
        """

        if self._ground_truth is None:
            self._ground_truth = set(map(lambda x: x.name, heldout_genes))

        if self._score_type is ScoreTypes.TOP_K:
            '''
                If the user passes in a list of tuples ensure we
                sort it by score, otherwise we assume it is passed
                in sorted order
            '''
            if all(isinstance(v, tuple) for v in predictions):
                # need to do this for type checking :/
                new_predictions: List[Tuple[Gene, float]
                                      ] = list(predictions)  # type: ignore
                sorted_predictions = sorted(
                    new_predictions,
                    key=lambda x: x[1],
                    reverse=True)  # type: ignore
                predictions = list(
                    map(lambda x: x[0], sorted_predictions))

            predictions = predictions[0:self._k]
            genes_in_ground_truth = list(filter(
                lambda x: x.name in self._ground_truth, predictions))

            self._scores.append({"top_k_predictions": genes_in_ground_truth,
                                 "score": len(genes_in_ground_truth)
                                 / len(self._ground_truth),
                                 "ground_truth": self._ground_truth,
                                 "predictions": predictions})
            return
        else:
            raise Exception(
                f"{self._score_type} is not a valid score function")

    def training_summary(self) -> str:
        if self._training and self._scores:
            self._training_scores.append(self._scores)
            self._scores = []

    def testing_summary(self) -> dict:
        """
            Raises:
                ZeroDivisionError: ...
        """
        # If we are in testing and have scores left
        # ensure that we include them.
        self._collate_testing_scores()

        merged_scores = list(
            itertools.chain.from_iterable(self._testing_scores))
        scores = list(map(lambda x: x['score'], merged_scores))

        predicted_genes = list(
            map(lambda x: (x['ground_truth'], x['predictions']),
                merged_scores))

        return {"k": self._k, "score": scores,
                "avg_score": sum(scores) / len(scores),
                "predicted_genes": predicted_genes}

    def _collate_testing_scores(self) -> None:
        if not self._training and self._scores:
            self._testing_scores.append(self._scores)
            self._scores = []

    def compute_roc(self) -> None:

        self._collate_testing_scores()
        merged_scores = list(
            itertools.chain.from_iterable(self._testing_scores))
        ranked_predictions = [testing_score['predictions']
                              for testing_score in merged_scores]
        target = [testing_score['ground_truth']
                  for testing_score in merged_scores]
        self._roc = roc_of_ranked_list(ranked_predictions, target)

    def auc_roc(self, threshold: Optional[float]) -> float:
        if not hasattr(self, '_roc'):
            raise Exception("ROC not computed")
        if self._roc is None:
            raise Exception("ROC not computed")
        tpr = list(map(lambda x: x[0], self._roc))
        fpr = list(map(lambda x: x[1], self._roc))
        if threshold is not None:
            fpr, tpr = get_fpr_tpr_for_thresh(fpr, tpr, threshold)
            return auc_from_fpr_tpr(fpr, tpr)
        else:
            return auc_from_fpr_tpr(fpr, tpr)


def roc_of_ranked_list(ranked_predictions: List[List[Gene]],
                       target: List[List[Gene]],
                       num_intervals: int = 20) -> List[Tuple[float, float]]:
    """
    Given an ordered list of ranked predictions and a target list,
    for each interval in the list, calculate the ROC score.

    Args:
        ranked_predictions: A list of ranked predictions
        target: A list of ground truth
        intervals: A list of intervals to calculate ROC scores
    """

    roc_scores = []
    num_predictions = len(ranked_predictions[0])
    intervals = list(range(1, num_predictions + 1,
                     num_predictions // num_intervals))
    for i, interval in enumerate(intervals):
        # Get the top i predictions
        top_i_predictions = [predictions[:interval]
                             for predictions in ranked_predictions]
        # Get the ground truth for those predictions
        # Calculate the ROC score
        roc_scores.append(roc_score(top_i_predictions, target,
                          end_range=len(ranked_predictions[0])))
    return roc_scores


def roc_score(ranked_predictions: List[List[Gene]],
              target: List[List[Gene]],
              end_range: int) -> Tuple[float, float]:
    """
    Given a ranked list of predictions and a target list,
    calculate the true positive rate and false positive rate.

    Args:
        ranked_predictions: A list of ranked predictions
        target: A list of ground truth
    Returns:
        Tuple of true positive rate and false positive rate
    """
    tp = 0
    for pred, tar in zip(ranked_predictions, target):
        pred = [gene.name for gene in pred]
        pred = set(pred)
        tar = set(tar)
        tp += len(pred & tar)

    return (tp / (len(target[0]) * len(target)),
            len(ranked_predictions[0]) / end_range)

# From: https://stackoverflow.com/questions/39537443/how-to-calculate-a-partial-area-under-the-curve-auc/53464588
# By: Ami Tavory


def auc_from_fpr_tpr(fpr: List[float], tpr: List[float],
                     trapezoid: bool = False) -> float:
    import numpy as np
    fpr = np.array(fpr)
    tpr = np.array(tpr)
    inds = [i for (i, (s, e)) in enumerate(
        zip(fpr[: -1], fpr[1:])) if s != e] + [len(fpr) - 1]
    fpr, tpr = fpr[inds], tpr[inds]
    area = 0
    ft = list(zip(fpr, tpr))
    for p0, p1 in zip(ft[: -1], ft[1:]):
        area += (p1[0] - p0[0]) * ((p1[1] + p0[1]) / 2 if trapezoid else p0[1])
    return area


def get_fpr_tpr_for_thresh(fpr: List[float],
                           tpr: List[float],
                           thresh: float) -> Tuple[List[float], List[float]]:
    import bisect
    p = bisect.bisect_left(fpr, thresh)
    assert p < len(fpr), f"Attempted to bisect_left on value {thresh=} in {fpr=} that yielded an index {p=} out of bounds"
    fpr = fpr.copy()
    fpr[p] = thresh
    return fpr[: p + 1], tpr[: p + 1]
