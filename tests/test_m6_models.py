import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from epicurus_neo.m6.models import fit_predict


def _toy():
    rng = np.random.default_rng(0)
    n = 60
    signal = rng.normal(size=n)
    frame = pd.DataFrame(
        {
            "feat_a": signal,
            "feat_b": rng.normal(size=n),
            "label": (signal > 0).astype(int),
        }
    )
    return frame.iloc[:40], frame.iloc[40:]


def test_prevalence_scores_are_constant_train_positive_rate():
    train, evaluation = _toy()
    scores = fit_predict("prevalence", train, evaluation, ["feat_a", "feat_b"])
    assert len(scores) == len(evaluation)
    assert np.allclose(scores, train.label.mean())


def test_learned_models_recover_a_separable_signal_and_are_deterministic():
    train, evaluation = _toy()
    for name in ("logistic", "boosting"):
        first = fit_predict(name, train, evaluation, ["feat_a", "feat_b"])
        second = fit_predict(name, train, evaluation, ["feat_a", "feat_b"])
        assert np.array_equal(first, second)  # determinism
        # feat_a separates the classes; ranking must beat coin-flip AUROC.
        assert roc_auc_score(evaluation.label, first) > 0.8
