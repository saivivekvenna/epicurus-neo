import numpy as np
import pytest

from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.event_a import load_event_a_frame
from epicurus_neo.m6.transfer import (
    TEACHER_SCORE_COLUMN,
    add_teacher_score,
    assert_no_event_a_leakage,
    train_frozen_teacher,
)


def test_frozen_teacher_scores_every_event_b_candidate_including_long_slps():
    teacher = train_frozen_teacher(load_event_a_frame(), seed=17)
    event_b = load_label_frame()
    scored = add_teacher_score(event_b, teacher)
    series = scored[TEACHER_SCORE_COLUMN]
    assert len(scored) == len(event_b)
    # Length-agnostic core features mean the teacher scores EVERY candidate, including
    # long class-II SLPs it never saw in training — no NOT_VIABLE compatibility collapse.
    assert np.isfinite(series.to_numpy()).all()
    assert (series >= 0).all() and (series <= 1).all()
    long_slps = scored[scored.mutant_peptide.map(len) > 15]
    assert len(long_slps) > 0
    assert np.isfinite(long_slps[TEACHER_SCORE_COLUMN].to_numpy()).all()


def test_teacher_is_deterministic():
    event_b = load_label_frame()
    first = add_teacher_score(event_b, train_frozen_teacher(load_event_a_frame(), seed=17))
    second = add_teacher_score(event_b, train_frozen_teacher(load_event_a_frame(), seed=17))
    assert np.array_equal(
        first[TEACHER_SCORE_COLUMN].to_numpy(), second[TEACHER_SCORE_COLUMN].to_numpy()
    )


def test_leakage_guard_fires_on_injected_overlap():
    event_a = load_event_a_frame()
    event_b = load_label_frame()
    poisoned = event_a.copy()
    poisoned.loc[poisoned.index[0], "mutant_peptide"] = event_b.mutant_peptide.iloc[0]
    poisoned.loc[poisoned.index[0], "hla_allele"] = event_b.hla_allele.iloc[0]
    with pytest.raises(AssertionError):
        assert_no_event_a_leakage(poisoned, event_b.head(1))


def test_real_corpora_have_no_event_a_event_b_leakage():
    # If this ever fails, the teacher has seen a test peptide and M6B is invalid.
    assert assert_no_event_a_leakage(load_event_a_frame(), load_label_frame()) is None
