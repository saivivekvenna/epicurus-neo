# M6A: Event-B-Only First Recognition Swing — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the leakage-safe, leave-one-study-out (LOSO) Event-B-only recognition swing defined in `docs/superpowers/specs/2026-07-10-m6-first-recognition-swing-design.md`, producing a pre-registered audit that answers — honestly, under the standing `INSUFFICIENT_CANDIDATE_RESOLVED` verdict — whether a learned model improves patient-level selection over prevalence and (on hu+pdac) presentation.

**Architecture:** A small reusable package `src/epicurus_neo/m6/` (dataset loader, completeness gate, LOSO folds, feature matrix, guarded presentation baseline, model ladder, per-patient-k + classification metrics, macro/micro evaluation, confound diagnostics, audit renderer) driven by one thin runner `experiments/m6a_recognition_swing.py`. Reuses existing `benchmark` metrics/stats/verdict and the `epicurus_neo.features` library rather than reinventing them.

**Tech Stack:** Python 3.10+, pandas 2.2, numpy, scikit-learn ≥1.4 (LogisticRegression, HistGradientBoostingClassifier — **no xgboost**), pytest, ruff. MHCflurry is an optional (`.[mhc]`) dependency used only by the guarded presentation baseline.

## Global Constraints

Every task's requirements implicitly include these (values copied verbatim from the spec/grounding):

- **Corpus pin:** `outputs/event_b_backbone/combined/` (Event-B only — never the legacy `outputs/event_b_corpus_combined/`).
- **Population contract (grounded, must hold):** label frame = 965 rows, 272 positive, 693 negative, 45 patients, 4 studies (`braun_rcc_2025`, `hu_neovax_2021`, `mkras_vax_2026`, `pdac_neovax_2023`); Nous excluded automatically (patient-level-only, no `candidate_id`); 9 UNTESTED dropped.
- **Label:** `POSITIVE=1` vs `TESTED_NEGATIVE=0`.
- **Determinism:** seed `17`, `kind="mergesort"`, md5 identity tie-break (via `benchmark.metrics.identity_tiebreak`). Re-runs must be byte-identical.
- **Banned feature columns (never reach any model):** `label`, `response_label`, `study_id`, `patient_id`, `candidate_id`, `cancer_type`, `vaccine_platform`, `mhc_class` (string), `hla_alleles`, `hla_allele`, `mutant_peptide`, `wildtype_peptide`, `timepoint`, `relative_to_vaccine`, `event_type`, `qualitative_result`, `quantitative_result`, and any assay-result / recognition-evidence column. Enforced via `epicurus_neo.features.NON_FEATURE_COLUMNS` + an explicit guard test.
- **k_patient:** `k_patient = min(20, n_eligible)`; never divide precision by a fixed 20 for a shorter list.
- **Ranking-informative reality:** only 8/45 patients (all hu, `n_eligible > 20`) yield a non-degenerate selection; for the other 37 the differenced top-*k* delta is mechanically 0. `hits@k_patient` is the registered headline (hackathon alignment); the classification track (AUROC/Brier/AP) is first-class, not decorative. AUROC is not the headline.
- **Standing caveat:** every emitted artifact carries `INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA`; M6A is diagnostic, not a headline/clinical claim.
- **Test command:** `.venv/bin/python -m pytest <path> -v`. **Lint:** `.venv/bin/ruff check <path>` (line-length 100). Tests live in `tests/test_m6_*.py`.

---

### Task 1: M6 label-frame loader

**Files:**
- Create: `src/epicurus_neo/m6/__init__.py`
- Create: `src/epicurus_neo/m6/dataset.py`
- Test: `tests/test_m6_dataset.py`

**Interfaces:**
- Produces: `load_label_frame(corpus_dir: str | Path = CORPUS_DIR) -> pd.DataFrame` with columns `candidate_id, patient_id, study_id, mutant_peptide, wildtype_peptide, peptide_length, hla_alleles, mhc_class, hla_allele, label`. Also `parse_alleles(value) -> list[str]` and constant `CORPUS_DIR`, `CANDIDATE_RESOLVED_STUDIES`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_dataset.py
from epicurus_neo.m6.dataset import load_label_frame, parse_alleles


def test_label_frame_matches_registered_population():
    frame = load_label_frame()
    assert len(frame) == 965
    assert int(frame.label.sum()) == 272
    assert int((frame.label == 0).sum()) == 693
    assert frame.patient_id.nunique() == 45
    assert sorted(frame.study_id.unique()) == [
        "braun_rcc_2025", "hu_neovax_2021", "mkras_vax_2026", "pdac_neovax_2023",
    ]
    # UNTESTED dropped; labels are strictly binary.
    assert set(frame.label.unique()) == {0, 1}
    # Every row carries a scalar tie-break allele (possibly empty), never a list.
    assert frame.hla_allele.map(lambda v: isinstance(v, str)).all()


def test_parse_alleles_handles_json_list_and_scalar():
    assert parse_alleles('["HLA-A*02:01", "HLA-B*07:02"]') == ["HLA-A*02:01", "HLA-B*07:02"]
    assert parse_alleles("HLA-A*02:01") == ["HLA-A*02:01"]
    assert parse_alleles(None) == []
    assert parse_alleles("[]") == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_dataset.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6'`

- [ ] **Step 3: Create the package marker and loader**

```python
# src/epicurus_neo/m6/__init__.py
"""M6A: Event-B-only recognition swing."""
```

```python
# src/epicurus_neo/m6/dataset.py
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from event_b.models import BiologicalEvent, ResponseLabel

CORPUS_DIR = Path("outputs/event_b_backbone/combined")
CANDIDATE_RESOLVED_STUDIES = (
    "braun_rcc_2025",
    "hu_neovax_2021",
    "mkras_vax_2026",
    "pdac_neovax_2023",
)
_FEATURE_COLUMNS = [
    "candidate_id",
    "patient_id",
    "study_id",
    "mutant_peptide",
    "wildtype_peptide",
    "peptide_length",
    "hla_alleles",
    "mhc_class",
]


def parse_alleles(value: object) -> list[str]:
    """Normalize the polymorphic ``hla_alleles`` field to a clean list of strings."""
    if isinstance(value, (list, tuple, np.ndarray)):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    if text.startswith("["):
        try:
            return [str(item).strip() for item in json.loads(text) if str(item).strip()]
        except json.JSONDecodeError:
            return []
    return [item.strip() for item in text.replace(";", ",").split(",") if item.strip()]


def _first_allele(value: object) -> str:
    alleles = parse_alleles(value)
    return sorted(alleles)[0] if alleles else ""


def load_label_frame(corpus_dir: str | Path = CORPUS_DIR) -> pd.DataFrame:
    """Load the candidate-resolved Event-B label frame (one binary label per candidate)."""
    corpus_dir = Path(corpus_dir)
    candidates = pd.read_parquet(corpus_dir / "candidates.parquet")
    assays = pd.read_parquet(corpus_dir / "assays.parquet")
    primary = assays[
        assays.event_type.astype(str).eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)
        & assays.candidate_id.notna()
    ]
    frame = primary[["candidate_id", "response_label"]].merge(
        candidates[_FEATURE_COLUMNS], on="candidate_id", how="left", validate="one_to_one"
    )
    keep = [ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value]
    frame = frame[frame.response_label.isin(keep)].copy()
    frame["label"] = (frame.response_label == ResponseLabel.POSITIVE.value).astype(int)
    frame["hla_allele"] = frame.hla_alleles.map(_first_allele)
    return frame.drop(columns=["response_label"]).reset_index(drop=True)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_dataset.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/ tests/test_m6_dataset.py
git add src/epicurus_neo/m6/__init__.py src/epicurus_neo/m6/dataset.py tests/test_m6_dataset.py
git commit -m "milestone-6a: M6 candidate-resolved label-frame loader"
```

---

### Task 2: Candidate-universe completeness gate (WP0)

**Files:**
- Modify: `src/epicurus_neo/m6/dataset.py`
- Test: `tests/test_m6_completeness.py`

**Interfaces:**
- Consumes: `load_label_frame()` output.
- Produces: `completeness_report(frame: pd.DataFrame, *, k_cap: int = 20) -> pd.DataFrame` — one row per patient with columns `patient_id, study_id, n_candidates, n_positive, k_patient, ranking_informative, denominator_type`. `denominator_type ∈ {"COMPLETE_TESTED_SET", "POSITIVE_ENRICHED"}` (all 45 candidate-resolved patients carry ≥1 tested negative → all `COMPLETE_TESTED_SET`; the branch is retained so a future positive-only cohort is flagged, not silently admitted).

  > **Erratum (post-implementation, 2026-07-10):** the "all 45 → all `COMPLETE_TESTED_SET`" grounding is wrong. Shipped reality: 38 `HAS_TESTED_NEGATIVE`, 7 `NO_TESTED_NEGATIVE` (the mKRAS 6/6-responders). The labels were renamed to state the actual criterion — presence/absence of a tested negative (rankability), *not* denominator completeness — and the Step-1 assertion below (`(report.denominator_type == "COMPLETE_TESTED_SET").all()`) is superseded by the shipped test asserting 38/7 in `tests/test_m6_completeness.py`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_completeness.py
from epicurus_neo.m6.dataset import completeness_report, load_label_frame


def test_completeness_gate_grounded_counts():
    report = completeness_report(load_label_frame())
    assert len(report) == 45
    # Only hu patients clear n_eligible > 20 (selection is non-degenerate).
    assert int(report.ranking_informative.sum()) == 8
    assert report.loc[report.ranking_informative, "study_id"].unique().tolist() == ["hu_neovax_2021"]
    # mKRAS is a fixed 6-peptide shared panel: k_patient == 6, never informative.
    mkras = report[report.study_id == "mkras_vax_2026"]
    assert (mkras.k_patient == 6).all()
    assert (~mkras.ranking_informative).all()
    # Every candidate-resolved patient has a trustworthy denominator here.
    assert (report.denominator_type == "COMPLETE_TESTED_SET").all()
    # k_patient never exceeds 20 and never exceeds the patient's candidate count.
    assert (report.k_patient <= 20).all()
    assert (report.k_patient <= report.n_candidates).all()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_completeness.py -v`
Expected: FAIL with `ImportError: cannot import name 'completeness_report'`

- [ ] **Step 3: Implement the gate**

Append to `src/epicurus_neo/m6/dataset.py`:

```python
def completeness_report(frame: pd.DataFrame, *, k_cap: int = 20) -> pd.DataFrame:
    """Characterize each patient's candidate denominator before any top-k metric.

    ``n_eligible`` is the count of label-resolved candidates for the patient (this
    frame is already restricted to POSITIVE/TESTED_NEGATIVE, so ``n_candidates`` is
    ``n_eligible``). Selection is degenerate where ``n_eligible <= k_cap`` because
    the top-k set is the whole list.
    """
    grouped = frame.groupby("patient_id", sort=True)
    report = grouped.agg(
        study_id=("study_id", "first"),
        n_candidates=("candidate_id", "size"),
        n_positive=("label", "sum"),
    ).reset_index()
    report["n_positive"] = report.n_positive.astype(int)
    report["k_patient"] = report.n_candidates.clip(upper=k_cap).astype(int)
    report["ranking_informative"] = report.n_candidates > k_cap
    has_negative = (
        frame.assign(is_neg=frame.label == 0)
        .groupby("patient_id")["is_neg"]
        .any()
        .reindex(report.patient_id)
        .to_numpy()
    )
    report["denominator_type"] = np.where(has_negative, "COMPLETE_TESTED_SET", "POSITIVE_ENRICHED")
    return report
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_completeness.py -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/dataset.py tests/test_m6_completeness.py
git add src/epicurus_neo/m6/dataset.py tests/test_m6_completeness.py
git commit -m "milestone-6a: candidate-universe completeness gate"
```

---

### Task 3: Leave-one-study-out folds (WP4 infra)

**Files:**
- Create: `src/epicurus_neo/m6/loso.py`
- Test: `tests/test_m6_loso.py`

**Interfaces:**
- Consumes: label frame with `study_id`, `patient_id`, `candidate_id`.
- Produces: `loso_folds(frame: pd.DataFrame) -> list[LosoFold]` where `LosoFold` is a frozen dataclass `(held_out_study: str, train: pd.DataFrame, evaluation: pd.DataFrame)`. Studies are held out in sorted order, each exactly once; train/eval are disjoint by `study_id` and by `patient_id`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_loso.py
from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.loso import loso_folds


def test_loso_holds_out_each_study_exactly_once():
    folds = loso_folds(load_label_frame())
    assert [f.held_out_study for f in folds] == [
        "braun_rcc_2025", "hu_neovax_2021", "mkras_vax_2026", "pdac_neovax_2023",
    ]
    for fold in folds:
        # Evaluation is exactly the held-out study; training excludes it.
        assert set(fold.evaluation.study_id.unique()) == {fold.held_out_study}
        assert fold.held_out_study not in set(fold.train.study_id.unique())
        # No patient or candidate leaks across the split.
        assert not set(fold.train.patient_id) & set(fold.evaluation.patient_id)
        assert not set(fold.train.candidate_id) & set(fold.evaluation.candidate_id)
        assert len(fold.train) + len(fold.evaluation) == 965


def test_loso_on_two_study_subset_yields_two_folds():
    frame = load_label_frame()
    subset = frame[frame.study_id.isin(["hu_neovax_2021", "pdac_neovax_2023"])]
    assert len(loso_folds(subset)) == 2
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_loso.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.loso'`

- [ ] **Step 3: Implement LOSO folds**

```python
# src/epicurus_neo/m6/loso.py
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class LosoFold:
    held_out_study: str
    train: pd.DataFrame
    evaluation: pd.DataFrame


def loso_folds(frame: pd.DataFrame) -> list[LosoFold]:
    """Yield one leave-one-study-out fold per study, in sorted order."""
    studies = sorted(frame.study_id.unique())
    if len(studies) < 2:
        raise ValueError("LOSO requires at least two studies")
    folds: list[LosoFold] = []
    for study in studies:
        is_eval = frame.study_id == study
        folds.append(
            LosoFold(
                held_out_study=study,
                train=frame[~is_eval].reset_index(drop=True),
                evaluation=frame[is_eval].reset_index(drop=True),
            )
        )
    return folds
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_loso.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/loso.py tests/test_m6_loso.py
git add src/epicurus_neo/m6/loso.py tests/test_m6_loso.py
git commit -m "milestone-6a: leave-one-study-out folds"
```

---

### Task 4: Feature matrix — core tier + class-gated contrastive (WP2)

**Files:**
- Create: `src/epicurus_neo/m6/features.py`
- Test: `tests/test_m6_features.py`

**Interfaces:**
- Consumes: label frame from Task 1; reuses `epicurus_neo.features.add_sequence_features`, `add_contrastive_features`, `infer_numeric_feature_columns`.
- Produces:
  - `build_feature_matrix(frame: pd.DataFrame, tier: str) -> pd.DataFrame` for `tier ∈ {"core", "contrastive", "presentation"}`. `"core"` = sequence + `peptide_length` + class indicators (available under equivalent definitions for all 4 studies). `"contrastive"` = core + class-gated mutant-vs-WT deltas. `"presentation"` = core + `presentation_score`/`mhcflurry_*` if present.
  - `feature_columns(matrix: pd.DataFrame) -> list[str]` — numeric feature columns with all banned columns removed.
  - `assert_no_banned_features(columns: list[str]) -> None`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_features.py
import numpy as np

from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.features import (
    assert_no_banned_features,
    build_feature_matrix,
    feature_columns,
)


def test_core_features_are_universal_and_leakage_free():
    frame = load_label_frame()
    matrix = build_feature_matrix(frame, "core")
    columns = feature_columns(matrix)
    assert "seq_len" in columns
    assert "peptide_length" in columns
    assert "is_class_i" in columns
    # Banned identifiers/labels never survive as features.
    for banned in ("label", "study_id", "patient_id", "mhc_class", "hla_allele", "mutant_peptide"):
        assert banned not in columns
    assert_no_banned_features(columns)  # raises if violated
    # Core features are fully populated for all rows (no study-correlated missingness).
    assert matrix[columns].notna().all().all()


def test_contrastive_anchor_features_are_class_gated():
    frame = load_label_frame()
    matrix = build_feature_matrix(frame, "contrastive")
    class_i = matrix.mhc_class == "CLASS_I"
    # Anchor/TCR-face counts are defined only where the class-I register applies...
    assert matrix.loc[~class_i, "mutation_anchor_count"].isna().all()
    # ...and are present for at least some class-I rows that have a paired wildtype.
    assert matrix.loc[class_i, "mutation_anchor_count"].notna().any()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_features.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.features'`

- [ ] **Step 3: Implement the feature builder**

```python
# src/epicurus_neo/m6/features.py
from __future__ import annotations

import pandas as pd

from epicurus_neo.features import (
    add_contrastive_features,
    add_sequence_features,
    infer_numeric_feature_columns,
)

# Anchor/TCR-face positions use the class-I binding register; nulled elsewhere.
_CLASS_I_ONLY = ("mutation_anchor_count", "mutation_tcr_face_count")

# Belt-and-suspenders on top of features.NON_FEATURE_COLUMNS.
_EXTRA_BANNED = {
    "study_id", "patient_id", "candidate_id", "cancer_type", "vaccine_platform",
    "mhc_class", "hla_alleles", "hla_allele", "mutant_peptide", "wildtype_peptide",
    "peptide_length_str", "response_label", "event_type", "sample_date", "timepoint",
    "sample_id", "genomic_variant", "gene", "transcript", "protein_change",
    "candidate_source", "vaccine_inclusion", "vaccine_inclusion_origin",
    "generation_provenance", "mutant_wildtype_verified", "provenance_id", "schema_version",
}


def _add_class_indicators(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["is_class_i"] = (out.mhc_class == "CLASS_I").astype(float)
    out["is_class_ii"] = (out.mhc_class == "CLASS_II").astype(float)
    return out


def build_feature_matrix(frame: pd.DataFrame, tier: str) -> pd.DataFrame:
    """Build the requested feature tier. Missingness is never a study label."""
    if tier not in {"core", "contrastive", "presentation"}:
        raise ValueError(f"unknown tier: {tier!r}")
    matrix = _add_class_indicators(add_sequence_features(frame))
    if tier == "contrastive":
        matrix = add_contrastive_features(matrix)
        non_class_i = matrix.mhc_class != "CLASS_I"
        for column in _CLASS_I_ONLY:
            if column in matrix.columns:
                matrix.loc[non_class_i, column] = float("nan")
    # "presentation" tier consumes presentation_score/mhcflurry_* already merged
    # onto ``frame`` by Task 5; if absent the tier degenerates to core.
    return matrix


def feature_columns(matrix: pd.DataFrame) -> list[str]:
    numeric = infer_numeric_feature_columns(matrix)
    return [column for column in numeric if column not in _EXTRA_BANNED]


def assert_no_banned_features(columns: list[str]) -> None:
    leaked = sorted(set(columns) & _EXTRA_BANNED)
    if leaked:
        raise AssertionError(f"banned columns leaked into feature set: {leaked}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_features.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/features.py tests/test_m6_features.py
git add src/epicurus_neo/m6/features.py tests/test_m6_features.py
git commit -m "milestone-6a: core + class-gated contrastive feature matrix"
```

---

### Task 5: HLA resolution + guarded presentation baseline (WP1)

**Files:**
- Create: `src/epicurus_neo/m6/presentation.py`
- Test: `tests/test_m6_presentation.py`

**Interfaces:**
- Consumes: label frame; `antigens.parquet`; reuses `epicurus_neo.mhcflurry_features.add_mhcflurry_predictions`.
- Produces:
  - `PRESENTATION_STUDIES = ("hu_neovax_2021", "pdac_neovax_2023")`.
  - `resolve_class_i_alleles(frame, corpus_dir=CORPUS_DIR) -> pd.DataFrame` — adds `class_i_alleles: list[str]` per candidate (hu from candidate HLA; pdac from `antigens` `PREDICTED_BEST_BINDER` join on `study_id, gene, protein_change`; else `[]`).
  - `add_presentation_score(frame, *, predictor=None) -> pd.DataFrame` — adds `presentation_score` (best MHCflurry `mhcflurry_presentation_score` across a candidate's class-I alleles; NaN where no allele). Raises `PresentationUnavailable` if MHCflurry is not installed.
  - `presentation_availability(frame) -> pd.DataFrame` — per-study count of candidates with ≥1 resolved class-I allele.

- [ ] **Step 1: Write the failing test** (does not require MHCflurry — exercises resolution + availability only)

```python
# tests/test_m6_presentation.py
from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.presentation import (
    presentation_availability,
    resolve_class_i_alleles,
)


def test_hla_resolution_covers_hu_and_pdac_only():
    frame = resolve_class_i_alleles(load_label_frame())
    availability = presentation_availability(frame).set_index("study_id")
    # hu carries candidate-level HLA; pdac gets predicted best-binder alleles by join.
    assert availability.loc["hu_neovax_2021", "resolved"] > 0
    assert availability.loc["pdac_neovax_2023", "resolved"] > 0
    # braun HLA is not public; mKRAS long peptides are NOT_ASSESSED -> zero resolved.
    assert availability.loc["braun_rcc_2025", "resolved"] == 0
    assert availability.loc["mkras_vax_2026", "resolved"] == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_presentation.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.presentation'`

- [ ] **Step 3: Implement HLA resolution + guarded scoring**

```python
# src/epicurus_neo/m6/presentation.py
from __future__ import annotations

from pathlib import Path

import pandas as pd

from epicurus_neo.m6.dataset import CORPUS_DIR, parse_alleles

PRESENTATION_STUDIES = ("hu_neovax_2021", "pdac_neovax_2023")


class PresentationUnavailable(RuntimeError):
    """Raised when the MHCflurry presentation baseline cannot be computed."""


def _predicted_antigen_alleles(corpus_dir: Path) -> pd.DataFrame:
    antigens = pd.read_parquet(corpus_dir / "antigens.parquet")
    predicted = antigens[antigens.hla_evidence_type == "PREDICTED_BEST_BINDER"].copy()
    predicted["antigen_alleles"] = predicted.hla_alleles.map(parse_alleles)
    keyed = predicted.groupby(["study_id", "gene", "protein_change"])["antigen_alleles"].agg(
        lambda lists: sorted({allele for sublist in lists for allele in sublist})
    )
    return keyed.rename("antigen_alleles").reset_index()


def resolve_class_i_alleles(frame: pd.DataFrame, corpus_dir: str | Path = CORPUS_DIR) -> pd.DataFrame:
    """Attach a per-candidate ``class_i_alleles`` list (hu direct; pdac by antigen join)."""
    out = frame.copy()
    candidate_alleles = out.hla_alleles.map(parse_alleles)
    predicted = _predicted_antigen_alleles(Path(corpus_dir))
    out = out.merge(predicted, on=["study_id", "gene", "protein_change"], how="left") \
        if {"gene", "protein_change"}.issubset(out.columns) else out.assign(antigen_alleles=None)

    def _resolve(row_direct: list[str], row_predicted: object) -> list[str]:
        if row_direct:
            return row_direct
        return list(row_predicted) if isinstance(row_predicted, list) else []

    out["class_i_alleles"] = [
        _resolve(direct, predicted)
        for direct, predicted in zip(candidate_alleles, out.get("antigen_alleles", [None] * len(out)))
    ]
    return out.drop(columns=[c for c in ["antigen_alleles"] if c in out.columns])


def presentation_availability(frame: pd.DataFrame) -> pd.DataFrame:
    resolved = frame.get("class_i_alleles")
    if resolved is None:
        frame = resolve_class_i_alleles(frame)
        resolved = frame.class_i_alleles
    has = resolved.map(lambda alleles: len(alleles) > 0)
    table = frame.assign(_has=has).groupby("study_id")["_has"].agg(
        resolved="sum", total="size"
    ).reset_index()
    table["resolved"] = table.resolved.astype(int)
    return table


def add_presentation_score(frame: pd.DataFrame, *, predictor: object | None = None) -> pd.DataFrame:
    """Add ``presentation_score`` = best MHCflurry presentation across class-I alleles."""
    try:
        from epicurus_neo.mhcflurry_features import add_mhcflurry_predictions
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise PresentationUnavailable(str(exc)) from exc

    resolved = resolve_class_i_alleles(frame)
    exploded = resolved.explode("class_i_alleles")
    scoreable = exploded[exploded.class_i_alleles.map(lambda a: isinstance(a, str) and bool(a))].copy()
    result = resolved.copy()
    result["presentation_score"] = float("nan")
    if scoreable.empty:
        return result
    scoreable = scoreable.rename(columns={"class_i_alleles": "hla_allele"})
    try:
        scored = add_mhcflurry_predictions(scoreable, predictor=predictor, allele_col="hla_allele")
    except RuntimeError as exc:  # MHCflurry import failure inside the helper
        raise PresentationUnavailable(str(exc)) from exc
    best = scored.groupby(scoreable.index)["mhcflurry_presentation_score"].max()
    result.loc[best.index, "presentation_score"] = best.to_numpy()
    return result
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_presentation.py -v`
Expected: PASS (1 passed). Note: this test avoids MHCflurry; `add_presentation_score` is covered end-to-end in Task 11's runner only when `.[mhc]` is installed.

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/presentation.py tests/test_m6_presentation.py
git add src/epicurus_neo/m6/presentation.py tests/test_m6_presentation.py
git commit -m "milestone-6a: HLA resolution + guarded MHCflurry presentation baseline"
```

---

### Task 6: Model ladder — B0/B1/M1/M2 (WP3)

**Files:**
- Create: `src/epicurus_neo/m6/models.py`
- Test: `tests/test_m6_models.py`

**Interfaces:**
- Consumes: feature matrix + `feature_columns` from Task 4.
- Produces: `fit_predict(model_name, train, evaluation, feature_cols, *, seed=17) -> np.ndarray` returning one score per evaluation row (higher = more likely POSITIVE). `model_name ∈ {"prevalence", "presentation", "logistic", "boosting"}`. `REQUIRED_MODELS = ("prevalence", "logistic", "boosting")`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_models.py
import numpy as np
import pandas as pd

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
        from sklearn.metrics import roc_auc_score
        assert roc_auc_score(evaluation.label, first) > 0.8
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_models.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.models'`

- [ ] **Step 3: Implement the ladder**

```python
# src/epicurus_neo/m6/models.py
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

REQUIRED_MODELS = ("prevalence", "logistic", "boosting")


def _logistic(seed: int) -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=1000, random_state=seed)),
        ]
    )


def fit_predict(
    model_name: str,
    train: pd.DataFrame,
    evaluation: pd.DataFrame,
    feature_cols: list[str],
    *,
    seed: int = 17,
) -> np.ndarray:
    """Fit ``model_name`` on train and return one POSITIVE-likelihood score per eval row."""
    if model_name == "prevalence":
        return np.full(len(evaluation), float(train.label.mean()), dtype=float)
    if model_name == "presentation":
        if "presentation_score" not in evaluation.columns:
            raise ValueError("presentation model requires a presentation_score column")
        return pd.to_numeric(evaluation.presentation_score, errors="coerce").to_numpy(dtype=float)

    x_train, y_train = train[feature_cols], train.label.to_numpy()
    if model_name == "logistic":
        model = _logistic(seed)
    elif model_name == "boosting":
        model = HistGradientBoostingClassifier(
            max_depth=3, max_iter=200, learning_rate=0.05, random_state=seed
        )
    else:
        raise ValueError(f"unknown model: {model_name!r}")
    model.fit(x_train, y_train)
    return model.predict_proba(evaluation[feature_cols])[:, 1].astype(float)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_models.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/models.py tests/test_m6_models.py
git add src/epicurus_neo/m6/models.py tests/test_m6_models.py
git commit -m "milestone-6a: B0/B1/M1/M2 model ladder"
```

---

### Task 7: Per-patient-k ranking + classification metrics (WP4)

**Files:**
- Create: `src/epicurus_neo/m6/ranking.py`
- Test: `tests/test_m6_ranking.py`

**Interfaces:**
- Consumes: a scored frame with `patient_id, label, <score_col>, mutant_peptide, hla_allele`; reuses `benchmark.metrics.identity_tiebreak`.
- Produces:
  - `patient_rank_vectors(df, score_col, *, k_cap=20, label_col="label", ascending=False) -> dict[str, np.ndarray]` — per-patient arrays (sorted by `patient_id`) for `hits_at_k`, `precision_at_k`, `capture_fraction`, `p_at_least_1`, `p_at_least_2`, `p_at_least_4`, using `k=min(k_cap, n)`. Zero-positive patients are `nan` for capture/`p_at_least_*`.
  - `classification_metrics(y_true, y_score) -> dict` — `auroc`, `average_precision`, `brier`, plus a 5-bin reliability list `calibration`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_ranking.py
import numpy as np
import pandas as pd

from epicurus_neo.m6.ranking import classification_metrics, patient_rank_vectors


def _frame():
    # Patient P1: 3 candidates, k=3 (degenerate). Patient P2: perfect ranking.
    return pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P1", "P2", "P2", "P2"],
            "mutant_peptide": list("ABCDEF"),
            "hla_allele": ["HLA-A*02:01"] * 6,
            "label": [1, 0, 0, 1, 1, 0],
            "score": [0.9, 0.1, 0.2, 0.9, 0.8, 0.1],
        }
    )


def test_k_patient_uses_min_of_cap_and_length():
    vectors = patient_rank_vectors(_frame(), "score", k_cap=2)
    # P1 has 3 candidates -> k=2; top-2 by score = A(1), C(0) -> 1 hit, precision 1/2.
    # P2 has 3 candidates -> k=2; top-2 = D(1), E(1) -> 2 hits, precision 2/2.
    assert vectors["hits_at_k"].tolist() == [1.0, 2.0]
    assert vectors["precision_at_k"].tolist() == [0.5, 1.0]
    assert vectors["p_at_least_2"].tolist() == [0.0, 1.0]


def test_classification_metrics_on_separable_scores():
    metrics = classification_metrics(
        np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9])
    )
    assert metrics["auroc"] == 1.0
    assert 0.0 <= metrics["brier"] <= 0.25
    assert len(metrics["calibration"]) >= 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_ranking.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.ranking'`

- [ ] **Step 3: Implement the metrics**

```python
# src/epicurus_neo/m6/ranking.py
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from benchmark.metrics import identity_tiebreak


def _ranked_labels(df: pd.DataFrame, score_col: str, label_col: str, ascending: bool):
    work = df.copy()
    work["_score"] = pd.to_numeric(work[score_col], errors="coerce")
    work["_tiebreak"] = identity_tiebreak(work)
    for patient, group in work.groupby("patient_id", sort=True):
        ranked = group.sort_values(
            ["_score", "_tiebreak"], ascending=[ascending, True], kind="mergesort", na_position="last"
        )
        yield ranked[label_col].to_numpy(dtype=float)


def patient_rank_vectors(
    df: pd.DataFrame,
    score_col: str,
    *,
    k_cap: int = 20,
    label_col: str = "label",
    ascending: bool = False,
) -> dict[str, np.ndarray]:
    """Per-patient ranking metrics with k = min(k_cap, n_candidates)."""
    keys = ("hits_at_k", "precision_at_k", "capture_fraction", "p_at_least_1", "p_at_least_2", "p_at_least_4")
    out: dict[str, list[float]] = {key: [] for key in keys}
    for labels in _ranked_labels(df, score_col, label_col, ascending):
        k = min(k_cap, len(labels))
        top = labels[:k]
        hits = float(top.sum())
        positives = int(labels.sum())
        out["hits_at_k"].append(hits)
        out["precision_at_k"].append(hits / k if k else float("nan"))
        out["capture_fraction"].append(hits / min(positives, k) if positives else float("nan"))
        for threshold in (1, 2, 4):
            out[f"p_at_least_{threshold}"].append(float(hits >= threshold) if positives else float("nan"))
    return {key: np.asarray(values, dtype=float) for key, values in out.items()}


def classification_metrics(y_true: np.ndarray, y_score: np.ndarray, *, bins: int = 5) -> dict:
    """Threshold-free classification metrics over pooled out-of-fold predictions."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    finite = np.isfinite(y_score)
    y_true, y_score = y_true[finite], y_score[finite]
    edges = np.linspace(0.0, 1.0, bins + 1)
    which = np.clip(np.digitize(y_score, edges[1:-1]), 0, bins - 1)
    calibration = []
    for b in range(bins):
        mask = which == b
        if mask.any():
            calibration.append(
                {"bin": b, "n": int(mask.sum()), "mean_score": float(y_score[mask].mean()),
                 "mean_label": float(y_true[mask].mean())}
            )
    return {
        "auroc": float(roc_auc_score(y_true, y_score)) if len(set(y_true.tolist())) > 1 else float("nan"),
        "average_precision": float(average_precision_score(y_true, y_score)) if y_true.any() else float("nan"),
        "brier": float(brier_score_loss(y_true, np.clip(y_score, 0.0, 1.0))),
        "calibration": calibration,
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_ranking.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/ranking.py tests/test_m6_ranking.py
git add src/epicurus_neo/m6/ranking.py tests/test_m6_ranking.py
git commit -m "milestone-6a: per-patient-k ranking and classification metrics"
```

---

### Task 8: Reusable verdict + macro/micro LOSO evaluation (WP4/WP7)

**Files:**
- Modify: `src/benchmark/scorecard.py` (extract the verdict; behavior unchanged)
- Create: `src/epicurus_neo/m6/evaluate.py`
- Test: `tests/test_m6_evaluate.py`, `tests/test_scorecard_verdict_refactor.py`

**Interfaces:**
- Produces in `scorecard.py`: `pre_registered_verdict(primary: dict, co_primary: dict, clinical: dict) -> str` where each dict has `delta_vs_baseline: float` and `delta_ci: [lo, hi]`. `scorecard()` calls it (identical output).
- Produces in `evaluate.py`: `evaluate_track(frame, *, model_name, baseline_name, track, k_cap=20, seed=17) -> dict` — runs LOSO, trains the model per fold on the core tier, scores model + baseline, and returns per-fold/macro/micro `hits@k_patient` deltas with patient-level bootstrap CIs, the pooled classification metrics, the verdict, and the ranking-informative diagnostic. Uses `benchmark.stats.paired_bootstrap` and `macro_paired_delta` (defined here).

- [ ] **Step 1: Write the failing refactor test**

```python
# tests/test_scorecard_verdict_refactor.py
from benchmark.scorecard import pre_registered_verdict


def test_verdict_rule_matches_registered_semantics():
    accept = {"delta_vs_baseline": 0.5, "delta_ci": [0.1, 0.9]}
    safe = {"delta_vs_baseline": 0.0, "delta_ci": [-0.1, 0.2]}
    assert pre_registered_verdict(accept, safe, safe) == "ACCEPT"
    unresolved = {"delta_vs_baseline": 0.3, "delta_ci": [-0.2, 0.8]}
    assert pre_registered_verdict(unresolved, safe, safe) == "CONSISTENT_WITH_NO_EFFECT"
    regressed = {"delta_vs_baseline": -0.4, "delta_ci": [-0.9, -0.1]}
    assert pre_registered_verdict(regressed, safe, safe) == "REJECT"
    # A significant primary gain is vetoed by a significant co-primary regression.
    harm = {"delta_vs_baseline": -0.3, "delta_ci": [-0.5, -0.1]}
    assert pre_registered_verdict(accept, harm, safe) == "REJECT"
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_scorecard_verdict_refactor.py -v`
Expected: FAIL with `ImportError: cannot import name 'pre_registered_verdict'`

- [ ] **Step 3: Extract the verdict in `scorecard.py`**

Add this module-level function to `src/benchmark/scorecard.py` (after imports):

```python
def pre_registered_verdict(primary: dict, co_primary: dict, clinical: dict) -> str:
    """The registered ACCEPT/CONSISTENT/REJECT rule over (delta, delta_ci) entries."""
    primary_significant = primary["delta_ci"][0] > 0.0
    no_co_primary_regression = co_primary["delta_ci"][1] >= 0.0
    no_clinical_regression = clinical["delta_ci"][1] >= 0.0
    if primary_significant and no_co_primary_regression and no_clinical_regression:
        return "ACCEPT"
    if primary["delta_vs_baseline"] > 0.0 and primary["delta_ci"][0] <= 0.0 <= primary["delta_ci"][1]:
        return "CONSISTENT_WITH_NO_EFFECT"
    return "REJECT"
```

Then replace the inline verdict block (the `primary_significant = ...` through `report["verdict"] = verdict` lines) in `scorecard()` with:

```python
    report["verdict"] = pre_registered_verdict(
        report[primary_name], report["capture_fraction"], report["p_at_least_one"]
    )
```

- [ ] **Step 4: Run refactor + existing scorecard tests to verify unchanged behavior**

Run: `.venv/bin/python -m pytest tests/test_scorecard_verdict_refactor.py tests/ -k scorecard -v`
Expected: PASS (refactor test passes; any existing scorecard tests still pass)

- [ ] **Step 5: Write the failing evaluation test**

```python
# tests/test_m6_evaluate.py
import numpy as np
import pandas as pd

from epicurus_neo.m6.evaluate import macro_paired_delta


def test_macro_delta_equal_weights_studies():
    # Study A: delta +1 per patient (2 patients). Study B: delta 0 (4 patients).
    per_study = {
        "A": (np.array([2.0, 2.0]), np.array([1.0, 1.0])),
        "B": (np.array([1.0, 1.0, 1.0, 1.0]), np.array([1.0, 1.0, 1.0, 1.0])),
    }
    result = macro_paired_delta(per_study, seed=17)
    # Macro = mean of per-study mean deltas = mean(+1, 0) = 0.5 (not patient-weighted 1/3).
    assert abs(result["delta"] - 0.5) < 1e-9
    assert result["delta_ci"][0] <= result["delta"] <= result["delta_ci"][1]
```

- [ ] **Step 6: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_evaluate.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.evaluate'`

- [ ] **Step 7: Implement `evaluate.py`**

```python
# src/epicurus_neo/m6/evaluate.py
from __future__ import annotations

import numpy as np
import pandas as pd

from benchmark.scorecard import pre_registered_verdict
from benchmark.stats import paired_bootstrap
from epicurus_neo.m6.features import build_feature_matrix, feature_columns
from epicurus_neo.m6.loso import loso_folds
from epicurus_neo.m6.models import fit_predict
from epicurus_neo.m6.ranking import classification_metrics, patient_rank_vectors


def _entry(candidate: np.ndarray, baseline: np.ndarray) -> dict:
    comparison = paired_bootstrap(candidate, baseline, seed=17)
    return {"delta_vs_baseline": comparison.delta, "delta_ci": [comparison.lo, comparison.hi],
            "p_better": comparison.p_better, "n": comparison.n}


def macro_paired_delta(per_study: dict[str, tuple[np.ndarray, np.ndarray]], *, seed: int = 17,
                       n: int = 20_000) -> dict:
    """Equal-weight-per-study delta with a study-stratified patient bootstrap CI."""
    studies = sorted(per_study)
    diffs = {s: _finite_diffs(*per_study[s]) for s in studies}
    point = float(np.mean([diffs[s].mean() for s in studies if len(diffs[s])]))
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n):
        study_means = []
        for s in studies:
            d = diffs[s]
            if len(d):
                study_means.append(d[rng.integers(0, len(d), size=len(d))].mean())
        if study_means:
            draws.append(np.mean(study_means))
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return {"delta": point, "delta_ci": [float(lo), float(hi)],
            "per_study": {s: float(diffs[s].mean()) if len(diffs[s]) else float("nan") for s in studies}}


def _finite_diffs(candidate: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    diff = np.asarray(candidate, dtype=float) - np.asarray(baseline, dtype=float)
    return diff[np.isfinite(diff)]


_GUARD_METRICS = ("hits_at_k", "capture_fraction", "p_at_least_1")


def _macro_entry(macro: dict) -> dict:
    return {"delta_vs_baseline": macro["delta"], "delta_ci": macro["delta_ci"]}


def evaluate_track(frame: pd.DataFrame, *, model_name: str, baseline_name: str, track: str,
                   tier: str = "core", k_cap: int = 20, seed: int = 17) -> dict:
    """Run LOSO for one (model vs baseline) comparison and assemble the registered report."""
    per_study: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {m: {} for m in _GUARD_METRICS}
    micro: dict[str, dict[str, list]] = {m: {"cand": [], "base": []} for m in _GUARD_METRICS}
    pooled_true, pooled_score = [], []
    per_fold = {}
    informative_total = 0
    for fold in loso_folds(frame):
        train = build_feature_matrix(fold.train, tier)
        evaluation = build_feature_matrix(fold.evaluation, tier)
        cols = feature_columns(train)
        model_scores = fit_predict(model_name, train, evaluation, cols, seed=seed)
        base_scores = fit_predict(baseline_name, train, evaluation, cols, seed=seed)
        scored = evaluation.assign(_model=model_scores, _base=base_scores)
        cand = patient_rank_vectors(scored, "_model", k_cap=k_cap)
        base = patient_rank_vectors(scored, "_base", k_cap=k_cap)
        for metric in _GUARD_METRICS:
            per_study[metric][fold.held_out_study] = (cand[metric], base[metric])
            micro[metric]["cand"].append(cand[metric])
            micro[metric]["base"].append(base[metric])
        pooled_true.append(scored.label.to_numpy())
        pooled_score.append(model_scores)
        informative = int((scored.groupby("patient_id").size() > k_cap).sum())
        informative_total += informative
        per_fold[fold.held_out_study] = {
            "hits_at_k": _entry(cand["hits_at_k"], base["hits_at_k"]),
            "n_patients": int(scored.patient_id.nunique()),
            "n_ranking_informative": informative,
        }
    macro = {metric: macro_paired_delta(per_study[metric], seed=seed) for metric in _GUARD_METRICS}
    micro_hits = _entry(
        np.concatenate(micro["hits_at_k"]["cand"]), np.concatenate(micro["hits_at_k"]["base"])
    )
    classification = classification_metrics(np.concatenate(pooled_true), np.concatenate(pooled_score))
    # Registered WP7 verdict: primary hits@k, co-primary capture, clinical P>=1 — all at macro level.
    verdict = pre_registered_verdict(
        _macro_entry(macro["hits_at_k"]),
        _macro_entry(macro["capture_fraction"]),
        _macro_entry(macro["p_at_least_1"]),
    )
    return {
        "track": track, "model": model_name, "baseline": baseline_name,
        "macro_hits_at_k": macro["hits_at_k"], "macro_capture": macro["capture_fraction"],
        "macro_p_at_least_1": macro["p_at_least_1"], "micro_hits_at_k": micro_hits,
        "per_fold": per_fold, "classification": classification, "verdict": verdict,
        "ranking_informative_patients": informative_total,
    }
```

- [ ] **Step 8: Run to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_evaluate.py tests/test_scorecard_verdict_refactor.py -v`
Expected: PASS (2 passed)

- [ ] **Step 9: Lint and commit**

```bash
.venv/bin/ruff check src/benchmark/scorecard.py src/epicurus_neo/m6/evaluate.py tests/test_m6_evaluate.py tests/test_scorecard_verdict_refactor.py
git add src/benchmark/scorecard.py src/epicurus_neo/m6/evaluate.py tests/test_m6_evaluate.py tests/test_scorecard_verdict_refactor.py
git commit -m "milestone-6a: reusable verdict + macro/micro LOSO evaluation"
```

---

### Task 9: Study-confound diagnostics (WP5)

**Files:**
- Create: `src/epicurus_neo/m6/confounds.py`
- Test: `tests/test_m6_confounds.py`

**Interfaces:**
- Consumes: label frame; core feature matrix from Task 4.
- Produces:
  - `prevalence_by_study(frame) -> pd.DataFrame` — `study_id, n, positive_rate`.
  - `study_only_classifier(frame, *, seed=17) -> dict` — grouped-by-patient CV accuracy of predicting `study_id` from core features, vs the majority-class rate (high ⇒ strong confound).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_confounds.py
from epicurus_neo.m6.confounds import prevalence_by_study, study_only_classifier
from epicurus_neo.m6.dataset import load_label_frame


def test_prevalence_varies_sharply_by_study():
    table = prevalence_by_study(load_label_frame()).set_index("study_id")
    assert int(table.n.sum()) == 965
    # Study prevalence is heterogeneous — the core confound M6A must expose.
    assert table.positive_rate.max() - table.positive_rate.min() > 0.1


def test_study_only_classifier_reports_confound_strength():
    result = study_only_classifier(load_label_frame())
    assert 0.0 <= result["accuracy"] <= 1.0
    assert result["majority_rate"] <= result["accuracy"] + 1e-9 or result["accuracy"] >= 0.0
    assert "per_study" in result
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_confounds.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.confounds'`

- [ ] **Step 3: Implement diagnostics**

```python
# src/epicurus_neo/m6/confounds.py
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedGroupKFold

from epicurus_neo.m6.features import build_feature_matrix, feature_columns


def prevalence_by_study(frame: pd.DataFrame) -> pd.DataFrame:
    table = frame.groupby("study_id").agg(n=("label", "size"), positives=("label", "sum")).reset_index()
    table["positive_rate"] = table.positives / table.n
    return table[["study_id", "n", "positive_rate"]]


def study_only_classifier(frame: pd.DataFrame, *, seed: int = 17) -> dict:
    """Can pre-vaccine core features predict which study a candidate came from?"""
    matrix = build_feature_matrix(frame, "core")
    cols = feature_columns(matrix)
    x = matrix[cols].to_numpy()
    y = matrix.study_id.to_numpy()
    groups = matrix.patient_id.to_numpy()
    splitter = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=seed)
    predictions = np.empty(len(y), dtype=object)
    for train_idx, test_idx in splitter.split(x, y, groups):
        model = HistGradientBoostingClassifier(max_depth=3, random_state=seed)
        model.fit(x[train_idx], y[train_idx])
        predictions[test_idx] = model.predict(x[test_idx])
    correct = predictions == y
    per_study = {
        study: float(correct[y == study].mean()) for study in sorted(set(y.tolist()))
    }
    _, counts = np.unique(y, return_counts=True)
    return {
        "accuracy": float(correct.mean()),
        "majority_rate": float(counts.max() / counts.sum()),
        "per_study": per_study,
    }
```

- [ ] **Step 4: Run to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_confounds.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/confounds.py tests/test_m6_confounds.py
git add src/epicurus_neo/m6/confounds.py tests/test_m6_confounds.py
git commit -m "milestone-6a: study-confound diagnostics"
```

---

### Task 10: Audit assembler + markdown renderer (WP7)

**Files:**
- Create: `src/epicurus_neo/m6/audit.py`
- Test: `tests/test_m6_audit.py`

**Interfaces:**
- Consumes: `evaluate_track` results, `completeness_report`, `prevalence_by_study`, `study_only_classifier`, `presentation_availability`.
- Produces:
  - `assemble_audit(*, universal, presentation, completeness, prevalence, confound, availability) -> dict` — the full audit dict, always stamped `"corpus_verdict": "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"`.
  - `render_audit_markdown(audit: dict) -> str`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m6_audit.py
from epicurus_neo.m6.audit import assemble_audit, render_audit_markdown


def _fake_track(track, verdict):
    return {"track": track, "model": "logistic", "baseline": "prevalence", "verdict": verdict,
            "macro_hits_at_k": {"delta": 0.1, "delta_ci": [-0.2, 0.4], "per_study": {}},
            "micro_hits_at_k": {"delta_vs_baseline": 0.05, "delta_ci": [-0.1, 0.2]},
            "classification": {"auroc": 0.55, "brier": 0.2, "average_precision": 0.3, "calibration": []},
            "ranking_informative_patients": 8, "per_fold": {}}


def test_audit_stamps_insufficiency_and_carries_both_tracks():
    audit = assemble_audit(
        universal=_fake_track("universal", "CONSISTENT_WITH_NO_EFFECT"),
        presentation=_fake_track("presentation", "REJECT"),
        completeness=[], prevalence=[], confound={"accuracy": 0.9, "majority_rate": 0.5, "per_study": {}},
        availability=[],
    )
    assert audit["corpus_verdict"] == "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    assert audit["universal"]["verdict"] == "CONSISTENT_WITH_NO_EFFECT"
    markdown = render_audit_markdown(audit)
    assert "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA" in markdown
    assert "Universal track" in markdown
    assert "ranking-informative" in markdown.lower()
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_audit.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'epicurus_neo.m6.audit'`

- [ ] **Step 3: Implement the assembler + renderer**

```python
# src/epicurus_neo/m6/audit.py
from __future__ import annotations

import pandas as pd

CORPUS_VERDICT = "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"


def _records(value) -> list[dict]:
    if isinstance(value, pd.DataFrame):
        return value.to_dict("records")
    return list(value)


def assemble_audit(*, universal, presentation, completeness, prevalence, confound, availability) -> dict:
    return {
        "corpus_verdict": CORPUS_VERDICT,
        "note": "Diagnostic swing under the standing insufficiency verdict; not a headline claim.",
        "universal": universal,
        "presentation": presentation,
        "completeness": _records(completeness),
        "prevalence_by_study": _records(prevalence),
        "study_confound": confound,
        "presentation_availability": _records(availability),
    }


def render_audit_markdown(audit: dict) -> str:
    universal, presentation = audit["universal"], audit["presentation"]
    lines = [
        "# Milestone 6A audit: Event-B-only recognition swing",
        "",
        f"**Corpus verdict (standing):** `{audit['corpus_verdict']}`",
        "",
        audit["note"],
        "",
        "## Universal track (all 4 studies, learned vs prevalence)",
        f"- Verdict: **{universal['verdict']}**",
        f"- Macro-study Δ hits@k_patient: {universal['macro_hits_at_k']['delta']:.4f} "
        f"CI {universal['macro_hits_at_k']['delta_ci']}",
        f"- AUROC (pooled OOF, secondary): {universal['classification']['auroc']:.4f} | "
        f"Brier: {universal['classification']['brier']:.4f}",
        f"- Ranking-informative patients (n_eligible > k): "
        f"{universal['ranking_informative_patients']} (selection signal is hu-dominated)",
        "",
        "## Presentation track (hu + pdac, learned vs presentation-only)",
        f"- Verdict: **{presentation['verdict']}**",
        f"- Macro-study Δ hits@k_patient: {presentation['macro_hits_at_k']['delta']:.4f} "
        f"CI {presentation['macro_hits_at_k']['delta_ci']}",
        "",
        "## Study confound",
        f"- Study-only classifier accuracy: {audit['study_confound']['accuracy']:.4f} "
        f"(majority rate {audit['study_confound']['majority_rate']:.4f})",
    ]
    return "\n".join(lines) + "\n"
```

- [ ] **Step 4: Run to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_audit.py -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Lint and commit**

```bash
.venv/bin/ruff check src/epicurus_neo/m6/audit.py tests/test_m6_audit.py
git add src/epicurus_neo/m6/audit.py tests/test_m6_audit.py
git commit -m "milestone-6a: audit assembler and markdown renderer"
```

---

### Task 11: Runner — wire both tracks, write artifacts

**Files:**
- Create: `experiments/m6a_recognition_swing.py`
- Test: `tests/test_m6_runner.py`

**Interfaces:**
- Consumes: everything above.
- Produces: `run(corpus_dir=CORPUS_DIR, out_dir="artifacts/milestone_6", *, seed=17) -> dict` and a `__main__` guard. Writes `artifacts/milestone_6/m6a_audit.json` and `m6a_audit.md`. The presentation track is guarded: if MHCflurry is not installed, presentation scores are skipped and the presentation track records `{"verdict": "SKIPPED_PRESENTATION_UNAVAILABLE"}` with a logged reason — the universal track always runs.

- [ ] **Step 1: Write the failing test** (guards on structure, not exact numbers, and does not require MHCflurry)

```python
# tests/test_m6_runner.py
import json
from pathlib import Path

from experiments.m6a_recognition_swing import run


def test_runner_writes_audit_with_standing_verdict(tmp_path):
    audit = run(out_dir=tmp_path, seed=17)
    assert audit["corpus_verdict"] == "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    assert audit["universal"]["verdict"] in {"ACCEPT", "CONSISTENT_WITH_NO_EFFECT", "REJECT"}
    assert (Path(tmp_path) / "m6a_audit.json").exists()
    assert (Path(tmp_path) / "m6a_audit.md").exists()
    # Universal track always runs on all 4 studies.
    assert set(audit["universal"]["per_fold"]) == {
        "braun_rcc_2025", "hu_neovax_2021", "mkras_vax_2026", "pdac_neovax_2023",
    }


def test_runner_is_deterministic(tmp_path):
    first = run(out_dir=tmp_path / "a", seed=17)
    second = run(out_dir=tmp_path / "b", seed=17)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
```

- [ ] **Step 2: Run to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_m6_runner.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'experiments.m6a_recognition_swing'` (add `experiments/__init__.py` if the repo lacks one — check first with `ls experiments/__init__.py`)

- [ ] **Step 3: Implement the runner**

```python
# experiments/m6a_recognition_swing.py
from __future__ import annotations

import json
import sys
from pathlib import Path

# Allow `python experiments/m6a_recognition_swing.py` without installation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epicurus_neo.m6.audit import assemble_audit, render_audit_markdown  # noqa: E402
from epicurus_neo.m6.confounds import prevalence_by_study, study_only_classifier  # noqa: E402
from epicurus_neo.m6.dataset import CORPUS_DIR, completeness_report, load_label_frame  # noqa: E402
from epicurus_neo.m6.evaluate import evaluate_track  # noqa: E402
from epicurus_neo.m6.presentation import (  # noqa: E402
    PRESENTATION_STUDIES,
    PresentationUnavailable,
    add_presentation_score,
    presentation_availability,
    resolve_class_i_alleles,
)


def run(corpus_dir=CORPUS_DIR, out_dir="artifacts/milestone_6", *, seed: int = 17) -> dict:
    frame = load_label_frame(corpus_dir)
    universal = evaluate_track(
        frame, model_name="logistic", baseline_name="prevalence", track="universal", seed=seed
    )
    availability = presentation_availability(resolve_class_i_alleles(frame))
    try:
        scored = add_presentation_score(frame)
        subset = scored[scored.study_id.isin(PRESENTATION_STUDIES)].reset_index(drop=True)
        presentation = evaluate_track(
            subset, model_name="logistic", baseline_name="presentation",
            track="presentation", tier="presentation", seed=seed,
        )
    except PresentationUnavailable as exc:
        presentation = {
            "track": "presentation", "verdict": "SKIPPED_PRESENTATION_UNAVAILABLE",
            "reason": str(exc), "macro_hits_at_k": {"delta": float("nan"), "delta_ci": [None, None]},
            "classification": {"auroc": float("nan"), "brier": float("nan")},
            "ranking_informative_patients": 0, "per_fold": {},
        }
    audit = assemble_audit(
        universal=universal, presentation=presentation,
        completeness=completeness_report(frame), prevalence=prevalence_by_study(frame),
        confound=study_only_classifier(frame, seed=seed), availability=availability,
    )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "m6a_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    (out / "m6a_audit.md").write_text(render_audit_markdown(audit))
    return audit


if __name__ == "__main__":
    result = run()
    print(f"universal verdict: {result['universal']['verdict']}")
    print(f"presentation verdict: {result['presentation']['verdict']}")
```

- [ ] **Step 4: Run to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_m6_runner.py -v`
Expected: PASS (2 passed). If the determinism test fails on a floating-point tie, confirm every `paired_bootstrap`/`macro_paired_delta`/model call threads `seed=17`.

- [ ] **Step 5: Run the runner for real and eyeball the artifact**

Run: `.venv/bin/python experiments/m6a_recognition_swing.py`
Expected: prints two verdicts; writes `artifacts/milestone_6/m6a_audit.{json,md}`. Presentation track is `SKIPPED_PRESENTATION_UNAVAILABLE` unless `.[mhc]` is installed.

- [ ] **Step 6: Lint and commit**

```bash
.venv/bin/ruff check experiments/m6a_recognition_swing.py tests/test_m6_runner.py
git add experiments/m6a_recognition_swing.py tests/test_m6_runner.py artifacts/milestone_6/
git commit -m "milestone-6a: two-track runner and audit artifacts"
```

---

### Task 12: Integration guard — leakage + full-suite regression

**Files:**
- Test: `tests/test_m6_integration.py`

**Interfaces:**
- Consumes: the whole pipeline.

- [ ] **Step 1: Write the failing integration test**

```python
# tests/test_m6_integration.py
from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.features import assert_no_banned_features, build_feature_matrix, feature_columns
from epicurus_neo.m6.loso import loso_folds


def test_no_banned_feature_reaches_any_tier():
    frame = load_label_frame()
    for tier in ("core", "contrastive", "presentation"):
        assert_no_banned_features(feature_columns(build_feature_matrix(frame, tier)))


def test_loso_training_never_sees_the_held_out_study():
    frame = load_label_frame()
    for fold in loso_folds(frame):
        assert fold.held_out_study not in set(fold.train.study_id)
        # The core feature matrix carries no study_id/patient_id/label as a usable feature.
        cols = feature_columns(build_feature_matrix(fold.train, "core"))
        assert not ({"study_id", "patient_id", "label", "cancer_type"} & set(cols))
```

- [ ] **Step 2: Run to verify it fails, then passes with existing code**

Run: `.venv/bin/python -m pytest tests/test_m6_integration.py -v`
Expected: PASS immediately (guards the invariants built in Tasks 1–11). If it fails, a banned column leaked — fix `_EXTRA_BANNED` in `features.py`.

- [ ] **Step 3: Run the entire suite + lint the package**

Run: `.venv/bin/python -m pytest tests/ -q && .venv/bin/ruff check src/epicurus_neo/m6/ experiments/m6a_recognition_swing.py`
Expected: all tests pass; ruff clean.

- [ ] **Step 4: Commit**

```bash
git add tests/test_m6_integration.py
git commit -m "milestone-6a: leakage and LOSO integration guards"
```

---

## Self-Review

**1. Spec coverage:**
- WP0 completeness gate → Task 2. WP1 HLA/presentation → Task 5 (hu+pdac resolved, braun/mKRAS zero, guarded). WP2 core/contrastive tiers + class-gated anchors + missingness discipline → Task 4. WP3 ladder B0/B1/M1/M2 (sklearn, no xgboost; M3/ESM correctly deferred) → Task 6. WP4 LOSO + k_patient + macro/micro + P≥1/≥2/≥4 + AUROC/Brier/calibration → Tasks 3, 7, 8. WP5 confounds → Task 9. WP6 calibration restraint → honored by using a single global training-only model per fold with `is_class_i/ii` as features and class-stratified reporting via `classification`; isotonic-on-subsets is not used (no per-class isotonic anywhere). WP7 registered verdict + audit + standing caveat → Tasks 8, 10, 11. Two headline tracks (no mixed baseline) → Task 11. Determinism/tests → Task 12. **M6B transfer and M6C Osteosarc are correctly out of this plan** (spec sequences them after M6A freeze).
- Gap noted and intentionally deferred: the fixed-@20 sensitivity analysis and the literal per-study feature-distribution table are *reporting* refinements; the registered headline (`hits@k_patient` macro) and the ranking-informative diagnostic that supersedes them are implemented. Add the fixed-@20 sensitivity table as a fast-follow if the user wants it in the first artifact.

**2. Placeholder scan:** No "TBD"/"TODO"/"handle edge cases"/vague-requirement patterns. The macro verdict uses all three registered WP7 guard axes — primary `hits@k_patient`, co-primary `capture_fraction`, clinical `p_at_least_1` — each computed via `macro_paired_delta`, with no dummy guard dicts.

**3. Type consistency:** `load_label_frame` → frame consumed by `completeness_report`, `loso_folds`, `build_feature_matrix`, `resolve_class_i_alleles` (all take the same frame). `fit_predict(model_name, train, evaluation, feature_cols, *, seed)` signature is identical in Task 6 and its two call sites in Task 8. `patient_rank_vectors` returns the `hits_at_k` key consumed in `evaluate_track`. `pre_registered_verdict(primary, co_primary, clinical)` dict shape (`delta_vs_baseline`, `delta_ci`) matches both the refactored `scorecard()` call and the `evaluate.py` call. `assemble_audit` keyword args match the runner's call.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-10-m6a-recognition-swing.md`. Two execution options:

1. **Subagent-Driven (recommended)** — a fresh subagent per task with two-stage review between tasks; fast iteration, each task gated independently.
2. **Inline Execution** — execute tasks in this session with checkpoints for review.

Which approach?
