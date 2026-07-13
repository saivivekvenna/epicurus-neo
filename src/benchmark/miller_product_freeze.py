"""Label-blind, patient-agnostic product selections for Miller universes."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from benchmark.end_to_end_product import assert_label_free
from benchmark.universal_portfolio import (
    select_evidence_lane_portfolio,
    select_rank_fusion_cap1,
)
from epicurus_neo.product import InferenceConfig, normalize_product_candidates, score_product_candidates


K = 20
POLICY_ID = "miller-universal-product-portfolios-v1"
ARM_IDS = (
    "shipped_epicurus_product",
    "epicurus_plain",
    "epicurus_mutation_cap1",
    "prime_plain",
    "prime_mutation_cap1",
    "rank_fusion_cap1",
    "evidence_lane_portfolio",
)
FUSION_COLUMNS = ("recognition_score", "epicurus_lower_evidence_score")
LANE_COLUMNS = (*FUSION_COLUMNS, "presented_evidence_score")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tie_key(frame: pd.DataFrame) -> pd.Series:
    raw = frame["mutant_peptide"].astype(str) + "|" + frame["hla_allele"].astype(str)
    return raw.map(lambda value: hashlib.md5(value.encode()).hexdigest())


def adapt_universe(raw: pd.DataFrame) -> pd.DataFrame:
    """Adapt a frozen lossless universe without inventing recognition evidence."""
    assert_label_free(raw)
    adapted = normalize_product_candidates(raw, source_name="miller_lossless_raw_reconstruction")
    # PRIME ranks in these universes come from the pinned genuine PRIME binary. Keep an explicit
    # availability bit so a product fallback score can never be mistaken for genuine PRIME evidence.
    adapted["genuine_prime_available"] = (
        pd.to_numeric(raw["prime_rank"], errors="coerce").notna().to_numpy()
        if "prime_rank" in raw
        else False
    )
    adapted["frozen_epicurus_score"] = (
        pd.to_numeric(raw["epicurus"], errors="coerce").to_numpy()
        if "epicurus" in raw
        else float("nan")
    )
    assert_label_free(adapted)
    return adapted


def _ordered(frame: pd.DataFrame, score_col: str, *, cap1: bool, k: int) -> pd.DataFrame:
    work = frame.copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work[work[score_col].notna()].copy()
    work["_tie"] = _tie_key(work)
    work = work.sort_values([score_col, "_tie"], ascending=[False, True], kind="mergesort")
    if cap1:
        work = work.drop_duplicates("mutation_id", keep="first")
    return work.head(k).drop(columns="_tie").reset_index(drop=True)


def build_selections(raw: pd.DataFrame, *, k: int = K) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Score once, then construct six pre-registered ordered portfolios without labels."""
    adapted = adapt_universe(raw)
    scored = score_product_candidates(adapted, InferenceConfig(k=k))
    assert_label_free(scored)
    valid = scored[scored["deterministic_gate_pass"].astype(bool)].copy()

    shipped = scored[scored["selected"].astype(bool)].sort_values("rank", kind="mergesort").head(k)
    # The generalization target is the shipped end-to-end product.  The legacy v0.1
    # research score remains in the frozen universe for provenance/diagnostics, but it
    # must not silently stand in for the product score in the preregistered Epicurus arms.
    epicurus_plain = _ordered(valid, "epicurus_lower_evidence_score", cap1=False, k=k)
    epicurus_cap1 = _ordered(valid, "epicurus_lower_evidence_score", cap1=True, k=k)

    prime_pool = valid[valid["genuine_prime_available"].astype(bool)].copy()
    prime_plain = _ordered(prime_pool, "recognition_score", cap1=False, k=k)
    prime_cap1 = _ordered(prime_pool, "recognition_score", cap1=True, k=k)

    portfolio_pool = valid.copy()
    portfolio_pool.loc[
        ~portfolio_pool["genuine_prime_available"].astype(bool), "recognition_score"
    ] = float("nan")
    portfolio_pool.loc[
        ~portfolio_pool["presented_evidence_available"].astype(bool), "presented_evidence_score"
    ] = float("nan")
    fusion = select_rank_fusion_cap1(portfolio_pool, FUSION_COLUMNS, k=k)
    lanes = select_evidence_lane_portfolio(portfolio_pool, LANE_COLUMNS, k=k)
    selections = {
        "shipped_epicurus_product": shipped.reset_index(drop=True),
        "epicurus_plain": epicurus_plain,
        "epicurus_mutation_cap1": epicurus_cap1,
        "prime_plain": prime_plain,
        "prime_mutation_cap1": prime_cap1,
        "rank_fusion_cap1": fusion,
        "evidence_lane_portfolio": lanes,
    }
    for selected in selections.values():
        assert_label_free(selected)
    return scored, selections


def freeze_product_selections(raw: pd.DataFrame, freeze_dir: Path, *, k: int = K) -> dict:
    """Write ordered product selections and return manifest metadata for the enclosing freeze."""
    if raw.empty:
        scored = pd.DataFrame(
            columns=(*LANE_COLUMNS, "frozen_epicurus_score", "genuine_prime_available",
                     "translated_evidence_available",
                     "presented_evidence_available", "recognized_evidence_available",
                     "coverage_evidence_available")
        )
        selections = {
            arm_id: pd.DataFrame(columns=("candidate_id", "mutation_id")) for arm_id in ARM_IDS
        }
    else:
        scored, selections = build_selections(raw, k=k)
    freeze_dir.mkdir(parents=True, exist_ok=True)
    arms: dict[str, dict] = {}
    hashes: dict[str, str] = {}
    for arm_id in ARM_IDS:
        selected = selections[arm_id].copy()
        selected.insert(0, "selection_rank", range(1, len(selected) + 1))
        filename = f"product_select_{arm_id}.csv"
        path = freeze_dir / filename
        selected.to_csv(path, index=False)
        hashes[filename] = _sha256(path)
        arms[arm_id] = {
            "selection_file": filename,
            "n_selected": int(len(selected)),
            "n_unique_mutations": int(selected["mutation_id"].astype(str).nunique()),
            "saturated": bool(len(selected) >= k),
            "ordered_candidate_ids": selected["candidate_id"].astype(str).tolist(),
            "ordered_mutation_ids": selected["mutation_id"].astype(str).tolist(),
        }

    availability = {
        name: int(scored[column].astype(bool).sum())
        for name, column in {
            "translated": "translated_evidence_available",
            "presented": "presented_evidence_available",
            "recognized": "recognized_evidence_available",
            "coverage": "coverage_evidence_available",
        }.items()
    }
    availability["genuine_prime_available"] = int(scored["genuine_prime_available"].astype(bool).sum())
    availability["frozen_epicurus_available"] = int(
        pd.to_numeric(scored["frozen_epicurus_score"], errors="coerce").notna().sum()
    )
    availability["shipped_epicurus_score_available"] = int(
        pd.to_numeric(scored["epicurus_lower_evidence_score"], errors="coerce").notna().sum()
    )
    return {
        "policy_id": POLICY_ID,
        "k": k,
        "labels_opened": False,
        "prime_provenance": "genuine PRIME percentile rank from frozen universe prime_rank; converted once to higher-is-better recognition_score",
        "epicurus_arm_score": "shipped epicurus_lower_evidence_score; legacy frozen_epicurus_score is diagnostic only",
        "shipped_product_arm": "shipped_epicurus_product",
        "preregistered_arm_ids": [arm for arm in ARM_IDS if arm != "shipped_epicurus_product"],
        "feature_availability_rows": availability,
        "arms": arms,
        "sha256": hashes,
    }
