"""NCI crosswalk — join the fragmented NCI representations into one source-grounded, MIL-aware view.

DATA AUDIT ONLY. No model fitting, no tuning, nothing touched for training here. This module reconciles three
NCI files that are the SAME cohort split across incompatible representations:

  * Gartner NmersTrainingSet.txt   — MUTATION-level (25mer) candidates, 70 patients, `Screening Status`
    (CD8 / '-' screened-negative / unscreened) + expression & VAF deciles. This is the ascertainment (the
    "working outcomes"), measured at the mutation/minigene level.
  * MullerNCItrain min             — minimal peptide-HLA INSTANCES for 56 of those patients, with presentation
    features (Score_EL, binding, stability, agretopicity, quantification) and an exact instance label
    (`VALIDATED`). It has NO mutation key, so it joins to Gartner only by patient + substring containment.
  * MmpsTestingSet_extract.tsv     — the exact held-out TESTING minimal peptide-HLA table for 26 patients, WITH
    a stable parent `key`, an instance `minkey`, parent `Screening Status`, and instance `epitope status`.

Multiple-instance (MIL) discipline — the point of this module:
  * A MUTATION is a BAG. `Screening Status` is measured once per mutation (per `key`), NOT per minimal
    peptide-HLA. A screened-negative parent therefore makes each child a NEGATIVE_BAG_CHILD (bag-level
    negative), NOT an independently tested-negative instance. We NEVER expand one mutation-level assay into
    many independent experiments.
  * `VALIDATED` / `epitope status` are INSTANCE-level (an exact confirmed minimal peptide-HLA).
  * Every row keeps BOTH `mutation_family_label` (bag) and `instance_label` (instance) + an `ambiguity` field.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

GARTNER_TRAIN = Path("data/raw/gartner_nci/NmersTrainingSet.txt")
GARTNER_TEST_RAW = Path("data/raw/gartner_nci/NmersTestingSet.txt")           # raw exact-key ascertainment
GARTNER_TEST_EXTRACT = Path("data/raw/gartner_nci/MmpsTestingSet_extract.tsv")  # minimal instances (propagated)
MULLER_MIN = Path("data/raw/neoranking_mirror/MullerNCItrain.train-data.min.tsv")
CACHE = Path("data/raw/gartner_nci")

# Gartner TRAIN mutation-level Screening Status vocabulary -> bag ascertainment state.
STATUS_MAP = {"CD8": "CD8_POS", "-": "SCREENED_NEG", "unscreened": "UNSCREENED"}
# Gartner TEST raw (NmersTestingSet) uses a DIFFERENT vocabulary for the SAME states: 1=CD8+, 0=screened-neg.
TEST_RAW_STATUS_MAP = {"1": "CD8_POS", "0": "SCREENED_NEG", "unscreened": "UNSCREENED"}
MIN_LEN, MAX_LEN = 8, 12   # minimal-epitope window lengths to index for containment

# The assay/screening unit is the TRANSCRIPT-specific 25mer `key` (patient_variant_25mer), NOT the genomic
# mutation. Source proof: within genomic variant `4359|11:17097026-17097026 C>T` one transcript 25mer was
# screened CD8+ (…RKDKDA) while its sibling (…RKRAGF) is `unscreened` in raw NmersTestingSet — a per-genomic
# screen would carry an identical status. So `candidate_parent_key`/`mutation_key` (transcript) is the assay
# (bag) unit for weighting; `genomic_mutation_family_id` (patient + raw Variant key) clusters transcript
# isoforms and is used ONLY for leakage-blocking folds, never as the assay unit.
def _genomic_family_id(patient_id: str, variant_key: str) -> str:
    return f"{patient_id}|{variant_key}"


@dataclass
class Crosswalk:
    instances: pd.DataFrame
    audit: dict


# --------------------------------------------------------------------------------------------------
# source loaders (normalized, no imputation)
# --------------------------------------------------------------------------------------------------
def load_gartner_train() -> pd.DataFrame:
    g = pd.read_csv(GARTNER_TRAIN, sep="\t", dtype=str)
    patient = g["ID"].astype(str)
    variant_key = g["Variant key"].astype(str).str.strip()
    out = pd.DataFrame({
        "patient_id": patient,
        "mutation_key": g["key"].astype(str),               # transcript 25mer key = assay (bag) unit
        "genomic_variant_key": variant_key,                  # raw genomic locus (shared across isoforms)
        "genomic_mutation_family_id": (patient + "|" + variant_key),  # leakage-blocking unit
        "mut_epitope": g["Mut Epitope"].astype(str).str.upper(),
        "screening_status": g["Screening Status"].astype(str),
        "bag_state": g["Screening Status"].map(STATUS_MAP).fillna("UNKNOWN"),
        "expr_decile": pd.to_numeric(g["Gene Expression Decile for this sample(1=lowest expression-10=highest expression)"], errors="coerce"),
        "vaf_decile": pd.to_numeric(g["Exome VAF Decile"], errors="coerce"),
    })
    return out


def load_gartner_test_raw() -> pd.DataFrame:
    """Raw NmersTestingSet — the exact-key ascertainment source of truth for the held-out TEST set. Its
    transcript-key `Screening Status` (1/0/unscreened) is authoritative; the MmpsTestingSet_extract propagates a
    genomic-family positive onto unscreened sibling transcripts, so we reconcile the extract against THIS."""
    t = pd.read_csv(GARTNER_TEST_RAW, sep="\t", dtype=str)
    patient = t["ID"].astype(str)
    variant_key = t["Variant key"].astype(str).str.strip()
    status = t["Screening Status"].astype(str).str.strip()
    return pd.DataFrame({
        "patient_id": patient,
        "mutation_key": t["key"].astype(str),
        "genomic_variant_key": variant_key,
        "genomic_mutation_family_id": (patient + "|" + variant_key),
        "raw_screening_status": status,
        "raw_bag_state": status.map(TEST_RAW_STATUS_MAP).fillna("UNKNOWN"),
    })


def load_muller_min() -> pd.DataFrame:
    m = pd.read_csv(MULLER_MIN, sep="\t", dtype=str)
    m["patient_id"] = m["PatientID"].str.replace("Muller-NCI-train_", "", regex=False)
    return pd.DataFrame({
        "patient_id": m["patient_id"].astype(str),
        "peptide": m["MT_pep_x"].astype(str).str.upper(),
        "hla_allele": m["HLA_type_x"].astype(str),
        "score_el": pd.to_numeric(m["Score_EL"], errors="coerce"),
        "bind_aff": pd.to_numeric(m["MT_BindAff"], errors="coerce"),
        "bind_stab": pd.to_numeric(m["BindStab"], errors="coerce"),
        "agretopicity": pd.to_numeric(m["Agretopicity"], errors="coerce"),
        "quantification": pd.to_numeric(m["Quantification"], errors="coerce"),
        "ln_num_tested": pd.to_numeric(m["ln_NumTested"], errors="coerce"),
        "validated": (m["VALIDATED"].astype(str) == "1").astype(int),
    })


# --------------------------------------------------------------------------------------------------
# TRAIN crosswalk: Gartner mutation ascertainment  x  Müller minimal instances (containment)
# --------------------------------------------------------------------------------------------------
def _substring_parent_index(train: pd.DataFrame) -> dict:
    """(patient, minimal-peptide) -> list of (transcript_key, genomic_family_id, bag_state). Built by
    enumerating all 8..12mer windows of each Gartner 25mer parent per patient. A minimal peptide may hit several
    parents -> we keep ALL candidate parent links (ambiguity is preserved, never silently resolved). We carry
    BOTH the transcript key (assay/bag unit) and the genomic mutation family (leakage-blocking unit)."""
    idx: dict = defaultdict(list)
    for pid, ep, key, gfam, state in zip(train["patient_id"], train["mut_epitope"], train["mutation_key"],
                                         train["genomic_mutation_family_id"], train["bag_state"]):
        ep = str(ep)
        seen = set()
        for k in range(MIN_LEN, MAX_LEN + 1):
            for s in range(len(ep) - k + 1):
                sub = ep[s:s + k]
                if sub in seen:
                    continue
                seen.add(sub)
                idx[(pid, sub)].append((key, gfam, state))
    return idx


def _instance_label(validated: int, states: set) -> tuple[str, str]:
    """Return (instance_label, ambiguity) per the MIL rules. `states` = set of bag states of candidate parents.

    Instance is confirmed positive (VALIDATED) -> POSITIVE_EXACT. Otherwise the label is BAG-derived and its
    confidence is bag-level, never instance-level:
      * only screened-negative parents  -> NEGATIVE_BAG_CHILD (bag-level negative; NOT an independent test)
      * any CD8-positive parent         -> AMBIGUOUS_POSITIVE_BAG (mutation is immunogenic; this instance was
                                            not the confirmed one -> cannot be called negative)
      * only unscreened parents         -> UNTESTED (never coerced to negative)
      * mixed screened/unscreened       -> MIXED_AMBIGUOUS
    """
    if validated == 1:
        return "POSITIVE_EXACT", "none"
    if states == {"SCREENED_NEG"}:
        return "NEGATIVE_BAG_CHILD", "bag_level_negative"
    if "CD8_POS" in states:
        return "AMBIGUOUS_POSITIVE_BAG", "positive_bag_unconfirmed_instance"
    if states == {"UNSCREENED"}:
        return "UNTESTED", "unscreened"
    return "MIXED_AMBIGUOUS", "mixed_parent_ascertainment"


def build_train_crosswalk(*, use_cache: bool = True) -> Crosswalk:
    cache = CACHE / "_nci_train_crosswalk.parquet"
    if use_cache and cache.exists():
        inst = pd.read_parquet(cache)
        return Crosswalk(inst, _train_audit(inst))

    train = load_gartner_train()
    mu = load_muller_min()
    idx = _substring_parent_index(train)

    (labels, ambig, family_labels, parent_keys, gfam_lists,
     n_tkeys, n_gfams, resolved_tkey, resolved_gfam) = ([] for _ in range(9))
    for pid, pep, val in zip(mu["patient_id"], mu["peptide"], mu["validated"]):
        parents = idx.get((pid, pep), [])
        states = {st for _, _, st in parents}
        tkeys = sorted({k for k, _, _ in parents})            # transcript (assay/bag) units
        gfams = sorted({gf for _, gf, _ in parents})          # genomic (leakage-blocking) units
        lab, amb = _instance_label(int(val), states)
        labels.append(lab)
        ambig.append(amb)
        family_labels.append("|".join(sorted(states)) if states else "NO_PARENT")
        parent_keys.append(";".join(tkeys))
        gfam_lists.append(";".join(gfams))
        n_tkeys.append(len(tkeys))
        n_gfams.append(len(gfams))
        # transcript key resolves only with exactly one candidate 25mer; genomic family can resolve even when
        # several transcript isoforms of the SAME locus match (isoforms share one genomic family).
        resolved_tkey.append(tkeys[0] if len(tkeys) == 1 else np.nan)
        resolved_gfam.append(gfams[0] if len(gfams) == 1 else np.nan)

    inst = mu.copy()
    inst["instance_label"] = labels
    inst["mutation_family_label"] = family_labels          # BAG-level ascertainment (may be multi-state)
    inst["ambiguity"] = ambig
    inst["n_candidate_parents"] = n_tkeys                   # transcript-key candidates (back-compat name)
    inst["n_candidate_genomic_families"] = n_gfams
    inst["candidate_parent_keys"] = parent_keys            # transcript 25mer keys (assay units)
    inst["candidate_genomic_families"] = gfam_lists
    # resolved single-parent ids: transcript for assay weighting, genomic family for leakage-blocking.
    inst["resolved_parent_key"] = resolved_tkey
    inst["resolved_genomic_family_id"] = resolved_gfam
    inst["cohort"] = "gartner_train_x_muller"
    if use_cache:
        inst.to_parquet(cache)
    return Crosswalk(inst, _train_audit(inst))


def _train_audit(inst: pd.DataFrame) -> dict:
    train = load_gartner_train()
    # transcript-key (assay/bag unit) ascertainment counts
    bag = train.groupby("bag_state")["mutation_key"].nunique().to_dict()
    # genomic-family (leakage-blocking unit) ascertainment counts (isoforms of one locus collapse)
    gfam_states = train.groupby("genomic_mutation_family_id")["bag_state"].agg(set)
    gfam_pos = int(gfam_states.apply(lambda s: "CD8_POS" in s).sum())
    gfam_neg = int(gfam_states.apply(lambda s: ("CD8_POS" not in s) and ("SCREENED_NEG" in s)).sum())
    gfam_unt = int(gfam_states.apply(lambda s: s <= {"UNSCREENED"}).sum())
    # effective independent negatives = distinct screened-negative BAGS that actually have Müller children
    neg_children = inst[inst["instance_label"] == "NEGATIVE_BAG_CHILD"]
    # a NEGATIVE_BAG_CHILD with exactly one candidate parent has a resolvable negative bag id
    neg_bags = neg_children.loc[neg_children["resolved_parent_key"].notna(), "resolved_parent_key"].nunique()
    neg_gfams = neg_children.loc[neg_children["resolved_genomic_family_id"].notna(),
                                 "resolved_genomic_family_id"].nunique()
    pos_patients = inst.loc[inst["instance_label"] == "POSITIVE_EXACT", "patient_id"].nunique()
    return {
        "instances_total": int(len(inst)),
        "muller_patients": int(inst["patient_id"].nunique()),
        "crosswalk_coverage_rows_with_a_parent": float((inst["n_candidate_parents"] > 0).mean()),
        "instance_label_counts": inst["instance_label"].value_counts().to_dict(),
        "positive_exact": int((inst["instance_label"] == "POSITIVE_EXACT").sum()),
        "positive_exact_map_to_CD8_parent": int(
            ((inst["instance_label"] == "POSITIVE_EXACT")
             & inst["mutation_family_label"].str.contains("CD8_POS")).sum()),
        "positive_patients": int(pos_patients),
        "ambiguity_counts": inst["ambiguity"].value_counts().to_dict(),
        "assay_unit": "transcript 25mer `key` (candidate_parent_key)",
        "leakage_block_unit": "genomic_mutation_family_id (patient + raw Variant key)",
        "gartner_train_bags_transcript_key": {"CD8_POS": int(bag.get("CD8_POS", 0)),
                                              "SCREENED_NEG": int(bag.get("SCREENED_NEG", 0)),
                                              "UNSCREENED": int(bag.get("UNSCREENED", 0))},
        "gartner_train_genomic_families": {"total": int(len(gfam_states)), "any_CD8_positive": gfam_pos,
                                           "any_screened_neg_no_pos": gfam_neg, "untested_only": gfam_unt},
        "EFFECTIVE_independent_negatives_note": "NEGATIVE_BAG_CHILD rows are children of mutation bags, NOT "
            "independent tested-negatives. Effective independent negative transcript BAGS (resolvable "
            f"single-parent) = {int(neg_bags)} (= {int(neg_gfams)} genomic families); do NOT count the "
            f"{int(len(neg_children))} child rows as independent experiments.",
        "effective_negative_bags_resolvable": int(neg_bags),
        "effective_negative_genomic_families_resolvable": int(neg_gfams),
    }


# --------------------------------------------------------------------------------------------------
# TEST crosswalk: exact, by `key` (NO fuzzy containment). Held out for model fitting.
# --------------------------------------------------------------------------------------------------
def build_test_crosswalk(*, use_cache: bool = True) -> Crosswalk:
    """Exact TEST crosswalk by transcript `key`. The MmpsTestingSet_extract PROPAGATES a genomic-family positive
    onto unscreened sibling transcripts (see the RPS13/4359 case), so its `Screening Status` is NOT a safe
    exact-key ascertainment source. We reconcile every key against the raw NmersTestingSet `Screening Status`
    (authoritative, transcript-key level): raw wins, the Mmps value is retained as `mmps_screening_status`, and a
    `source_conflict` flag marks keys where the two disagree (the extra sibling is kept UNTESTED, not positive)."""
    cache = CACHE / "_nci_test_crosswalk.parquet"
    if use_cache and cache.exists():
        inst = pd.read_parquet(cache)
        return Crosswalk(inst, _test_audit(inst))

    raw = load_gartner_test_raw().drop_duplicates("mutation_key").set_index("mutation_key")
    raw_bag = raw["raw_bag_state"]
    raw_status = raw["raw_screening_status"]
    raw_gfam = raw["genomic_mutation_family_id"]

    cols = ["ID", "Mut Epitope", "HLA", "Peptide Mutant", "key", "minkey", "Screening Status", "epitope status"]
    parts = []
    for ch in pd.read_csv(GARTNER_TEST_EXTRACT, sep="\t", usecols=cols, dtype=str, chunksize=400000):
        ch = ch.rename(columns={"ID": "patient_id", "HLA": "hla_allele", "Peptide Mutant": "peptide",
                                "key": "mutation_key", "minkey": "instance_key",
                                "Screening Status": "mmps_screening_status", "epitope status": "epitope_status"})
        # Mmps-propagated bag state (recorded for the conflict diagnostic only) ...
        ch["mmps_bag_state"] = ch["mmps_screening_status"].map(STATUS_MAP).fillna("UNKNOWN")
        # ... RECONCILED to the authoritative raw exact-key ascertainment.
        ch["bag_state"] = ch["mutation_key"].map(raw_bag).fillna("UNKNOWN")
        ch["raw_screening_status"] = ch["mutation_key"].map(raw_status)
        ch["genomic_mutation_family_id"] = ch["mutation_key"].map(raw_gfam)
        ch["source_conflict"] = (ch["bag_state"] != ch["mmps_bag_state"])
        ch["validated"] = (ch["epitope_status"].astype(str) == "1.0").astype(int)
        parts.append(ch[["patient_id", "peptide", "hla_allele", "mutation_key", "instance_key",
                         "genomic_mutation_family_id", "bag_state", "mmps_bag_state",
                         "raw_screening_status", "mmps_screening_status", "source_conflict",
                         "epitope_status", "validated"]])
    inst = pd.concat(parts, ignore_index=True)
    # instance label uses the RECONCILED (raw) bag state
    inst["instance_label"] = [_test_instance_label(v, s, e) for v, s, e
                              in zip(inst["validated"], inst["bag_state"], inst["epitope_status"])]
    inst["cohort"] = "gartner_test"
    if use_cache:
        inst.to_parquet(cache)
    return Crosswalk(inst, _test_audit(inst))


def _test_instance_label(validated: int, bag_state: str, epitope_status: str) -> str:
    """Exact table: instance label from `epitope status` (1 -> POSITIVE_EXACT; 0 -> screened but this instance
    not immunogenic; NA -> unscreened instance). Bag state from parent `Screening Status`."""
    if validated == 1:
        return "POSITIVE_EXACT"
    if str(epitope_status) == "0.0":
        return "NEGATIVE_BAG_CHILD" if bag_state == "SCREENED_NEG" else "AMBIGUOUS_POSITIVE_BAG"
    return "UNTESTED"     # epitope status NA


def _test_audit(inst: pd.DataFrame) -> dict:
    muts = inst.drop_duplicates("mutation_key")
    tkey_bags = muts["bag_state"].value_counts().to_dict()
    # genomic-family (leakage-blocking unit) counts
    gfam_states = inst.groupby("genomic_mutation_family_id")["bag_state"].agg(set)
    gfam_pos = int(gfam_states.apply(lambda s: "CD8_POS" in s).sum())
    gfam_neg = int(gfam_states.apply(lambda s: ("CD8_POS" not in s) and ("SCREENED_NEG" in s)).sum())
    gfam_unt = int(gfam_states.apply(lambda s: s <= {"UNSCREENED"}).sum())
    # source-conflict diagnostic: transcript keys where Mmps propagation disagrees with raw exact-key
    conflict_keys = muts.loc[muts["source_conflict"], "mutation_key"]
    example = None
    if len(conflict_keys):
        k = conflict_keys.iloc[0]
        row = muts[muts["mutation_key"] == k].iloc[0]
        example = {"key": k, "raw_screening_status": row["raw_screening_status"],
                   "raw_bag_state": row["bag_state"], "mmps_bag_state": row["mmps_bag_state"]}
    return {
        "instances_total": int(len(inst)),
        "patients": int(inst["patient_id"].nunique()),
        "parent_mutations_transcript_keys": int(inst["mutation_key"].nunique()),
        "genomic_mutation_families": int(inst["genomic_mutation_family_id"].nunique()),
        "instance_label_counts": inst["instance_label"].value_counts().to_dict(),
        "positive_exact": int((inst["instance_label"] == "POSITIVE_EXACT").sum()),
        "assay_unit": "transcript 25mer `key`",
        "leakage_block_unit": "genomic_mutation_family_id (patient + raw Variant key)",
        # transcript-key (assay unit) ascertainment — reconciled to raw NmersTestingSet
        "bag_state_counts_transcript_key": {"CD8_POS": int(tkey_bags.get("CD8_POS", 0)),
                                            "SCREENED_NEG": int(tkey_bags.get("SCREENED_NEG", 0)),
                                            "UNSCREENED": int(tkey_bags.get("UNSCREENED", 0))},
        # genomic-family (leakage unit) ascertainment
        "genomic_family_counts": {"total": int(len(gfam_states)), "any_CD8_positive": gfam_pos,
                                  "any_screened_neg_no_pos": gfam_neg, "untested_only": gfam_unt},
        "mmps_vs_raw_source_conflicts": {
            "n_conflicting_transcript_keys": int(len(conflict_keys)),
            "example": example,
            "resolution": "raw NmersTestingSet exact-key Screening Status wins; the Mmps-propagated positive on "
                          "the unscreened sibling is kept UNTESTED (not a CD8 bag). Mmps `Screening Status` is a "
                          "GENOMIC-family-propagated label, not a safe transcript-key ascertainment.",
        },
        "role": "HOLDOUT (semi-consumed: aggregate metrics already observed); never fit/tuned",
    }
