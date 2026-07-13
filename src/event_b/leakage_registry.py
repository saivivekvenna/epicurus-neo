"""Cross-source leakage-control registry for a PRIME augmentation experiment.

Given the incumbent's training sources (PRIME/public — a documented blocker), CEDAR, Zhao,
and the current Event-B backbone, decide which CEDAR rows are GENUINELY NOVEL and safe to add
to a challenger's training set, and how to group evaluation splits so no peptide leaks between
train and test.

Rules:
    * Exact leakage: a CEDAR peptide identical (canonical uppercase) to any Zhao / backbone
      peptide is removed from the novel training pool (it is neither novel nor test-safe).
    * Near-duplicate leakage: a CEDAR peptide with normalized edit similarity >= threshold to a
      Zhao / backbone peptide is removed (k-mer-blocked to stay tractable).
    * PRIME-training leakage: RESERVED and BLOCKED — PRIME 2.0's published training table is not
      available locally, so CEDAR rows that are in PRIME's own training set cannot yet be excluded.
      This is recorded as an explicit blocker, not silently ignored.
    * Contradictions: CEDAR (peptide, HLA) pairs with both POSITIVE and TESTED_NEGATIVE across
      the literature are AMBIGUOUS. They are excluded from the primary novel training set and kept
      only for an explicit sensitivity arm.
    * Grouping: evaluation/CV groups are keyed by PMID (study) and by peptide-similarity cluster
      so neither a study nor a near-duplicate peptide spans train and test.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from event_b.splits import peptide_cluster_ids, peptide_similarity


def canonical_peptide(value: object) -> str:
    return "".join(ch for ch in str(value).upper() if ch.isalpha())


def canonical_hla(value: object) -> str:
    """Normalize HLA spellings to LOCUS*FIELD1:FIELD2 (e.g. HLA-A*02:01 / HLA-A02:01 -> A*02:01).

    Generic restrictions ('HLA class I/II', blanks) map to '' so they never form a spurious key.
    """
    text = str(value).strip().upper()
    if not text or text in {"NAN", "<NA>"} or "CLASS" in text:
        return ""
    text = text.replace("HLA-", "").replace("HLA", "").replace(" ", "").replace("*", "")
    if len(text) >= 5 and text[0] in "ABCEG":
        rest = text[1:].replace(":", "")
        if len(rest) >= 4 and rest[:4].isalnum():
            return f"{text[0]}*{rest[:2]}:{rest[2:4]}"
    return text


def _kmer_index(peptides: set[str], k: int = 5) -> dict[str, list[str]]:
    index: dict[str, list[str]] = {}
    for peptide in peptides:
        grams = {peptide[i : i + k] for i in range(max(1, len(peptide) - k + 1))} or {peptide}
        for gram in grams:
            index.setdefault(gram, []).append(peptide)
    return index


def near_duplicate(
    query: str, reference_index: dict[str, list[str]], *, threshold: float, k: int = 5
) -> str | None:
    """Return a reference peptide with similarity >= threshold to ``query``, else None.

    k-mer blocked, and length-pruned: a pair whose lengths differ by more than
    ``(1 - threshold)`` of the longer length can never reach the similarity threshold, so it is
    skipped without computing the edit distance.
    """
    grams = {query[i : i + k] for i in range(max(1, len(query) - k + 1))} or {query}
    seen: set[str] = set()
    ql = len(query)
    for gram in grams:
        for candidate in reference_index.get(gram, ()):
            if candidate in seen:
                continue
            seen.add(candidate)
            longer = max(ql, len(candidate))
            if longer and abs(ql - len(candidate)) / longer > (1.0 - threshold):
                continue
            if peptide_similarity(query, candidate) >= threshold:
                return candidate
    return None


class _UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


# k-mer buckets larger than this are unioned wholesale instead of verified pairwise. Wholesale
# union is leakage-CONSERVATIVE (it can only merge groups, never split a similar pair apart) and
# it bounds the cost so clustering never degrades to O(n^2) on a common k-mer.
_BUCKET_UNION_CAP = 400


def peptide_similarity_clusters(
    peptides, *, threshold: float = 0.8, k: int = 4
) -> dict[str, int]:
    """Map each peptide to a similarity-cluster id via k-mer-blocked union-find.

    Two peptides are unioned when they share a k-mer AND their normalized edit similarity is
    >= threshold. Candidate pairs come only from shared k-mer buckets (never all-pairs), and a
    length prune skips pairs that cannot reach the threshold, so the cost stays near-linear.
    """
    unique = sorted({canonical_peptide(p) for p in peptides} - {""})
    index = {p: i for i, p in enumerate(unique)}
    uf = _UnionFind(len(unique))
    buckets: dict[str, list[int]] = {}
    for i, pep in enumerate(unique):
        grams = {pep[j : j + k] for j in range(max(1, len(pep) - k + 1))} or {pep}
        for gram in grams:
            buckets.setdefault(gram, []).append(i)
    for members in buckets.values():
        if len(members) > _BUCKET_UNION_CAP:
            base = members[0]
            for other in members[1:]:
                uf.union(base, other)
            continue
        for a in range(len(members)):
            ia = members[a]
            pa = unique[ia]
            for b in range(a + 1, len(members)):
                ib = members[b]
                if uf.find(ia) == uf.find(ib):
                    continue
                pb = unique[ib]
                longer = max(len(pa), len(pb))
                if longer and abs(len(pa) - len(pb)) / longer > (1.0 - threshold):
                    continue
                if peptide_similarity(pa, pb) >= threshold:
                    uf.union(ia, ib)
    return {pep: uf.find(index[pep]) for pep in unique}


def leakage_safe_component_groups(
    pmids, peptides, *, threshold: float = 0.8, k: int = 4
) -> list[str]:
    """Connected components of the bipartite PMID <-> peptide-similarity-cluster graph.

    Rows in the same component share a PMID or a peptide-similarity cluster (transitively), so a
    GroupKFold on these ids can never place the same study OR two similar peptides on both sides
    of a train/test split.
    """
    pmids = list(pmids)
    peptides = [canonical_peptide(p) for p in peptides]
    n = len(pmids)
    uf = _UnionFind(n)
    first_pmid: dict[str, int] = {}
    for i, pm in enumerate(pmids):
        key = str(pm)
        if key in first_pmid:
            uf.union(first_pmid[key], i)
        else:
            first_pmid[key] = i
    clusters = peptide_similarity_clusters(peptides, threshold=threshold, k=k)
    first_cluster: dict[int, int] = {}
    for i, pep in enumerate(peptides):
        root = clusters.get(pep)
        if root is None:
            continue
        if root in first_cluster:
            uf.union(first_cluster[root], i)
        else:
            first_cluster[root] = i
    return [f"component:{uf.find(i)}" for i in range(n)]


def assert_group_leakage_safe(
    frame: pd.DataFrame,
    fold: pd.Series,
    *,
    pmid_col: str = "pmid",
    peptide_col: str = "peptide",
    threshold: float = 0.8,
) -> None:
    """Assert no PMID, similarity cluster, or near-duplicate pair spans two folds."""
    work = frame.copy()
    work["_fold"] = list(fold)
    if work.groupby(pmid_col)["_fold"].nunique().gt(1).any():
        raise AssertionError("a PMID spans multiple folds")
    clusters = peptide_similarity_clusters(work[peptide_col], threshold=threshold)
    work["_cluster"] = work[peptide_col].map(lambda p: clusters.get(canonical_peptide(p)))
    if work.groupby("_cluster")["_fold"].nunique().gt(1).any():
        raise AssertionError("a peptide-similarity cluster spans multiple folds")


@dataclass(frozen=True)
class LeakageRegistry:
    novel: pd.DataFrame  # CEDAR rows safe to add to training (novel, non-contradictory)
    ambiguous: pd.DataFrame  # contradictory CEDAR rows (sensitivity-only)
    report: dict


def build_leakage_registry(
    cedar: pd.DataFrame,
    *,
    zhao_peptides: set[str],
    backbone_peptides: set[str],
    similarity_threshold: float = 0.8,
    prime_training_peptides: set[str] | None = None,
) -> LeakageRegistry:
    """Partition CEDAR into a novel training pool, an ambiguous pool, and a leakage report."""
    work = cedar.copy()
    work["pep"] = work["peptide"].map(canonical_peptide)
    work["hla_key"] = work["mhc_allele"].map(canonical_hla)

    protected = {canonical_peptide(p) for p in (zhao_peptides | backbone_peptides)}
    protected.discard("")
    exact_mask = work["pep"].isin(protected)

    zhao_canon = {canonical_peptide(p) for p in zhao_peptides} - {""}
    backbone_canon = {canonical_peptide(p) for p in backbone_peptides} - {""}
    reference_index = _kmer_index(zhao_canon | backbone_canon)
    # Test near-duplication once per UNIQUE non-exact peptide, then broadcast to rows.
    unique_peptides = set(work.loc[~exact_mask, "pep"]) - {""}
    near_peptides = {
        pep
        for pep in unique_peptides
        if near_duplicate(pep, reference_index, threshold=similarity_threshold) is not None
    }
    near_mask = work["pep"].isin(near_peptides) & ~exact_mask

    prime_leak_mask = pd.Series(False, index=work.index)
    prime_blocked = prime_training_peptides is None
    if prime_training_peptides:
        prime_canon = {canonical_peptide(p) for p in prime_training_peptides} - {""}
        prime_leak_mask = work["pep"].isin(prime_canon)

    # Contradictions: (peptide, hla) pairs with both POSITIVE and TESTED_NEGATIVE.
    pair = work.groupby(["pep", "hla_key"], dropna=False)["response_label"]
    contradictory = pair.transform("nunique").gt(1)

    leaked = exact_mask | near_mask | prime_leak_mask
    novel_mask = ~leaked & ~contradictory
    ambiguous_mask = ~leaked & contradictory

    novel = work.loc[novel_mask].copy()
    ambiguous = work.loc[ambiguous_mask].copy()

    # Study grouping key for study-held-out splits (PMID). Peptide-similarity clustering is
    # available via peptide_cluster_ids but is O(n^2) and not needed for PMID-grouped CV, so it
    # is computed lazily by callers that request it, not eagerly over the whole novel pool.
    if not novel.empty:
        novel["study_group"] = novel["pmid"].fillna("<no_pmid>").astype(str)

    report = {
        "cedar_rows_in": int(len(work)),
        "excluded_exact_overlap": int(exact_mask.sum()),
        "excluded_near_duplicate": int((near_mask & ~exact_mask).sum()),
        "excluded_prime_training_overlap": int(prime_leak_mask.sum()),
        "excluded_contradictory_ambiguous": int(ambiguous_mask.sum()),
        "novel_rows": int(novel_mask.sum()),
        "novel_unique_peptides": int(novel["pep"].nunique()) if not novel.empty else 0,
        "novel_unique_pmids": int(novel["pmid"].nunique()) if not novel.empty else 0,
        "novel_label_counts": (
            novel["response_label"].value_counts().to_dict() if not novel.empty else {}
        ),
        "similarity_threshold": similarity_threshold,
        "prime_training_exclusion": (
            "BLOCKED: PRIME 2.0 published training table not available locally; CEDAR rows in "
            "PRIME's own training set cannot yet be excluded. Any real-PRIME augmentation result "
            "is invalid until this exclusion is applied."
            if prime_blocked
            else f"applied ({int(prime_leak_mask.sum())} rows removed)"
        ),
    }
    return LeakageRegistry(novel, ambiguous, report)
