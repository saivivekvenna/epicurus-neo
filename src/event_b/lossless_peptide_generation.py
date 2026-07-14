"""Input-only lossless variant-to-peptide generator (frozen ``lossless-peptide-generation-1.0.0``).

Frozen exploratory protocol:
``docs/superpowers/specs/2026-07-12-osteosarc-lossless-peptide-recovery-exploratory-protocol.md``
(paired copy ``artifacts/milestone_7_decision/peptide_recovery/EXPLORATORY_PROTOCOL.md``).

Given a raw GRCh38 allele ``(chrom, pos, ref, alt)`` and the Ensembl reference only, this module
derives a genomic HGVS, resolves the MANE Select / canonical protein-coding transcript via the
Ensembl VEP REST endpoint, fetches the reference protein / CDS, and enumerates every standard-amino-
acid class-I 8-14mer window that spans the mutated residue (missense) or a novel-frame residue
(frameshift). Windows are crossed with a patient class-I HLA panel to form ``(peptide, HLA)``
candidates that a downstream genuine-PRIME step can score.

Design boundary (hard): this module reads ONLY the raw variant allele, the Ensembl reference (VEP /
sequence), and — for the runner-supplied HLA panel — a class-I allele list. It never reads or imports
any functional-assay, therapeutic, or measured-outcome table, and the frozen import/input-hygiene
test enforces that. Genuine-PRIME scoring and any outcome-label join happen strictly downstream, in
the runner, after candidate generation is complete.

All Ensembl responses are cached to a gitignored raw cache with their URL and a SHA-256 of the bytes,
so an offline rerun reproduces the online run exactly and fails closed on a cache miss. Frozen
per-variant transcript / HGVS / residue expectations are verified before any window is emitted; any
drift aborts rather than silently falling back to a convenient isoform.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

POLICY_ID = "lossless-peptide-generation-1.0.0"

STD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 8, 14

ENSEMBL_BASE_URL = "https://rest.ensembl.org"

# Standard genetic code (DNA codons). Stop codons map to ``*``.
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ---------------------------------------------------------------------------
# 1. Raw GRCh38 allele -> genomic HGVS
# ---------------------------------------------------------------------------
def genomic_hgvs(chrom: object, pos: int, ref: str, alt: str) -> str:
    """Convert a raw ``(chrom, pos, ref, alt)`` GRCh38 allele to a genomic HGVS ``g.`` string.

    Supports SNVs, equal-length multi-nucleotide substitutions, and simple left-anchored
    deletions. Other classes fail closed — this generator never guesses insertion or complex-indel
    HGVS.
    """
    contig = str(chrom).strip()
    if contig.lower().startswith("chr"):
        contig = contig[3:]
    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()
    pos = int(pos)

    if len(ref) == 1 and len(alt) == 1:
        return f"{contig}:g.{pos}{ref}>{alt}"

    # Equal-length block substitution (MNV). HGVS represents this as a genomic delins; keeping
    # the block intact is essential because atomizing it into independent SNVs can invent peptides
    # that do not exist on the observed haplotype.
    if len(ref) == len(alt) and len(ref) > 1:
        end = pos + len(ref) - 1
        return f"{contig}:g.{pos}_{end}delins{alt}"

    # Left-anchored deletion: alt is the shared prefix of ref (VCF convention).
    if len(ref) > len(alt) and ref.startswith(alt) and len(alt) >= 1:
        del_start = pos + len(alt)
        del_end = pos + len(ref) - 1
        if del_start == del_end:
            return f"{contig}:g.{del_start}del"
        return f"{contig}:g.{del_start}_{del_end}del"

    raise ValueError(
        f"unsupported variant class for genomic HGVS: {contig}:{pos} {ref}>{alt} "
        "(only SNV and simple left-anchored deletion are handled; aborting rather than guessing)"
    )


# ---------------------------------------------------------------------------
# 2. Ensembl REST client with URL/SHA offline cache (fail-closed)
# ---------------------------------------------------------------------------
class CacheMiss(RuntimeError):
    """Raised when an offline Ensembl client is asked for an uncached URL."""


def _default_fetcher(url: str) -> bytes:
    """Fetch an Ensembl response, retrying only transient transport/server failures.

    The retry policy cannot alter a successful response and every accepted response is still
    content-addressed by :class:`EnsemblClient`. Permanent client errors fail closed immediately.
    """
    request = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    max_attempts = 6
    for attempt in range(max_attempts):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code != 429 and not 500 <= exc.code < 600:
                raise
            error: BaseException = exc
            retry_after = exc.headers.get("Retry-After") if exc.headers else None
        except (TimeoutError, urllib.error.URLError) as exc:
            error = exc
            retry_after = None

        if attempt == max_attempts - 1:
            raise error
        try:
            server_delay = float(retry_after) if retry_after is not None else 0.0
        except ValueError:
            server_delay = 0.0
        time.sleep(min(60.0, max(server_delay, float(2 ** (attempt + 1)))))

    raise RuntimeError("unreachable Ensembl retry state")


class EnsemblClient:
    """Fetch + cache Ensembl VEP / sequence responses, keyed by URL with a SHA-256 of the bytes.

    Online, responses are written under ``cache_dir`` (a gitignored raw path) and recorded in
    ``manifest.json``. Offline, cached bytes are served without network access and an uncached URL
    raises :class:`CacheMiss` — so an offline rerun reproduces the online run exactly or fails closed.
    Full third-party sequences live only in this gitignored cache; they are never committed.
    """

    def __init__(
        self,
        cache_dir: str | Path,
        *,
        offline: bool = False,
        fetcher: Callable[[str], bytes] | None = None,
        base_url: str = ENSEMBL_BASE_URL,
    ) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.fetcher = fetcher
        self.base_url = base_url.rstrip("/")
        self.manifest_path = self.cache_dir / "manifest.json"
        self.manifest: dict[str, dict] = (
            json.loads(self.manifest_path.read_text()) if self.manifest_path.exists() else {}
        )
        # exact per-run access log: every URL actually CONSUMED this run (cache hit OR network fetch),
        # with the cache filename + content SHA. Additive provenance; does not affect caching behavior.
        self.accessed: dict[str, dict] = {}

    def _get(self, url: str) -> tuple[bytes, str]:
        entry = self.manifest.get(url)
        if entry is not None:
            data = (self.cache_dir / entry["file"]).read_bytes()
            self.accessed[url] = {"file": entry["file"], "sha256": entry["sha256"]}
            return data, entry["sha256"]
        if self.offline:
            raise CacheMiss(f"offline Ensembl client has no cache for {url}")
        fetch = self.fetcher or _default_fetcher
        data = fetch(url)
        sha = hashlib.sha256(data).hexdigest()
        filename = f"{sha}.json"
        (self.cache_dir / filename).write_bytes(data)
        self.manifest[url] = {"file": filename, "sha256": sha}
        self.manifest_path.write_text(json.dumps(self.manifest, indent=2, sort_keys=True) + "\n")
        self.accessed[url] = {"file": filename, "sha256": sha}
        return data, sha

    def vep_hgvs(self, hgvs: str) -> dict:
        url = f"{self.base_url}/vep/human/hgvs/{hgvs}?mane=1;canonical=1;hgvs=1"
        data, sha = self._get(url)
        return {"url": url, "sha256": sha, "json": json.loads(data)}

    def sequence(self, transcript_id: str, seq_type: str) -> dict:
        url = f"{self.base_url}/sequence/id/{transcript_id}?type={seq_type}"
        data, sha = self._get(url)
        payload = json.loads(data)
        return {"url": url, "sha256": sha, "seq": str(payload["seq"]).upper()}


# ---------------------------------------------------------------------------
# 3. MANE/VEP transcript selection + fail-closed verification
# ---------------------------------------------------------------------------
def _hgvs_suffix(value: object) -> str:
    """``ENST00000367409.9:c.6535G>A`` -> ``c.6535G>A`` (drop the versioned transcript prefix)."""
    text = str(value or "")
    return text.split(":", 1)[1] if ":" in text else text


def select_transcript(vep_json: list, *, expected_consequence: str,
                      require_mane_refseq: bool = True) -> dict:  # noqa: D417
    """Select the MANE Select protein-coding transcript whose consequence matches; else canonical.

    Fails closed (``ValueError``) when neither a MANE Select nor a canonical protein-coding
    transcript with the expected consequence is present — never a best-effort isoform.

    ``require_mane_refseq`` (default True = the frozen strict behavior) aborts if the chosen transcript
    is canonical but carries no MANE Select RefSeq. Setting it False is a LABEL-BLIND relaxation used by
    the complete-denominator benchmark: Ensembl's canonical protein-coding transcript is accepted when no
    MANE Select exists (some genes have no MANE entry). This never consults recognition labels and still
    refuses non-canonical best-effort isoforms.
    """
    if not vep_json:
        raise ValueError("empty VEP response")
    consequences = vep_json[0].get("transcript_consequences") or []

    def _match(pool: list[dict]) -> dict | None:
        for tc in pool:
            if (
                tc.get("biotype") == "protein_coding"
                and expected_consequence in (tc.get("consequence_terms") or [])
            ):
                return tc
        return None

    chosen = _match([tc for tc in consequences if tc.get("mane_select")])
    if chosen is None:
        chosen = _match([tc for tc in consequences if tc.get("canonical")])
    if chosen is None:
        raise ValueError(
            f"no MANE Select / canonical protein-coding transcript with consequence "
            f"'{expected_consequence}' in VEP response"
        )
    if require_mane_refseq and not chosen.get("mane_select"):
        raise ValueError("selected transcript is canonical but lacks a MANE Select RefSeq")

    return {
        "transcript_id": str(chosen["transcript_id"]),
        "mane_refseq": str(chosen.get("mane_select") or ""),
        "protein_start": int(chosen["protein_start"]),
        "amino_acids": str(chosen["amino_acids"]),
        "hgvsc": _hgvs_suffix(chosen.get("hgvsc")),
        "hgvsp": _hgvs_suffix(chosen.get("hgvsp")),
        "gene_id": str(chosen.get("gene_id", "")),
        "gene_symbol": str(chosen.get("gene_symbol", "")),
        "consequence_terms": list(chosen.get("consequence_terms") or []),
    }


def verify_transcript(selected: dict, expected: dict) -> None:
    """Assert the selected transcript matches the frozen per-variant expectation; abort on any drift."""
    mismatches = {
        key: (selected.get(key), value)
        for key, value in expected.items()
        if selected.get(key) != value
    }
    if mismatches:
        raise ValueError(f"frozen transcript expectation mismatch (aborting): {mismatches}")


# ---------------------------------------------------------------------------
# 4. Window enumeration
# ---------------------------------------------------------------------------
def enumerate_windows_covering(
    seq: str,
    positions_1based: Iterable[int],
    *,
    min_len: int = MIN_LEN,
    max_len: int = MAX_LEN,
) -> list[str]:
    """Every ``min_len..max_len`` window fully inside ``seq`` that spans at least one given position.

    Positions are 1-based within ``seq``. Windows containing a non-standard amino acid are dropped
    (so a frozen count that assumes standard AAs fails closed if a residue is unexpected). Order is
    deterministic: ascending length, then ascending start.
    """
    targets = {int(p) for p in positions_1based}
    if not targets:
        return []
    lo, hi = min(targets), max(targets)
    windows: list[str] = []
    length = len(seq)
    for window_len in range(min_len, max_len + 1):
        for start in range(0, length - window_len + 1):
            first, last = start + 1, start + window_len  # 1-based inclusive span
            if last < lo or first > hi:
                continue
            if not any(first <= p <= last for p in targets):
                continue
            sub = seq[start : start + window_len]
            if set(sub).issubset(STD_AA):
                windows.append(sub)
    return windows


def missense_windows(protein_seq: str, protein_pos_1based: int, wt_aa: str, mut_aa: str) -> list[str]:
    """All windows spanning a missense-substituted residue. Fail-closed if the reference AA differs."""
    protein_seq = protein_seq.upper()
    index = protein_pos_1based - 1
    if not (0 <= index < len(protein_seq)):
        raise ValueError(f"protein position {protein_pos_1based} out of range for length {len(protein_seq)}")
    observed = protein_seq[index]
    if observed != wt_aa.upper():
        raise ValueError(
            f"reference residue mismatch at protein position {protein_pos_1based}: "
            f"expected {wt_aa!r}, Ensembl protein has {observed!r} (aborting)"
        )
    mutated = protein_seq[:index] + mut_aa.upper() + protein_seq[index + 1 :]
    return enumerate_windows_covering(mutated, {protein_pos_1based})


def substitution_windows(
    protein_seq: str,
    protein_start_1based: int,
    wt_segment: str,
    mutant_segment: str,
) -> list[str]:
    """All windows spanning changed residues in an equal-length protein substitution.

    VEP can represent one genomic MNV as a multi-residue ``amino_acids`` block (for example
    ``HQ/HE``). The complete reference block is verified, the block is replaced atomically, and
    only positions whose amino acid actually changes anchor emitted windows. This preserves the
    observed haplotype and avoids treating adjacent bases as independent variants.
    """
    protein_seq = protein_seq.upper()
    wt_segment = wt_segment.upper()
    mutant_segment = mutant_segment.upper()
    if not wt_segment or len(wt_segment) != len(mutant_segment):
        raise ValueError("protein substitution requires non-empty equal-length WT and mutant segments")
    if not (set(wt_segment) | set(mutant_segment)).issubset(STD_AA):
        raise ValueError("protein substitution contains a non-standard amino acid")
    start = protein_start_1based - 1
    if not (0 <= start and start + len(wt_segment) <= len(protein_seq)):
        raise ValueError(
            f"protein substitution {protein_start_1based}+{len(wt_segment)} exceeds "
            f"reference length {len(protein_seq)}"
        )
    observed = protein_seq[start : start + len(wt_segment)]
    if observed != wt_segment:
        raise ValueError(
            f"reference segment mismatch at protein position {protein_start_1based}: "
            f"expected {wt_segment!r}, Ensembl protein has {observed!r} (aborting)"
        )
    changed = {
        protein_start_1based + offset
        for offset, (wt, mutant) in enumerate(zip(wt_segment, mutant_segment, strict=True))
        if wt != mutant
    }
    if not changed:
        raise ValueError("protein substitution has no changed amino-acid positions")
    mutated = protein_seq[:start] + mutant_segment + protein_seq[start + len(wt_segment) :]
    return enumerate_windows_covering(mutated, changed)


def inframe_windows(protein_seq: str, protein_start_1based: int, wt_seg: str, mut_seg: str) -> list[str]:
    """All windows spanning an IN-FRAME deletion/insertion/delins junction. VEP ``amino_acids`` is
    ``WT/MUT`` (``MUT`` == ``-`` or empty for a pure deletion). The reference segment is verified against
    the fetched protein (fail-closed); the junction residues (first kept residue on each side) anchor the
    windows so every emitted peptide contains novel context created by the edit. Downstream is wild-type
    (in-frame), so only the junction neighbourhood is novel."""
    protein_seq = protein_seq.upper()
    wt_seg = "" if wt_seg.strip() in {"-", ""} else wt_seg.upper()
    mut_seg = "" if mut_seg.strip() in {"-", ""} else mut_seg.upper()
    start = protein_start_1based - 1
    if not (0 <= start <= len(protein_seq)):
        raise ValueError(f"protein start {protein_start_1based} out of range for length {len(protein_seq)}")
    observed = protein_seq[start : start + len(wt_seg)]
    if observed != wt_seg:
        raise ValueError(
            f"reference segment mismatch at protein position {protein_start_1based}: "
            f"expected {wt_seg!r}, Ensembl protein has {observed!r} (aborting)"
        )
    mutated = protein_seq[:start] + mut_seg + protein_seq[start + len(wt_seg):]
    # novel junction positions in the MUTATED protein: the inserted/delins residues, plus the residue
    # immediately flanking the edit on each side (their adjacency is new).
    lo = max(1, start)                       # 1-based position just before the edit (flank)
    hi = start + max(len(mut_seg), 1)        # through the last inserted residue / first kept-after residue
    positions = set(range(lo, min(hi, len(mutated)) + 1))
    return enumerate_windows_covering(mutated, positions)


def translate_to_stop(nt: str) -> str:
    """Translate a nucleotide sequence from its start to the first stop; the stop is excluded."""
    nt = nt.upper()
    residues: list[str] = []
    for i in range(0, len(nt) - 2, 3):
        aa = CODON_TABLE.get(nt[i : i + 3])
        if aa is None:  # non-standard codon (e.g. contains N) — stop translating, never guess
            break
        if aa == "*":
            break
        residues.append(aa)
    return "".join(residues)


def frameshift_novel_protein(
    cds: str, cds_start_1based: int, del_c_start: int, del_c_end: int
) -> str:
    """Apply a 1-based inclusive CDS deletion, then translate to the first stop.

    ``cds_start_1based`` is the CDS position (c.) of the first nucleotide of ``cds`` (1 for a full
    CDS; a codon-aligned offset for a slice), so the deletion is addressed in c. coordinates.
    """
    if (cds_start_1based - 1) % 3 != 0:
        raise ValueError(f"cds_start_1based {cds_start_1based} is not codon-aligned")
    i0 = del_c_start - cds_start_1based
    j0 = del_c_end - cds_start_1based
    if not (0 <= i0 <= j0 < len(cds)):
        raise ValueError(
            f"deletion c.{del_c_start}_{del_c_end} out of range for CDS slice starting at "
            f"c.{cds_start_1based} (length {len(cds)})"
        )
    edited = cds[:i0] + cds[j0 + 1 :]
    return translate_to_stop(edited)


def frameshift_windows(
    cds: str,
    cds_start_1based: int,
    del_c_start: int,
    del_c_end: int,
    *,
    protein_start: int,
) -> tuple[str, list[str]]:
    """Novel protein + every window spanning >=1 novel-frame residue (position >= ``protein_start``).

    No purely pre-junction (reference-frame) window is emitted; the stop is never included.
    """
    novel = frameshift_novel_protein(cds, cds_start_1based, del_c_start, del_c_end)
    offset_protein = (cds_start_1based - 1) // 3 + 1
    novel_local_start = protein_start - offset_protein + 1
    if novel_local_start < 1:
        raise ValueError(
            f"protein_start {protein_start} precedes the CDS slice offset {offset_protein}"
        )
    positions = range(novel_local_start, len(novel) + 1)
    return novel, enumerate_windows_covering(novel, positions)


# ---------------------------------------------------------------------------
# 5. Patient class-I HLA panel + candidate assembly
# ---------------------------------------------------------------------------
def read_hla_panel(pvac_tsv_path: str | Path, *, column: str = "HLA Allele") -> list[str]:
    """Read the patient's class-I HLA panel programmatically from the pVAC candidate table."""
    frame = pd.read_csv(pvac_tsv_path, sep="\t", usecols=[column], low_memory=False)
    return sorted({str(a).strip() for a in frame[column] if str(a).strip()})


def _parse_frameshift_hgvsc(hgvsc: str) -> tuple[int, int]:
    """``c.2603_2630del`` -> ``(2603, 2630)``; ``c.100del`` -> ``(100, 100)``."""
    match = re.fullmatch(r"c\.(\d+)(?:_(\d+))?del[A-Z]*", hgvsc)
    if not match:
        raise ValueError(f"cannot parse frameshift CDS deletion from hgvsc {hgvsc!r}")
    start = int(match.group(1))
    end = int(match.group(2)) if match.group(2) else start
    return start, end


def generate_variant_candidates(
    variant: dict,
    client: EnsemblClient,
    hla_panel: list[str],
    *,
    expected: dict | None = None,
    require_mane_refseq: bool = True,
) -> dict:
    """Generate ``(peptide, HLA)`` candidates for one variant (network via ``client``).

    ``variant`` carries ``chrom, pos, ref, alt, gene, source_variant_type`` — ``missense``/``snv``,
    ``frameshift``/``deletion``, or ``inframe``/``inframe_deletion``/``inframe_insertion``. Returns the
    candidate rows plus a provenance record (Ensembl URLs + SHAs, transcript fields, window counts).

    ``expected`` is OPTIONAL: when a frozen expected-transcript dict is supplied it is verified before any
    window is emitted (target-conditioned reconstruction). When ``expected is None`` the generator runs
    LABEL-BLIND — it trusts VEP's MANE/canonical transcript with the variant's own consequence, which is
    what an end-to-end benchmark over the complete variant universe requires (no per-target verification).
    ``require_mane_refseq`` is threaded to ``select_transcript`` (see there for the label-blind relaxation).
    """
    kind = str(variant["source_variant_type"]).lower()
    is_frameshift = kind in {"frameshift", "deletion", "frameshift_variant"}
    is_inframe = kind in {"inframe", "inframe_deletion", "inframe_insertion"}
    consequence = ("frameshift_variant" if is_frameshift
                   else kind if is_inframe else "missense_variant")

    hgvs = genomic_hgvs(variant["chrom"], variant["pos"], variant["ref"], variant["alt"])
    vep = client.vep_hgvs(hgvs)
    selected = select_transcript(vep["json"], expected_consequence=consequence,
                                 require_mane_refseq=require_mane_refseq)
    if expected is not None:
        verify_transcript(selected, expected)

    provenance: dict = {
        "variant": {k: variant[k] for k in ("chrom", "pos", "ref", "alt", "gene")},
        "genomic_hgvs": hgvs,
        "consequence": consequence,
        "transcript_id": selected["transcript_id"],
        "mane_refseq": selected["mane_refseq"],
        "protein_start": selected["protein_start"],
        "amino_acids": selected["amino_acids"],
        "hgvsc": selected["hgvsc"],
        "hgvsp": selected["hgvsp"],
        "gene_id": selected["gene_id"],
        "gene_symbol": selected.get("gene_symbol", ""),   # exposed for downstream schema (additive)
        "ensembl": {"vep": {"url": vep["url"], "sha256": vep["sha256"]}},
    }

    if is_frameshift:
        cds = client.sequence(selected["transcript_id"], "cds")
        protein = client.sequence(selected["transcript_id"], "protein")
        del_start, del_end = _parse_frameshift_hgvsc(selected["hgvsc"])
        # Fail-closed: the unedited CDS must translate to the fetched reference protein prefix.
        reference_translation = translate_to_stop(cds["seq"])
        prefix_len = selected["protein_start"] - 1
        if reference_translation[:prefix_len] != protein["seq"][:prefix_len]:
            raise ValueError("CDS translation does not match the reference protein prefix (aborting)")
        novel, windows = frameshift_windows(
            cds["seq"], 1, del_start, del_end, protein_start=selected["protein_start"]
        )
        junction_context = novel[max(0, selected["protein_start"] - 8) :][:44]
        provenance["ensembl"]["cds"] = {"url": cds["url"], "sha256": cds["sha256"]}
        provenance["ensembl"]["protein"] = {"url": protein["url"], "sha256": protein["sha256"]}
        provenance["novel_junction_context"] = junction_context
    elif is_inframe:
        protein = client.sequence(selected["transcript_id"], "protein")
        wt_seg, mut_seg = (selected["amino_acids"].split("/", 1) + [""])[:2]
        windows = inframe_windows(protein["seq"], selected["protein_start"], wt_seg, mut_seg)
        provenance["ensembl"]["protein"] = {"url": protein["url"], "sha256": protein["sha256"]}
    else:
        protein = client.sequence(selected["transcript_id"], "protein")
        wt_aa, mut_aa = selected["amino_acids"].split("/")
        if len(wt_aa) == len(mut_aa) == 1:
            windows = missense_windows(protein["seq"], selected["protein_start"], wt_aa, mut_aa)
        else:
            windows = substitution_windows(
                protein["seq"], selected["protein_start"], wt_aa, mut_aa
            )
        provenance["ensembl"]["protein"] = {"url": protein["url"], "sha256": protein["sha256"]}

    provenance["n_windows"] = len(windows)
    provenance["n_unique_peptides"] = len(set(windows))

    unique_peptides = sorted(set(windows))
    rows = []
    for peptide in unique_peptides:
        for hla in hla_panel:
            rows.append(
                {
                    "patient_id": "sid",
                    "mutation_id": f"{variant['gene']}-{variant['chrom']}-{variant['pos']}",
                    "gene_symbol": variant["gene"],
                    "chrom": variant["chrom"],
                    "pos": int(variant["pos"]),
                    "ref": variant["ref"],
                    "alt": variant["alt"],
                    "source_variant_type": "SNV" if not is_frameshift else "frameshift",
                    "mhc_class": "I",
                    "mutant_peptide": peptide,
                    "hla_allele": hla,
                    "candidate_source": "lossless_recovery",
                }
            )
    candidates = pd.DataFrame(rows)
    return {"candidates": candidates, "provenance": provenance, "windows": windows}


# ---------------------------------------------------------------------------
# 6. Candidate union (stable genomic candidate identity)
# ---------------------------------------------------------------------------
def union_candidates(
    frames: Iterable[pd.DataFrame],
    *,
    identity_cols: tuple[str, ...] = ("patient_id", "mutation_id", "mutant_peptide", "hla_allele"),
) -> pd.DataFrame:
    """Concatenate candidate frames and dedup on the stable candidate identity (keep first).

    ``keep='first'`` makes the earlier frame (the incumbent pVAC set) the representative for any
    duplicate ``(patient, mutation, peptide, HLA)`` route; distinct routes are never collapsed.
    """
    parts = [f for f in frames if f is not None and len(f)]
    if not parts:
        return pd.DataFrame()
    combined = pd.concat(parts, ignore_index=True, sort=False)
    keys = [c for c in identity_cols if c in combined.columns]
    return combined.drop_duplicates(subset=keys, keep="first").reset_index(drop=True)
