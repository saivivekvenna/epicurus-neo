"""Public reconstruction of the osteosarc.com (Sid Sijbrandij) recognition record.

Reconstructs the *full* public osteosarc.com evidence graph into an evidence-graded ledger,
per the frozen preregistration
``docs/superpowers/specs/2026-07-12-osteosarc-sid-reconstruction-preregistration.md``
(paired protocol under ``artifacts/milestone_7_decision/osteosarc_sid_reconstruction/``).

This SUPERSEDES the ``dd3efd1`` diagnostic (``scripts/osteosarc_rank.py``), which used the
21-mutation curated set as the universe and treated the other 20 as assumed negatives. Here we
parse the site structurally (stdlib ``html.parser`` mini-DOM — no new dependency), keep the two
label streams (site ELISPOT vs Hudson IFNγ/TCR) separate, never invent negatives from absence,
and never propagate pool-positivity to member peptides.

Nothing here fits, tunes, or compares a model. The frozen Epicurus config is not touched.

Runtime inputs are fetched with a stable UA + retry/backoff and cached (gitignored) under
``data/raw/osteosarc/site_cache/`` with URL+content SHA256 recorded, so a rerun is network-free.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path

# --------------------------------------------------------------------------------------------
# Paths / constants
# --------------------------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/osteosarc"
SITE_CACHE = RAW / "site_cache"
ART = ROOT / "artifacts/milestone_7_decision/osteosarc_sid_reconstruction"

BASE = "https://osteosarc.com"
INDEX_URL = f"{BASE}/variants/"
VAFS_LONG_URL = f"{BASE}/variants/variant_vafs_long.tsv"
VAFS_COLS_URL = f"{BASE}/variants/variant_vafs_long.columns.tsv"
USER_AGENT = (
    "epicurus-neo-research/1.0 "
    "(+contact: vihaan.sharma@mail.utoronto.ca; academic recognition-benchmark reconstruction)"
)

# Frozen site headline invariants (spec §4).
EXPECT_VARIANTS = 182
EXPECT_VACCINE = 44
EXPECT_ELISPOT = 14

LABEL_STATES = {
    "POSITIVE_STRONG", "POSITIVE_WEAK", "POSITIVE", "NEGATIVE", "AMBIGUOUS", "UNTESTED",
}
RESOLUTION_STATES = {
    "INDIVIDUAL_PEPTIDE", "MUTATION_LONG_PEPTIDE", "POOL", "MUTATION_TCR", "UNKNOWN",
}

# Stable Ensembl gene ids for the Hudson-recognized genes (so expression can be reported even
# when the gene is absent from the pVACtools candidate table — that absence is the finding).
RECOGNIZED_ENSG = {
    "ASPM": "ENSG00000066279", "DYNC1H1": "ENSG00000197102", "MAP2": "ENSG00000078018",
}

_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


# --------------------------------------------------------------------------------------------
# Minimal, robust stdlib DOM (html.parser) — NOT a one-off grep
# --------------------------------------------------------------------------------------------

_VOID = {"br", "img", "input", "meta", "link", "hr", "source", "area", "base", "col", "embed",
         "param", "track", "wbr"}


class _El:
    __slots__ = ("tag", "attrs", "children", "parent")

    def __init__(self, tag: str, attrs, parent):
        self.tag = tag
        self.attrs = {k: (v or "") for k, v in attrs}
        self.children: list = []
        self.parent = parent

    def cls(self) -> list[str]:
        return self.attrs.get("class", "").split()

    def has_class(self, name: str) -> bool:
        return name in self.cls()

    def text(self) -> str:
        out: list[str] = []
        stack = [self]
        # depth-first, preserving order
        def rec(n: "_El"):
            for c in n.children:
                if isinstance(c, str):
                    out.append(c)
                else:
                    rec(c)
        rec(self)
        return "".join(out).replace("\xa0", " ")

    def iter_desc(self):
        for c in self.children:
            if isinstance(c, _El):
                yield c
                yield from c.iter_desc()

    def find_all(self, tag: str | None = None, cls: str | None = None) -> list["_El"]:
        res = []
        for d in self.iter_desc():
            if tag is not None and d.tag != tag:
                continue
            if cls is not None and not d.has_class(cls):
                continue
            res.append(d)
        return res

    def find(self, tag: str | None = None, cls: str | None = None) -> "_El | None":
        for d in self.iter_desc():
            if tag is not None and d.tag != tag:
                continue
            if cls is not None and not d.has_class(cls):
                continue
            return d
        return None


class _TreeBuilder(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root = _El("#root", [], None)
        self.cur = self.root

    def handle_starttag(self, tag, attrs):
        el = _El(tag, attrs, self.cur)
        self.cur.children.append(el)
        if tag not in _VOID:
            self.cur = el

    def handle_startendtag(self, tag, attrs):
        self.cur.children.append(_El(tag, attrs, self.cur))

    def handle_endtag(self, tag):
        node = self.cur
        while node is not self.root and node.tag != tag:
            node = node.parent
        if node is not self.root and node.tag == tag:
            self.cur = node.parent

    def handle_data(self, data):
        self.cur.children.append(data)


def parse_html(text: str) -> _El:
    b = _TreeBuilder()
    b.feed(text)
    b.close()
    return b.root


def _ws(s: str) -> str:
    return " ".join((s or "").split())


# --------------------------------------------------------------------------------------------
# Protein-change normalization (3-letter <-> 1-letter) for coordinate/protein reachability joins
# --------------------------------------------------------------------------------------------

def normalize_protein(p: str) -> str:
    """'p.Gly2179Arg' -> 'G2179R'; 'p.V314I' -> 'V314I'; 'p.Leu867fs' -> 'L867FS' (upper, no 'p.')."""
    s = _ws(p or "")
    if not s:
        return ""
    s = re.sub(r"^p\.", "", s, flags=re.I)
    # expand 3-letter codes
    s = re.sub(r"[A-Z][a-z]{2}", lambda m: _AA3TO1.get(m.group(0), m.group(0)), s)
    return s.upper().replace(" ", "")


def protein_fs_position(p: str) -> int | None:
    """Frameshift residue position, e.g. p.Gly868fs -> 868 (used to disambiguate MAP2)."""
    m = re.search(r"(\d+)\s*fs", normalize_protein(p), flags=re.I)
    return int(m.group(1)) if m else None


# --------------------------------------------------------------------------------------------
# Fetching with cache + provenance
# --------------------------------------------------------------------------------------------

def _cache_name(url: str) -> str:
    u = url.split("#")[0]
    if u.rstrip("/").endswith("/variants"):
        return "variants_index.html"
    if "/variant/" in u:
        slug = u.rstrip("/").split("/variant/")[1]
        return f"variant__{slug}.html"
    return u.rstrip("/").split("/")[-1]


@dataclass
class Fetcher:
    cache_dir: Path = SITE_CACHE
    offline: bool = False
    refresh: bool = False
    delay: float = 0.05
    retries: int = 4
    provenance: dict = field(default_factory=dict)

    def __post_init__(self):
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get(self, url: str) -> str:
        path = self.cache_dir / _cache_name(url)
        if path.exists() and not self.refresh:
            text = path.read_text()
            self._record(url, path, text, from_cache=True)
            return text
        if self.offline:
            raise FileNotFoundError(f"offline and not cached: {url} -> {path}")
        text = self._download(url)
        path.write_text(text)
        self._record(url, path, text, from_cache=False, http_status=200)
        if self.delay:
            time.sleep(self.delay)
        return text

    def _download(self, url: str) -> str:
        last = None
        for attempt in range(self.retries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=45) as resp:
                    return resp.read().decode("utf-8", "replace")
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
                last = e
                time.sleep(0.5 * (2 ** attempt))
        raise RuntimeError(f"failed to fetch {url} after {self.retries} attempts: {last}")

    def _record(self, url, path, text, *, from_cache, http_status=None):
        self.provenance[url] = {
            "url": url,
            "cache_file": str(path.relative_to(ROOT)) if str(path).startswith(str(ROOT)) else str(path),
            "bytes": len(text.encode("utf-8")),
            "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            "http_status": http_status if http_status is not None else "cache",
            "from_cache": from_cache,
        }


# --------------------------------------------------------------------------------------------
# Index parsing
# --------------------------------------------------------------------------------------------

def parse_index(html: str) -> list[dict]:
    root = parse_html(html)
    rows = [tr for tr in root.find_all("tr") if "data-gene" in tr.attrs]
    out = []
    for tr in rows:
        a = tr.find("a")
        href = a.attrs.get("href", "") if a else ""
        slug = href.rstrip("/").split("/variant/")[-1] if "/variant/" in href else ""
        loc = tr.attrs.get("data-location-key", "")
        out.append({
            "slug": slug,
            "gene": tr.attrs.get("data-gene", "").upper(),
            "chr": tr.attrs.get("data-chr", ""),
            "location_key": loc,
            "protein_attr": _ws(tr.attrs.get("data-protein", "")),
            "n_vaccines": _to_int(tr.attrs.get("data-vaccines", "0")),
            "n_pipelines": _to_int(tr.attrs.get("data-pipelines", "0")),
            "elispot_flag": _to_int(tr.attrs.get("data-elispot", "0")),
            "tumor_vaf": _to_float(tr.attrs.get("data-tumor", "")),
        })
    return out


def _to_int(s):
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return 0


def _to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------------------------
# Variant page parsing
# --------------------------------------------------------------------------------------------

def _section_after_h2(root: _El, title_substr: str) -> _El | None:
    """Return the <section> whose <h2> contains title_substr."""
    for sec in root.find_all("section"):
        h2 = sec.find("h2")
        if h2 and title_substr.lower() in _ws(h2.text()).lower():
            return sec
    return None


def parse_variant_page(html: str, variant_id: str) -> dict:
    root = parse_html(html)
    h1 = root.find("h1")
    gene = _ws(h1.text()) if h1 else ""

    # --- Genomic context (dl of dt/dd) ---
    genomic = _parse_dl(_section_after_h2(root, "Genomic context"))
    protein_sec = _section_after_h2(root, "Protein context")
    protein_change = ""
    if protein_sec:
        dd = protein_sec.find("dd")
        # first dd is usually the HGVS protein; fall back to any dd containing 'p.'
        for d in protein_sec.find_all("dd"):
            t = _ws(d.text())
            if t.startswith("p.") or "fs" in t:
                protein_change = t
                break
        if not protein_change and dd:
            protein_change = _ws(dd.text())

    # --- Detection pills ---
    detected, notdetected = [], []
    det = _section_after_h2(root, "Detection")
    if det:
        for pill in det.find_all("span", cls="pill"):
            label = _ws(pill.text()).lstrip("✓").strip()
            if pill.has_class("pill-on"):
                detected.append(label)
            elif pill.has_class("pill-off"):
                notdetected.append(label)

    # --- Vaccines targeting ---
    vaccines, elispot_summary = [], ""
    vax = _section_after_h2(root, "Vaccines targeting")
    if vax:
        for pill in vax.find_all("span", cls="pill-vax"):
            vaccines.append(_ws(pill.text()))
        es = vax.find("span", cls="pill-elispot")
        if es:
            elispot_summary = _ws(es.text()).replace("ELISPOT:", "").strip()

    # --- Peptide blocks + experiments ---
    blocks = _parse_peptide_blocks(root, variant_id, gene)

    # --- Predicted neoantigen peptides (captured for audit only; not a frozen table) ---
    predicted = _parse_predicted(root)

    return {
        "variant_id": variant_id,
        "gene": gene,
        "genomic": genomic,
        "protein_change": protein_change,
        "detected_pipelines": detected,
        "notdetected_pipelines": notdetected,
        "vaccines_targeting": vaccines,
        "elispot_summary": elispot_summary,
        "peptide_blocks": blocks,
        "predicted_peptides": predicted,
    }


def _parse_dl(sec: _El | None) -> dict:
    out: dict[str, str] = {}
    if sec is None:
        return out
    dl = sec.find("dl")
    if dl is None:
        return out
    key = None
    for child in dl.children:
        if isinstance(child, _El) and child.tag == "dt":
            key = _ws(child.text())
        elif isinstance(child, _El) and child.tag == "dd" and key is not None:
            out[key] = _ws(child.text())
            key = None
    return out


def _parse_peptide_blocks(root: _El, variant_id: str, gene: str) -> list[dict]:
    sec = _section_after_h2(root, "Vaccine peptides")
    if sec is None:
        return []
    blocks = []
    for bi, blk in enumerate(sec.find_all("section", cls="peptide-block")):
        header = blk.find(cls="peptide-header")
        sources = [_ws(s.text()) for s in blk.find_all("span", cls="pep-source")]
        seq_el = blk.find(cls="peptide-seq")
        peptide_seq, minimal = "", ""
        if seq_el is not None:
            for sp in seq_el.find_all("span"):
                t = _ws(sp.text())
                peptide_seq += t
                if sp.has_class("ep"):
                    minimal += t
        aux = blk.find(cls="peptide-aux")
        aux_txt = _ws(aux.text()) if aux else ""
        aa = _first_int(aux_txt, r"(\d+)\s*aa")
        declared = _first_int(aux_txt, r"(\d+)\s*experiments?")
        exps = _parse_exp_table(blk, variant_id, bi, gene, peptide_seq, minimal)
        blocks.append({
            "variant_id": variant_id,
            "gene": gene,
            "block_index": bi,
            "sources": sources,
            "peptide_seq": peptide_seq,
            "minimal_epitope": minimal,
            "peptide_len": len(peptide_seq),
            "aa_declared": aa,
            "declared_experiment_count": declared if declared is not None else 0,
            "experiments": exps,
        })
    return blocks


def _parse_exp_table(blk: _El, variant_id: str, block_index: int, gene: str,
                     peptide_seq: str, minimal: str) -> list[dict]:
    table = blk.find("table", cls="exp-table")
    if table is None:
        return []
    is_minimal_only = bool(minimal) and peptide_seq == minimal
    out = []
    tbody = table.find("tbody") or table
    for tr in tbody.find_all("tr"):
        cells = [td for td in tr.children if isinstance(td, _El) and td.tag == "td"]
        if len(cells) < 6:
            continue
        date = _ws(cells[0].text())
        experiment = _ws(cells[1].text())
        jlf_id = _ws(cells[2].text())
        pool_raw = _ws(cells[3].text())
        res_span = cells[4].find("span", cls="res")
        result_raw = _ws(res_span.text()) if res_span else _ws(cells[4].text())
        result_class = next((c for c in (res_span.cls() if res_span else []) if c.startswith("res-")), "")
        notes = " ¦ ".join(_ws(d.text()) for d in cells[5].find_all("div")) or _ws(cells[5].text())
        label = label_state_from_text(result_raw)
        resolution = resolution_state(pool_raw, is_minimal_only)
        out.append({
            "variant_id": variant_id, "gene": gene, "block_index": block_index,
            "peptide_seq": peptide_seq, "minimal_epitope": minimal,
            "jlf_peptide_id": jlf_id, "exp_date": date, "experiment_name": experiment,
            "pool_raw": pool_raw, "result_raw": result_raw, "result_class": result_class,
            "label_state": label, "resolution_state": resolution, "notes_raw": notes,
        })
    return out


def _parse_predicted(root: _El) -> list[dict]:
    sec = _section_after_h2(root, "Predicted neoantigen")
    if sec is None:
        return []
    table = sec.find("table", cls="peptides-table")
    if table is None:
        return []
    out = []
    tbody = table.find("tbody") or table
    for tr in tbody.find_all("tr"):
        cells = [td for td in tr.children if isinstance(td, _El) and td.tag == "td"]
        if len(cells) < 4:
            continue
        out.append({
            "source": _ws(cells[0].text()), "peptide": _ws(cells[1].text()),
            "hla_allele": _ws(cells[2].text()), "ic50_mut": _ws(cells[3].text()),
        })
    return out


def _first_int(text: str, pattern: str) -> int | None:
    m = re.search(pattern, text)
    return int(m.group(1)) if m else None


# --------------------------------------------------------------------------------------------
# Label / resolution state machines (frozen — spec §3)
# --------------------------------------------------------------------------------------------

def label_state_from_text(result_text: str) -> str:
    t = _ws(result_text).lower()
    if not t:
        return "AMBIGUOUS"
    if t.startswith("positive (strong"):
        return "POSITIVE_STRONG"
    if t.startswith("positive (weak"):
        return "POSITIVE_WEAK"
    if t == "positive":
        return "POSITIVE"
    if t == "negative":
        return "NEGATIVE"
    return "AMBIGUOUS"


_BLANK_POOL = {"", "—", "-", "na", "n/a", "none"}


def resolution_state(pool_raw: str, is_minimal_only: bool) -> str:
    p = _ws(pool_raw).lower()
    if p and p not in _BLANK_POOL:
        return "POOL"
    # blank/NA pool: never treated as proof of individual testing.
    if is_minimal_only:
        return "INDIVIDUAL_PEPTIDE"
    return "MUTATION_LONG_PEPTIDE"


def experiment_key(row: dict) -> str:
    raw = "|".join(str(row.get(k, "")) for k in (
        "variant_id", "block_index", "jlf_peptide_id", "exp_date",
        "experiment_name", "pool_raw", "result_raw"))
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:16]


# --------------------------------------------------------------------------------------------
# Local TSV loaders
# --------------------------------------------------------------------------------------------

def load_vafs_long(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_hudson(paths: list[tuple[str, str, Path]]) -> list[dict]:
    """paths: list of (timepoint, pool_kind, file). Emits per-clonotype rows; label_state at
    (timepoint,mutation) granularity = POSITIVE iff any mutation-specific clonotype, else UNTESTED."""
    rows = []
    for timepoint, pool_kind, fp in paths:
        if not fp.exists():
            continue
        sha = hashlib.sha256(fp.read_bytes()).hexdigest()
        with fp.open() as fh:
            data = list(csv.DictReader(fh, delimiter="\t"))
        # rollup: which (mutation) has any mutation-specific clonotype at this file
        specific = {r["mutation"] for r in data
                    if str(r.get("is_mutation_specific", "")).upper() == "TRUE"
                    and r["mutation"] not in ("", "NA")}
        for r in data:
            mut = r.get("mutation", "")
            if mut in ("", "NA"):
                continue
            gene = mut.split(".")[0]
            protein = mut.split(".", 1)[1] if "." in mut else ""
            is_spec = str(r.get("is_mutation_specific", "")).upper() == "TRUE"
            label = "POSITIVE" if mut in specific else "UNTESTED"
            rows.append({
                "timepoint": timepoint, "mutation_label": mut, "gene": gene,
                "protein_change": protein, "is_mutation_specific": is_spec,
                "trb": r.get("TRB", ""),
                "log2fc_umi": r.get("log2FC_umi", ""),
                "fold_expansion": r.get("fold_expansion", ""),
                "baseline_pct": r.get("baseline_pct", ""),
                "poststim_pct": r.get("poststim_pct", ""),
                "pool_kind": pool_kind, "label_state": label,
                "resolution_state": "MUTATION_TCR",
                "source_file": fp.name, "source_sha256": sha,
            })
    return rows


def load_pvactools_genes(path: Path) -> tuple[set[tuple[str, str]], dict[str, str]]:
    """Return {(gene, normalized_protein)} present, and {gene: ensembl_id} map."""
    present, ensg = set(), {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = (r.get("Gene Name") or "").upper()
            hgvsp = r.get("HGVSp") or ""
            prot = hgvsp.split(":")[-1] if ":" in hgvsp else hgvsp
            present.add((g, normalize_protein(prot)))
            eid = (r.get("Ensembl Gene ID") or "").split(".")[0]
            if g and eid:
                ensg.setdefault(g, eid)
    return present, ensg


def load_curated(path: Path) -> set[tuple[str, str]]:
    present = set()
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            present.add(((r.get("Gene") or "").upper(), normalize_protein(r.get("AA Change") or "")))
    return present


def load_rsem_tpm(path: Path) -> dict[str, float]:
    out = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            eid = (r.get("gene_id") or "").split(".")[0]
            try:
                out[eid] = float(r.get("TPM"))
            except (TypeError, ValueError):
                pass
    return out


# --------------------------------------------------------------------------------------------
# vafs_long variant index (coordinate/protein reachability)
# --------------------------------------------------------------------------------------------

def index_vafs_variants(vafs_rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in vafs_rows:
        vid = r["variant_id"]
        v = out.setdefault(vid, {
            "variant_id": vid, "gene": (r.get("gene") or "").upper(),
            "chrom": r.get("chrom", ""), "pos": r.get("pos", ""),
            "ref": r.get("ref", ""), "alt": r.get("alt", ""),
            "protein_change": r.get("protein_change", ""),
            "consequence": r.get("consequence", ""),
            "on_variants_page": r.get("on_variants_page", ""),
            "pipelines": set(),
        })
        if r.get("pipeline"):
            v["pipelines"].add(r["pipeline"])
        if not v["protein_change"] and r.get("protein_change"):
            v["protein_change"] = r["protein_change"]
    return out


def match_vafs_variant(gene: str, protein: str, vafs_idx: dict[str, dict]) -> dict | None:
    """Coordinate/protein match — never gene alone. Missense by normalized protein; frameshift by
    (gene, fs-position). Returns the vafs variant dict or None."""
    gene = gene.upper()
    npro = normalize_protein(protein)
    fs = protein_fs_position(protein)
    cands = [v for v in vafs_idx.values() if v["gene"] == gene]
    for v in cands:
        if npro and normalize_protein(v["protein_change"]) == npro:
            return v
    if fs is not None:
        fs_hits = [v for v in cands if protein_fs_position(v["protein_change"]) == fs]
        if len(fs_hits) == 1:
            return fs_hits[0]
        if fs_hits:
            return fs_hits[0]  # ambiguity flagged by caller via adjudication
    return None


# --------------------------------------------------------------------------------------------
# Table assembly
# --------------------------------------------------------------------------------------------

def build_variant_catalog(index_rows, page_map, vafs_idx) -> tuple[list[dict], list[str]]:
    catalog, mismatches = [], []
    for ir in sorted(index_rows, key=lambda r: r["slug"]):
        slug = ir["slug"]
        page = page_map.get(slug, {})
        vaf = vafs_idx.get(slug, {})
        gen = page.get("genomic", {})
        chrom = gen.get("Chromosome", "") or vaf.get("chrom", "")
        pos = _digits(gen.get("Position (hg38)", "")) or vaf.get("pos", "")
        refalt = _split_refalt(gen.get("Ref > Alt", ""))
        page_vax = page.get("vaccines_targeting", [])
        # cross-check index counts vs page (gene identity + vaccine count)
        if page.get("gene") and page["gene"].upper() != ir["gene"]:
            mismatches.append(f"gene mismatch {slug}: index={ir['gene']} page={page['gene']}")
        if ir["n_vaccines"] != len(page_vax):
            mismatches.append(
                f"vaccine-count mismatch {slug}: index={ir['n_vaccines']} page={len(page_vax)}")
        catalog.append({
            "variant_id": slug,
            "gene": ir["gene"],
            "chrom": chrom,
            "pos": pos,
            "chr_pos": f"{chrom}:{pos}" if chrom and pos else "",
            "protein_change_raw": page.get("protein_change", "") or vaf.get("protein_change", ""),
            "consequence": vaf.get("consequence", ""),
            "variant_type": _variant_type(refalt),
            "ref": refalt[0] or vaf.get("ref", ""),
            "alt": refalt[1] or vaf.get("alt", ""),
            "change": gen.get("cDNA change", ""),
            "n_vaccines": ir["n_vaccines"],
            "n_pipelines": ir["n_pipelines"],
            "elispot_positive_flag": ir["elispot_flag"],
            "tumor_vaf": ir["tumor_vaf"],
            "detected_pipelines": "; ".join(page.get("detected_pipelines", [])),
            "notdetected_pipelines": "; ".join(page.get("notdetected_pipelines", [])),
            "vaccines_targeting": "; ".join(page_vax),
            "elispot_summary": page.get("elispot_summary", ""),
            "source_url": f"{BASE}/variant/{slug}/",
        })
    return catalog, mismatches


def build_peptide_inventory(page_map) -> tuple[list[dict], list[str]]:
    rows, length_mismatches = [], []
    for slug in sorted(page_map):
        for blk in page_map[slug].get("peptide_blocks", []):
            if blk["aa_declared"] is not None and blk["aa_declared"] != blk["peptide_len"]:
                length_mismatches.append(
                    f"length mismatch {slug}#{blk['block_index']}: "
                    f"aa_declared={blk['aa_declared']} seq_len={blk['peptide_len']}")
            rows.append({
                "variant_id": blk["variant_id"], "gene": blk["gene"],
                "block_index": blk["block_index"], "sources": "; ".join(blk["sources"]),
                "peptide_seq": blk["peptide_seq"], "minimal_epitope": blk["minimal_epitope"],
                "peptide_len": blk["peptide_len"],
                "declared_experiment_count": blk["declared_experiment_count"],
                "parsed_experiment_count": len(blk["experiments"]),
                "source_url": f"{BASE}/variant/{slug}/",
            })
    return rows, length_mismatches


def build_assay_ledger(page_map) -> list[dict]:
    rows = []
    for slug in sorted(page_map):
        for blk in page_map[slug].get("peptide_blocks", []):
            if not blk["experiments"]:
                # vaccine-included peptide with zero experiments -> UNTESTED (never NEGATIVE)
                rows.append({
                    "experiment_key": experiment_key({
                        "variant_id": slug, "block_index": blk["block_index"],
                        "jlf_peptide_id": "", "exp_date": "", "experiment_name": "",
                        "pool_raw": "", "result_raw": "__UNTESTED__"}),
                    "variant_id": slug, "gene": blk["gene"], "block_index": blk["block_index"],
                    "peptide_seq": blk["peptide_seq"], "minimal_epitope": blk["minimal_epitope"],
                    "jlf_peptide_id": "", "exp_date": "", "experiment_name": "",
                    "pool_raw": "", "result_raw": "", "result_class": "",
                    "label_state": "UNTESTED", "resolution_state": "UNKNOWN", "notes_raw": "",
                    "source_url": f"{BASE}/variant/{slug}/",
                })
                continue
            for ex in blk["experiments"]:
                row = dict(ex)
                row["experiment_key"] = experiment_key(ex)
                row["source_url"] = f"{BASE}/variant/{slug}/"
                rows.append(row)
    # canonical column order
    cols = ["experiment_key", "variant_id", "gene", "block_index", "peptide_seq",
            "minimal_epitope", "jlf_peptide_id", "exp_date", "experiment_name", "pool_raw",
            "result_raw", "result_class", "label_state", "resolution_state", "notes_raw",
            "source_url"]
    return [{c: r.get(c, "") for c in cols} for r in rows]


def build_reachability(hudson_rows, catalog, vafs_idx, pv_present, pv_ensg, curated, tpm) -> tuple[list[dict], dict]:
    cat_by_id = {c["variant_id"]: c for c in catalog}
    pv_genes = {g for g, _ in pv_present}
    # Hudson-recognized (gene, protein) at any timepoint
    recognized = sorted({(r["gene"].upper(), r["protein_change"]) for r in hudson_rows
                         if r["is_mutation_specific"]})
    site_positive = {c["variant_id"] for c in catalog if c["elispot_positive_flag"]}

    def gene_tpm(gene, vid):
        eid = RECOGNIZED_ENSG.get(gene) or pv_ensg.get(gene)
        return tpm.get(eid) if eid else None

    def funnel_row(gene, protein, source_tag, recognized_flag):
        gene = gene.upper()
        npro = normalize_protein(protein)
        v = match_vafs_variant(gene, protein, vafs_idx)
        vid = v["variant_id"] if v else ""
        cat = cat_by_id.get(vid, {})
        in_vafs = v is not None
        detected = sorted(v["pipelines"]) if v else []
        in_pv = (gene, npro) in pv_present
        in_curated = (gene, npro) in curated
        in_vax = bool(cat.get("n_vaccines", 0))
        has_elispot = bool(cat.get("elispot_positive_flag", 0)) or (vid in site_positive)
        # automated Epicurus funnel gates, in order
        gates = [("vafs_long_detected", in_vafs and bool(detected)),
                 ("pvactools_2025_candidate", in_pv),
                 ("curated_21_shortlist", in_curated)]
        first_fail = next((name for name, ok in gates if not ok), "reached_shortlist")
        return {
            "target_id": vid or f"{gene}-{npro}",
            "gene": gene, "protein_change": protein,
            "chrom": v["chrom"] if v else "", "pos": v["pos"] if v else "",
            "ref": v["ref"] if v else "", "alt": v["alt"] if v else "",
            "recognized_by": source_tag,
            "in_vafs_long": in_vafs, "detected_pipelines": "; ".join(detected),
            "in_pvactools_2025": in_pv, "in_curated_21": in_curated,
            "in_vaccine": in_vax, "has_site_elispot": has_elispot,
            "site_elispot_best": cat.get("elispot_summary", ""),
            "hudson_recognized": recognized_flag,
            "gene_tpm": gene_tpm(gene, vid),
            "first_failure_stage": first_fail,
        }

    rows, seen = [], {}
    # 1) every Hudson-recognized mutation (required, fully traced)
    for gene, protein in recognized:
        r = funnel_row(gene, protein, "hudson_ifng_tcr", True)
        rows.append(r)
        seen[r["target_id"]] = r
    # 2) every site ELISPOT-positive variant
    for c in catalog:
        if not c["elispot_positive_flag"]:
            continue
        if c["variant_id"] in seen:
            continue
        r = funnel_row(c["gene"], c["protein_change_raw"], "site_elispot_positive", False)
        r["target_id"] = c["variant_id"]  # anchor to site coordinate
        r["has_site_elispot"] = True
        rows.append(r)
        seen[r["target_id"]] = r

    # explicit adjudications required by the spec
    adjudication = _adjudicate(hudson_rows, vafs_idx, pv_genes, cat_by_id)
    for r in rows:
        key = (r["gene"], normalize_protein(r["protein_change"]))
        if key in adjudication:
            r["adjudication"] = adjudication[key]
        else:
            r["adjudication"] = ""
    rows.sort(key=lambda r: (r["recognized_by"], r["gene"], r["target_id"]))
    return rows, adjudication


def _adjudicate(hudson_rows, vafs_idx, pv_genes, cat_by_id) -> dict:
    out = {}
    # ASPM p.G2179R
    aspm = match_vafs_variant("ASPM", "p.G2179R", vafs_idx)
    out[("ASPM", "G2179R")] = (
        f"ASPM p.G2179R (Hudson-recognized May+Aug) = site variant "
        f"{aspm['variant_id'] if aspm else '?'} (p.Gly2179Arg), called by "
        f"{'/'.join(sorted(aspm['pipelines'])) if aspm else '?'}. The page lists an mRNA peptide construct "
        f"but zero vaccines targeting it and zero site-ELISPOT experiments; the Hudson recognition stream "
        f"is the positive evidence. It is ABSENT from the pVACtools 2025.01 candidate universe (no ASPM). "
        f"Lost at candidate generation, NOT variant calling — corrects dd3efd1's 'off-callset' claim.")
    # DYNC1H1 p.V314I
    dyn = match_vafs_variant("DYNC1H1", "p.V314I", vafs_idx)
    out[("DYNC1H1", "V314I")] = (
        f"DYNC1H1 p.V314I (recognized May+Aug) = site variant {dyn['variant_id'] if dyn else '?'} "
        f"(p.Val314Ile); present in pVACtools 2025 AND curated-21 -> reaches the shortlist (the only "
        f"recognized target Epicurus could rank).")
    # MAP2 frameshift
    map2 = [v for v in vafs_idx.values() if v["gene"] == "MAP2" and protein_fs_position(v["protein_change"])]
    labels = sorted({f"{v['variant_id']}({v['protein_change']})" for v in map2})
    out[("MAP2", "GYCVFNKYTV868FS")] = (
        f"MAP2 frameshift (Hudson label p.GYCVFNKYTV868fs, recognized May): the site carries two "
        f"overlapping frameshift annotations {labels} 4bp apart. They have different genomic alleles and "
        f"different neo-frame sequences, so they remain distinct records. The "
        f"vaccine (JLF V1/V2/V3) + ELISPOT-strong record is on MAP2-chr2-209694768 (p.Leu867fs); "
        f"the Hudson '868fs' label position-matches MAP2-chr2-209694772 (p.Gly868fs, no vaccine/experiments). "
        f"The Hudson neo-frame GYCVFNKYTV differs from the Leu867fs vaccine neo-frame (RVVPFTKAL), "
        f"supporting the Gly868fs mapping. BOTH are called by DRAGEN/Sarek/oncoanalyser and absent from pVACtools "
        f"2025 candidate universe -> lost at candidate generation. The prior 'low-TPM expression-filter drop' "
        f"attribution is not confirmable from these files (MAP2 is simply absent from the pVACtools output).")
    return out


# --------------------------------------------------------------------------------------------
# helpers for catalog parsing
# --------------------------------------------------------------------------------------------

def _digits(s: str) -> str:
    return re.sub(r"[^0-9]", "", s or "")


def _split_refalt(s: str) -> tuple[str, str]:
    m = re.search(r"([ACGT\-]+)\s*>\s*([ACGT\-]+)", (s or "").replace("\xa0", " ").upper())
    return (m.group(1), m.group(2)) if m else ("", "")


def _variant_type(refalt: tuple[str, str]) -> str:
    ref, alt = refalt
    if not ref or not alt:
        return ""
    if ref in ("-",) or (len(alt) > len(ref)):
        return "insertion"
    if alt in ("-",) or (len(ref) > len(alt)):
        return "deletion"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    return "mnv"


# --------------------------------------------------------------------------------------------
# Integrity checks (fail-fast — spec §4)
# --------------------------------------------------------------------------------------------

class IntegrityError(AssertionError):
    pass


def run_integrity_checks(catalog, inventory, ledger, hudson, funnel, mismatches):
    errs = []
    if len(catalog) != EXPECT_VARIANTS:
        errs.append(f"variant count {len(catalog)} != {EXPECT_VARIANTS}")
    nvax = sum(1 for c in catalog if c["n_vaccines"] > 0)
    if nvax != EXPECT_VACCINE:
        errs.append(f"vaccine-targeted count {nvax} != {EXPECT_VACCINE}")
    neli = sum(1 for c in catalog if c["elispot_positive_flag"])
    if neli != EXPECT_ELISPOT:
        errs.append(f"ELISPOT-positive count {neli} != {EXPECT_ELISPOT}")
    for m in mismatches:
        errs.append(f"cross-check: {m}")
    # per-block experiment counts. Robust rule (see spec §9 registered deviation): the parser must
    # never MISS a declared experiment (parsed >= declared), and long-vaccine-peptide blocks
    # (declared>0) must match exactly. The site's peptide-aux counter under-counts by omitting the
    # short individual P-series peptide rows (declared=0, parsed=1); those are reported, not failed.
    for r in inventory:
        d, p = r["declared_experiment_count"], r["parsed_experiment_count"]
        if p < d:
            errs.append(f"parser MISS {r['variant_id']}#{r['block_index']}: declared={d} parsed={p}")
        elif p > d and d != 0:
            errs.append(f"unexpected over-parse {r['variant_id']}#{r['block_index']}: "
                        f"declared={d} parsed={p}")
    # dedup uniqueness
    keys = [r["experiment_key"] for r in ledger]
    if len(keys) != len(set(keys)):
        errs.append(f"duplicate experiment_key(s): {len(keys) - len(set(keys))}")
    hkeys = [(r["timepoint"], r["pool_kind"], r["trb"], r["mutation_label"]) for r in hudson]
    if len(hkeys) != len(set(hkeys)):
        errs.append(f"duplicate hudson key(s): {len(hkeys) - len(set(hkeys))}")
    # enum membership
    for r in ledger:
        if r["label_state"] not in LABEL_STATES:
            errs.append(f"bad label_state {r['label_state']}")
        if r["resolution_state"] not in RESOLUTION_STATES:
            errs.append(f"bad resolution_state {r['resolution_state']}")
    # all 3 Hudson positives traced with a non-null first_failure_stage
    recog_genes = {r["gene"] for r in hudson if r["is_mutation_specific"]}
    for gene in ("ASPM", "DYNC1H1", "MAP2"):
        if gene not in recog_genes:
            errs.append(f"expected Hudson-recognized gene missing: {gene}")
        hits = [r for r in funnel if r["gene"] == gene and r["hudson_recognized"]]
        if not hits or any(not r["first_failure_stage"] for r in hits):
            errs.append(f"reachability trace missing/blank for recognized {gene}")
    if errs:
        raise IntegrityError("; ".join(errs))


# --------------------------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------------------------

CATALOG_COLS = ["variant_id", "gene", "chrom", "pos", "chr_pos", "protein_change_raw",
                "consequence", "variant_type", "ref", "alt", "change", "n_vaccines",
                "n_pipelines", "elispot_positive_flag", "tumor_vaf", "detected_pipelines",
                "notdetected_pipelines", "vaccines_targeting", "elispot_summary", "source_url"]
INVENTORY_COLS = ["variant_id", "gene", "block_index", "sources", "peptide_seq",
                  "minimal_epitope", "peptide_len", "declared_experiment_count",
                  "parsed_experiment_count", "source_url"]
LEDGER_COLS = ["experiment_key", "variant_id", "gene", "block_index", "peptide_seq",
               "minimal_epitope", "jlf_peptide_id", "exp_date", "experiment_name", "pool_raw",
               "result_raw", "result_class", "label_state", "resolution_state", "notes_raw",
               "source_url"]
HUDSON_COLS = ["timepoint", "mutation_label", "gene", "protein_change", "is_mutation_specific",
               "trb", "log2fc_umi", "fold_expansion", "baseline_pct", "poststim_pct",
               "pool_kind", "label_state", "resolution_state", "source_file", "source_sha256"]
FUNNEL_COLS = ["target_id", "gene", "protein_change", "chrom", "pos", "ref", "alt",
               "recognized_by", "in_vafs_long", "detected_pipelines", "in_pvactools_2025",
               "in_curated_21", "in_vaccine", "has_site_elispot", "site_elispot_best",
               "hudson_recognized", "gene_tpm", "first_failure_stage", "adjudication"]


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fieldnames)
        for r in rows:
            w.writerow([_fmt(r.get(c)) for c in fieldnames])


# --------------------------------------------------------------------------------------------
# Summary / contradictions / audit
# --------------------------------------------------------------------------------------------

def _find_contradictions(ledger: list[dict]) -> list[dict]:
    """Same peptide with BOTH a positive and a negative result — reported, never collapsed."""
    pos = {"POSITIVE_STRONG", "POSITIVE_WEAK", "POSITIVE"}
    by_pep: dict[tuple, list[dict]] = {}
    for r in ledger:
        if r["label_state"] == "UNTESTED":
            continue
        key = (r["variant_id"], r["peptide_seq"] or r["minimal_epitope"])
        by_pep.setdefault(key, []).append(r)
    out = []
    for (vid, pep), rows in sorted(by_pep.items()):
        states = {r["label_state"] for r in rows}
        if states & pos and "NEGATIVE" in states:
            out.append({
                "variant_id": vid, "peptide": pep,
                "states": sorted(states),
                "evidence": [{"date": r["exp_date"], "experiment": r["experiment_name"],
                              "pool": r["pool_raw"], "result": r["result_raw"],
                              "resolution": r["resolution_state"]} for r in rows],
            })
    return out


def summarize(catalog, inventory, ledger, hudson, funnel, contradictions) -> dict:
    exp_rows = [r for r in ledger if r["label_state"] != "UNTESTED"]
    resolved_rows = [r for r in exp_rows if r["resolution_state"] in
                     ("INDIVIDUAL_PEPTIDE", "MUTATION_LONG_PEPTIDE")]
    resolved_units: dict[tuple[str, str], list[dict]] = {}
    for row in resolved_rows:
        resolved_units.setdefault((row["variant_id"], row["peptide_seq"]), []).append(row)
    positive_states = {"POSITIVE_STRONG", "POSITIVE_WEAK", "POSITIVE"}
    resolved_positive_units = sum(
        any(row["label_state"] in positive_states for row in rows)
        for rows in resolved_units.values())
    resolved_negative_units = sum(
        any(row["label_state"] == "NEGATIVE" for row in rows)
        for rows in resolved_units.values())
    resolved_contradictory_units = sum(
        any(row["label_state"] in positive_states for row in rows)
        and any(row["label_state"] == "NEGATIVE" for row in rows)
        for rows in resolved_units.values())
    def n_res(state):
        return sum(1 for r in exp_rows if r["resolution_state"] == state)
    def n_lab(state):
        return sum(1 for r in exp_rows if r["label_state"] == state)
    defensible_neg = sum(1 for r in exp_rows if r["label_state"] == "NEGATIVE"
                         and r["resolution_state"] in ("INDIVIDUAL_PEPTIDE", "MUTATION_LONG_PEPTIDE"))
    pool_neg = sum(1 for r in exp_rows if r["label_state"] == "NEGATIVE" and r["resolution_state"] == "POOL")
    site_pos_variants = {c["variant_id"] for c in catalog if c["elispot_positive_flag"]}
    hudson_recog = sorted({r["gene"] for r in hudson if r["is_mutation_specific"]})
    hudson_recog_variants = {r["target_id"] for r in funnel if r["hudson_recognized"]}
    overlap = sorted(site_pos_variants & hudson_recog_variants)
    hudson_tests = sorted({(r["timepoint"], r["mutation_label"]) for r in hudson})
    return {
        "unique_variants": len(catalog),
        "vaccine_targeted_variants": sum(1 for c in catalog if c["n_vaccines"] > 0),
        "site_elispot_positive_variants": len(site_pos_variants),
        "peptide_blocks": len(inventory),
        "assay_ledger_rows": len(ledger),
        "site_experiment_rows_tested": len(exp_rows),
        "individual_peptide_tests": n_res("INDIVIDUAL_PEPTIDE"),
        "long_peptide_tests": n_res("MUTATION_LONG_PEPTIDE"),
        "pool_tests": n_res("POOL"),
        "hudson_mutation_tcr_rows": len(hudson),
        "hudson_distinct_timepoint_mutation_tests": len(hudson_tests),
        "positives_strong": n_lab("POSITIVE_STRONG"),
        "positives_weak": n_lab("POSITIVE_WEAK"),
        "positives_unqualified": n_lab("POSITIVE"),
        "negatives_total": n_lab("NEGATIVE"),
        "defensible_negatives_individual_or_longpeptide": defensible_neg,
        "pool_negatives_not_perpeptide_defensible": pool_neg,
        "resolved_nonpool_rows": len(resolved_rows),
        "resolved_nonpool_positive_rows": sum(
            r["label_state"] in positive_states for r in resolved_rows),
        "resolved_nonpool_negative_rows": sum(
            r["label_state"] == "NEGATIVE" for r in resolved_rows),
        "resolved_unique_peptide_units": len(resolved_units),
        "resolved_unique_positive_units": resolved_positive_units,
        "resolved_unique_negative_units": resolved_negative_units,
        "resolved_unique_contradictory_units": resolved_contradictory_units,
        "ambiguous": n_lab("AMBIGUOUS"),
        "untested_vaccine_peptide_blocks": sum(1 for r in ledger if r["label_state"] == "UNTESTED"),
        "contradictions": len(contradictions),
        "hudson_recognized_genes": hudson_recog,
        "site_vs_hudson_overlap_variants": overlap,
        "recognized_first_failure_stage": {
            r["gene"]: r["first_failure_stage"] for r in funnel if r["hudson_recognized"]},
    }


def _render_report(summary, funnel, adjudication, provenance) -> str:
    L = []
    A = L.append
    A("# osteosarc.com (Sid) — public reconstruction REPORT\n")
    A("> Reconstructed from the public osteosarc.com site (182 variant pages + VAF TSVs) plus the "
      "local public pVACtools/RSEM/Hudson inputs, per the frozen preregistration "
      "(`docs/superpowers/specs/2026-07-12-osteosarc-sid-reconstruction-preregistration.md`). No model "
      "is fit, tuned, or compared here. Supersedes the `dd3efd1` diagnostic (assumed-negatives + "
      "single-positive AUROC — descriptive only).\n")
    s = summary
    A("## 1. Evidence-graded counts\n")
    A(f"- Unique site variants: **{s['unique_variants']}**  (vaccine-targeted **{s['vaccine_targeted_variants']}**, "
      f"site-ELISPOT-positive **{s['site_elispot_positive_variants']}**).")
    A(f"- Peptide blocks (long vaccine peptides): **{s['peptide_blocks']}**;  assay-ledger rows: "
      f"**{s['assay_ledger_rows']}** ({s['site_experiment_rows_tested']} real experiment rows + "
      f"{s['untested_vaccine_peptide_blocks']} UNTESTED vaccine peptides).")
    A(f"- Site ELISPOT tests by resolution: individual-peptide **{s['individual_peptide_tests']}**, "
      f"long-peptide **{s['long_peptide_tests']}**, pool **{s['pool_tests']}**.")
    A(f"- Site positives: strong **{s['positives_strong']}**, weak **{s['positives_weak']}**, "
      f"unqualified **{s['positives_unqualified']}**;  negatives total **{s['negatives_total']}** "
      f"(defensible individual/long-peptide **{s['defensible_negatives_individual_or_longpeptide']}**, "
      f"pool-only **{s['pool_negatives_not_perpeptide_defensible']}**);  ambiguous **{s['ambiguous']}**.")
    A(f"- Hudson IFNγ/TCR stream (SEPARATE modality): **{s['hudson_mutation_tcr_rows']}** clonotype rows "
      f"across **{s['hudson_distinct_timepoint_mutation_tests']}** (timepoint,mutation) tests; "
      f"mutation-specific recognized genes **{s['hudson_recognized_genes']}**.")
    A(f"- Contradictions (same peptide positive *and* negative across protocol/timepoint): "
      f"**{s['contradictions']}** — reported, never collapsed (see AUDIT.json).\n")
    A("## 2. Do the 14 site-ELISPOT-positive variants overlap the Hudson positives?\n")
    A(f"Overlap (by site coordinate) = **{s['site_vs_hudson_overlap_variants'] or 'none forced'}**. "
      "The two streams are different assays (peptide ELISPOT vs IFNγ/TCR expansion); overlap is reported "
      "without asserting equivalence.\n")
    A("## 3. True evidence-supported denominator for ranking today\n")
    A("The defensible within-patient evidence set is NOT the 21 curated mutations with the other 20 "
      "called negative. Excluding pool-only and untested rows leaves "
      f"**{s['resolved_nonpool_rows']} assay rows** "
      f"({s['resolved_nonpool_positive_rows']} positive, {s['resolved_nonpool_negative_rows']} negative) "
      f"across **{s['resolved_unique_peptide_units']} unique (variant, peptide) units**. At the unit level, "
      f"{s['resolved_unique_positive_units']} are ever-positive and {s['resolved_unique_negative_units']} "
      f"are ever-negative, with **{s['resolved_unique_contradictory_units']} in both groups** across "
      "protocol/timepoint. Therefore this is an evidence ledger, not a clean binary reranker denominator; "
      "ordinary AUROC requires a separately frozen longitudinal label policy.\n")
    A("## 4. Where recognized targets are lost (reachability)\n")
    A("| recognized gene | first failure stage in the automated funnel |")
    A("|---|---|")
    for gene, stage in sorted(s["recognized_first_failure_stage"].items()):
        A(f"| {gene} | `{stage}` |")
    A("\nAdjudications:\n")
    for key, txt in sorted(adjudication.items()):
        A(f"- **{key[0]} {key[1]}** — {txt}")
    A("\n## 5. Changes justified for Epicurus (proposed, NOT fit to this patient)\n")
    A("1. **Multi-caller / longitudinal union at candidate generation.** 2 of 3 Hudson-recognized "
      "neoantigens (ASPM, MAP2 Gly868fs) were *called by DRAGEN/Sarek/"
      "oncoanalyser* yet dropped by the single pVACtools 2025.01 candidate step. The recoverable loss is "
      "candidate **recall**, upstream of any ranker.")
    A("2. **No hard TPM/tier drop; carry low-evidence candidates with a flag** rather than filtering them "
      "out before ranking.")
    A("3. **Evidence tiers + honest abstention + diversity/uncertainty in the top-20**, since the "
      "recognition axis is unobserved for most candidates.")
    A("\nEach must be validated on the independent labeled cohorts (multimer/Gartner/IMPROVE/CheckMate), "
      "never on this single patient.\n")
    A("## Provenance\n")
    A(f"- {len(provenance)} fetched/hashed URLs (see PROVENANCE.json); rerun is network-free from "
      "`data/raw/osteosarc/site_cache/`.\n")
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------------------------
# Orchestrator
# --------------------------------------------------------------------------------------------

def build(*, offline: bool = False, refresh: bool = False, out_dir: Path = ART,
          raw_dir: Path = RAW, cache_dir: Path = SITE_CACHE) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    fetcher = Fetcher(cache_dir=cache_dir, offline=offline, refresh=refresh)

    # --- runtime site data ---
    index_rows = parse_index(fetcher.get(INDEX_URL))
    fetcher.get(VAFS_COLS_URL)
    vafs_rows = _parse_tsv(fetcher.get(VAFS_LONG_URL))
    vafs_idx = index_vafs_variants(vafs_rows)

    page_map = {}
    for ir in index_rows:
        slug = ir["slug"]
        page_map[slug] = parse_variant_page(fetcher.get(f"{BASE}/variant/{slug}/"), slug)

    # --- local public inputs ---
    hudson_paths = [
        ("May", "tumor", raw_dir / "May_all_expanders.tsv"),
        ("Aug", "tumor", raw_dir / "Aug_all_expanders.tsv"),
        ("May", "bystander", raw_dir / "May_bys_all_expanders.tsv"),
        ("Aug", "bystander", raw_dir / "Aug_bys_all_expanders.tsv"),
    ]
    hudson = load_hudson(hudson_paths)
    pv_present, pv_ensg = load_pvactools_genes(raw_dir / "pvactools_all_epitopes.tsv")
    curated = load_curated(raw_dir / "pvactools_curated_aggregated.tsv")
    tpm = load_rsem_tpm(raw_dir / "rsem.2025.01.genes.results")
    local_prov = _hash_local_inputs(raw_dir, hudson_paths)

    # --- assemble tables ---
    catalog, cat_mismatch = build_variant_catalog(index_rows, page_map, vafs_idx)
    inventory, len_mismatch = build_peptide_inventory(page_map)
    ledger = build_assay_ledger(page_map)
    funnel, adjudication = build_reachability(
        hudson, catalog, vafs_idx, pv_present, pv_ensg, curated, tpm)

    # --- integrity (fail-fast; no partial write) ---
    run_integrity_checks(catalog, inventory, ledger, hudson, funnel,
                         cat_mismatch + len_mismatch)

    # --- derive + write ---
    contradictions = _find_contradictions(ledger)
    exp_count_discrepancies = [
        {"variant_id": r["variant_id"], "block_index": r["block_index"],
         "peptide_seq": r["peptide_seq"], "peptide_len": r["peptide_len"],
         "declared": r["declared_experiment_count"], "parsed": r["parsed_experiment_count"],
         "note": "site peptide-aux under-counts: short individual P-series peptide row is rendered "
                 "in the block table but omitted from the block's declared experiment total"}
        for r in inventory if r["parsed_experiment_count"] > r["declared_experiment_count"]]
    summary = summarize(catalog, inventory, ledger, hudson, funnel, contradictions)

    write_csv(out_dir / "variant_catalog.csv", CATALOG_COLS,
              sorted(catalog, key=lambda r: r["variant_id"]))
    write_csv(out_dir / "peptide_inventory.csv", INVENTORY_COLS,
              sorted(inventory, key=lambda r: (r["variant_id"], r["block_index"])))
    write_csv(out_dir / "assay_ledger.csv", LEDGER_COLS,
              sorted(ledger, key=lambda r: r["experiment_key"]))
    write_csv(out_dir / "hudson_tcr_labels.csv", HUDSON_COLS,
              sorted(hudson, key=lambda r: (r["timepoint"], r["pool_kind"], r["mutation_label"], r["trb"])))
    write_csv(out_dir / "reachability_funnel.csv", FUNNEL_COLS, funnel)

    audit = {
        "role": "PUBLIC recognition-record reconstruction (evidence-graded ledger). Not a model, not a "
                "gate. Supersedes dd3efd1 (assumed-negatives + single-positive AUROC).",
        "invariants": {"variants": EXPECT_VARIANTS, "vaccine_targeted": EXPECT_VACCINE,
                       "site_elispot_positive": EXPECT_ELISPOT, "all_passed": True},
        "summary": summary,
        "contradictions": contradictions,
        "site_declared_experiment_count_discrepancies": exp_count_discrepancies,
        "adjudications": {f"{k[0]} {k[1]}": v for k, v in sorted(adjudication.items())},
        "label_states": sorted(LABEL_STATES),
        "resolution_states": sorted(RESOLUTION_STATES),
    }
    (out_dir / "AUDIT.json").write_text(json.dumps(audit, indent=2, default=str) + "\n")
    provenance = {"runtime_urls": fetcher.provenance, "local_inputs": local_prov}
    (out_dir / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2, default=str) + "\n")
    (out_dir / "REPORT.md").write_text(
        _render_report(summary, funnel, adjudication, fetcher.provenance))

    return {"summary": summary, "contradictions": contradictions, "out_dir": str(out_dir),
            "n_provenance": len(fetcher.provenance)}


def _parse_tsv(text: str) -> list[dict]:
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def _hash_local_inputs(raw_dir: Path, hudson_paths) -> dict:
    out = {}
    files = ["pvactools_all_epitopes.tsv", "pvactools_curated_aggregated.tsv",
             "rsem.2025.01.genes.results"] + [p.name for _, _, p in hudson_paths]
    for name in files:
        fp = raw_dir / name
        if fp.exists():
            out[name] = {"bytes": fp.stat().st_size,
                         "sha256": hashlib.sha256(fp.read_bytes()).hexdigest()}
    return out
