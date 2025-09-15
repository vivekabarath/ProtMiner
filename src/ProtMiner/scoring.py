import re
import pandas as pd

def _ensure_cols(d: pd.DataFrame, cols):
    for c in cols:
        if c not in d.columns:
            d[c] = pd.NA
    return d

def _bool(x):
    return x.fillna(False).astype(bool)

def _tmhmm_ranges(topology: str):
    if not isinstance(topology, str) or not topology:
        return []
    toks = topology.strip().split()
    helices = []
    for i in range(0, len(toks)-1, 2):
        label = toks[i]
        rng = toks[i+1]
        if label.upper().startswith("TMHELIX"):
            m = re.match(r"(\d+)-(\d+)", rng)
            if m:
                s, e = map(int, m.groups())
                helices.append((s, e))
    return helices

def _evidence_sources_row(row):
    parts = []
    if row.get("ev_HMM", False): parts.append("HMM")
    if row.get("ev_IPR", False): parts.append("IPR")
    if row.get("ev_DIAM", False): parts.append("DIAMOND")
    if row.get("ev_MOTIF", False): parts.append("MOTIF")
    if row.get("ev_LEN", False): parts.append("LEN")
    if row.get("ev_TM", False): parts.append("TM")
    if row.get("ev_SP", False): parts.append("SP")
    return ";".join(parts)

def score_and_filter(df: pd.DataFrame,
                     min_hmm_cov: float, max_hmm_ie: float,
                     min_pident: float, min_qcov: float, min_scov: float,
                     require_two_evidences: bool,
                     allow_pfam=set(), allow_ipr=set(),
                     len_min=None, len_max=None,
                     tm_min=None, tm_max=None, tm_exact=None,
                     tm_nterm_window=None, tm_cterm_window=None,
                     tm_first_helix_nterm=False,
                     sp_required=False, sp_forbidden=False,
                     topology_metadata_only: bool=False,
                     scoring_mode: str="inclusive") -> pd.DataFrame:

    d = df.copy()
    d = _ensure_cols(d, ["HMM_cov","HMM_ie","Signatures","InterProAccs","Motif_hits",
                         "pident","qcovhsp","scovhsp","Length","TM_pred","Topology","SignalP"])

    # Cast numeric
    for col in ["HMM_cov","HMM_ie","pident","qcovhsp","scovhsp","Length","TM_pred"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    # Evidence definitions
    d["ev_HMM"] = (d["HMM_cov"].fillna(0) >= min_hmm_cov) & (d["HMM_ie"].fillna(1) <= max_hmm_ie)

    have_diam = all(c in d.columns for c in ["pident","qcovhsp","scovhsp"])
    if have_diam:
        d["ev_DIAM"] = (d["pident"].fillna(0) >= min_pident) &                        (d["qcovhsp"].fillna(0) >= min_qcov) &                        (d["scovhsp"].fillna(0) >= min_scov)
    else:
        d["ev_DIAM"] = False

    def contains_any(series: pd.Series, keys: set):
        if not keys:
            return pd.Series(False, index=series.index)
        pat = "|".join(sorted(keys))
        return series.fillna("").str.contains(pat)

    if allow_pfam or allow_ipr:
        d["ev_IPR"] = _bool(contains_any(d["Signatures"], allow_pfam) | contains_any(d["InterProAccs"], allow_ipr))
    else:
        d["ev_IPR"] = d["Signatures"].fillna("").ne("") | d["InterProAccs"].fillna("").ne("")

    d["ev_MOTIF"] = d["Motif_hits"].fillna("").ne("")

    if len_min is not None or len_max is not None:
        ok_len = pd.Series(True, index=d.index)
        if len_min is not None: ok_len &= (d["Length"].fillna(0) >= len_min)
        if len_max is not None: ok_len &= (d["Length"].fillna(10**9) <= len_max)
        d["ev_LEN"] = _bool(ok_len)
    else:
        d["ev_LEN"] = False

    if topology_metadata_only:
        d["ev_TM"] = False
        d["ev_SP"] = False
    else:
        d["ev_TM"] = False
        d["ev_SP"] = False

    evidence_cols = ["ev_HMM","ev_DIAM","ev_IPR","ev_MOTIF","ev_LEN","ev_TM","ev_SP"]
    d["evidence_count"] = d[evidence_cols].sum(axis=1)

    if scoring_mode == "inclusive":
        keep = d[["ev_HMM","ev_DIAM","ev_IPR","ev_MOTIF"]].any(axis=1)
    else:
        if require_two_evidences:
            keep = d["ev_HMM"] & (d["evidence_count"] >= 2)
        else:
            keep = d[evidence_cols].any(axis=1)

    d["keep_candidate"] = keep
    d["evidence_sources"] = d.apply(_evidence_sources_row, axis=1)

    # Ensure uniqueness
    d = d.drop_duplicates(subset=["Protein_ID"], keep="first")
    return d
