from pathlib import Path
import pandas as pd

def parse_hmmer_domtbl(domtbl: Path) -> pd.DataFrame:
    rows = []
    if not domtbl.exists() or domtbl.stat().st_size == 0:
        return pd.DataFrame(columns=[
            "Protein_ID","query","qlen","tlen","i_evalue","dom_score",
            "hmm_from","hmm_to","ali_from","ali_to","env_from","env_to","dom_cov_hmm","ali_len"
        ])
    with domtbl.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            try:
                target = parts[0]; tlen = int(parts[2]); query = parts[3]; qlen = int(parts[5])
                i_eval = float(parts[12]); dom_score = float(parts[13])
                hmm_from = int(parts[15]); hmm_to = int(parts[16])
                ali_from = int(parts[17]); ali_to = int(parts[18])
                env_from = int(parts[19]); env_to = int(parts[20])
                dom_cov = (hmm_to - hmm_from + 1) / qlen if qlen > 0 else 0.0
                ali_len = (ali_to - ali_from + 1)
                rows.append({
                    "Protein_ID": target, "query": query, "qlen": qlen, "tlen": tlen,
                    "i_evalue": i_eval, "dom_score": dom_score,
                    "hmm_from": hmm_from, "hmm_to": hmm_to,
                    "ali_from": ali_from, "ali_to": ali_to,
                    "env_from": env_from, "env_to": env_to,
                    "dom_cov_hmm": dom_cov, "ali_len": ali_len
                })
            except Exception:
                continue
    df = pd.DataFrame(rows)
    if df.empty: return df
    df = df.sort_values(["Protein_ID","dom_cov_hmm","dom_score"], ascending=[True, False, False])\
           .drop_duplicates(subset=["Protein_ID"], keep="first")
    return df

def parse_diamond_tsv(tsv: Path) -> pd.DataFrame:
    cols = ["Protein_ID","Subject","pident","aln_len","evalue","bitscore","qcovhsp","scovhsp","stitle"]
    if not tsv.exists() or tsv.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(tsv, sep="\t", names=cols)
    df = df.sort_values(["Protein_ID","bitscore","qcovhsp","pident"], ascending=[True, False, False, False])\
           .drop_duplicates("Protein_ID")
    for num in ["pident","aln_len","evalue","bitscore","qcovhsp","scovhsp"]:
        df[num] = pd.to_numeric(df[num], errors="coerce")
    return df

def parse_interpro_tsv(tsv: Path) -> pd.DataFrame:
    if not tsv.exists() or tsv.stat().st_size == 0:
        return pd.DataFrame(columns=["Protein_ID","Signatures","SignatureDescs","InterProAccs","InterProDescs","GO_terms"])
    cols = ["Protein_ID","MD5","Length","Analysis","SignatureAcc","SignatureDesc",
            "Start","End","Score","Status","Date","InterProAcc","InterProDesc","GO","Pathways"]
    df = pd.read_csv(tsv, sep="\t", header=None, names=cols, dtype=str)
    agg = df.groupby("Protein_ID").agg({
        "SignatureAcc": lambda x: ";".join(sorted(set(v for v in x if isinstance(v, str)))),
        "SignatureDesc": lambda x: ";".join(sorted(set(v for v in x if isinstance(v, str)))),
        "InterProAcc": lambda x: ";".join(sorted(set(v for v in x if isinstance(v, str)))),
        "InterProDesc": lambda x: ";".join(sorted(set(v for v in x if isinstance(v, str)))),
        "GO": lambda x: ";".join(sorted(set(";".join(v for v in x if isinstance(v, str)).split("|"))))
    }).reset_index()
    return agg.rename(columns={
        "SignatureAcc":"Signatures",
        "SignatureDesc":"SignatureDescs",
        "InterProAcc":"InterProAccs",
        "InterProDesc":"InterProDescs",
        "GO":"GO_terms"
    })
