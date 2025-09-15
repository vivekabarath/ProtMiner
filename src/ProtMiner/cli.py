#!/usr/bin/env python3
import argparse, re, sys, shutil, subprocess as sp
from pathlib import Path
import pandas as pd
from Bio import SeqIO

from .utils import have, which, run_logged, ensure_dir
from .install import ensure_python_modules, try_install
from .parsers import parse_hmmer_domtbl, parse_diamond_tsv, parse_interpro_tsv
from .scoring import score_and_filter

PFAM_RE = re.compile(r"^PF\d{5}$")
IPR_RE  = re.compile(r"^IPR\d{6}$")

def lengths_of_fasta(fasta: Path) -> pd.DataFrame:
    return pd.DataFrame([{"Protein_ID": r.id, "Length": len(r.seq)} for r in SeqIO.parse(str(fasta), "fasta")])

def write_fasta_subset(all_faa: Path, ids, out_faa: Path):
    idset = set(map(str, ids))
    with out_faa.open("w") as w:
        for rec in SeqIO.parse(str(all_faa), "fasta"):
            if rec.id in idset:
                SeqIO.write(rec, w, "fasta")

def scan_motif(fasta: Path, motif: str):
    rx = re.compile(motif)
    hits = {}
    for r in SeqIO.parse(str(fasta), "fasta"):
        s = str(r.seq)
        for m in rx.finditer(s):
            hits.setdefault(r.id, []).append((m.start()+1, m.end(), m.group(0)))
    return {k: ";".join([f"{a}-{b}:{s}" for (a,b,s) in v]) for k,v in hits.items()}

def parse_id_list(s: str, pat: re.Pattern, label: str):
    if not s: return set()
    vals = {x.strip() for x in s.split(",") if x.strip()}
    bad  = [x for x in vals if not pat.match(x)]
    if bad:
        sys.stderr.write(f"[WARN] Some {label} look malformed and will be ignored: {', '.join(bad)}\n")
        vals = {x for x in vals if pat.match(x)}
    return vals

def check_env():
    tools = ["hmmsearch","diamond","interproscan.sh","seqkit","tmhmm","deeptmhmm","signalp"]
    print("Python modules:")
    try:
        import pandas; print("  [OK] pandas")
    except Exception as e:
        print("  [MISSING] pandas", e)
    try:
        import Bio.SeqIO; print("  [OK] biopython")
    except Exception as e:
        print("  [MISSING] biopython", e)
    print("\nExternal tools:")
    for t in tools:
        p = which(t)
        status = "[OK]" if p else "[MISSING]"
        print(f"  {status} {t} {('-> '+p) if p else ''}")

def main():
    ap = argparse.ArgumentParser(description="ProtScout — inclusive class-aware protein discovery")
    ap.add_argument("--check-env", action="store_true", help="Print environment diagnostics and exit")

    # core
    ap.add_argument("--fasta", help="Input protein FASTA")
    ap.add_argument("--outdir", help="Output directory")
    ap.add_argument("--hmm", nargs="+", help="One or more HMM files for target domain(s)")
    ap.add_argument("--label", default="Target proteins", help="Label for the run")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ref-fasta", help="Curated reference FASTA for DIAMOND DB")
    ap.add_argument("--motif", help="Regex motif to check (e.g., 'H[ST]GTT[ST]')")
    ap.add_argument("--auto-install-py", action="store_true", help="pip-install pandas/biopython if missing")
    ap.add_argument("--auto-install-tools", action="store_true", help="Try installing external tools via conda/apt/brew")

    # thresholds
    ap.add_argument("--min-hmm-cov", type=float, default=0.50)
    ap.add_argument("--max-hmm-ie",  type=float, default=1e-5)
    ap.add_argument("--min-diam-pident", type=float, default=25.0)
    ap.add_argument("--min-diam-qcov",   type=float, default=60.0)
    ap.add_argument("--min-diam-scov",   type=float, default=60.0)
    ap.add_argument("--require-two-evidences", action="store_true")

    # class-aware allow lists
    ap.add_argument("--allow-pfam", help="Comma-separated Pfam IDs (e.g., PF00201)")
    ap.add_argument("--allow-ipr",  help="Comma-separated InterPro IDs (e.g., IPR001296)")

    # length options
    ap.add_argument("--len-min", type=int, help="Minimum AA length (optional)")
    ap.add_argument("--len-max", type=int, help="Maximum AA length (optional)")

    # topology/signal (metadata by default; disabled in inclusive scoring)
    ap.add_argument("--tm-min", type=int, help="Min TM helices")
    ap.add_argument("--tm-max", type=int, help="Max TM helices")
    ap.add_argument("--tm-exact", type=int, help="Exact TM helices")
    ap.add_argument("--tm-nterm-window", type=int, help="A helix starts within first N aa (TMHMM topology required)")
    ap.add_argument("--tm-cterm-window", type=int, help="A helix ends within last N aa (TMHMM topology required)")
    ap.add_argument("--tm-first-helix-nterm", action="store_true", help="First helix must be N-terminal (TMHMM)")
    ap.add_argument("--tm-predictor", choices=["auto","tmhmm","deeptmhmm"], default="auto", help="Which TM predictor to use")

    spg = ap.add_mutually_exclusive_group()
    spg.add_argument("--sp-required", action="store_true", help="Require a SignalP hit")
    spg.add_argument("--sp-forbidden", action="store_true", help="Require absence of SignalP hit")

    ap.add_argument("--topology-metadata-only", action="store_true",
                    help="Always treat TM/SignalP as metadata only (never as evidence)")

    ap.add_argument("--scoring-mode", choices=["inclusive","thresholded"], default="inclusive",
                    help="inclusive: keep any positive (default). thresholded: classic rules.")

    args = ap.parse_args()

    if args.check_env:
        check_env()
        return

    # Required for a run
    for req in ["fasta","outdir","hmm"]:
        if getattr(args, req) is None:
            ap.error(f"--{req} is required unless --check-env is used")

    # Python deps
    ensure_python_modules(args.auto_install_py)

    outdir = Path(args.outdir); ensure_dir(outdir)
    log = outdir / "pipeline.log"
    with log.open("w") as L:
        L.write(f"# ProtScout — {args.label}\n")

    fasta = Path(args.fasta)
    if not fasta.exists():
        sys.stderr.write(f"[ERROR] FASTA not found: {fasta}\n"); sys.exit(2)
    for hmm in args.hmm:
        if not Path(hmm).exists():
            sys.stderr.write(f"[ERROR] HMM not found: {hmm}\n"); sys.exit(2)

    # External tools: check/try-install
    for tool in ["hmmsearch","diamond","interproscan.sh","seqkit","tmhmm","signalp","deeptmhmm"]:
        try_install(tool, args.auto_install_tools)

    have_hmm   = have("hmmsearch")
    have_diam  = have("diamond")
    have_ipr   = have("interproscan.sh")
    have_seq   = have("seqkit")
    have_tmh   = have("tmhmm")
    have_dtmh  = shutil.which("deeptmhmm") is not None
    have_sigp  = have("signalp")

    # Deduplicate
    dedup_faa = outdir / "proteins.dedup.faa"
    if have_seq:
        run_logged(["seqkit","rmdup","-s", str(fasta)], log)
        sp.run(f"seqkit rmdup -s {fasta} > {dedup_faa}", shell=True)
    else:
        dedup_faa.write_text(fasta.read_text())

    # HMMER
    if not have_hmm:
        sys.stderr.write("[ERROR] hmmsearch is required for domain discovery.\n"); sys.exit(4)

    all_hmm = []
    for hmm in args.hmm:
        domtbl = outdir / (Path(hmm).stem + ".domtblout")
        run_logged(["hmmsearch", "--cut_ga", "--domtblout", str(domtbl), str(Path(hmm).resolve()), str(dedup_faa)], log)
        df = parse_hmmer_domtbl(domtbl)
        if not df.empty:
            df = df.rename(columns={"dom_cov_hmm":"HMM_cov","i_evalue":"HMM_ie","dom_score":"HMM_bits","query":"HMM_model"})
            all_hmm.append(df[["Protein_ID","tlen","qlen","HMM_model","HMM_cov","HMM_ie","HMM_bits"]])
    hmm_df = pd.concat(all_hmm, ignore_index=True) if all_hmm else pd.DataFrame(
        columns=["Protein_ID","HMM_model","HMM_cov","HMM_ie","HMM_bits","tlen","qlen"]
    )
    hmm_df = hmm_df.sort_values(["Protein_ID","HMM_cov","HMM_bits"], ascending=[True, False, False])\
                   .drop_duplicates("Protein_ID")

    # DIAMOND (optional)
    diam_df = pd.DataFrame()
    if args.ref_fasta and have_diam:
        db = outdir / "ref.dmnd"
        run_logged(["diamond","makedb","--in", str(Path(args.ref_fasta).resolve()), "-d", str(db)], log)
        diam_tsv = outdir / "diamond.tsv"
        run_logged([
            "diamond","blastp",
            "-q", str(dedup_faa),
            "-d", str(db),
            "-o", str(diam_tsv),
            "--outfmt","6","qseqid","sseqid","pident","length","evalue","bitscore","qcovhsp","scovhsp","stitle",
            "--max-target-seqs","5","--evalue","1e-5","--threads",str(args.threads)
        ], log)
        diam_df = parse_diamond_tsv(diam_tsv)
    elif args.ref_fasta and not have_diam:
        sys.stderr.write("[WARN] --ref-fasta given but diamond missing; skipping DIAMOND.\n")

    # InterProScan (optional)
    ipr_df = pd.DataFrame()
    if have_ipr:
        ipr_tsv = outdir / "interpro.tsv"
        run_logged(["interproscan.sh","-i",str(dedup_faa),"-f","tsv","-appl","Pfam,SMART,CDD","--goterms","--iprlookup","-o",str(ipr_tsv)], log)
        ipr_df = parse_interpro_tsv(ipr_tsv)

    # TM predictor (metadata unless thresholded with filters)
    tm_df = pd.DataFrame()
    if have_tmh:
        tm_out = outdir / "tmhmm.out"
        sp.run(f"tmhmm {dedup_faa} > {tm_out}", shell=True)
        rows = []
        with tm_out.open() as f:
            for line in f:
                if not line.strip() or line.startswith("#"): continue
                parts = line.strip().split()
                prot = parts[0]; tm_pred=None; topo=None
                for tok in parts:
                    if tok.startswith("PredHel="): tm_pred = int(tok.split("=")[1])
                    if tok.startswith("Topology="): topo = tok.split("=",1)[1]
                rows.append({"Protein_ID": prot, "TM_pred": tm_pred, "Topology": topo})
        tm_df = pd.DataFrame(rows)
    elif have_dtmh:
        sp.run(f"deeptmhmm --fasta {dedup_faa} --outdir {outdir} --format tsv", shell=True)
        tsv = outdir / "deeptmhmm.tsv"
        if tsv.exists():
            with tsv.open() as f:
                header = f.readline().strip().split("\t")
                id_idx = header.index("protein_id") if "protein_id" in header else 0
                hel_idx = header.index("number_of_tmh") if "number_of_tmh" in header else 2
                rows = []
                for line in f:
                    if not line.strip(): continue
                    parts = line.strip().split("\t")
                    prot = parts[id_idx]
                    val = parts[hel_idx]
                    tm_pred = int(val) if val.isdigit() else None
                    rows.append({"Protein_ID": prot, "TM_pred": tm_pred, "Topology": pd.NA})
            tm_df = pd.DataFrame(rows)

    # SignalP (optional; metadata)
    sp_df = pd.DataFrame()
    if have_sigp:
        prefix = outdir / "signalp"
        run_logged(["signalp","-fasta", str(dedup_faa), "-org","euk","-format","short","-verbose","-prefix", str(prefix)], log)
        short_files = list(outdir.glob("signalp*summary*")) + list(outdir.glob("signalp*.short*"))
        rows = []
        for p in short_files:
            with p.open() as f:
                for line in f:
                    if not line.strip() or line.startswith("#"): continue
                    parts = line.strip().split()
                    rows.append({"Protein_ID": parts[0], "SignalP": "Yes"})
        sp_df = pd.DataFrame(rows if rows else [], columns=["Protein_ID","SignalP"])

    # Motif
    motif_series = None
    if args.motif:
        rx = re.compile(args.motif)
        motif_hits = {}
        for r in SeqIO.parse(str(dedup_faa), "fasta"):
            s = str(r.seq)
            for m in rx.finditer(s):
                motif_hits.setdefault(r.id, []).append((m.start()+1, m.end(), m.group(0)))
        motif_series = {k: ";".join([f"{a}-{b}:{s}" for (a,b,s) in v]) for k,v in motif_hits.items()}

    # Merge
    len_df = lengths_of_fasta(dedup_faa)
    master = len_df.copy()
    for d in [hmm_df, diam_df, ipr_df, tm_df, sp_df]:
        if not d.empty:
            master = master.merge(d, on="Protein_ID", how="left")
    if motif_series:
        ms = pd.Series(motif_series, name="Motif_hits"); ms.index.name = "Protein_ID"
        master = master.merge(ms.reset_index(), on="Protein_ID", how="left")

    # Ensure numeric columns exist and cast
    for col in ["HMM_cov","HMM_ie","HMM_bits","pident","qcovhsp","scovhsp","bitscore","TM_pred","Length"]:
        if col not in master.columns: master[col] = pd.NA
        master[col] = pd.to_numeric(master[col], errors="coerce")

    # Parse allow lists
    allow_pfam = parse_id_list(args.allow_pfam or "", re.compile(r"^PF\d{5}$"), "Pfam IDs")
    allow_ipr  = parse_id_list(args.allow_ipr  or "", re.compile(r"^IPR\d{6}$"), "InterPro IDs")

    # Score
    scored = score_and_filter(
        df=master,
        min_hmm_cov=args.min_hmm_cov, max_hmm_ie=args.max_hmm_ie,
        min_pident=args.min_diam_pident, min_qcov=args.min_diam_qcov, min_scov=args.min_diam_scov,
        require_two_evidences=args.require_two_evidences,
        allow_pfam=allow_pfam, allow_ipr=allow_ipr,
        len_min=args.len_min, len_max=args.len_max,
        tm_min=args.tm_min, tm_max=args.tm_max, tm_exact=args.tm_exact,
        tm_nterm_window=args.tm_nterm_window, tm_cterm_window=args.tm_cterm_window,
        tm_first_helix_nterm=args.tm_first_helix_nterm,
        sp_required=args.sp_required, sp_forbidden=args.sp_forbidden,
        topology_metadata_only=args.topology_metadata_only,
        scoring_mode=args.scoring_mode
    )

    # Outputs
    all_tsv = Path(outdir / "all_proteins_annotation.tsv"); scored.to_csv(all_tsv, sep="\t", index=False)
    cand = scored[scored["keep_candidate"] == True].copy()
    cand = cand.sort_values(["HMM_cov","HMM_bits","pident"], ascending=[False, False, False])
    cand_tsv = Path(outdir / "candidates.tsv"); cand.to_csv(cand_tsv, sep="\t", index=False)

    keep_ids = cand["Protein_ID"].astype(str).tolist()
    if keep_ids:
        write_fasta_subset(dedup_faa, keep_ids, outdir / "candidates.faa")

    total = len(master)
    print(f"[DONE] {args.label}: {total} sequences scanned; {len(keep_ids)} candidate(s).")
    print(f"- Candidates table: {cand_tsv}")
    print(f"- Candidates FASTA: {outdir/'candidates.faa'}")
    print(f"- Full annotations: {all_tsv}")
    print(f"- Log: {log}")

if __name__ == "__main__":
    main()
