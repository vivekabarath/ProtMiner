# README.md

# ProtMiner 🧬🧭

> **Inclusive protein family discovery** using HMMER, DIAMOND, InterPro, and motif evidence — with optional topology/signal peptide metadata.

&#x20;&#x20;

---

## Table of Contents

- [Overview](#overview)
- [Installation](#-installation)
- [Input Requirements](#-input-requirements)
- [Usage Examples](#-usage-examples)
  - [Inclusive Mode (Default)](#1-inclusive-mode-default)
  - [Thresholded Mode (Stricter)](#2-thresholded-mode-stricter)
  - [Add DIAMOND Evidence](#3-add-diamond-evidence)
  - [Add Motif Evidence](#4-add-motif-evidence)
  - [Apply Size Restriction](#5-apply-size-restriction)
  - [Topology Filters](#6-topology-filters)
  - [Signal Peptides](#7-signal-peptides)
  - [Maximum Parameters Example](#8-maximum-parameters-all-features-enabled)
- [Outputs](#-output-files)
- [Example Workflows](#-example-workflows)
- [Best Practices](#-best-practices)
- [Code Structure & Extensibility](#-code-structure--extensibility)
- [Roadmap](#-roadmap)

---

## Overview

ProtMiner identifies proteins belonging to a given family by combining multiple evidence sources:

- **HMMER** (Pfam / custom HMMs)
- **DIAMOND** (reference similarity search)
- **InterProScan** (Pfam/SMART/CDD/IPR signatures)
- **Motif search** (regex)

Optional metadata filters:

- **Transmembrane topology** (TMHMM / DeepTMHMM)
- **Signal peptides** (SignalP)
- **Protein length ranges**

---

## 🚀 Installation

### Clone and install

```bash
git clone https://github.com/<yourusername>/ProtMiner.git
cd ProtMiner

# Create environment
conda env create -f environment.yml
conda activate protminer

# Install package
pip install -e .
```

### Verify environment

```bash
find-proteins --check-env
```

This prints which Python modules and external tools (`hmmsearch`, `diamond`, `interproscan.sh`, `seqkit`, `tmhmm`, `deeptmhmm`, `signalp`) are detected.

> **Note:** If `deeptmhmm` is unavailable for your conda mirror, create the env without it and add later (or install via `pip`).

---

## 📂 Input Requirements

- **Protein FASTA** (`--fasta`)
- **HMM profile(s)** (`--hmm`)
- *(Optional)* **Reference FASTA** for DIAMOND (`--ref-fasta`)
- *(Optional)* **Regex motif** (`--motif`)
- *(Optional)* Allowed **Pfam** (`--allow-pfam`) and **InterPro** (`--allow-ipr`) IDs

**Tip:** InterPro evidence only counts for IDs you list via `--allow-ipr` (and Pfam via `--allow-pfam`). This keeps unrelated IPR/Pfam hits from producing false positives.

---

## 🛠 Usage Examples

### 1. Inclusive Mode (Default)

Keep any protein with **HMM, InterPro, DIAMOND, or motif** evidence.

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_ugt \
  --hmm PF00201.hmm \
  --allow-pfam PF00201 \
  --allow-ipr IPR001296,IPR035595 \
  --scoring-mode inclusive
```

### 2. Thresholded Mode (Stricter)

Require HMM hit plus ≥2 evidence sources if `--require-two-evidences` is set.

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_strict \
  --hmm PF00067.hmm \
  --allow-pfam PF00067 \
  --allow-ipr IPR002401 \
  --require-two-evidences \
  --scoring-mode thresholded
```

### 3. Add DIAMOND Evidence

Use curated reference proteins to increase confidence.

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_diamond \
  --hmm PF00201.hmm \
  --ref-fasta curated_references.faa \
  --allow-pfam PF00201 \
  --allow-ipr IPR001296 \
  --scoring-mode inclusive
```

### 4. Add Motif Evidence

Regex pattern (not PROSITE syntax).

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_motif \
  --hmm PF00201.hmm \
  --allow-pfam PF00201 \
  --allow-ipr IPR001296 \
  --motif '[FW].{2}[QL].{2}[LIVMYA][LIMV].{4,6}[LVGAC][LVFYAHM][LIVMF][STAGCM][HNQ][STAGC]G.{2}[STAG].{3}[STAGL][LIVMFA].{4,5}[PQR][LIVMTA].{3}[PA].{2,3}[DES][QEHNR]' \
  --scoring-mode inclusive
```

### 5. Apply Size Restriction

Filter candidates by length.

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_sized \
  --hmm PF00201.hmm \
  --len-min 200 --len-max 600 \
  --scoring-mode inclusive
```

### 6. Topology Filters

**Soluble families (no TM domains):**

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_soluble \
  --hmm PF00201.hmm \
  --tm-max 0 \
  --scoring-mode thresholded
```

**ER-anchored CYPs (single TM helix):**

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_cyp \
  --hmm PF00067.hmm \
  --tm-exact 1 \
  --scoring-mode thresholded
```

### 7. Signal Peptides

Require or forbid N-terminal signal peptide.

```bash
# Require SP
find-proteins \
  --fasta proteome.faa \
  --outdir out_sp \
  --hmm PF00201.hmm \
  --sp-required

# Forbid SP
find-proteins \
  --fasta proteome.faa \
  --outdir out_no_sp \
  --hmm PF00201.hmm \
  --sp-forbidden
```

### 8. Maximum Parameters (All Features Enabled)

A comprehensive run combining **HMM, Pfam, InterPro, DIAMOND, motif, length restriction, TM topology, and SignalP**:

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_full \
  --hmm PF00201.hmm \
  --allow-pfam PF00201 \
  --allow-ipr IPR001296,IPR035595,IPR050271,IPR002213 \
  --ref-fasta curated_references.faa \
  --motif '[FW].{2}[QL].{2}[LIVMYA][LIMV].{4,6}[LVGAC][LVFYAHM][LIVMF][STAGCM][HNQ][STAGC]G.{2}[STAG].{3}[STAGL][LIVMFA].{4,5}[PQR][LIVMTA].{3}[PA].{2,3}[DES][QEHNR]' \
  --len-min 200 --len-max 600 \
  --tm-exact 1 \
  --sp-required \
  --require-two-evidences \
  --scoring-mode thresholded
```

This setup:

- Combines evidence sources (HMM, IPR/Pfam, DIAMOND, motif)
- Restricts to expected length range
- Requires exactly 1 TM helix and a signal peptide
- Uses stricter thresholding (≥2 evidences alongside HMM)

---

## 📊 Output Files

- `all_proteins_annotation.tsv` — all evidence merged
- `candidates.tsv` — final positives (deduplicated)
- `candidates.faa` — FASTA of candidates
- `pipeline.log` — commands and logs

Each row includes:

- `Protein_ID`
- Evidence flags (`ev_HMM`, `ev_DIAM`, `ev_IPR`, `ev_MOTIF`)
- `evidence_sources` (summary of which sources triggered inclusion)
- Metadata: `Length`, `TM_pred`, `SignalP`

---

## 🧪 Example Workflows

**UGTs (inclusive evidence-only, with motif + size restriction)**

```bash
find-proteins \
  --fasta plant_proteome.faa \
  --outdir out_ugt \
  --hmm PF00201.hmm \
  --allow-pfam PF00201 \
  --allow-ipr IPR001296,IPR035595,IPR050271 \
  --motif '[FW].{2}[QL].{2}[LIVMYA]...' \
  --len-min 380 --len-max 520 \
  --scoring-mode inclusive
```

**CYP450s (thresholded with topology)**

```bash
find-proteins \
  --fasta proteome.faa \
  --outdir out_cyp \
  --hmm PF00067.hmm \
  --allow-pfam PF00067 \
  --allow-ipr IPR002401 \
  --tm-exact 1 \
  --require-two-evidences \
  --scoring-mode thresholded
```

---

## ⚠️ Best Practices

- Always include at least one **HMM profile**.
- Use `--allow-pfam` / `--allow-ipr` to restrict InterPro evidence to your family.
- Motifs must be **regex** (convert PROSITE manually).
- TMHMM requires registration; prefer **DeepTMHMM** via bioconda.
- SignalP is license-restricted; enable only if installed.
- Add `--topology-metadata-only` to treat TM/SignalP as annotations only.
- Start with **inclusive mode**; move to `thresholded` for stricter filtering.

---

## 🧱 Code Structure & Extensibility

### Repository layout

```
ProtMiner/
├─ README.md
├─ environment.yml
├─ src/
│  └─ protminer/
│     ├─ __init__.py
│     ├─ cli.py            # command‑line interface & pipeline orchestrator
│     ├─ scoring.py        # evidence logic (inclusive/thresholded) & final filtering
│     ├─ parsers.py        # readers for HMMER/DIAMOND/InterPro outputs
│     ├─ utils.py          # helpers: which(), have(), run_logged(), ensure_dir()
│     └─ install.py        # dependency checks & (optional) auto‑install hints
├─ examples/
│  ├─ proteins.demo.faa
│  └─ PF00201.demo.hmm
└─ docs/ (optional site)
```

### `cli.py` — pipeline orchestrator

**What it does:**

- Parses CLI flags (input FASTA, HMMs, motif, DIAMOND reference, IPR/Pfam allow‑lists, topology/signal/length options, scoring mode).
- Runs environment checks and *optional* tool discovery/installation hints (`install.try_install`).
- Executes stages:
  1. **Deduplicate FASTA** (via `seqkit rmdup` if available).
  2. **HMMER** (`hmmsearch --domtblout`) → parse with `parsers.parse_hmmer_domtbl`.
  3. **DIAMOND** (if `--ref-fasta`) → make DB + `diamond blastp` → parse with `parsers.parse_diamond_tsv`.
  4. **InterProScan** (if available) → TSV → `parsers.parse_interpro_tsv`.
  5. **TM predictors / SignalP** (metadata unless you choose thresholded filters).
  6. **Motif** regex scan (inline), producing `Motif_hits`.
  7. **Merge** per‑protein evidence into a master `DataFrame`.
  8. **Score & filter** via `scoring.score_and_filter`.
  9. **Emit outputs** (`all_proteins_annotation.tsv`, `candidates.tsv`, `candidates.faa`, `pipeline.log`).

**Key functions & I/O:**

- `lengths_of_fasta(fasta)->DataFrame`: `Protein_ID`, `Length`.
- `write_fasta_subset(all_faa, ids, out_faa)`: exports candidate sequences.
- `scan_motif(fasta, motif)->dict`: `{Protein_ID: "start-end:motif;..."}`.

**Extend it:**

- To add a new analysis stage, compute a `DataFrame` with a `Protein_ID` column and merge it into `master` before calling `score_and_filter`.

### `scoring.py` — evidence logic

**Core entry point:** `score_and_filter(df, ..., scoring_mode="inclusive") -> DataFrame`

**Inputs (selected):**

- Evidence columns expected (if present): `HMM_cov`, `HMM_ie`, `pident`, `qcovhsp`, `scovhsp`, `Signatures`, `InterProAccs`, `Motif_hits`, `TM_pred`, `Topology`, `SignalP`, `Length`.
- Allow‑lists: `allow_pfam` (e.g., `{"PF00201"}`), `allow_ipr` (e.g., `{"IPR001296"}`).
- Topology/SignalP/length options (used mainly in thresholded mode unless you override).

**What it computes:**

- Boolean evidence flags per protein: `ev_HMM`, `ev_DIAM`, `ev_IPR`, `ev_MOTIF`, plus metadata‑driven `ev_LEN`, `ev_TM`, `ev_SP`.
- `evidence_count` and \`\` according to mode:
  - **inclusive (default):** keep if any of `HMM|DIAMOND|IPR|MOTIF` is true.
  - **thresholded:** optionally require `--require-two-evidences` in addition to HMM (stricter).
- `evidence_sources`: semicolon list (e.g., `HMM;IPR;MOTIF`).
- Deduplicates by `Protein_ID`.

**Extend it:** add a new evidence source by:

1. Ensuring your `cli.py` stage merges a column (e.g., `NewTool_score`) onto `master`.
2. In `scoring.py`, set `d["ev_NEWTOOL"] = (your threshold logic)`.
3. Add `"ev_NEWTOOL"` to the evidence list and to the `evidence_sources` aggregator.

### `parsers.py` — file readers

- `parse_hmmer_domtbl(path)->DataFrame`
  - Reads `--domtblout` lines, keeps **best domain per protein**, emits columns: `Protein_ID`, `HMM_model`, `HMM_cov`, `HMM_ie`, `HMM_bits`, `qlen`, `tlen`, etc.
  - `HMM_cov` computed as domain length ÷ HMM model length.
- `parse_diamond_tsv(path)->DataFrame`
  - Expects outfmt 6 with columns: `qseqid sseqid pident length evalue bitscore qcovhsp scovhsp stitle`.
  - Keeps best hit per protein and casts numeric fields.
- `parse_interpro_tsv(path)->DataFrame`
  - Aggregates per protein: `Signatures` (Pfam/SMART/CDD IDs), `SignatureDescs`, `InterProAccs`, `InterProDescs`, `GO_terms`.

**Extend it:**

- Add a `parse_*` for any new tool and return a `DataFrame` keyed by `Protein_ID`.

### `utils.py` — helpers

- `have(tool)->bool`, `which(tool)->str|None`: CLI discovery.
- `run_logged(cmd, log)->int`: runs a command and appends stdout/stderr to `pipeline.log` with a `# CMD:` header.
- `ensure_dir(path)`: safe `mkdir -p`.

### \`install.py

