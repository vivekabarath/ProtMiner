# README.md

# ProtScout 🧬🧭

> **Inclusive protein family discovery** using HMMER, DIAMOND, InterPro, and motif evidence — with optional topology/signal peptide metadata.

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Conda](https://img.shields.io/badge/conda-bioconda%20%7C%20conda--forge-blue)](#-installation)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](#-installation)

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
- [Roadmap](#-roadmap)

---

## Overview
ProtScout identifies proteins belonging to a given family by combining multiple evidence sources:
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
git clone https://github.com/<yourusername>/ProtScout.git
cd ProtScout

# Create environment
conda env create -f environment.yml
conda activate protscout

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

## 🔮 Roadmap
- Preset configs (`presets/ugt.yaml`, `presets/cyp.yaml`)
- PROSITE → regex conversion (`--motif-prosite`)
- Containerized distribution (Docker/Singularity)

