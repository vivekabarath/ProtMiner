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

## 📊 Output Files — Detailed

ProtMiner produces several key outputs in the chosen `--outdir`:

### 1. `all_proteins_annotation.tsv`

- Comprehensive table of all proteins scanned with evidence and metadata.
- Columns include: see explanation below.

**Column explanations:**

- `Protein_ID`: Identifier of the protein from FASTA.
- `Length`: Sequence length in amino acids.
- `HMM_model`: Pfam/custom HMM model ID matched.
- `HMM_cov`: Fractional coverage of the HMM domain.
- `HMM_ie`: Independent E-value from HMMER.
- `HMM_bits`: Bitscore from HMMER.
- `pident`: Percent identity to best DIAMOND hit.
- `qcovhsp`/`scovhsp`: Query and subject coverage of HSP in DIAMOND.
- `bitscore`: Alignment bitscore in DIAMOND.
- `stitle`: Subject title (annotation) from DIAMOND DB.
- `Signatures`: Raw domain signatures (Pfam/SMART/CDD).
- `InterProAccs`: InterPro accession numbers.
- `GO_terms`: Gene Ontology annotations from InterPro.
- `Motif_hits`: Regex motif matches with positions.
- `TM_pred`: Number of transmembrane segments predicted.
- `Topology`: Predicted topology string (optional).
- `SignalP`: Signal peptide prediction.
- `ev_*`: Boolean evidence flags (1 if present, else 0).
- `evidence_count`: Total number of evidence sources supporting the protein.
- `evidence_sources`: Semicolon list of which sources triggered inclusion.
- `keep_candidate`: True if selected as final positive.

**Example snippet:**

```tsv
Protein_ID	Length	HMM_model	HMM_cov	pident	Signatures	InterProAccs	Motif_hits	TM_pred	SignalP	ev_HMM	ev_DIAM	ev_IPR	ev_MOTIF	evidence_count	evidence_sources	keep_candidate
protA	482	PF00201	0.95	62.1	PF00201;PF03106	IPR001296;IPR035595	45-72:FWx2QL...	0	None	1	1	1	1	4	HMM;DIAM;IPR;MOTIF	True
protB	615	-	-	25.7	-	-	-	1	SP	0	1	0	0	1	DIAM	True
```

### 2. `candidates.tsv`

- Subset of `all_proteins_annotation.tsv` containing only the final positives (`keep_candidate=True`).
- Deduplicated by `Protein_ID`.

### 3. `candidates.faa`

- FASTA file of the sequences from `candidates.tsv`.
- Headers usually include `Protein_ID` and evidence summary.

**Example snippet:**

```fasta
>protA | HMM;DIAM;IPR;MOTIF
MASSKTFVLLVLLAVAA...WQLEQLLL...
>protB | DIAM
MTDSQVKDLLELAAKRV...
```

### 4. `pipeline.log`

- Full log of the run, including environment checks, executed commands, and warnings.
- Useful for debugging.

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

## 📊 Output Files

- `all_proteins_annotation.tsv` — all evidence merged
- `candidates.tsv` — final positives (deduplicated)
- `candidates.faa` — FASTA of candidates
- `pipeline.log` — commands and logs

---

## 📊 Output Files — Detailed

ProtMiner produces several key outputs in the chosen `--outdir`:

### 1. `all_proteins_annotation.tsv`

- Comprehensive table of all proteins scanned with evidence and metadata.
- Columns include: see explanation below.

**Column explanations:**

- `Protein_ID`: Identifier of the protein from FASTA.
- `Length`: Sequence length in amino acids.
- `HMM_model`: Pfam/custom HMM model ID matched.
- `HMM_cov`: Fractional coverage of the HMM domain.
- `HMM_ie`: Independent E-value from HMMER.
- `HMM_bits`: Bitscore from HMMER.
- `pident`: Percent identity to best DIAMOND hit.
- `qcovhsp`/`scovhsp`: Query and subject coverage of HSP in DIAMOND.
- `bitscore`: Alignment bitscore in DIAMOND.
- `stitle`: Subject title (annotation) from DIAMOND DB.
- `Signatures`: Raw domain signatures (Pfam/SMART/CDD).
- `InterProAccs`: InterPro accession numbers.
- `GO_terms`: Gene Ontology annotations from InterPro.
- `Motif_hits`: Regex motif matches with positions.
- `TM_pred`: Number of transmembrane segments predicted.
- `Topology`: Predicted topology string (optional).
- `SignalP`: Signal peptide prediction.
- `ev_*`: Boolean evidence flags (1 if present, else 0).
- `evidence_count`: Total number of evidence sources supporting the protein.
- `evidence_sources`: Semicolon list of which sources triggered inclusion.
- `keep_candidate`: True if selected as final positive.

**Example snippet:**

```tsv
Protein_ID	Length	HMM_model	HMM_cov	pident	Signatures	InterProAccs	Motif_hits	TM_pred	SignalP	ev_HMM	ev_DIAM	ev_IPR	ev_MOTIF	evidence_count	evidence_sources	keep_candidate
protA	482	PF00201	0.95	62.1	PF00201;PF03106	IPR001296;IPR035595	45-72:FWx2QL...	0	None	1	1	1	1	4	HMM;DIAM;IPR;MOTIF	True
protB	615	-	-	25.7	-	-	-	1	SP	0	1	0	0	1	DIAM	True
```

### 2. `candidates.tsv`

- Subset of `all_proteins_annotation.tsv` containing only the final positives (`keep_candidate=True`).
- Deduplicated by `Protein_ID`.

### 3. `candidates.faa`

- FASTA file of the sequences from `candidates.tsv`.
- Headers usually include `Protein_ID` and evidence summary.

**Example snippet:**

```fasta
>protA | HMM;DIAM;IPR;MOTIF
MASSKTFVLLVLLAVAA...WQLEQLLL...
>protB | DIAM
MTDSQVKDLLELAAKRV...
```

### 4. `pipeline.log`

- Full log of the run, including environment checks, executed commands, and warnings.
- Useful for debugging.

---

## 🧱 Code Structure & Extensibility

### Repository layout

```
ProtMiner/
├─ README.md
├─ environment.yml
├─ src/
   └─ protminer/
      ├─ __init__.py
      ├─ cli.py            # command-line interface & pipeline orchestrator
      ├─ scoring.py        # evidence logic (inclusive/thresholded) & final filtering
      ├─ parsers.py        # readers for HMMER/DIAMOND/InterPro outputs
      ├─ utils.py          # helpers: which(), have(), run_logged(), ensure_dir()
      └─ install.py        # dependency checks & (optional) auto-install hints

```

### `cli.py` — pipeline orchestrator

- Parses CLI flags (input FASTA, HMMs, motif, DIAMOND reference, IPR/Pfam allow-lists, topology/signal/length options, scoring mode).
- Runs environment checks and optional tool discovery/installation hints.
- Executes stages: deduplication, HMMER, DIAMOND, InterProScan, TM/SignalP, motif scan, merge evidence, scoring/filtering, output.

### `scoring.py` — evidence logic

- Core function: `score_and_filter(df, ..., scoring_mode="inclusive")`.
- Sets evidence flags, counts, and selects candidates according to inclusive/thresholded rules.
- Extendable: add new evidence flags by merging new columns and adjusting logic.

### `parsers.py` — file readers

- `parse_hmmer_domtbl`, `parse_diamond_tsv`, `parse_interpro_tsv`.
- Each returns a DataFrame keyed by `Protein_ID`.

### `utils.py` — helpers

- `have()`, `which()`, `run_logged()`, `ensure_dir()`.

### `install.py` — dependency checks

- Ensures Python modules and external tools exist; prints installation hints.

---

## 🔮 Roadmap

- Add visualization utilities (plots, summaries).
- Extend to additional domain databases.
- Provide Docker image for consistent environments.
- Improve motif scanning with PROSITE/ELM conversion.
- Add benchmarking datasets and tests.

