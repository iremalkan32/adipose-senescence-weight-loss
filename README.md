# Adipose tissue senescence signatures after surgical and dietary weight loss

Analysis code for:

> Alkan İ, Korkmaz Ülüfer D. *Adipose tissue senescence signatures following surgical and dietary weight loss: an integrative analysis of four longitudinal transcriptomic datasets.* Submitted to *Biogerontology*.

The study re-analyses four public Gene Expression Omnibus (GEO) datasets of human abdominal subcutaneous adipose tissue sampled before and after weight loss. It scores senescence with the SenMayo gene set, splits SenMayo into immune-enriched and stromal/adipocyte-enriched components using Human Protein Atlas single-cell data, and adjusts for tissue composition.

| Dataset | Intervention | Follow-up | Platform |
|---|---|---|---|
| GSE65540 | Bariatric surgery | 3 months | RNA-seq |
| GSE59034 | Roux-en-Y gastric bypass | 2 years | Affymetrix HuGene 1.1 ST (GPL11532) |
| GSE199063 | Roux-en-Y gastric bypass | 2 and 5 years | Affymetrix Clariom D (GPL23126) |
| GSE77962 | Low- or very-low-calorie diet | End of weight loss and maintenance | Affymetrix HuGene 1.1 ST (GPL11532) |

GSE59034 and GSE199063 derive from the same Swedish source cohort; `R/03_participant_fingerprint.R` quantifies their participant overlap.

## Repository contents

| Path | Purpose | Manuscript section |
|---|---|---|
| `R/01_main_analysis.R` | Data download from GEO, probe-to-gene mapping, ssGSEA scoring, HPA-based classification of SenMayo genes, paired tests, comparisons with lean controls, random-gene-set specificity, composition adjustment, limma, sensitivity analyses (classification thresholds, leave-one-gene-out) | Methods; Results |
| `R/02_additional_analyses.R` | Linear mixed models and drop-out analysis in GSE199063, TOST equivalence tests, expression-matched random gene sets, SASP-free SenMayo subset, Fridman-specific genes, Holm and Benjamini–Hochberg corrections | Methods (Statistics); Results |
| `R/03_participant_fingerprint.R` | Expression fingerprint to check GSE59034 before/after pairing and participant overlap with GSE199063 | Methods (Cohort independence); Results (Overlap) |
| `R/04_geo_metadata_comparison.R` | Prints and compares GEO series and sample metadata of GSE59034 and GSE199063 | Methods (Datasets) |
| `R/05_publication_figures.R` | Draws Figs. 1–4 from the summary statistics produced by scripts 01–02 | Figures |
| `gene_sets/` | SenMayo (MSigDB SAUL_SEN_MAYO, 125 genes) and Fridman senescence-up (MSigDB FRIDMAN_SENESCENCE_UP, 77 genes) | Methods (Gene sets) |
| `notebooks/` | Google Colab notebooks that run the scripts end to end (outputs removed) | — |
| `results/` | GSE59034 before/after pairing check and GSE59034–GSE199063 participant matches | Results (Overlap) |

Code comments and console messages are in Turkish; function and file names are described above.

## How to reproduce

The analyses were run on Google Colaboratory with an R runtime (R 4.6.1, Bioconductor 3.24).

1. Open `notebooks/01_main_and_additional_analyses.ipynb` in Colab and choose **Runtime › Run all** (about 45–90 minutes). Cell 1 writes the gene-set files, cell 2 runs `R/01_main_analysis.R`, and the last cell runs `R/02_additional_analyses.R`. Outputs are written to `senesans_pilot_cikti/` and `ek_analiz_cikti/`.
2. Run `notebooks/02_normalisation_sensitivity.ipynb` for the TMM and probe-averaging sensitivity analyses.
3. Run `notebooks/03_participant_fingerprint.ipynb` for the participant-overlap analysis.
4. Run `notebooks/04_geo_metadata_comparison.ipynb` to print the GEO metadata of the two Swedish datasets.

All data are downloaded directly from the NCBI GEO FTP server; no data files are stored in this repository. The random seed is fixed (`set.seed(2026)`).

Main R packages: GEOquery, limma, Biobase, nlme, ggplot2, edgeR (sensitivity analysis only) and org.Hs.eg.db (fallback gene annotation).

## Data sources

- Gene Expression Omnibus: GSE65540, GSE59034, GSE199063, GSE77962.
- MSigDB: SAUL_SEN_MAYO, FRIDMAN_SENESCENCE_UP.
- Human Protein Atlas version 23, single-cell type RNA data (`rna_single_cell_type.tsv`).

## Use of generative AI

Generative AI (Claude, Anthropic) was used to help write the analysis code. The authors checked the code and results and take full responsibility for them.

## Licence

MIT (see `LICENSE`).
