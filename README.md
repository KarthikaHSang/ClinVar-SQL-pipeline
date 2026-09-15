# ClinVar Hereditary Cancer Variant Analysis

A bioinformatics pipeline analyzing ClinVar genetic variant data for 9
hereditary cancer genes (BRCA1, BRCA2, TP53, MLH1, MSH2, MSH6, PMS2, APC,
PTEN), combining data engineering, SQL database design, sequence validation
against NCBI, statistical hypothesis testing, and machine learning.

## What this project does

1. **ETL pipeline** (`Sql-pipeline.py`) - downloads ClinVar's variant summary
   file, filters to the 9-gene panel, deduplicates by genome assembly and
   variant ID, and loads the result into a normalized SQLite database.
2. **Normalized SQL schema** - 4 tables (`genes`, `variants`, `phenotypes`,
   `variant_phenotypes`) with enforced primary and foreign key constraints,
   queried through 5 externalized `.sql` files in `sql/`.
3. **Protein-change parsing** - extracts amino acid change and position from
   HGVS notation (e.g. `p.Arg248Gln`) using regex.
4. **Sequence validation** - fetches real protein sequences from NCBI via
   Biopython and cross-checks parsed positions against them (99.6% match
   rate for TP53; all discrepancies explained as stop-loss variants).
5. **Domain-pathogenicity analysis** (`all_genes_domain_analysis.py`) - tests
   whether pathogenic variants cluster within known functional domains,
   using a chi-squared test across all 9 genes.
6. **Machine learning** (`build_features.py`, `train_classifier.py`) - trains
   a Random Forest classifier to predict pathogenicity from variant features.

Full methods and results are in [`domain_pathogenicity_report.md`](domain_pathogenicity_report.md).

## Key findings

- TP53 and PTEN show the strongest statistical clustering of pathogenic
  variants within their core functional domains (p < 1e-8).
- BRCA1, MLH1, and PMS2 show no significant clustering for the single domain
  tested, most likely because these genes have multiple clinically important
  domains not fully captured by a single-domain analysis.
- A Random Forest classifier achieves 81% accuracy predicting pathogenicity,
  though feature importance reveals that raw position and submission count
  outweigh domain membership - an interesting tension with the statistical
  findings, discussed in the report.

## Tech stack

Python, pandas, SQLite, Biopython (Entrez/SeqIO), scipy (statistical testing),
scikit-learn (Random Forest classification)

## Setup


Note: the first run of `Sql-pipeline.py` downloads ClinVar's full variant
summary file (~300MB) if no local cache is present.