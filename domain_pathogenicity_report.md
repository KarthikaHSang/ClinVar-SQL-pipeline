# ClinVar Variant Domain Analysis: Pathogenicity Clustering in Functional Protein Domains

## Overview

This analysis examines whether pathogenic variants in 9 hereditary cancer genes
(BRCA1, BRCA2, TP53, MLH1, MSH2, MSH6, PMS2, APC, PTEN) cluster within known
functional protein domains, using ClinVar variant data cross-validated against
real NCBI protein sequences.

## Data Pipeline

1. Downloaded and filtered ClinVar's `variant_summary.txt` to the 9-gene panel
2. Deduplicated by genome assembly (GRCh38 only) and by VariationID
3. Normalized into a 4-table SQLite schema (genes, variants, phenotypes,
   variant_phenotypes) with enforced primary/foreign key constraints
4. Parsed HGVS protein notation (e.g. `p.Arg248Gln`) from variant names using
   regex, extracting reference amino acid, position, and alternate amino acid
   for ~65% of variants (the remainder are frameshift or splice-site variants,
   which don't have a single-position protein notation)
5. Fetched real protein sequences for all 9 genes from NCBI via Biopython's
   Entrez interface, and validated parsed positions against them

## Validation

For TP53, 2,143 of 2,152 parsed variants (99.6%) matched the actual amino acid
at their claimed position in the real NCBI sequence. All 9 mismatches were
stop-codon (Ter) variants at position 394 - one past the end of the 393-residue
protein - representing stop-loss/read-through mutations rather than parsing
errors.

## Domain Selection

For each gene, a single functional domain was chosen from NCBI's annotated
protein features, based on known biological significance:

| Gene  | Domain                          | Range      |
|-------|----------------------------------|------------|
| TP53  | DNA-binding domain (P53)         | 109-288    |
| PTEN  | Phosphatase domain (PTP_PTEN)     | 24-181     |
| APC   | Beta-catenin regulation region    | 960-1337   |
| BRCA1 | Tandem BRCT domains               | 1650-1855  |
| BRCA2 | DNA-binding domain (BRCA2DBD)      | 2804-3054  |
| MLH1  | MutL N-terminal (ATPase) domain   | 6-315      |
| MSH2  | MutS domain                       | 175-853    |
| MSH6  | MutS domain                       | 407-1360   |
| PMS2  | MutL domain                        | 14-853     |

**Limitation**: this is a simplification. Several genes (notably BRCA1, which
also has a well-known pathogenic RING domain at positions 7-99) have more than
one clinically relevant domain, and only one was selected per gene for this
analysis. A more complete analysis would incorporate all annotated domains,
ideally sourced from UniProt's curated domain boundaries rather than NCBI's
mixed region/site annotations.

## Results

Pathogenic rate inside vs. outside the selected domain, with chi-squared test
of independence:

| Gene  | Domain      | n (in) | n (out) | Rate (in) | Rate (out) | Chi2   | p-value    |
|-------|-------------|--------|---------|-----------|------------|--------|------------|
| TP53  | 109-288     | 1162   | 990     | 22.98%    | 5.15%      | 133.48 | 7.1e-31    |
| PTEN  | 24-181      | 848    | 968     | 21.93%    | 11.88%     | 32.32  | 1.3e-08    |
| APC   | 960-1337    | 1624   | 10601   | 8.44%     | 4.56%      | 43.23  | 4.9e-11    |
| MSH2  | 175-853     | 3486   | 1411    | 9.75%     | 6.02%      | 17.16  | 3.4e-05    |
| BRCA2 | 2804-3054   | 1234   | 13596   | 10.29%    | 7.78%      | 9.36   | 2.2e-03    |
| MSH6  | 407-1360    | 5309   | 2252    | 7.38%     | 5.77%      | 6.14   | 1.3e-02    |
| PMS2  | 14-853      | 4207   | 141     | 5.16%     | 9.22%      | 3.72   | 5.4e-02    |
| MLH1  | 6-315       | 1513   | 2004    | 11.17%    | 10.13%     | 0.88   | 3.5e-01    |
| BRCA1 | 1650-1855   | 1845   | 7404    | 8.02%     | 8.94%      | 1.45   | 2.3e-01    |

## Discussion

- TP53 and PTEN show the strongest domain clustering, both roughly doubling
  or more the pathogenic rate inside vs. outside the selected domain. This
  aligns with known biology: TP53's DNA-binding domain and PTEN's catalytic
  phosphatase domain are essential to each protein's tumor-suppressor function.
- APC, MSH2, BRCA2, and MSH6 show statistically significant but more modest
  clustering - real, but with smaller effect sizes.
- MLH1, PMS2, and BRCA1 show no significant clustering for the domain
  selected here. This does not mean these genes lack pathogenic hotspots -
  it more likely reflects the limitation of selecting only one domain per gene.
  BRCA1 in particular has other well-documented pathogenic regions (e.g. the
  N-terminal RING domain) not included in this analysis.

## Next Steps

- Incorporate multiple domains per gene, ideally from UniProt
- Extend the amino-acid-level analysis to a predictive model (see companion
  scikit-learn analysis)