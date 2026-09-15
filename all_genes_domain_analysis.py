import sqlite3
import pandas as pd
from scipy.stats import chi2_contingency

domains = {
    "TP53":  (109, 288),
    "APC":   (960, 1337),
    "BRCA1": (1650, 1855),
    "BRCA2": (2804, 3054),
    "MLH1":  (6, 315),
    "MSH2":  (175, 853),
    "MSH6":  (407, 1360),
    "PMS2":  (14, 853),
    "PTEN":  (24, 181),
}

conn = sqlite3.connect("clinvar.db")

results = []

for gene, (start, end) in domains.items():
    query = """
    SELECT v.position, v.clinical_significance
    FROM variants v
    JOIN genes g ON v.gene_id = g.gene_id
    WHERE g.gene_symbol = ? AND v.position IS NOT NULL
    """
    df = pd.read_sql_query(query, conn, params=(gene,))

    df["in_domain"] = df["position"].apply(lambda p: start <= p <= end)
    df["is_pathogenic"] = df["clinical_significance"].str.contains("Pathogenic", na=False)

    contingency = pd.crosstab(df["in_domain"], df["is_pathogenic"])

    if contingency.shape == (2, 2):
        chi2, p_value, dof, expected = chi2_contingency(contingency)
    else:
        chi2, p_value = None, None

    in_domain_df = df[df["in_domain"]]
    out_domain_df = df[~df["in_domain"]]

    in_rate = in_domain_df["is_pathogenic"].mean() if len(in_domain_df) else None
    out_rate = out_domain_df["is_pathogenic"].mean() if len(out_domain_df) else None

    results.append({
        "gene": gene,
        "domain_range": f"{start}-{end}",
        "n_in_domain": len(in_domain_df),
        "n_out_domain": len(out_domain_df),
        "pathogenic_rate_in": round(in_rate, 4) if in_rate is not None else None,
        "pathogenic_rate_out": round(out_rate, 4) if out_rate is not None else None,
        "chi2": round(chi2, 2) if chi2 is not None else None,
        "p_value": p_value,
    })

conn.close()

results_df = pd.DataFrame(results)
pd.set_option('display.width', None)
print(results_df)

results_df.to_csv("domain_pathogenicity_all_genes.csv", index=False)
print("\nSaved to domain_pathogenicity_all_genes.csv")

code domain_pathogenicity_report.md