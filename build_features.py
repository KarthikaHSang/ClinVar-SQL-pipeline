import sqlite3
import pandas as pd

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
query = """
SELECT v.variation_id, v.position, v.type, v.review_status,
       v.number_submitters, v.clinical_significance, g.gene_symbol
FROM variants v
JOIN genes g ON v.gene_id = g.gene_id
WHERE v.position IS NOT NULL
"""
df = pd.read_sql_query(query, conn)
conn.close()

print("Total variants with a position:", len(df))
print(df["clinical_significance"].value_counts().head(10))

# Keep only clearly Pathogenic or clearly Benign variants (drop the ambiguous
# middle ground - Uncertain significance, conflicting classifications, etc.
# - since those aren't a clean label for a classifier to learn from)
df = df[df["clinical_significance"].isin(["Pathogenic", "Benign"])].copy()
df["label"] = (df["clinical_significance"] == "Pathogenic").astype(int)

print("\nAfter filtering to clean Pathogenic/Benign labels:", len(df))
print(df["label"].value_counts())

# Add the in_domain feature using the same domain boundaries as before
def in_domain(row):
    if row["gene_symbol"] not in domains:
        return False
    start, end = domains[row["gene_symbol"]]
    return start <= row["position"] <= end

df["in_domain"] = df.apply(in_domain, axis=1)

print("\nFeature dataframe preview:")
print(df[["variation_id", "gene_symbol", "position", "in_domain", "type",
          "review_status", "number_submitters", "label"]].head(10))

df.to_csv("ml_features.csv", index=False)
print("\nSaved to ml_features.csv, shape:", df.shape)