print("SCRIPT STARTED")

import pandas as pd
import os

gene_panel = ["BRCA1", "BRCA2", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "APC", "PTEN"]

if os.path.exists("clinvar_filtered.csv"):
    print("Loading filtered subset...")
    df_filtered = pd.read_csv("clinvar_filtered.csv")
else:
    if os.path.exists("variant_summary.txt.gz"):
        print("Loading local full file...")
        clinvar = pd.read_csv("variant_summary.txt.gz", sep="\t", compression="gzip", low_memory=False)
    else:
        print("Downloading full file (one-time, ~300MB, may take a few minutes)...")
        url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        clinvar = pd.read_csv(url, sep="\t", compression="gzip", low_memory=False)
        clinvar.to_csv("variant_summary_full_cache.csv", index=False)

    print(clinvar.shape)
    df_filtered = clinvar[clinvar["GeneSymbol"].isin(gene_panel)].copy()
    df_filtered.to_csv("clinvar_filtered.csv", index=False)

print(df_filtered.shape)
print(df_filtered["ClinicalSignificance"].value_counts())
print(df_filtered["PhenotypeList"].nunique())
print(df_filtered.isnull().sum())

print(df_filtered["Assembly"].value_counts())

df_filtered = df_filtered[df_filtered["Assembly"] == "GRCh38"].copy()
print(df_filtered.shape)

df_filtered.to_csv("clinvar_filtered_grch38.csv", index=False)

genes_df = df_filtered[["GeneID", "GeneSymbol", "HGNC_ID"]].drop_duplicates(subset="GeneID").reset_index(drop=True)
genes_df = genes_df.rename(columns={"GeneID": "gene_id", "GeneSymbol": "gene_symbol", "HGNC_ID": "hgnc_id"})
print(genes_df.shape)
print(genes_df)

variants_df = df_filtered[[
    "VariationID", "#AlleleID", "GeneID", "Name", "Type",
    "Chromosome", "Start", "Stop", "ReferenceAllele", "AlternateAllele",
    "ClinicalSignificance", "ReviewStatus", "NumberSubmitters", "LastEvaluated"
]].drop_duplicates(subset="VariationID").reset_index(drop=True)

variants_df = variants_df.rename(columns={
    "VariationID": "variation_id",
    "#AlleleID": "allele_id",
    "GeneID": "gene_id",
    "Name": "name",
    "Type": "type",
    "Chromosome": "chromosome",
    "Start": "start",
    "Stop": "stop",
    "ReferenceAllele": "reference_allele",
    "AlternateAllele": "alternate_allele",
    "ClinicalSignificance": "clinical_significance",
    "ReviewStatus": "review_status",
    "NumberSubmitters": "number_submitters",
    "LastEvaluated": "last_evaluated"
})

print(variants_df.shape)
print(variants_df.head())

phenotype_rows = df_filtered[["VariationID", "PhenotypeList"]].drop_duplicates(subset="VariationID").copy()
phenotype_rows["PhenotypeList"] = phenotype_rows["PhenotypeList"].str.split(r"[|;]")
phenotype_rows = phenotype_rows.explode("PhenotypeList")
phenotype_rows["PhenotypeList"] = phenotype_rows["PhenotypeList"].str.strip()
phenotype_rows = phenotype_rows[~phenotype_rows["PhenotypeList"].isin(["not provided", "not specified"])]
phenotype_rows = phenotype_rows.drop_duplicates(subset=["VariationID", "PhenotypeList"])

print(phenotype_rows.shape)
print(phenotype_rows.head(10))

phenotypes_df = phenotype_rows["PhenotypeList"].drop_duplicates().reset_index(drop=True).to_frame()
phenotypes_df["phenotype_id"] = phenotypes_df.index + 1
phenotypes_df = phenotypes_df.rename(columns={"PhenotypeList": "phenotype_name"})
phenotypes_df = phenotypes_df[["phenotype_id", "phenotype_name"]]

print(phenotypes_df.shape)
print(phenotypes_df.head())

variant_phenotypes_df = phenotype_rows.merge(
    phenotypes_df, left_on="PhenotypeList", right_on="phenotype_name"
)[["VariationID", "phenotype_id"]]

variant_phenotypes_df = variant_phenotypes_df.rename(columns={"VariationID": "variation_id"})

print(variant_phenotypes_df.shape)
print(variant_phenotypes_df.head())

import re

sample_names = variants_df["name"].head(10)
pattern = r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|=|Ter|\*)"

for n in sample_names:
    match = re.search(pattern, n)
    if match:
        print(n, "->", match.groups())
    else:
        print(n, "-> NO MATCH")

def parse_protein_change(name):
    match = re.search(pattern, name)
    if match:
        return pd.Series(match.groups())
    else:
        return pd.Series([None, None, None])

variants_df[["ref_aa", "position", "alt_aa"]] = variants_df["name"].apply(parse_protein_change)

print(variants_df[["variation_id", "name", "ref_aa", "position", "alt_aa"]].head(10))
print(variants_df["ref_aa"].isnull().sum(), "variants had no parseable protein change")

variants_df["position"] = pd.to_numeric(variants_df["position"], errors="coerce")
print(variants_df["position"].dtype)

import sqlite3
conn = sqlite3.connect("clinvar.db")
cursor = conn.cursor()

cursor.execute("DROP TABLE IF EXISTS genes")
cursor.execute("""
    CREATE TABLE genes (
        gene_id INTEGER PRIMARY KEY,
        gene_symbol TEXT NOT NULL,
        hgnc_id TEXT
    )
""")
genes_df.to_sql("genes", conn, if_exists="append", index=False)
print("genes table created with PRIMARY KEY constraint")

cursor.execute("DROP TABLE IF EXISTS variants")
cursor.execute("""
    CREATE TABLE variants (
        variation_id INTEGER PRIMARY KEY,
        allele_id INTEGER,
        gene_id INTEGER,
        name TEXT,
        type TEXT,
        chromosome TEXT,
        start INTEGER,
        stop INTEGER,
        reference_allele TEXT,
        alternate_allele TEXT,
        clinical_significance TEXT,
        review_status TEXT,
        number_submitters INTEGER,
        last_evaluated TEXT,
        ref_aa TEXT,
        position REAL,
        alt_aa TEXT,
        FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
    )
""")
variants_df.to_sql("variants", conn, if_exists="append", index=False)
print("variants table created with PRIMARY KEY and FOREIGN KEY constraint")

cursor.execute("DROP TABLE IF EXISTS phenotypes")
cursor.execute("""
    CREATE TABLE phenotypes (
        phenotype_id INTEGER PRIMARY KEY,
        phenotype_name TEXT NOT NULL
    )
""")
phenotypes_df.to_sql("phenotypes", conn, if_exists="append", index=False)
print("phenotypes table created with PRIMARY KEY constraint")

cursor.execute("DROP TABLE IF EXISTS variant_phenotypes")
cursor.execute("""
    CREATE TABLE variant_phenotypes (
        variation_id INTEGER,
        phenotype_id INTEGER,
        PRIMARY KEY (variation_id, phenotype_id),
        FOREIGN KEY (variation_id) REFERENCES variants(variation_id),
        FOREIGN KEY (phenotype_id) REFERENCES phenotypes(phenotype_id)
    )
""")
variant_phenotypes_df.to_sql("variant_phenotypes", conn, if_exists="append", index=False)
print("variant_phenotypes table created with composite PRIMARY KEY and 2 FOREIGN KEYs")

conn.close()
print("Loaded 4 normalized tables into clinvar.db")

conn = sqlite3.connect("clinvar.db")
conn.execute("PRAGMA foreign_keys = ON")
for table in ["genes", "variants", "phenotypes", "variant_phenotypes"]:
    check = pd.read_sql(f"SELECT COUNT(*) FROM {table}", conn)
    print(table, check.iloc[0, 0])
conn.close()

conn = sqlite3.connect("clinvar.db")
with open("sql/variants_with_genes.sql") as f:
    query1 = f.read()
result1 = pd.read_sql_query(query1, conn)
pd.set_option('display.width', None)
print(result1)
conn.close()

conn = sqlite3.connect("clinvar.db")
with open("sql/pathogenic_by_gene.sql") as f:
    query2 = f.read()
result2 = pd.read_sql_query(query2, conn)
print(result2)
conn.close()

conn = sqlite3.connect("clinvar.db")
with open("sql/phenotypes_by_pathogenic.sql") as f:
    query3 = f.read()
result3 = pd.read_sql_query(query3, conn)
print(result3)
conn.close()

conn = sqlite3.connect("clinvar.db")
gene_of_interest = "BRCA1"
with open("sql/gene_phenotype_breakdown.sql") as f:
    query4 = f.read()
result4 = pd.read_sql_query(query4, conn, params=(gene_of_interest,))
print(result4)
conn.close()

conn = sqlite3.connect("clinvar.db")
with open("sql/variants_most_phenotypes.sql") as f:
    query5 = f.read()
result5 = pd.read_sql_query(query5, conn)
print(result5)
conn.close()



