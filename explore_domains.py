from Bio import Entrez, SeqIO

Entrez.email = "Karthikasang@gmail.com"

gene_protein_ids = {
    "APC": "NP_000029.2",
    "BRCA1": "NP_009225.1",
    "BRCA2": "NP_000050.3",
    "MLH1": "NP_000240.1",
    "MSH2": "NP_000242.1",
    "MSH6": "NP_000170.1",
    "PMS2": "NP_000526.2",
    "PTEN": "NP_000305.3",
}

for gene, protein_id in gene_protein_ids.items():
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    print(f"\n=== {gene} ({protein_id}), length {len(record.seq)} ===")
    for feature in record.features:
        if feature.type == "Region":
            name = feature.qualifiers.get("region_name", [""])[0]
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            print(f"  [{start}-{end}] {name}")
