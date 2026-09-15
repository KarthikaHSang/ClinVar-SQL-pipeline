SELECT v.variation_id, v.name, g.gene_symbol, v.clinical_significance
FROM variants v
JOIN genes g ON v.gene_id = g.gene_id
LIMIT 10;
