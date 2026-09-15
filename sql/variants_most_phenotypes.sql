SELECT v.variation_id, v.name, g.gene_symbol, COUNT(*) AS phenotype_count
FROM variants v
JOIN genes g ON v.gene_id = g.gene_id
JOIN variant_phenotypes vp ON v.variation_id = vp.variation_id
GROUP BY v.variation_id, v.name, g.gene_symbol
ORDER BY phenotype_count DESC
LIMIT 15;
