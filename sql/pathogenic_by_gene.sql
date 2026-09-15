SELECT g.gene_symbol, COUNT(*) AS pathogenic_count
FROM variants v
JOIN genes g ON v.gene_id = g.gene_id
WHERE v.clinical_significance = 'Pathogenic'
GROUP BY g.gene_symbol
ORDER BY pathogenic_count DESC;
