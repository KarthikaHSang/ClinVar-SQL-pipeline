SELECT p.phenotype_name, COUNT(*) AS variant_count
FROM variants v
JOIN variant_phenotypes vp ON v.variation_id = vp.variation_id
JOIN phenotypes p ON vp.phenotype_id = p.phenotype_id
WHERE v.clinical_significance = 'Pathogenic'
GROUP BY p.phenotype_name
ORDER BY variant_count DESC
LIMIT 15;
