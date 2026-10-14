# `calc_gsea_enrichment_on_state_specific_genes`

Calculates Gene Set Enrichment Analysis (GSEA) enrichment scores for state-specific genes based on their pattern activity scores.

```r
calc_gsea_enrichment_on_state_specific_genes(
  gene_patterns_within_state_graph,
  gene_df,
  sig_thresh = 0.1,
  seed_key = NULL
)
```

## Arguments

- **gene_patterns_within_state_graph**  
  *data.frame*  
  A data frame containing gene pattern activity scores.

- **gene_df**  
  *data.frame*  
  A data frame containing gene information with columns `gene_short_name` and `gs_name`.

- **sig_thresh**  
  *numeric*  
  Significance threshold for adjusted p-values. Default is `0.1`.

- **seed_key**  
  *character*  
  Optional string identifying this call, used to derive a deterministic RNG seed (`fgseaMultilevel` is Monte Carlo). Default `NULL`, in which case the seed is derived from the gene ranking and gene sets themselves, so results are reproducible per input without supplying anything.
## Value

A tibble containing the GSEA results filtered by the specified significance threshold.

## Examples

```r
gene_patterns_within_state_graph <- data.frame(
  pattern_activity_score = matrix(runif(100), ncol = 1),
  gene_short_name = sample(letters, 100, replace = TRUE)
)
gene_df <- data.frame(
  gene_short_name = sample(letters, 100, replace = TRUE),
  gs_name = sample(LETTERS[1:5], 100, replace = TRUE)
)
result <- calc_gsea_enrichment_on_state_specific_genes(gene_patterns_within_state_graph, gene_df)
```
