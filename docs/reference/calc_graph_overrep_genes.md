# `calc_graph_overrep_genes`

Calculates overrepresented pathways for genes in a given dataset, grouped by cell state. Processes the input data, simplifies gene interpretation, and computes overrepresented pathways for each cell state.

```r
calc_graph_overrep_genes(
  ref_degs,
  pathway_type = c("GO:BP", "GO:MF"),
  species = "Danio rerio"
)
```

## Arguments

- **ref_degs**  
  *data.frame*  
  A data frame containing reference differentially expressed genes (DEGs). Should include columns such as `cell_state`, `gene_short_name`, `gene_id`, and `interpretation`. If the data frame contains a `gene_class_scores` column, it will be unnested.

- **pathway_type**  
  *character vector*  
  The type of pathways to analyze: `"GO:BP"` (Biological Process) or `"GO:MF"` (Molecular Function). Defaults to `c("GO:BP", "GO:MF")`.

- **species**  
  *character*  
  The species for pathway analysis. Defaults to `"Danio rerio"`.

## Value

A nested data frame where each row corresponds to a cell state, and the `overrep_pathways` column contains the overrepresented pathways for that cell state.

## Details

Simplifies the `interpretation` column into categories (`"Up"`, `"Down"`, `"Maintained"`, or retains the original value if none match), groups by `cell_state`, and calculates overrepresented pathways for each group using `calc_overrep_genes()`.

## Examples

```r
ref_degs <- data.frame(
  cell_state = c("State1", "State2"),
  gene_short_name = c("GeneA", "GeneB"),
  gene_id = c("ID1", "ID2"),
  interpretation = c("Upregulated", "Downregulated")
)
result <- calc_graph_overrep_genes(ref_degs, pathway_type = "GO:BP", species = "Danio rerio")
```
