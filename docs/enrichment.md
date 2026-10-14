# Functional enrichment in Platt

Once you have per-gene pattern calls — from [Reference DEGs](https://cole-trapnell-lab.github.io/platt/wt_degs/) or [Perturbation DEGs](https://cole-trapnell-lab.github.io/platt/perturbation_degs/) — the next question is often at the level of programs rather than individual genes: is a cell state's differentially-expressed gene set enriched for a particular biological process, like proliferation, apoptosis, or a stress response? Platt wraps [`fgsea`](https://bioconductor.org/packages/fgsea/) for this, at two granularities: rank-based GSEA within a single cell state, or over-representation analysis (ORA) across the whole graph at once.

## Rank-based GSEA within a single cell state

Any table with a `pattern_activity_score` column — such as `condensate_genes` from the [Reference DEGs](https://cole-trapnell-lab.github.io/platt/wt_degs/#comparing-genes-in-a-single-cell-state) page — can be handed to `calc_gsea_enrichment_on_state_specific_genes()`, which uses that score to rank the state's genes and runs `fgsea::fgsea()` against a gene set table you supply:

* `gene_patterns_within_state_graph` - a table with `gene_short_name` and `pattern_activity_score` columns, e.g. `condensate_genes`
* `gene_df` - a gene set table with `gene_short_name` and `gs_name` columns; [`msigdbr`](https://cran.r-project.org/package=msigdbr) is a convenient source for zebrafish GO/pathway gene sets, but any curated list works the same way
* `sig_thresh` - `padj` cutoff for the returned results, default `0.1`

```r
library(msigdbr)

gene_sets = msigdbr(species = "Danio rerio", subcategory = "GO:BP") %>%
  dplyr::distinct(gs_name, gene_short_name = gene_symbol)

condensate_gsea = calc_gsea_enrichment_on_state_specific_genes(
  condensate_genes,
  gene_sets,
  sig_thresh = 0.1
)
```

The result is a standard `fgsea` output tibble — `pathway`, `pval`, `padj`, `ES`, `NES`, `size`, and the leading-edge gene list — already filtered down to `padj < sig_thresh`. A pathway with a large positive `NES` is enriched among the genes with the highest `pattern_activity_score` at this state (e.g. its most strongly activated/upregulated genes); a large negative `NES` is enriched among its most strongly deactivated/downregulated genes.

_`gene_df` doesn't have to come from `msigdbr` — a handful of curated proliferation/apoptosis/stress/senescence gene sets, tested the same way, is what typically informs the `fitness_label` calls on the [Phenotyping](https://cole-trapnell-lab.github.io/platt/phenotyping/) page._

## Over-representation across the whole graph

For a graph-wide pass, `calc_graph_overrep_genes()` takes the same kind of nested DEG table `compare_genes_over_graph()` returns (see [Reference DEGs](https://cole-trapnell-lab.github.io/platt/wt_degs/#running-degs-over-a-graph)) and computes ORA independently at every cell state:

* `ref_degs` - a DEG table with `cell_state`, `gene_short_name`, `gene_id`, and `interpretation` columns — `pf_graph_degs` from the Reference DEGs page works directly (its nested `gene_class_scores` column is unnested automatically)
* `pathway_type` - `"GO:BP"` or `"GO:MF"`, forwarded to `msigdbr` to build the gene set table
* `species` - forwarded to `msigdbr`, default `"Danio rerio"`

```r
pf_graph_overrep = calc_graph_overrep_genes(pf_graph_degs, pathway_type = "GO:BP")
```

Under the hood, this collapses each gene's finer `interpretation` call (the same category names from [Gene expression patterns](https://cole-trapnell-lab.github.io/platt/wt_degs/#gene-expression-patterns)) down to `"Up"`, `"Down"`, or `"Maintained"`, then runs ORA (`fgsea::fora()`) separately for the up- and down-regulated gene sets at each state — `"Maintained"` genes aren't tested, since there's no direction of change to look for. The result is a nested tibble, one row per cell state:

```r
pf_graph_overrep %>%
  tidyr::unnest(overrep_pathways) %>%
  dplyr::select(cell_state, interpretation_simple, pathway, padj)
```

| cell_state                  | interpretation_simple | pathway                          | padj   |
|------------------------------|------------------------|-----------------------------------|--------|
| pectoral fin condensate      | Up                     | GOBP_SKELETAL_SYSTEM_DEVELOPMENT | 0.004  |
| pectoral fin condensate      | Up                     | GOBP_COLLAGEN_FIBRIL_ORGANIZATION | 0.017  |
| pectoral fin bud progenitor  | Down                   | GOBP_CELL_DIFFERENTIATION        | 0.031  |

`padj` is already filtered to the `q_val_thresh` default (`0.05`) inside `calc_overrep_genes()`, the per-state worker `calc_graph_overrep_genes()` calls internally. If you already have a plain gene list for one state rather than a full graph-wide DEG table, `calc_overrep_genes()` and `calc_pathway_enrichment_on_state_specific_genes()` (a thinner wrapper around `clusterProfiler::enricher()`) are the lower-level building blocks — most users will only need the graph-wide wrapper above.

## Visualizing enrichment results

There's no dedicated `plot_*` wrapper for enrichment results yet, but the unnested tibble drops straight into `ggplot2` — a dot plot with point size mapped to significance and color mapped to direction is a common way to scan results across many cell states at once:

```r
pf_graph_overrep %>%
  tidyr::unnest(overrep_pathways) %>%
  dplyr::group_by(cell_state) %>%
  dplyr::slice_min(padj, n = 5) %>%
  ggplot(aes(x = cell_state, y = pathway, size = -log10(padj), color = interpretation_simple)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

_For more on turning a handful of enrichment results into the fitness axes (proliferation, apoptosis, stress, senescence) used on a phenotype call, see the [Phenotyping](https://cole-trapnell-lab.github.io/platt/phenotyping/) page — Platt doesn't yet package that specific mapping for you, so going from enrichment results to a single `fitness_label` is still a judgment call you make yourself._
