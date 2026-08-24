# `summarize_impact_in_lineage_context`

Builds a lineage-wide impact table — one row per cell type — from your abundance/fitness/identity calls, DEGs, and pathway enrichment, ready for `impact_to_phenos()`. Optionally has an LLM write a per-cell-type mechanistic summary; passing `llm_fun = NULL` (the default) skips that step and returns the same tabular impact calls with no LLM involved.

```r
summarize_impact_in_lineage_context(
  perturbation_description = NULL,
  target_gene_expression,
  dact_results,
  degs,
  ref_expression,
  ontologies,
  combined_psg,
  ai_notes_path = NULL,
  sig_p_val_thresh = 0.05,
  genes_of_interest = NULL,
  empirical_p_thresh = 1,
  llm_fun = NULL,
  max_lineage_depth = Inf,
  cell_types = NULL,
  primary_impact_summary = "",
  excluded_cell_types = NULL,
  abundance_phenotypes = NULL,
  fitness_phenotypes = NULL,
  identity_phenotypes = NULL,
  dact_expectations = NULL,
  effect_type = NULL,
  pre_cited_gene_claims = NULL,
  verbose = FALSE,
  sig_pathways_by_cell_type = NULL,
  sig_pathways_path = NULL,
  drop_llm_columns = FALSE,
  ...
)
```

## Arguments

- **perturbation_description**  
  Optional free-text description of the perturbation, used only to build context text for `llm_fun`.

- **target_gene_expression**  
  Data frame of target-gene expression calls by cell group; used to seed which cell types are analyzed when `cell_types` is not supplied.

- **dact_results**  
  Differential abundance result table (one row per cell type, e.g. from `dacts_when_abundant()`).

- **degs**  
  Differential expression table for the perturbation.

- **ref_expression**  
  Reference expression table.

- **ontologies**  
  Named list of ontology gene-set tables (`gs_name`/`gene_short_name` columns) used for per-cell-type pathway enrichment — the same shape `msigdbr` tables take on the [Functional Enrichment](https://cole-trapnell-lab.github.io/platt/enrichment/) page.

- **combined_psg**  
  A `cell_state_graph` object (or `igraph`) used for lineage traversal.

- **ai_notes_path**  
  Optional directory of precomputed per-cell-type markdown background files, forwarded to `llm_fun`'s context; ignored when `llm_fun` is `NULL`.

- **sig_p_val_thresh**  
  Significance cutoff used for DEG/pathway calls.

- **genes_of_interest**  
  Optional character vector used to highlight DEGs (e.g. a regulator/TF panel) when building the genes-of-interest line.

- **empirical_p_thresh**  
  Optional artifact-aware significance cutoff (`empirical_p`) applied to the genes-of-interest foreground; `1` disables it and falls back to nominal `sig_p_val_thresh`.

- **llm_fun**  
  Optional function used to generate a per-cell-type narrative summary. Defaults to `NULL`, which builds the impact table without calling an LLM — see Details.

- **max_lineage_depth**  
  Optional maximum graph depth from roots for included cell types.

- **cell_types**  
  Optional explicit cell types to seed lineage traversal, instead of inferring them from `target_gene_expression`.

- **primary_impact_summary**  
  Optional text summary of a prior (e.g. cell-autonomous) pass, forwarded to `llm_fun` as context.

- **excluded_cell_types**  
  Optional cell types to exclude from the result.

- **abundance_phenotypes**  
  Optional precomputed abundance phenotype table (`cell_group`, `abundance_code`, `abundance_severity`) — the same fields you'd otherwise type by hand into an impact table, see [Phenotyping](https://cole-trapnell-lab.github.io/platt/phenotyping/).

- **fitness_phenotypes**  
  Optional precomputed fitness phenotype table (`cell_group`, `fitness_label`, ...).

- **identity_phenotypes**  
  Optional precomputed identity phenotype table (`cell_group`, `identity_label`, ...).

- **dact_expectations**  
  Optional per-cell-type literature-expectation table, forwarded to `llm_fun` as context; ignored when `llm_fun` is `NULL`.

- **effect_type**  
  Optional label (e.g. `"Cell-autonomous"`) recorded for this pass and forwarded to `llm_fun` as context.

- **pre_cited_gene_claims**  
  Optional pre-cited gene-background text forwarded to `llm_fun`; ignored when `llm_fun` is `NULL`.

- **verbose**  
  Logical; if `TRUE`, prints debug messages instead of a progress bar.

- **sig_pathways_by_cell_type**  
  Optional named list of precomputed marker-specificity GSEA results, keyed by cell type.

- **sig_pathways_path**  
  Optional directory of precomputed `sig_pathways.rds` files, one per cell type.

- **drop_llm_columns**  
  If `TRUE`, drop the LLM-derived columns (`llm_summary`, `lineage_context_summary`, `llm_disrupted_pathways`, `llm_other_dysregulated_genes`) from the result — they are always empty when `llm_fun` is `NULL`.

- **...**  
  Additional arguments forwarded to `llm_fun`.

## Value

A tibble with one row per cell type, including `abundance_code`, `abundance_severity`, `identity_label`, `fitness_label`, per-cell-type DEGs, and pathway enrichment — and, unless `drop_llm_columns = TRUE`, the LLM narrative columns (always `NA`/empty when `llm_fun` is `NULL`).

## Details

An optional `llm_fun` can be supplied to additionally have an LLM write a short mechanistic summary for each cell type with a phenotype, genes of interest, or pathway enrichment of its own. Leaving `llm_fun` at its default of `NULL` skips that step entirely — no LLM server, Python, or `reticulate` setup required — and returns the same tabular impact calls (abundance/identity/fitness plus deterministic pathway enrichment) without any generated prose. Set `drop_llm_columns = TRUE` in that case to omit the always-empty LLM columns.

## Examples

```r
tbx16_impact_table <- summarize_impact_in_lineage_context(
  target_gene_expression = tbx16_target_gene_expression,
  dact_results           = tbx16_dact_summary_tbl,
  degs                   = tbx16_degs,
  ref_expression         = ref_expression,
  ontologies             = list(BP = gene_set_BP, MF = gene_set_MF, CC = gene_set_CC),
  combined_psg           = muscle_state_graph,
  genes_of_interest      = regulators$gene_short_name,
  abundance_phenotypes   = tbx16_abundance_phenotypes,
  fitness_phenotypes     = tbx16_fitness_phenotypes,
  identity_phenotypes    = tbx16_identity_phenotypes,
  llm_fun                = NULL,
  drop_llm_columns       = TRUE
)
```
