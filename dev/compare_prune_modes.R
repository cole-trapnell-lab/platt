#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(platt)
  library(tibble)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: dev/compare_prune_modes.R <ref_ccs.rds> <timeseries_graph.rds> <perturbation_ccm_tbl.rds> [output_prefix]")
}

ref_ccs <- readRDS(args[[1]])
timeseries_graph <- readRDS(args[[2]])
perturbation_ccm_tbl <- readRDS(args[[3]])
output_prefix <- if (length(args) >= 4) args[[4]] else NULL

report <- platt:::compare_discordant_pruning_modes(
  ref_ccs = ref_ccs,
  timeseries_graph = timeseries_graph,
  perturbation_ccm_tbl = perturbation_ccm_tbl
)

cat("\nEdge changes\n")
print(tibble::tibble(
  removed = report$edge_changes$n_removed,
  added = report$edge_changes$n_added
))

cat("\nDiscordant reachability violations\n")
print(report$discordant_violations)

cat("\nSupport distribution summary\n")
print(report$support_summary)

if (!is.null(output_prefix)) {
  readr::write_csv(report$edge_changes$removed, paste0(output_prefix, "_edges_removed.csv"))
  readr::write_csv(report$edge_changes$added, paste0(output_prefix, "_edges_added.csv"))
  readr::write_csv(report$discordant_violations, paste0(output_prefix, "_discordant_violations.csv"))
  readr::write_csv(report$support_summary, paste0(output_prefix, "_support_summary.csv"))
}
