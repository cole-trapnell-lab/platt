#' Annotate a within-state DEG table with an empirical FDR for the
#' low-expression control-sampling artifact
#'
#' In control-heavy perturbation designs the deep control arm detects lowly
#' expressed genes that the shallow perturbation arm misses, producing a large
#' excess of low-expression, "on-in-control" (down) DEG calls that are sampling
#' artifacts rather than biology. This helper estimates, per
#' (cell_group x expression-bin x direction), the rate of calls attributable to
#' that artifact and returns the input rows annotated with
#' \code{empirical_null_rate}, \code{observed_rate}, and \code{empirical_fdr}.
#'
#' The null is read directly off the perturbation's own \emph{unaffected} cell
#' types (perturbation biology only \emph{adds} calls, so cell types with no real
#' effect are self-null), matched to each cell type by abundance so that shallow
#' cell types are compared against shallow ones. No permutation, model fitting, or
#' re-run is required, so existing DEG tables can be decorated in place.
#'
#' @param deg_tbl A data.frame / tibble of within-state DEG results. Must contain
#'   \code{cell_group}, \code{log_mean_expression}, \code{perturb_to_ctrl_p_value}
#'   (the ashr local false sign rate used as the significance quantity), and
#'   \code{perturb_to_ctrl_shrunken_lfc}.
#' @param abundances Optional. Per-cell-type abundance used to stratify the null.
#'   Either a named numeric vector (\code{cell_group} -> total cells) or a
#'   data.frame with columns \code{cell_group} and \code{abund}. If \code{NULL},
#'   a depth proxy is derived from \code{mean_log_sf} (mean cells per pseudobulk),
#'   which must then be present in \code{deg_tbl}. Supplying real cell counts is
#'   more accurate and recommended.
#' @param unaffected_cell_types Optional character vector of \code{cell_group}
#'   names known to have no real perturbation effect (e.g. abundance "A0" and no
#'   fitness/identity phenotype), used as the self-null training set. If
#'   \code{NULL}, a lower-envelope fallback is used: within each abundance bin the
#'   cell types in the bottom tertile of overall call rate are treated as
#'   unaffected.
#' @param sig_thresh Significance threshold defining a call (default 0.05).
#' @param n_abund_bins Number of abundance strata (default 4).
#' @param expr_breaks Numeric breakpoints used to bin \code{log_mean_expression}.
#' @param lower_envelope_frac Fraction of each abundance bin's cell types treated
#'   as the self-null when \code{unaffected_cell_types} is \code{NULL} (default
#'   1/3).
#'
#' @return \code{deg_tbl} with three added columns joined per
#'   (\code{cell_group}, expression-bin, direction): \code{empirical_null_rate}
#'   (percent of admitted genes called in the matched unaffected cell types),
#'   \code{observed_rate} (percent called in this cell type), and
#'   \code{empirical_fdr} (\code{min(1, empirical_null_rate / observed_rate)}).
#'   Rows in strata with no unaffected estimate receive \code{empirical_fdr = 0}
#'   (i.e. not flagged), the conservative-for-discovery default.
#'
#' @details Estimation is direction-specific because the artifact is
#'   directional: in control-heavy designs it inflates down (on-in-control)
#'   calls, while up calls at high expression are largely genuine. Consumers can
#'   threshold \code{empirical_fdr} (e.g. drop calls with \code{empirical_fdr}
#'   above some cutoff) to suppress the artifact while retaining real signal.
#'
#' @export
annotate_empirical_fdr <- function(deg_tbl,
                                   abundances = NULL,
                                   unaffected_cell_types = NULL,
                                   sig_thresh = 0.05,
                                   n_abund_bins = 4L,
                                   expr_breaks = c(-Inf, -5, -4, -3.5, -3, -2.5, -2, -1, Inf),
                                   lower_envelope_frac = 1 / 3) {
  req <- c("cell_group", "log_mean_expression",
           "perturb_to_ctrl_p_value", "perturb_to_ctrl_shrunken_lfc")
  missing_cols <- setdiff(req, colnames(deg_tbl))
  if (length(missing_cols) > 0) {
    stop("deg_tbl is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  d <- as.data.frame(deg_tbl, stringsAsFactors = FALSE)
  d$cell_group <- as.character(d$cell_group)
  d$.row_id <- seq_len(nrow(d))

  finite_row <- is.finite(d$log_mean_expression) &
    is.finite(d$perturb_to_ctrl_p_value)

  # --- per-cell-type abundance ---
  if (is.null(abundances)) {
    if (!"mean_log_sf" %in% colnames(d)) {
      stop("Provide `abundances`, or a deg_tbl containing `mean_log_sf` for the ",
           "depth-proxy fallback.")
    }
    ab_agg <- stats::aggregate(mean_log_sf ~ cell_group, data = d,
                               FUN = function(x) exp(stats::median(x, na.rm = TRUE)))
    ab <- data.frame(cell_group = as.character(ab_agg$cell_group),
                     abund = ab_agg$mean_log_sf, stringsAsFactors = FALSE)
  } else if (is.data.frame(abundances)) {
    ab <- data.frame(cell_group = as.character(abundances$cell_group),
                     abund = as.numeric(abundances$abund), stringsAsFactors = FALSE)
  } else {
    ab <- data.frame(cell_group = names(abundances),
                     abund = as.numeric(abundances), stringsAsFactors = FALSE)
  }
  ab <- ab[is.finite(ab$abund) & ab$abund > 0, , drop = FALSE]

  # abundance bin (on log scale, quantile cuts)
  qs <- unique(stats::quantile(log(ab$abund),
                               probs = seq(0, 1, length.out = n_abund_bins + 1L),
                               na.rm = TRUE))
  ab$abund_bin <- if (length(qs) < 2) 1L else
    as.integer(cut(log(ab$abund), breaks = qs, include.lowest = TRUE))

  d <- merge(d, ab[, c("cell_group", "abund_bin")], by = "cell_group",
             all.x = TRUE, sort = FALSE)
  d <- d[order(d$.row_id), , drop = FALSE]

  d$direction <- ifelse(d$perturb_to_ctrl_shrunken_lfc < 0, "down", "up")
  d$called <- d$perturb_to_ctrl_p_value < sig_thresh
  d$expr_bin <- cut(d$log_mean_expression, breaks = expr_breaks)
  # Direction-specific call indicators over ALL genes in a (cell_group, expr_bin):
  # the artifact/biology "rate" is the fraction of admitted genes called in each
  # direction, not a within-direction detection probability.
  d$dn_call <- as.numeric(d$called & d$perturb_to_ctrl_shrunken_lfc < 0)
  d$up_call <- as.numeric(d$called & d$perturb_to_ctrl_shrunken_lfc > 0)

  usable <- finite_row & !is.na(d$abund_bin) & !is.na(d$expr_bin)

  # --- unaffected (self-null) cell types ---
  if (is.null(unaffected_cell_types)) {
    ct_rate <- stats::aggregate(cbind(dn_call, up_call) ~ cell_group + abund_bin,
                                data = d[usable, , drop = FALSE], FUN = mean)
    ct_rate$anycall <- ct_rate$dn_call + ct_rate$up_call
    ct_rate$thr <- stats::ave(ct_rate$anycall, ct_rate$abund_bin,
                              FUN = function(x) stats::quantile(x, lower_envelope_frac,
                                                               na.rm = TRUE))
    unaffected_cell_types <- ct_rate$cell_group[ct_rate$anycall <= ct_rate$thr]
  }
  d$unaff <- d$cell_group %in% unaffected_cell_types

  # --- null down/up rate per (abund_bin x expr_bin) from unaffected cells ---
  null_tbl <- stats::aggregate(cbind(dn_call, up_call) ~ abund_bin + expr_bin,
                               data = d[usable & d$unaff, , drop = FALSE],
                               FUN = function(x) 100 * mean(x))
  names(null_tbl)[names(null_tbl) %in% c("dn_call", "up_call")] <-
    c("null_down", "null_up")

  # --- observed down/up rate per (cell_group x abund_bin x expr_bin) ---
  obs_tbl <- stats::aggregate(cbind(dn_call, up_call) ~ cell_group + abund_bin + expr_bin,
                              data = d[usable, , drop = FALSE],
                              FUN = function(x) 100 * mean(x))
  names(obs_tbl)[names(obs_tbl) %in% c("dn_call", "up_call")] <-
    c("obs_down", "obs_up")

  ann <- merge(obs_tbl, null_tbl, by = c("abund_bin", "expr_bin"),
               all.x = TRUE, sort = FALSE)
  ann$null_down[is.na(ann$null_down)] <- 0
  ann$null_up[is.na(ann$null_up)] <- 0
  eps <- .Machine$double.eps
  ann$fdr_down <- pmin(1, ann$null_down / pmax(ann$obs_down, eps))
  ann$fdr_up <- pmin(1, ann$null_up / pmax(ann$obs_up, eps))

  d <- merge(d,
             ann[, c("cell_group", "abund_bin", "expr_bin",
                     "null_down", "null_up", "obs_down", "obs_up",
                     "fdr_down", "fdr_up")],
             by = c("cell_group", "abund_bin", "expr_bin"),
             all.x = TRUE, sort = FALSE)
  d <- d[order(d$.row_id), , drop = FALSE]
  is_down <- d$direction == "down"

  out <- as.data.frame(deg_tbl, stringsAsFactors = FALSE)
  out$empirical_null_rate <- ifelse(is_down, d$null_down, d$null_up)
  out$observed_rate <- ifelse(is_down, d$obs_down, d$obs_up)
  out$empirical_fdr <- ifelse(is_down, d$fdr_down, d$fdr_up)
  out
}
