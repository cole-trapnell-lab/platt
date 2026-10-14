# `plot_cell_type_control_kinetics`

Generates a kinetic plot of cell type control data over time.

```r
plot_cell_type_control_kinetics(
  control_ccm,
  cell_groups = NULL,
  start_time = NULL,
  stop_time = NULL,
  interval_step = 1,
  log_abund_detection_thresh = -3,
  delta_log_abund_loss_thresh = 0,
  interval_col = "timepoint",
  q_val = 0.01,
  min_log_abund = -5,
  log_scale = TRUE,
  group_nodes_by = "cell_type",
  nrow = 1,
  newdata = tibble::tibble(),
  color_points_by = NULL,
  size = 0.5,
  alpha = 0.5,
  linewidth = 1,
  raw_counts = FALSE,
  marginalize_over = "auto"
)
```

## Arguments

- **control_ccm**  
  *object*  
  A control cell count matrix object.

- **cell_groups**  
  *vector*  
  Cell groups to include in the plot. Default is `NULL`.

- **start_time**  
  *numeric*  
  Start time for the plot. If `NULL`, uses minimum timepoint in the data.

- **stop_time**  
  *numeric*  
  Stop time for the plot. If `NULL`, uses maximum timepoint in the data.

- **interval_step**  
  *numeric*  
  Step size for time intervals. Default is `1`.

- **log_abund_detection_thresh**  
  *numeric*  
  Log abundance detection threshold. Default is `-3`.

- **delta_log_abund_loss_thresh**  
  *numeric*  
  Delta log abundance loss threshold. Default is `0`.

- **interval_col**  
  *character*  
  Column name for time intervals. Default is `"timepoint"`.

- **q_val**  
  *numeric*  
  Q-value threshold for significance. Default is `0.01`.

- **min_log_abund**  
  *numeric*  
  Minimum log abundance. Default is `-5`.

- **log_scale**  
  *logical*  
  Whether to use a log scale for the y-axis. Default is `TRUE`.

- **group_nodes_by**  
  *character*  
  Column name to group nodes by. Default is `"cell_type"`.

- **nrow**  
  *integer*  
  Number of rows in the facet wrap. Default is `1`.

- **newdata**  
  *tibble*  
  New data to use in the plot. Default is an empty tibble.

- **color_points_by**  
  *character*  
  Column name to color points by. Default is `NULL`.

- **size**  
  *numeric*  
  Size of the points in the plot. Default is `0.5`.

- **alpha**  
  *numeric*  
  Alpha transparency of points. Default is `0.5`.

- **linewidth**  
  *numeric*  
  Line width of the fitted abundance curves. Default is `1`.

- **raw_counts**  
  *logical*  
  Whether to use raw counts. Default is `FALSE`.

- **marginalize_over**  
  *character or NULL*  
  `"auto"` (default): the drawn curve and the cell-type ordering are the model's marginal over its own nuisance factor (`estimate_abundances_marginal()`), or a plain prediction when the model has no such factor. A column name forces that factor; `NULL` restores the pre-2026-09 prediction at the factor's first level.

## Value

A `ggplot` object representing the kinetic plot.

## Examples

```r
# Example usage:
plot <- plot_cell_type_control_kinetics(control_ccm)
print(plot)
```
