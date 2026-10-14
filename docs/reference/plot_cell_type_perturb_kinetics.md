# `plot_cell_type_perturb_kinetics`

Generates a plot to visualize the kinetics of cell type perturbations over time.

```r
plot_cell_type_perturb_kinetics(
  perturbation_ccm,
  cell_groups = NULL,
  start_time = NULL,
  stop_time = NULL,
  interval_step = 1,
  log_abund_detection_thresh = -3,
  delta_log_abund_loss_thresh = 0,
  interval_col = "timepoint",
  q_val = 0.01,
  log_scale = TRUE,
  control_ccm = perturbation_ccm,
  control_start_time = start_time,
  control_stop_time = stop_time,
  group_nodes_by = "cell_type",
  newdata = tibble::tibble(),
  nrow = 1,
  size = 0.75,
  alpha = 0.5,
  raw_counts = FALSE
)
```

## Arguments

- **perturbation_ccm**  
  *CCM object*  
  A perturbation cell count matrix object.

- **cell_groups**  
  *vector*  
  Cell groups to include. If `NULL`, all groups are included.

- **start_time**  
  *numeric*  
  Start time. If `NULL`, uses minimum timepoint in data.

- **stop_time**  
  *numeric*  
  Stop time. If `NULL`, uses maximum timepoint in data.

- **interval_step**  
  *numeric*  
  Step size for time intervals.

- **log_abund_detection_thresh**  
  *numeric*  
  Log abundance detection threshold.

- **delta_log_abund_loss_thresh**  
  *numeric*  
  Change in log abundance loss threshold.

- **interval_col**  
  *character*  
  Column name for time intervals.

- **q_val**  
  *numeric*  
  Q-value threshold for significance.

- **log_scale**  
  *logical*  
  Use log scale for y-axis.

- **control_ccm**  
  *CCM object*  
  Control cell count matrix. Defaults to `perturbation_ccm`.

- **control_start_time**  
  *numeric*  
  Start time for control interval. Defaults to `start_time`.

- **control_stop_time**  
  *numeric*  
  Stop time for control interval. Defaults to `stop_time`.

- **group_nodes_by**  
  *character*  
  Column to group nodes by.

- **newdata**  
  *tibble*  
  New data for predictions.

- **nrow**  
  *integer*  
  Number of rows in facet wrap.

- **size**  
  *numeric*  
  Size of points in plot.

- **alpha**  
  *numeric*  
  Alpha transparency for difference shading.

- **raw_counts**  
  *logical*  
  Use raw counts.

## Value

A `ggplot` object representing the kinetics of cell type perturbations.

## Examples

```r
# Example usage:
plot <- plot_cell_type_perturb_kinetics(perturbation_ccm)
print(plot)
```
