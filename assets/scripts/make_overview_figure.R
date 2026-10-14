#!/usr/bin/env Rscript
#
# Generates docs/assets/platt_background.png, the homepage overview figure,
# entirely from code. All six panels are illustrative/synthetic (not real
# pipeline output) but drawn with the actual visual grammar Platt/Hooke docs
# use elsewhere: small directed trees with colored nodes for graph-shaped
# results, viridis/diverging fills for continuous values, and minimal
# scatter+fit panels for kinetics.
#
# Run from the platt_tutorial repo root:
#   Rscript docs/assets/scripts/make_overview_figure.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

title_theme <- theme(
  plot.title = element_text(size = 21, hjust = 0, margin = margin(b = 10))
)

# patchwork::plot_annotation(title=) only renders on the OUTERMOST composed
# object; once a titled sub-assembly is nested inside a further `/` or `|`,
# its title is silently dropped. So each row title is its own tiny ggplot,
# stacked above its content with plot_layout(), which survives nesting fine.
panel_title <- function(text) {
  ggplot() +
    labs(title = text) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(size = 21, hjust = 0, margin = margin(t = 16, b = 2)))
}

titled_panel <- function(text, content, title_height = 1.3, content_height = 7) {
  panel_title(text) / content + plot_layout(heights = c(title_height, content_height))
}

# ---------------------------------------------------------------------------
# shared helper: small directed tree of colored nodes
# ---------------------------------------------------------------------------

draw_tree <- function(edges, nodes, fill_scale = scale_fill_identity(),
                       show_legend = FALSE, node_size = 7, arrow_color = "grey65") {
  p <- ggplot()
  if (!is.null(edges) && nrow(edges) > 0) {
    p <- p + geom_segment(
      data = edges, aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.07, "in"), type = "closed"),
      color = arrow_color, linewidth = 0.5
    )
  }
  p <- p +
    geom_point(
      data = nodes, aes(x = x, y = y, fill = fill),
      shape = 21, size = node_size, color = "grey20", stroke = 0.4
    ) +
    fill_scale +
    coord_equal(clip = "off") +
    theme_void(base_size = 12) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(4, 4, 4, 4)
    )
  p
}

qual_pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
              "#E6AB02", "#A6761D", "#666666")

# ---------------------------------------------------------------------------
# panel 1: wild-type & perturbation kinetic modeling
# ---------------------------------------------------------------------------

mini_axis_theme <- theme_classic(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 6, color = "grey40"),
    axis.line = element_line(color = "grey50", linewidth = 0.3),
    axis.ticks = element_line(color = "grey50", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(4, 8, 2, 2)
  )

plot_kinetic_decline <- function() {
  x <- seq(15, 48, length.out = 150)
  y <- 12 / (1 + exp((x - 32) / 4))
  pts <- data.frame(
    x = x,
    y = pmax(y + rnorm(length(x), 0, 0.35 * y + 0.4), 0.05),
    grp = sample(1:3, length(x), TRUE)
  )
  ggplot(pts, aes(x, y)) +
    geom_point(aes(color = factor(grp)), size = 0.7, alpha = 0.55) +
    geom_line(data = data.frame(x = x, y = y), aes(x, y), color = "#B2182B", linewidth = 0.9) +
    scale_color_manual(values = c("#D6604D", "#4393C3", "#2166AC")) +
    scale_y_log10() +
    mini_axis_theme
}

plot_kinetic_rise <- function() {
  x <- seq(15, 48, length.out = 150)
  y <- 3 + 27 * exp(-((x - 33)^2) / (2 * 9^2))
  pts <- data.frame(
    x = x,
    y = pmax(y + rnorm(length(x), 0, 0.3 * y + 0.6), 0.1),
    grp = sample(1:3, length(x), TRUE)
  )
  ggplot(pts, aes(x, y)) +
    geom_point(aes(color = factor(grp)), size = 0.7, alpha = 0.55) +
    geom_line(data = data.frame(x = x, y = y), aes(x, y), color = "#B2182B", linewidth = 0.9) +
    scale_color_manual(values = c("#D6604D", "#4393C3", "#2166AC")) +
    scale_y_log10() +
    mini_axis_theme
}

plot_kinetic_ribbon <- function() {
  x <- seq(15, 38, length.out = 100)
  hi <- 60 * exp(-(x - 15) / 12)
  lo <- 60 * exp(-(x - 15) / 7)
  band <- data.frame(x = x, ymin = pmin(hi, lo), ymax = pmax(hi, lo))
  pts <- data.frame(x = sample(x, 60), y = NA)
  pts$y <- approx(x, (hi + lo) / 2, xout = pts$x)$y * exp(rnorm(60, 0, 0.3))

  ggplot() +
    geom_ribbon(data = band, aes(x = x, ymin = pmax(ymin, 0.05), ymax = ymax),
                fill = "#D6604D", alpha = 0.35) +
    geom_point(data = pts, aes(x, pmax(y, 0.05)), size = 0.6, color = "grey60", alpha = 0.6) +
    geom_line(data = band, aes(x, ymax), linewidth = 0.4, color = "#B2182B") +
    geom_line(data = band, aes(x, pmax(ymin, 0.05)), linewidth = 0.4, color = "#B2182B", linetype = "22") +
    scale_y_log10() +
    mini_axis_theme
}

kinetic_panel <- titled_panel(
  "wild-type & perturbation kinetic modeling",
  plot_kinetic_decline() | plot_kinetic_rise() | plot_kinetic_ribbon()
)

# ---------------------------------------------------------------------------
# panel 2: graph construction
# ---------------------------------------------------------------------------

tree_small <- draw_tree(
  edges = data.frame(x = c(1, 1, 1), y = c(2, 2, 2), xend = c(0.4, 1, 1.6), yend = c(1, 1, 1)),
  nodes = data.frame(x = c(1, 0.4, 1, 1.6), y = c(2, 1, 1, 1), fill = qual_pal[c(1, 2, 3, 4)])
)

tree_medium <- draw_tree(
  edges = data.frame(
    x = c(1, 1, 1, 0.55, 1.45),
    y = c(2.4, 2.4, 2.4, 1.6, 1.6),
    xend = c(0.55, 1, 1.45, 0.55, 1.45),
    yend = c(1.6, 1.6, 1.6, 0.8, 0.8)
  ),
  nodes = data.frame(
    x = c(1, 0.55, 1, 1.45, 0.55, 1.45),
    y = c(2.4, 1.6, 1.6, 1.6, 0.8, 0.8),
    fill = qual_pal[c(2, 5, 6, 1, 7, 8)]
  )
)

tree_fan <- draw_tree(
  edges = data.frame(
    x = rep(1, 6), y = rep(2.2, 6),
    xend = c(1, 0.35, 0.68, 1.32, 1.65, 1), yend = c(1.6, 1, 1, 1, 1, 1)
  ),
  nodes = data.frame(
    x = c(1, 1, 0.35, 0.68, 1.32, 1.65, 1),
    y = c(2.2, 1.6, 1, 1, 1, 1, 1),
    fill = c(qual_pal[5], qual_pal[2], qual_pal[6], qual_pal[3], qual_pal[7], qual_pal[8], qual_pal[4])
  )
)

graph_panel <- titled_panel("graph construction", tree_small | tree_medium | tree_fan)

# ---------------------------------------------------------------------------
# panel 3: identification of genetic requirements
# ---------------------------------------------------------------------------

viridis_fill <- scale_fill_gradientn(colors = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
                                      limits = c(0, 1.2), name = "gene expr")

req_branch <- draw_tree(
  edges = data.frame(
    x = c(1, 1, 0.6), y = c(2.6, 2.6, 1.8),
    xend = c(0.6, 1.4, 0.6), yend = c(1.8, 1.8, 1.0)
  ),
  nodes = data.frame(
    x = c(1, 0.6, 1.4, 0.6), y = c(2.6, 1.8, 1.8, 1.0),
    fill = c(0.5, 1.15, 0.15, 0.9)
  ),
  fill_scale = viridis_fill
)

req_star <- draw_tree(
  edges = data.frame(
    x = rep(1, 5), y = rep(2, 5),
    xend = c(0.3, 0.65, 1, 1.35, 1.7), yend = rep(1, 5)
  ),
  nodes = data.frame(
    x = c(1, 0.3, 0.65, 1, 1.35, 1.7), y = c(2, 1, 1, 1, 1, 1),
    fill = c(0.5, 0.1, 1.15, 0.2, 0.15, 0.1)
  ),
  fill_scale = viridis_fill, show_legend = TRUE
)

requirements_panel <- titled_panel("identification of genetic requirements", req_branch | req_star)

# ---------------------------------------------------------------------------
# panel 4: differential expression testing
# ---------------------------------------------------------------------------

fc_fill <- scale_fill_gradient2(low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0,
                                 limits = c(-1.5, 1.5), oob = scales::squish, name = "log2(fc)")

deg_simple <- draw_tree(
  edges = data.frame(x = 1, y = 2, xend = 1, yend = 1.2),
  nodes = data.frame(x = c(1, 1), y = c(2, 1.2), fill = c(0.4, -1.1)),
  fill_scale = fc_fill
)

deg_branch <- draw_tree(
  edges = data.frame(x = c(1, 1, 1), y = c(2.2, 1.4, 1.4), xend = c(1, 0.6, 1.4), yend = c(1.4, 0.6, 0.6)),
  nodes = data.frame(x = c(1, 1, 0.6, 1.4), y = c(2.2, 1.4, 0.6, 0.6), fill = c(-0.2, 1.2, -1.3, 0.1)),
  fill_scale = fc_fill, show_legend = TRUE
)

deg_panel <- titled_panel("differential expression testing", deg_simple | deg_branch)

# ---------------------------------------------------------------------------
# panel 5: phenotyping
# ---------------------------------------------------------------------------

grey_fill <- "grey82"

plot_transition <- function(title, subtitle, delta, label) {
  ctrl <- data.frame(x = 1, y = 1)
  outcome <- data.frame(x = 2.15, y = 1, delta = delta, size = 12 + min(abs(delta), 2.5) * 6)
  text_col <- if (abs(delta) >= 1) "white" else "grey20"

  ggplot() +
    geom_segment(x = 1.32, xend = 1.85, y = 1, yend = 1,
                 arrow = arrow(length = unit(0.09, "inches"), type = "closed"),
                 color = "grey60", linewidth = 0.6) +
    geom_point(data = ctrl, aes(x, y), size = 12, shape = 21, fill = grey_fill, color = "grey20", stroke = 0.5) +
    geom_text(data = ctrl, aes(x, y), label = "ctrl", size = 3, color = "grey20", fontface = "bold") +
    geom_point(data = outcome, aes(x, y, size = size, fill = delta), shape = 21, color = "grey20", stroke = 0.5) +
    geom_text(data = outcome, aes(x, y), label = label, size = 3, color = text_col, fontface = "bold") +
    scale_size_identity() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
                          limits = c(-2.5, 2.5), oob = scales::squish, name = "delta log\nabundance") +
    coord_cartesian(xlim = c(0.6, 2.8), ylim = c(0.5, 1.5), clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 7.5, color = "grey35", hjust = 0.5, margin = margin(b = 6)),
      plot.margin = margin(4, 8, 4, 8)
    )
}

phenotyping_content <- (
  plot_transition("No change", "|delta| < 0.5 or q > 0.1", 0, "+0.0") |
  plot_transition("Expansion", "delta >= 0.5, q < 0.1", 1.5, "+1.5") |
  plot_transition("Depletion", "delta <= -0.5, q < 0.1", -1.0, "-1.0")
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", legend.key.height = unit(0.3, "in"), legend.title = element_text(size = 7.5))

phenotyping_panel <- titled_panel("phenotyping", phenotyping_content)

# ---------------------------------------------------------------------------
# panel 6: gene pattern discovery
# ---------------------------------------------------------------------------

plot_pattern <- function(from_fill, to_fill, label) {
  draw_tree(
    edges = data.frame(x = 1, y = 2, xend = 1, yend = 1.2),
    nodes = data.frame(x = c(1, 1), y = c(2, 1), fill = c(from_fill, to_fill)),
    node_size = 9, arrow_color = "grey15"
  ) +
    annotate("text", x = 1, y = 0.55, label = label, fontface = "bold", size = 3.6) +
    coord_equal(xlim = c(0.6, 1.4), ylim = c(0.3, 2.3), clip = "off")
}

patterns_panel <- titled_panel(
  "gene pattern discovery",
  plot_pattern(grey_fill, "#21908C", "activated") |
    plot_pattern("#21908C", grey_fill, "deactivated") |
    plot_pattern("#440154", "#5DC863", "upregulated")
)

# ---------------------------------------------------------------------------
# assemble
# ---------------------------------------------------------------------------

full <- (kinetic_panel | graph_panel) /
  (requirements_panel | deg_panel) /
  (phenotyping_panel | patterns_panel)

dir.create("docs/assets", showWarnings = FALSE)
ggsave("docs/assets/platt_background.png", full, width = 15, height = 10.5, dpi = 150, bg = "white")
cat("wrote docs/assets/platt_background.png\n")
