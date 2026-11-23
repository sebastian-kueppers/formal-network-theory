# helper: net_panel()
# requires: ggplot2, png, grid

# exported network PNG -> ggplot panel (for patchwork figures)
net_panel <- function(path) {
  img <- readPNG(path)
  ggplot() +
    annotation_custom(rasterGrob(img, interpolate = TRUE),
                      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(aspect.ratio = nrow(img) / ncol(img),
          plot.margin  = margin(t = 5.5, r = 11, b = 5.5, l = 11))
}
