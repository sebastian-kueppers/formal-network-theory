# helper: pair_figure()
# requires: ggplot2, patchwork; net_panel(), scatter_states()

# per-pair figure: network (A) over scatter plot (B), printed and saved to file
pair_figure <- function(net.png, p, file) {
  fig <- net_panel(net.png) /
    scatter_states(p, colour = FALSE, title = NULL,
                   xlab = "stateE1", ylab = "stateE2") +
    plot_annotation(
      tag_levels = list(c("A", "B")),
      theme = theme(plot.tag = element_text(face = "bold", size = 16))
    )
  print(fig)
  ggsave(file, plot = fig, width = 6, height = 11, dpi = 300)
  invisible(fig)
}
