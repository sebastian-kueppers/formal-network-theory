# helper: scatter_states()
# requires: ggplot2

# scatter plot of the two emotion states (prep_pair() output); axes run from 0
# to the number of nodes of each emotion
#   colour = TRUE : points coloured by kind of goal activation
#   colour = FALSE: black points (panel figures, see pair_figure())
#   emotions      : names of E1 and E2, used for the axis titles and the legend,
#                   e.g. c("love", "sadness")
scatter_states <- function(p, colour = TRUE, title = p$label,
                           emotions = c("E1", "E2"),
                           xlab = paste("Intensity of", emotions[1]),
                           ylab = paste("Intensity of", emotions[2]),
                           check = TRUE) {
  d <- p$data
  
  if (check) {
    if (!isTRUE(all.equal(d$stateE1, unname(rowSums(d[p$names.E1])))) ||
        !isTRUE(all.equal(d$stateE2, unname(rowSums(d[p$names.E2])))))
      warning("stateE1/stateE2 are not the row sums of names.E1/names.E2; ",
              "max.E1/max.E2 may be wrong.")
  }
  if (anyNA(d$kind)) warning("Unexpected values in `kind` were set to NA and will be dropped.")
  
  g <- if (colour) ggplot(d, aes(stateE1, stateE2, colour = kind))
       else        ggplot(d, aes(stateE1, stateE2))
  g <- g +
    geom_point(alpha = 0.3, size = 0.8) +
    scale_x_continuous(limits = c(0, p$max.E1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, p$max.E2), expand = c(0, 0)) +
    # labs(x = xlab, y = ylab, title = title) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(aspect.ratio = 1,
          panel.grid   = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey70"),
          plot.margin  = margin(t = 5.5, r = 11, b = 5.5, l = 11))
  if (colour)
    g <- g +
      scale_colour_manual(
        values = c(E1 = "#E69F00", E2 = "#009E73", E1E2 = "black"),
        labels = c(E1 = emotions[1], E2 = emotions[2],
                   E1E2 = if (identical(emotions, c("E1", "E2"))) "E1E2"
                          else paste(emotions, collapse = " + ")),
        drop = FALSE, name = "Activated") +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
  g
}
