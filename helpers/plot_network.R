# helper: plot_network()
# requires: qgraph

# plot the network of one pair (prep_pair() output); nodes coloured
# E1-unique (orange) / shared (yellow-green) / E2-unique (green). Colours are
# matched by node name, so they do not depend on the node order.
#   labels : node labels in the order of union(names.E1, names.E2)
#            (default: node names)
#   file   : if given, the plot is written to this PNG instead of the screen
#   seed   : the spring layout is random; fixed per call for reproducible figures
#   ...    : further arguments to qgraph(), e.g. vsize
plot_network <- function(p, labels = NULL, title = NULL, file = NULL,
                         layout = "spring", seed = 1, ...) {
  keep <- union(p$names.E1, p$names.E2)
  W <- p$W
  if (!all(keep %in% colnames(W)))
    stop("Nodes missing from the network: ",
         paste(setdiff(keep, colnames(W)), collapse = ", "))
  W <- W[keep, keep, drop = FALSE]
  
  cols <- ifelse(keep %in% p$shared, "#B5C400",
                 ifelse(keep %in% p$names.E1, "#E69F00", "#009E73"))
  
  if (is.null(labels)) labels <- keep
  stopifnot(length(labels) == length(keep))
  
  set.seed(seed)
  if (!is.null(file)) png(file, width = 7, height = 7, units = "in", res = 300)
  g <- qgraph(W, weighted = TRUE, theme = "colorblind", layout = layout,
              edge.width = 0.5, color = cols, labels = labels,
              title = title, ...)
  if (!is.null(file)) dev.off()
  invisible(g)
}
