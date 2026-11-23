#layout for display only: nodes cluster by emotion, shared nodes between the clusters
#  within : attraction among nodes of the same emotion
#  between: weak attraction between the emotions (keeps the clusters close)
cluster_layout <- function(p, within = 1, between = 0.02, seed = 1) {
  nodes <- union(p$names.E1, p$names.E2)       #same order as in plot_network()
  in.E1 <- nodes %in% p$names.E1
  in.E2 <- nodes %in% p$names.E2
  L <- (outer(in.E1, in.E1) + outer(in.E2, in.E2)) * within
  L[L == 0] <- between
  L <- L + abs(p$W[nodes, nodes])              #real edges add a little structure
  diag(L) <- 0
  set.seed(seed)
  qgraph(L, layout = "spring", DoNotPlot = TRUE)$layout
}
