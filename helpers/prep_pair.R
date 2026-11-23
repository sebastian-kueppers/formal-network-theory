# helper: prep_pair()

# metadata + MIN for one simulated pair; the returned object p is the input of
# plot_network(), scatter_states(), pair_figure() and test_lshape()
#   names.E1 / names.E2: node names of the two emotions (empirical networks).
#   If NULL (synthetic networks, e.g. toy example), they are derived from
#   sim$settings: the model orders the nodes E1-unique, shared, E2-unique and
#   names them node1, node2, ...
prep_pair <- function(sim, names.E1 = NULL, names.E2 = NULL, label = sim$cell) {
  if (is.null(names.E1) || is.null(names.E2)) {
    s <- sim$settings
    if (is.null(s))
      stop("names.E1/names.E2 not given and sim$settings not available.")
    n.E1 <- s$n.nodes.in.network.E1
    n.E2 <- s$n.nodes.in.network.E2
    n.sh <- s$n.shared.nodes
    names.E1 <- paste0("node", seq_len(n.E1))
    names.E2 <- paste0("node", (n.E1 - n.sh + 1):(n.E1 + n.E2 - n.sh))
  }
  
  d <- sim$data
  d$kind <- factor(d$kind, levels = c("E1", "E2", "E1E2"))
  d$MIN  <- pmin(d$stateE1, d$stateE2)
  
  # synthetic networks come without dimnames; name them like the data columns
  W <- sim$weights.matrix
  if (is.null(colnames(W))) {
    nodes <- union(names.E1, names.E2)
    dimnames(W) <- list(nodes, nodes)
  }
  
  shared <- intersect(names.E1, names.E2)
  list(data = d, label = label,
       names.E1 = names.E1, names.E2 = names.E2, shared = shared,
       max.E1 = length(names.E1), max.E2 = length(names.E2),
       max.min = min(length(names.E1), length(names.E2)),
       n.shared = length(shared), W = W)
}
