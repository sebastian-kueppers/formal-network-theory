#' Create an Ising network object for cooccurrence.model()
#'
#' Bundles what cooccurrence.model() needs from an estimated Ising network: the edge weights,
#' the node thresholds (intercepts) and the coding of the data the network was estimated on.
#' The result can be passed as \code{existing.network}, whatever software produced the estimates
#' (e.g., glmnet, IsingFit, IsingSampler, bgms).
#'
#' @param weights Symmetric numeric matrix of edge weights with identical row and column names
#' (the node names) and a zero diagonal.
#' @param thresholds Numeric vector with one threshold (intercept) per node, named by node, or unnamed
#' in the column order of \code{weights}. May be \code{NULL} (default) if cooccurrence.model() is called with
#' fixed thresholds (\code{threshold.E1} and \code{threshold.E2}); it then stops if they are not supplied.
#' @param coding Coding of the data the estimates refer to: "0/1" (e.g., IsingFit, IsingSampler, glmnet
#' logistic regressions) or "-1/1". Must be given explicitly. Estimates on "0/1" are converted to the
#' "-1/1" scale used by the Glauber dynamics inside cooccurrence.model() (Haslbeck, Epskamp, Marsman &
#' Waldorp, 2021, Table 2); do not convert them beforehand.
#'
#' @return An object of class "ising_network": a list with elements \code{weights}, \code{thresholds}
#' (named, ordered like the columns of \code{weights}; \code{NULL} if none were given) and \code{coding}.
#'
#' @export
ising_network <- function(weights, thresholds = NULL, coding = c("0/1", "-1/1")) {
  if (missing(coding)) {
    stop("'coding' must be given explicitly: \"0/1\" or \"-1/1\" (the coding of the data the network was estimated on).")
  }
  coding <- match.arg(coding)
  
  weights <- as.matrix(weights)
  attr(weights, "intercepts") <- NULL
  if (!is.numeric(weights)) stop("'weights' must be a numeric matrix.")
  if (is.null(colnames(weights)) || !identical(rownames(weights), colnames(weights))) {
    stop("'weights' needs identical row and column names (the node names).")
  }
  if (!isTRUE(all.equal(unname(weights), unname(t(weights)), tolerance = 1e-8))) {
    stop("'weights' must be symmetric.")
  }
  if (any(abs(diag(weights)) > 1e-8)) stop("'weights' must have a zero diagonal.")
  
  if (!is.null(thresholds)) {
    if (is.matrix(thresholds)) {
      if (ncol(thresholds) != 1) stop("'thresholds' must be a vector or a one-column matrix.")
      tn <- rownames(thresholds)
      thresholds <- as.numeric(thresholds[, 1])
      names(thresholds) <- tn
    }
    if (!is.numeric(thresholds) || anyNA(thresholds)) stop("'thresholds' must be numeric without NA.")
    if (is.null(names(thresholds))) {
      if (length(thresholds) != ncol(weights)) {
        stop("Unnamed 'thresholds' need one value per node (", ncol(weights), "), but ",
             length(thresholds), " were given.")
      }
      names(thresholds) <- colnames(weights)
    }
    if (!all(colnames(weights) %in% names(thresholds))) {
      stop("No threshold for these nodes: ",
           paste(setdiff(colnames(weights), names(thresholds)), collapse = ", "))
    }
    thresholds <- thresholds[colnames(weights)]   # matched by name, never by position
  }
  
  structure(list(weights = weights,
                 thresholds = thresholds,   # NULL allowed (fixed thresholds in cooccurrence.model())
                 coding = coding),
            class = "ising_network")
}

#' Simulate Co-occurrence of Two Emotion Networks
#'
#' This function simulates the dynamics of two emotion networks based on goal activations
#' and Glauber dynamics. Networks can be generated as small-world, preferential attachment,
#' random, fully connected, or provided by the user.
#'
#' @param emotion.model Character. Type of network generation ("small.world.network",
#' "preferential.attachment.network", "random.network", "fully.connected.network",
#' "existing.network").
#' @param existing.network Optional. User-supplied network. Preferred: an Ising network object created with
#' \code{ising_network(weights, thresholds, coding)}, or a plain list with these three elements. Also accepted:
#' a bootnetResult from an Ising estimator (IsingFit, IsingSampler), or a weight matrix (then thresholds are read
#' from attribute "intercepts" or from \code{network.thresholds}, and \code{network.coding} must be set). qgraph
#' objects are accepted only together with \code{network.thresholds} and \code{network.coding}.
#' @param names.E1 Character vector. Node names for Emotion 1 (required if existing.network is provided).
#' @param names.E2 Character vector. Node names for Emotion 2 (required if existing.network is provided).
#' @param network.thresholds Optional. Numeric vector (named by node, or in the column order of
#' existing.network) or one-column matrix with the thresholds (intercepts) of existing.network. Only needed
#' if they are not stored in existing.network (bootnetResult$intercepts, or attribute "intercepts" of a matrix).
#' Mutually exclusive with fixed thresholds (\code{threshold.E1} and \code{threshold.E2}).
#' @param network.coding Character. Coding of the data existing.network was estimated on: "0/1", "-1/1",
#' or "auto" (default; taken from the 'coding' element of an ising_network, or inferred from bootnetResult$data;
#' not possible for matrices and qgraph objects). If both are given they must agree.
#' Weights and thresholds fitted on {0,1} data are converted to the {-1,1} scale (with fixed thresholds,
#' only the weights are converted; the fixed thresholds are already on the {-1,1} scale).
#' The conversion is exact only for unbiased estimates; for regularized fits (e.g.,
#' IsingFit) it is an approximation.
#' @param weight.E1 Numeric. Edge weights for Network 1. Default 0.1.
#' @param weight.E2 Numeric. Edge weights for Network 2. Default 0.1.
#' @param weight.pos.E1.E2 Numeric. Edge weight for positive cross-network edges. Default 0.08.
#' @param weight.neg.E1.E2 Numeric. Edge weight for negative cross-network edges. Default -0.08.
#' @param threshold.E1 Numeric. Thresholds for Network 1 nodes. Default -0.4. For emotion.model = "existing.network":
#' if both threshold.E1 and threshold.E2 are explicitly supplied, they replace the thresholds of the network
#' (fixed thresholds, {-1,1} scale; nodes only in E1 get threshold.E1, nodes only in E2 get threshold.E2,
#' shared nodes get their mean). Supplying only one of the two is an error.
#' @param threshold.E2 Numeric. Thresholds for Network 2 nodes. Default -0.4. See threshold.E1.
#' @param beta Numeric. Temperature parameter for Glauber dynamics. Default 1.
#' @param n.nodes.in.network.E1 Integer. Number of nodes for Network 1. Default 15.
#' @param n.nodes.in.network.E2 Integer. Number of nodes for Network 2. Default 10.
#' @param r.m.components.times.x Integer. Parameter for random network. Default 2.
#' @param sw.nei Integer. Small-world neighborhood size. Default 2.
#' @param sw.p Numeric. Small-world rewiring probability. Default 0.1.
#' @param pa.power Numeric. Preferential attachment power. Default 1.
#' @param pa.m Integer. Preferential attachment m parameter. Default 2.
#' @param n.shared.nodes Integer. Number of shared nodes between networks. Default 3.
#' @param n.cross.pos.edges Integer. Number of positive cross-network edges. Default 0.
#' @param n.cross.neg.edges Integer. Number of negative cross-network edges. Default 8.
#' @param cross.edge.type Character. "central" or "random" edge assignment. Default "random".
#' @param dist.goal.activation Numeric vector of length 3. Distribution of goal activations. Default c(0.45, 0.45, 0.1).
#' @param n.goal.activation Integer. Number of goal activations. Default 10000.
#' @param range.goal.activation Numeric. Range of possible goal activation intensities. Default 0.8.
#' @param n.Glauber Integer. Number of Glauber iterations. Default 1000.
#' @param n.ignored.Glauber Integer. Number of initial Glauber iterations to ignore. Default 50.
#'
#' @return A list containing:
#' \describe{
#'   \item{weights.matrix}{The symmetric weight matrix W of the network, on the {-1,1} scale used by the dynamics.}
#'   \item{initial.thresholds}{Vector of initial thresholds for all nodes (named if emotion.model = "existing.network").}
#'   \item{network.coding}{Only for existing networks: coding the supplied network was estimated on ("0/1" or "-1/1").}
#'   \item{data}{Data frame with final states after goal activations and Glauber dynamics, including states of shared nodes and E1/E2 networks.}
#'   \item{data.thresholds}{Matrix of thresholds at each step.}
#' }
#'
#' @importFrom igraph sample_pa sample_smallworld sample_gnm as_adj_list graph_from_adjacency_matrix
#' @importFrom qgraph getWmat
#' @importFrom Matrix forceSymmetric
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
cooccurrence.model <- function(
    emotion.model = c("small.world.network",
                      "preferential.attachment.network",
                      "random.network",
                      "fully.connected.network",
                      "existing.network"),
    existing.network = NULL,
    names.E1 = NULL,
    names.E2 = NULL,
    network.thresholds = NULL,
    network.coding = c("auto", "0/1", "-1/1"),
    weight.E1 = 0.1, weight.E2 = 0.1,
    weight.pos.E1.E2 = 0.08, weight.neg.E1.E2 = -0.08,
    threshold.E1 = -0.4, threshold.E2 = -0.4,
    beta = 1,
    n.nodes.in.network.E1 = 15,
    n.nodes.in.network.E2 = 10,
    r.m.components.times.x = 2,
    sw.nei = 2,
    sw.p = 0.1,
    pa.power = 1,
    pa.m = 2,
    n.shared.nodes = 3,
    n.cross.pos.edges = 0,
    n.cross.neg.edges = 8,
    cross.edge.type = c("random", "central"),
    dist.goal.activation = c(0.45, 0.45, 0.1),
    n.goal.activation = 10000,
    range.goal.activation = 0.8,
    n.Glauber = 1000,
    n.ignored.Glauber = 50
) {
  
  # Set emotion.model and cross.edge.type to the specified argument or set it to the first value ('small.world.network') as default
  emotion.model <- match.arg(emotion.model)
  cross.edge.type <- match.arg(cross.edge.type)
  network.coding <- match.arg(network.coding)
  
  supplied.args <- names(match.call())[-1]
  
  if (emotion.model != "existing.network" && !is.null(existing.network)) {
    emotion.model <- "existing.network"
    message("An existing network was provided via parameter 'existing.network', ",
            "but 'emotion.model' was not set to 'existing.network'. ",
            "The provided network will be used, and 'emotion.model' has been automatically set to 'existing.network'.")
  }
  
  # Fixed thresholds for an existing network: triggered if threshold.E1 and threshold.E2 are both
  # explicitly supplied. Then no thresholds are looked up in (or converted from) the network.
  fixed.thresholds <- FALSE
  if (emotion.model == "existing.network") {
    given <- c("threshold.E1", "threshold.E2") %in% supplied.args
    if (any(given) && !all(given)) {
      stop("For existing.network, supply either both threshold.E1 and threshold.E2 ",
           "(fixed thresholds) or neither (thresholds taken from the network).")
    }
    fixed.thresholds <- all(given)
    if (fixed.thresholds && !is.null(network.thresholds)) {
      stop("threshold.E1/threshold.E2 and network.thresholds were both supplied. Use only one.")
    }
  }
  
  # make sure weight.pos.E1.E2 is a positive and weight.neg.E1.E2 is a negative values
  if (weight.pos.E1.E2 < 0) {
    warning("weight.pos.E1.E2 must be a positive value. Changed weight.pos.E1.E2 to abs(weight.pos.E1.E2)")
  }
  weight.pos.E1.E2 <- abs(weight.pos.E1.E2)
  
  if (weight.neg.E1.E2 > 0) {
    warning("weight.neg.E1.E2 must be a negative value. Changed weight.neg.E1.E2 to -abs(weight.neg.E1.E2)")
  }
  weight.neg.E1.E2 <- -abs(weight.neg.E1.E2)
  
  ##################################################
  #Error Messages
  ##################################################
  
  # error if integer parameters are no whole integers
  if (!(n.nodes.in.network.E1 %% 1 == 0 && n.nodes.in.network.E2 %% 1 == 0)) {
    stop("n.nodes.in.network must be a whole number integer.")
  }
  
  if (!(r.m.components.times.x %% 1 == 0)) {
    stop("r.m.components.times.x must be a whole number integer.")
  }
  
  if (!(sw.nei %% 1 == 0)) {
    stop("sw.nei must be a whole number integer.")
  }
  
  if (!(n.shared.nodes %% 1 == 0)) {
    stop("n.shared.nodes must be a whole number integer.")
  }
  
  if (!(n.cross.pos.edges %% 1 == 0)) {
    stop("n.cross.pos.edges must be a whole number integer.")
  }
  
  if (!(n.cross.neg.edges %% 1 == 0)) {
    stop("n.cross.neg.edges must be a whole number integer.")
  }
  
  if (!(n.goal.activation %% 1 == 0)) {
    stop("n.goal.activation must be a whole number integer.")
  }
  
  if (!(n.Glauber %% 1 == 0)) {
    stop("n.Glauber must be a whole number integer.")
  }
  
  if (!(n.ignored.Glauber %% 1 == 0)) {
    stop("n.ignored.Glauber must be a whole number integer.")
  }
  
  # errors for exsiting network
  if (emotion.model == "existing.network") {
    # 1. Check that a network object was actually provided
    if (is.null(existing.network)) {
      stop("You selected emotion.model = 'existing.network', but no network was provided in existing.network.")
    }
    
    # 2. Check that the network is a valid format
    is.ising.list <- is.list(existing.network) && !is.data.frame(existing.network) &&
      !inherits(existing.network, c("bootnetResult", "qgraph"))
    if (!inherits(existing.network, c("bootnetResult", "qgraph")) && !is.matrix(existing.network) &&
        !is.ising.list) {
      stop("existing.network must be an ising_network (see ising_network()), an object estimated with ",
           "bootnet::estimateNetwork(), a qgraph object, or a matrix.")
    }
    if (is.ising.list && !all(c("weights", "thresholds") %in% names(existing.network))) {
      stop("A list passed as existing.network needs the elements 'weights' and 'thresholds' ",
           "(and 'coding'); see ising_network().")
    }
    
    # 3. check that variable names are provided and part of existing.network
    # Extract node labels (variable names)
    if (inherits(existing.network, "bootnetResult")) {
      # bootnet objects keep variable names in the data or colnames of $graph
      network_vars <- colnames(existing.network$graph)
      
    } else if (inherits(existing.network, "qgraph")) {
      # qgraph stores node labels in the graphAttributes
      network_vars <- existing.network$graphAttributes$Nodes$labels
      
    } else if (is.ising.list) {
      network_vars <- colnames(as.matrix(existing.network$weights))
      
    } else if (is.matrix(existing.network)) {
      # For matrices, we assume colnames correspond to variable names
      if (!is.null(colnames(existing.network))) {
        network_vars <- colnames(existing.network)
      } else {
        network_vars <- paste0("Node", seq_len(ncol(existing.network)))
        warning("existing.network is a matrix without colnames; generic node names 'Node1', 'Node2', ... were assigned.")
      }
    }
    
    if (is.null(names.E1) || length(names.E1) == 0) {
      stop("You selected emotion.model = 'existing.network', but did not provide names.E1.
         Please provide a character vector with variable names for Emotion E1.")
    }
    if (is.null(names.E2) || length(names.E2) == 0) {
      stop("You selected emotion.model = 'existing.network', but did not provide names.E2.
         Please provide a character vector with variable names for Emotion E2.")
    }
    if (!all(names.E1 %in% network_vars)) {
      stop("Some variables in names.E1 are not present in the columns of existing.network.
         Please check spelling or provide correct variable names.")
    }
    if (!all(names.E2 %in% network_vars)) {
      stop("Some variables in names.E2 are not present in the columns of existing.network.
         Please check spelling or provide correct variable names.")
    }
    
    
    # 4. check that the network is an Ising network, extract and check its thresholds,
    #    determine the coding it was estimated on, and convert weights and thresholds
    #    to the {-1,1} parameterization of the Glauber dynamics
    net.nodes <- union(names.E1, names.E2)
    
    # 4a. weights, stored thresholds, and (if available) the data the network was fitted on
    net.dat <- NULL
    net.tau <- NULL
    list.coding <- NULL
    if (inherits(existing.network, "bootnetResult")) {
      non.ising <- c("EBICglasso", "pcor", "cor", "huge", "ggmModSelect",
                     "TMFG", "LoGo", "GGMncv", "adalasso", "relimp")
      if (isTRUE(existing.network$default %in% non.ising)) {
        stop("existing.network was estimated with default = '", existing.network$default,
             "', which is not an Ising model. The Glauber dynamics require an Ising network ",
             "(e.g., bootnet default = 'IsingFit' or 'IsingSampler').")
      }
      W.existing <- existing.network$graph
      net.tau    <- existing.network$intercepts
      net.dat    <- existing.network$data
    } else if (inherits(existing.network, "qgraph")) {
      W.existing <- qgraph::getWmat(existing.network)
    } else if (is.ising.list) {
      W.existing  <- existing.network$weights
      net.tau     <- existing.network$thresholds
      list.coding <- existing.network$coding
    } else {
      W.existing <- existing.network
      net.tau    <- attr(existing.network, "intercepts")
    }
    if (!is.null(network.thresholds)) net.tau <- network.thresholds
    
    # 4b. weights: numeric, named, symmetric, zero diagonal
    W.existing <- as.matrix(W.existing)
    if (!is.numeric(W.existing)) stop("existing.network must contain numeric edge weights.")
    if (is.null(colnames(W.existing)) || is.null(rownames(W.existing)) ||
        !identical(rownames(W.existing), colnames(W.existing))) {
      stop("existing.network needs identical row and column names (the node names).")
    }
    if (!isTRUE(all.equal(unname(W.existing), unname(t(W.existing)), tolerance = 1e-8))) {
      stop("existing.network must be a symmetric weight matrix.")
    }
    if (any(abs(diag(W.existing)) > 1e-8)) stop("existing.network must have a zero diagonal.")
    diag(W.existing) <- 0
    
    # 4c. thresholds: present, numeric, named like the nodes (skipped for fixed thresholds)
    if (!fixed.thresholds) {
      if (is.null(net.tau)) {
        stop("No thresholds (intercepts) found in existing.network. Either it is not an Ising ",
             "network, or the intercepts were not stored. Supply them via ising_network() or network.thresholds.")
      }
      if (is.matrix(net.tau)) {
        if (ncol(net.tau) != 1) stop("network.thresholds must be a vector or a one-column matrix.")
        tau.names <- rownames(net.tau)
        net.tau <- as.numeric(net.tau[, 1])
        names(net.tau) <- tau.names
      }
      if (!is.numeric(net.tau) || anyNA(net.tau)) stop("Thresholds must be numeric and contain no NA.")
      if (is.null(names(net.tau))) {
        if (length(net.tau) != ncol(W.existing)) {
          stop("Unnamed thresholds must have one value per node of existing.network (",
               ncol(W.existing), "), but ", length(net.tau), " were found.")
        }
        names(net.tau) <- colnames(W.existing)
      }
      miss.tau <- setdiff(colnames(W.existing), names(net.tau))
      if (length(miss.tau) > 0) {
        stop("No threshold found for these nodes of existing.network: ",
             paste(miss.tau, collapse = ", "))
      }
      net.tau <- net.tau[colnames(W.existing)]   # match by name, never by position
    }
    
    # 4d. coding of the data the network was fitted on
    coding.used <- network.coding
    if (!is.null(list.coding)) {
      if (length(list.coding) != 1 || !list.coding %in% c("0/1", "-1/1")) {
        stop("The 'coding' element of existing.network must be \"0/1\" or \"-1/1\".")
      }
      if (coding.used == "auto") {
        coding.used <- list.coding
      } else if (coding.used != list.coding) {
        stop("network.coding = '", network.coding, "' contradicts the coding stored in existing.network ('",
             list.coding, "'). Remove one of them.")
      }
    }
    if (coding.used == "auto") {
      if (is.null(net.dat)) {
        stop("network.coding = 'auto' needs the coding stored in an ising_network or the fitted data ",
             "(bootnetResult$data). Set network.coding to '0/1' or '-1/1', or use ising_network(..., coding = ).")
      }
      net.v <- unique(as.vector(as.matrix(net.dat)))
      net.v <- net.v[!is.na(net.v)]
      coding.used <- if (all(net.v %in% c(0, 1))) "0/1"
      else if (all(net.v %in% c(-1, 1))) "-1/1"
      else stop("Cannot infer the coding: the fitted data are neither {0,1} nor {-1,1}. ",
                "Set network.coding explicitly.")
    }
    
    # 4e. convert {0,1} parameters to the {-1,1} parameterization used by the model
    # then keep only nodes assigned to E1 / E2
    if (coding.used == "0/1") {
      if (!fixed.thresholds) net.tau <- net.tau / 2 + rowSums(W.existing) / 4
      W.existing <- W.existing / 4
    }
    W.existing <- W.existing[net.nodes, net.nodes, drop = FALSE]
    
    if (fixed.thresholds) {
      # fixed thresholds are given on the {-1,1} scale and are not converted
      net.tau <- rep(NA_real_, length(net.nodes)); names(net.tau) <- net.nodes
      net.tau[setdiff(names.E1, names.E2)] <- threshold.E1
      net.tau[setdiff(names.E2, names.E1)] <- threshold.E2
      net.tau[intersect(names.E1, names.E2)] <- mean(c(threshold.E1, threshold.E2))
    } else {
      net.tau <- net.tau[net.nodes]
    }
    thresholds.existing <- matrix(net.tau, ncol = 1, dimnames = list(net.nodes, NULL))
    
  }
  
  #error if edge weights are unrealistically large/small
  if (abs(weight.E1) > 0.8) {
    stop("weight.E1 should not be larger 0.8 or smaller -0.8.")
  }
  if (abs(weight.E2) > 0.8) {
    stop("weight.E2 should not be larger 0.8 or smaller -0.8.")
  }
  
  #errors for cross emotion edge weights
  if (weight.pos.E1.E2 > 0.8) {
    stop("weight.pos.E1.E2 should not be larger 0.8.")
  }
  if (abs(weight.neg.E1.E2) > 0.8) {
    stop("abs(weight.neg.E1.E2) should not be larger 0.8.")
  }
  
  #error if thresholds are unrealistically large/small
  if (abs(threshold.E1) > 1) {
    stop("threshold.E1 should not be larger 1 or smaller -1.")
  }
  if (abs(threshold.E2) > 1) {
    stop("threshold.E2 should not be larger 1 or smaller -1.")
  }
  
  #error if beta is unrealistically large/small
  if (abs(beta) > 3) {
    stop("beta should not be larger 3 or smaller -3.")
  }
  
  #errors for dist.goal.activation
  if (length(dist.goal.activation) != 3) {
    stop("dist.goal.activation must include three values (i.e., fractions),
         one for the fraction of activations of Network 1,
         one for the fraction of activations of Network 2, and
         one for the fraction of activations of both networks simultaneously (in this order).")
  }
  if (!isTRUE(all.equal(sum(dist.goal.activation), 1))) {
    stop("The three values in dist.goal.activation must sum up to 1 (i.e., 100%)")
  }
  
  #error if n.goal.activation is too small to make sense
  if (n.goal.activation < 100) {
    stop("n.goal.activation is too small to lead to reasonable results. Set it to at least 100.")
  }
  
  #error if range.goal.activation is too small/large
  if (range.goal.activation < 0.1 | range.goal.activation > 2) {
    stop("range.goal.activation should be between 0.1 and 2.")
  }
  
  #errors for Glauber iteration arguments
  if (n.Glauber < 500) {
    stop("n.Glauber is too small to produce a stable distribution of possible states.
         Set it to at least 500.")
  }
  if (n.ignored.Glauber < 50 | n.ignored.Glauber > 300) {
    stop("n.ignored.Glauber should be between 50 and 300.
         Otherwise you exclude too few or too many iterations, respectively.")
  }
  
  # error for too many shared nodes
  if (n.shared.nodes >= n.nodes.in.network.E1 || n.shared.nodes >= n.nodes.in.network.E2) {
    stop("Too many shared nodes. n.nodes.in.network must be larger than n.shared.nodes.")
  }
  
  # error for too many cross edges
  max_cross_edges <- (n.nodes.in.network.E1 - n.shared.nodes) * (n.nodes.in.network.E2 - n.shared.nodes)
  
  if ((n.cross.pos.edges + n.cross.neg.edges) > max_cross_edges & cross.edge.type == "random") {
    stop(
      paste0("Too many shared edges. For this node configuration, at most a total of ",
             max_cross_edges,
             " cross edges is possible.")
    )
  } 
  
  if (emotion.model != "existing.network") {
    if ((n.cross.pos.edges + n.cross.neg.edges) > 0 & cross.edge.type == "central" & (n.cross.pos.edges + n.cross.neg.edges) > (min(n.nodes.in.network.E1, n.nodes.in.network.E2) - n.shared.nodes)) {
      stop(
        paste0(
          "Too many central shared edges requested. At most ",
          min(n.nodes.in.network.E1, n.nodes.in.network.E2) - n.shared.nodes,
          " central shared edges can be created given ",
          n.shared.nodes,
          " shared node(s). Please check your input."
        )
      )
    }
  }
  
  
  ##################################################
  #Warning Messages
  ##################################################
  
  if (emotion.model == "existing.network") {
    if (any(c("weight.E1", "weight.E2", "weight.pos.E1.E2", "weight.neg.E1.E2") %in% supplied.args)) {
      warning("When emotion.model = 'existing.network', weight.E1, weight.E2, weight.pos.E1.E2, and weight.neg.E1.E2 are ignored.")
    }
    
    if (!missing(n.nodes.in.network.E1) || !missing(n.nodes.in.network.E2)) {
      warning("When emotion.model = 'existing.network', n.nodes.in.network is ignored.")
    }
    
    if (!missing(n.shared.nodes)) {
      warning("When emotion.model = 'existing.network', n.shared.nodes is ignored.")
    }
    
    if (!missing(n.cross.pos.edges)) {
      warning("When emotion.model = 'existing.network', n.cross.pos.edges is ignored.")
    }
    
    if (!missing(n.cross.neg.edges)) {
      warning("When emotion.model = 'existing.network', n.cross.neg.edges is ignored.")
    }
    
    max.weight <- max(abs(W.existing[upper.tri(W.existing)]))
    if (max.weight > 0.8) {
      warning("existing.network has weights up to ", round(max.weight, 2),
              " on the {-1,1} scale, outside the range the model was designed for (|weight| <= 0.8).")
    }
    if (any(abs(thresholds.existing) > 1)) {
      warning("Thresholds of existing.network range from ",
              round(min(thresholds.existing), 2), " to ", round(max(thresholds.existing), 2),
              " on the {-1,1} scale, outside the range suggested for threshold.E1/threshold.E2 (|threshold| <= 1).")
    }
    
    vars_not_assigned <- setdiff(network_vars,
                                 union(names.E1, names.E2))
    
    if (length(vars_not_assigned) > 0) {
      warning("The following variables are not assigned to either E1 or E2 and will be ignored: ",
              paste(vars_not_assigned, collapse = ", "))
    }
  }
  
  if (!(emotion.model == "small.world.network")) {
    if (!missing(sw.nei) || !missing(sw.p)) {
      warning("When emotion.model is NOT small.world.network, sw.nei and sw.p are ignored.")
    }
  }
  
  if (!(emotion.model == "preferential.attachment.network")) {
    if (!missing(pa.power) || !missing(pa.m)) {
      warning("When emotion.model is NOT preferential.attachment.network, pa.power and pa.m are ignored.")
    }
  }
  
  if (!(emotion.model == "random.network")) {
    if (!missing(r.m.components.times.x)) {
      warning("When emotion.model is NOT random.network, r.m.components.times.x is ignored.")
    }
  }
  
  ##################################################
  #Initialize Network
  ##################################################
  
  # If existing network is passed, get network matrix W and thresholds from it
  if (emotion.model == "existing.network") {
    
    # weights and thresholds were extracted, checked, converted to the {-1,1} scale,
    # restricted to the nodes in names.E1 / names.E2 and ordered by union(names.E1, names.E2)
    # in the existing-network checks at the beginning of this function
    W <- W.existing
    thresholds <- thresholds.existing
    
    n.shared.nodes <- length(intersect(names.E1, names.E2))
  }
  
  n.nodes.unique.E1 <- n.nodes.in.network.E1 - n.shared.nodes
  n.nodes.unique.E2 <- n.nodes.in.network.E2 - n.shared.nodes
  
  if (!emotion.model == "existing.network") {
    # --- shared-node feasibility (new invariant under full-size generation) ---
    if (n.shared.nodes >= n.nodes.in.network.E1 || n.shared.nodes >= n.nodes.in.network.E2) {
      stop("n.shared.nodes (", n.shared.nodes, ") must be strictly smaller than both ",
           "n.nodes.in.network.E1 (", n.nodes.in.network.E1, ") and ",
           "n.nodes.in.network.E2 (", n.nodes.in.network.E2, "), ",
           "so each network keeps at least one unique node.")
    }
    
    # --- 'central' cross-edges: needs one unique node per edge in each network ---
    if (cross.edge.type == "central" &&
        (n.cross.pos.edges + n.cross.neg.edges) > min(n.nodes.unique.E1, n.nodes.unique.E2)) {
      stop("For cross.edge.type = 'central', the number of cross-edges (",
           n.cross.pos.edges + n.cross.neg.edges, ") cannot exceed the smaller number ",
           "of unique nodes (", min(n.nodes.unique.E1, n.nodes.unique.E2), ").")
    }
    
    # --- random network: requested edges cannot exceed the complete-graph maximum ---
    if (emotion.model == "random.network") {
      if (r.m.components.times.x > (n.nodes.in.network.E1 - 1) / 2 ||
          r.m.components.times.x > (n.nodes.in.network.E2 - 1) / 2) {
        stop("Requested random edges (n * r.m.components.times.x) exceed the maximum ",
             "possible n*(n-1)/2. Reduce r.m.components.times.x.")
      }
    }
    
    # --- small-world: ring needs enough nodes for the neighbourhood wiring ---
    if (emotion.model == "small.world.network" &&
        (n.nodes.in.network.E1 <= 2 * sw.nei || n.nodes.in.network.E2 <= 2 * sw.nei)) {
      warning("n.nodes.in.network is small relative to sw.nei; the small-world ring ",
              "may be (near-)complete or degenerate. Consider more nodes or smaller sw.nei.")
    }
  }
  
  
  ### Create networks
  
  if (emotion.model == "small.world.network") {
    network1 <- igraph::sample_smallworld(dim = 1,
                                          size = n.nodes.in.network.E1,
                                          nei = sw.nei,
                                          p = sw.p)
    
    network2 <- igraph::sample_smallworld(dim = 1,
                                          size = n.nodes.in.network.E2,
                                          nei = sw.nei,
                                          p = sw.p)
    
  } else if (emotion.model == "preferential.attachment.network") {
    network1 <- igraph::sample_pa(n = n.nodes.in.network.E1,
                                  power = pa.power,
                                  m = pa.m,
                                  directed = FALSE,
                                  out.pref = TRUE)
    
    network2 <- igraph::sample_pa(n = n.nodes.in.network.E2,
                                  power = pa.power,
                                  m = pa.m,
                                  directed = FALSE,
                                  out.pref = TRUE)
    
  } else if (emotion.model == "random.network") {
    network1 <- igraph::sample_gnm(n = n.nodes.in.network.E1,
                                   m = (n.nodes.in.network.E1 * r.m.components.times.x))
    
    network2 <- igraph::sample_gnm(n = n.nodes.in.network.E2,
                                   m = (n.nodes.in.network.E2 * r.m.components.times.x))
    
  } else if (emotion.model == "fully.connected.network") {
    network1 <- matrix(weight.E1, n.nodes.in.network.E1, n.nodes.in.network.E1)
    network1 <- igraph::graph_from_adjacency_matrix(network1, mode = 'undirected', diag = FALSE, weighted = TRUE)
    
    network2 <- matrix(weight.E2, n.nodes.in.network.E2, n.nodes.in.network.E2)
    network2 <- igraph::graph_from_adjacency_matrix(network2, mode = 'undirected', diag = FALSE, weighted = TRUE)
  }
  
  if (!emotion.model == "existing.network") {
    
    # Full-size weighted adjacencies (topology × scalar weight); keep symmetric
    A1 <- as.matrix(qgraph::getWmat(network1)); A1[A1 != 0] <- weight.E1   # n.E1 × n.E1
    A2 <- as.matrix(qgraph::getWmat(network2)); A2[A2 != 0] <- weight.E2   # n.E2 × n.E2
    
    # Total number of nodes
    total.nodes <- n.nodes.in.network.E1 + n.nodes.in.network.E2 - n.shared.nodes
    
    # Define indexes (combined-matrix layout: [uniqueE1 | shared | uniqueE2])
    non_shared_idx1 <- 1:n.nodes.unique.E1
    shared_idx <- if (n.shared.nodes > 0) (n.nodes.unique.E1 + 1):n.nodes.in.network.E1 else integer(0)
    non_shared_idx2 <- (n.nodes.in.network.E1 + 1):total.nodes
    
    idx1 <- 1:n.nodes.in.network.E1
    idx2 <- (n.nodes.unique.E1 + 1):total.nodes
    
    # Node partitions WITHIN each generated network (unique first, shared last)
    n1_u <- seq_len(n.nodes.unique.E1); n1_s <- n.nodes.unique.E1 + seq_len(n.shared.nodes)
    n2_u <- seq_len(n.nodes.unique.E2); n2_s <- n.nodes.unique.E2 + seq_len(n.shared.nodes)
    
    # --- Assemble combined matrix by block ---
    w <- matrix(0, total.nodes, total.nodes)
    
    # E1 side
    w[non_shared_idx1, non_shared_idx1] <- A1[n1_u, n1_u]   # uniqueE1 – uniqueE1
    w[non_shared_idx1, shared_idx]      <- A1[n1_u, n1_s]   # uniqueE1 – shared
    w[shared_idx,      non_shared_idx1] <- A1[n1_s, n1_u]
    
    # E2 side
    w[non_shared_idx2, non_shared_idx2] <- A2[n2_u, n2_u]   # uniqueE2 – uniqueE2
    w[non_shared_idx2, shared_idx]      <- A2[n2_u, n2_s]   # uniqueE2 – shared
    w[shared_idx,      non_shared_idx2] <- A2[n2_s, n2_u]
    
    # shared – shared: union of E1 and E2 edges, mean weight where both present
    if (n.shared.nodes > 0) {
      ss1 <- A1[n1_s, n1_s]; ss2 <- A2[n2_s, n2_s]
      ss  <- ss1 + ss2
      both <- (ss1 != 0) & (ss2 != 0)
      ss[both] <- ss[both] / 2
      w[shared_idx, shared_idx] <- ss
    }
    
    ### include shared edges (random or connecting the most central nodes in each network)
    if ((n.cross.pos.edges + n.cross.neg.edges) > 0 & cross.edge.type == 'random') {
      
      total.cross.edges <- n.cross.pos.edges + n.cross.neg.edges
      
      # Make sure there are enough 0-edges
      available_edges <- which(w[non_shared_idx1, non_shared_idx2] == 0)
      
      if (length(available_edges) < total.cross.edges) {
        stop("Not enough available edges to assign the requested number of shared edges.")
      }
      
      cross.edge.sample <- sample(available_edges, size = total.cross.edges)
      
      edge_weights <- c(rep(weight.pos.E1.E2, n.cross.pos.edges), rep(weight.neg.E1.E2, n.cross.neg.edges))
      edge_weights <- sample(edge_weights)
      
      w_block <- w[non_shared_idx1, non_shared_idx2]
      w_block[cross.edge.sample] <- edge_weights
      w[non_shared_idx1, non_shared_idx2] <- w_block
    }
    
    if ((n.cross.pos.edges + n.cross.neg.edges) > 0 & cross.edge.type == 'central') {
      
      deg1 <- lengths(igraph::as_adj_list(network1))
      deg2 <- lengths(igraph::as_adj_list(network2))
      
      deg1_subset <- deg1[non_shared_idx1]
      deg2_subset <- deg2[non_shared_idx2 - n.nodes.in.network.E1]
      
      edges1 <- cbind(deg1_subset, seq_along(deg1_subset))
      edges2 <- cbind(deg2_subset, seq_along(deg2_subset))
      
      edges1.sorted <- edges1[order(edges1[,1], decreasing = TRUE), , drop = FALSE]
      edges2.sorted <- edges2[order(edges2[,1], decreasing = TRUE), , drop = FALSE]
      
      total.cross.edges <- n.cross.pos.edges + n.cross.neg.edges
      
      edge_weights <- c(rep(weight.pos.E1.E2, n.cross.pos.edges), rep(weight.neg.E1.E2, n.cross.neg.edges))
      edge_weights <- sample(edge_weights)
      
      for (i in seq_len(total.cross.edges)) {
        
        row_index <- non_shared_idx1[ edges1.sorted[i, 2] ]
        col_index <- non_shared_idx2[ edges2.sorted[i, 2] ]
        
        w[row_index, col_index] <- edge_weights[i]
        w[col_index, row_index] <- edge_weights[i]
        
      }
    }
    
    # make weights matrix diagonal and then symmetrical
    w[lower.tri(w)] <- 0
    diag(w) <- 0
    W <- as.matrix(Matrix::forceSymmetric(w))
    
    ###initialize thresholds
    
    # create matrix with as many rows as nodes and their thresholds
    # I average the thresholds for shared nodes (if threshold.E1 and threshold.E2
    # are the same, this is equal to either threshold)
    
    thresholds <- matrix(c(
      rep(threshold.E1, n.nodes.unique.E1),                           # network 1 unique
      if(n.shared.nodes > 0) rep(mean(c(threshold.E1, threshold.E2)), n.shared.nodes) else NULL,  # shared nodes
      rep(threshold.E2, n.nodes.unique.E2)                            # network 2 unique
    ), nrow = nrow(W), ncol = 1, byrow = TRUE)
  }
  
  
  ##################################################
  #Emotion Activation
  ##################################################
  
  #create vector of kinds of goal activation
  
  vec.goal <- c(rep("E1", dist.goal.activation[1] * n.goal.activation),
                rep("E2", dist.goal.activation[2] * n.goal.activation),
                rep("E1E2", dist.goal.activation[3] * n.goal.activation))
  
  
  #create vector of intensities of goal activation from uniform distribution according to range.goal.activation
  vec.act <- runif(length(vec.goal), min = 0, max = range.goal.activation)
  
  # create indeces for E1, E2, and shared nodes
  if (emotion.model == "existing.network") {
    idx.E1 <- match(names.E1, colnames(W))
    idx.E1 <- idx.E1[!is.na(idx.E1)]
    
    idx.E2 <- match(names.E2, colnames(W))
    idx.E2 <- idx.E2[!is.na(idx.E2)]
    
    idx.shared <- intersect(idx.E1, idx.E2)
    n.shared.nodes <- length(idx.shared)
  } else {
    idx.E1 <- idx1
    idx.E2 <- idx2
    idx.shared <- shared_idx
  }
  
  ##################################################
  #Glauber dynamics, vectorised over all goal activations
  ##################################################
  # Each row of S is one independent chain (one goal activation). In each step,
  # every chain updates one randomly chosen node (random-site Glauber dynamics),
  # exactly as in the former per-activation loop. Equivalent in distribution to
  # versions <= 1.1.2, not bit-identical (different order of random draws).
  
  N <- nrow(W)                  # number of nodes
  K <- length(vec.goal)         # number of goal activations (= chains)
  rows <- seq_len(K)
  
  #thresholds per chain: initial thresholds + goal activation on the targeted nodes
  #idx.E1 and idx.E2 include shared nodes; E1E2 targets all nodes
  target.mask <- rbind(E1   = seq_len(N) %in% idx.E1,
                       E2   = seq_len(N) %in% idx.E2,
                       E1E2 = rep(TRUE, N))
  M   <- target.mask[match(vec.goal, rownames(target.mask)), , drop = FALSE]  # K x N, logical
  Tau <- matrix(as.vector(thresholds), K, N, byrow = TRUE) + M * vec.act       # K x N
  
  #starting states: all nodes off (-1) in every chain
  S <- matrix(-1, K, N)
  
  #averaging window: 50 adjacent iterations per chain, starting at a random iteration
  #after the ignored burn-in (same range as before: n.ignored.Glauber + 1 ... n.Glauber - 49)
  win.len   <- 50
  win.start <- n.ignored.Glauber +
    sample.int(n.Glauber - n.ignored.Glauber - win.len + 1, K, replace = TRUE)
  win.end   <- win.start + win.len - 1
  n.steps   <- max(win.end)     # iterations after the last window are not needed
  
  #running sum of node states within each chain's window
  acc <- matrix(0, K, N)
  
  #progress bar
  progress_bar <- utils::txtProgressBar(min = 0, max = n.steps, style = 3)
  
  for (t in seq_len(n.steps)) {
    
    #sample node that gets updated in each chain
    idx <- sample.int(N, K, replace = TRUE)
    ii  <- cbind(rows, idx)
    s.n <- S[ii]
    
    #local field of the chosen node: h_n = tau_n + sum_j W_nj s_j
    #flipping s_n changes the Hamiltonian by Delta H = 2 * s_n * h_n (zero diagonal),
    #so the flip probability equals 1 / (1 + exp(beta * Delta H)) as before
    h <- rowSums(W[idx, , drop = FALSE] * S) + Tau[ii]
    p <- 1 / (1 + exp(2 * beta * s.n * h))
    
    #flip chosen node, if probability is high enough
    flip <- runif(K) < p
    S[ii[flip, , drop = FALSE]] <- -s.n[flip]
    
    #accumulate states of chains whose window includes iteration t
    in.win <- t >= win.start & t <= win.end
    if (any(in.win)) acc[in.win, ] <- acc[in.win, ] + S[in.win, , drop = FALSE]
    
    #progress bar update (every 10 iterations)
    if (t %% 10 == 0 || t == n.steps) utils::setTxtProgressBar(progress_bar, t)
  }
  
  #close progress bar
  close(progress_bar)
  
  #final state: mean over the 50 window iterations, transformed from -1/1 to 0/1
  data <- 0.5 + 0.5 * acc / win.len
  
  #thresholds used in each goal activation
  data.thresholds <- Tau
  
  #calculate states of network
  if (n.shared.nodes > 0) {
    state.shared.nodes <- rowSums(data[, idx.shared, drop = FALSE])
  } else {
    state.shared.nodes <- rep(NA, nrow(data))
  }
  
  states.E1 <- rowSums(data[, idx.E1, drop = FALSE])
  states.E2 <- rowSums(data[, idx.E2, drop = FALSE])
  
  if (emotion.model == "existing.network") {
    node_names <- union(names.E1, names.E2)
  } else {
    node_names <- paste("node", 1:nrow(W), sep ="")
  }
  
  #combine all vectors and turn into final data frame
  data <- data.frame(cbind(vec.goal, vec.act, data, states.E1, states.E2, state.shared.nodes))
  names(data) <- c("kind", "intensity", node_names, "stateE1", "stateE2", "stateSN")
  data[c("intensity", node_names, "stateE1", "stateE2", "stateSN")] <-
    lapply(data[c("intensity", node_names, "stateE1", "stateE2", "stateSN")],
           as.numeric)
  
  #set colnames of data.thresholds
  colnames(data.thresholds) <- node_names
  
  #define returned objects
  initial.thresholds <- as.vector(thresholds)
  names(initial.thresholds) <- if (emotion.model == "existing.network") node_names else NULL
  
  out <- list(weights.matrix = W,
              initial.thresholds = initial.thresholds,
              data = data,
              data.thresholds = data.thresholds)
  if (emotion.model == "existing.network") out$network.coding <- coding.used
  return(out)
  
}