##################################################
##################################################
#Toy example (manuscript section "Toy Example"; Tables 3 and 5; Figure "Combined
#Network and Scatter Plots from the Toy Example")
#
#Three simulation cells: mixed emotions, positive blends, negative blends.
#
#Model: model/model_1_1_3.R
##################################################
##################################################

#all paths are relative to the project root: open Analysis.Rproj before running this script

source("model/model_1_1_3.R") #function for simulating the model

library(igraph)
library(qgraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(patchwork) #needed by helpers/pair_figure.R
library(png)       #needed by helpers/net_panel.R
library(grid)      #needed by helpers/net_panel.R

#shared helper functions (one file per function, also used by empirical example/simulation/simulation.R)
for (f in list.files("helpers", pattern = "\\.R$", full.names = TRUE)) source(f)

#output directory (relative to the project root)
out.dir <- "toy example/output"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

#TRUE: run the three simulations (takes a while) and save them
#FALSE: load the saved simulations from out.dir
run.sim <- TRUE

##################################################
##################################################
#Parameter settings
##################################################
##################################################

#Number of nodes per emotion
n.nodes <- 10


#Parameters fixed across cells (Table 3)
fixed.par <- list(
  emotion.model = "small.world.network",
  weight.E1 = 0.1, weight.E2 = 0.1,
  weight.pos.E1.E2 = 0.08, weight.neg.E1.E2 = -0.08,
  threshold.E1 = -0.4, threshold.E2 = -0.4,
  n.nodes.in.network.E1 = n.nodes, n.nodes.in.network.E2 = n.nodes,
  sw.nei = 2, sw.p = 0.1,
  cross.edge.type = "random",
  beta = 1,
  dist.goal.activation = c(0.45, 0.45, 0.1),
  n.goal.activation = 10000,
  range.goal.activation = 0.8,
  n.Glauber = 10000,
  n.ignored.Glauber = 50
)

#Cell-specific parameters (Table 3, scaled to 10 nodes per emotion)
#-> 5 nodes: mixed 0 shared / 8 neg; positive 1 shared / 8 pos; negative 1 shared / 2 pos
cell.par <- list(
  mixed    = list(n.shared.nodes = 0, n.cross.pos.edges = 0,n.cross.neg.edges = 16),
  positive = list(n.shared.nodes = 2, n.cross.pos.edges = 16, n.cross.neg.edges = 0),
  negative = list(n.shared.nodes = 2, n.cross.pos.edges = 4, n.cross.neg.edges = 0)
)

seed <- 4186

##################################################
##################################################
#Simulation functions
##################################################
##################################################

#simulate one cell; further arguments in ... override the settings above
simulate_toy <- function(cell, seed = 4186, ...) {
  args <- modifyList(c(fixed.par, cell.par[[cell]]), list(...))
  set.seed(seed)
  sim <- do.call(cooccurrence.model, args)

  sim$data$kind <- factor(sim$data$kind, levels = c("E1", "E2", "E1E2"))
  #MIN score: intensity of co-occurrence
  sim$data$MIN <- pmin(sim$data$stateE1, sim$data$stateE2)

  sim$cell <- cell
  sim$settings <- args
  sim
}

sim_mixed    <- function(...) simulate_toy("mixed", ...)
sim_positive <- function(...) simulate_toy("positive", ...)
sim_negative <- function(...) simulate_toy("negative", ...)

#run or load
if (run.sim) {
  toy <- list(mixed    = sim_mixed(seed = seed),
              positive = sim_positive(seed = seed),
              negative = sim_negative(seed = seed))
  saveRDS(toy, file.path(out.dir, "toy_sims_1-1-3.rds"))
} else {
  toy <- readRDS(file.path(out.dir, "toy_sims_1-1-3.rds"))
}

data.mixed    <- toy$mixed$data
data.positive <- toy$positive$data
data.negative <- toy$negative$data

#pair objects for the helpers (node names derived from the settings)
toy.p <- lapply(toy, prep_pair)

###check, whether settings are correct
lapply(toy, function(s) s$settings[c("emotion.model", "n.nodes.in.network.E1", "n.nodes.in.network.E2",
                                     "n.shared.nodes", "n.cross.pos.edges", "n.cross.neg.edges")])
lapply(toy, function(s) table(sign(s$weights.matrix[upper.tri(s$weights.matrix)])))
lapply(toy, function(s) table(s$data$kind))

##################################################
##################################################
#Simulation for mixed emotions
#a) emotions co-occur, especially in multi-goal situations
#b) mixed emotions correlate negatively
#c) the negative correlation in b) is L-shaped
##################################################
##################################################

###a) MIN is higher in E1E2 than in E1 or E2
min.mixed <- mean_by_kind(data.mixed, "MIN")
min.mixed
# A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.10  1.09 0.0163  4500
# 2 E2     2.08  1.07 0.0160  4500
# 3 E1E2   3.94  1.72 0.0543  1000


# --> A confirmed.


###b) overall negative correlation; weaker/positive within E1E2
cor.mixed <- cor_by_kind(data.mixed)
cor.mixed
# set     n      r   ci.l   ci.u            p
# 1 overall 10000 -0.474 -0.489 -0.459 0.000000e+00
# 2      E1  4500 -0.247 -0.274 -0.219 1.852084e-63
# 3      E2  4500 -0.258 -0.285 -0.231 1.748503e-69
# 4    E1E2  1000  0.285  0.228  0.341 3.286506e-20

#intensity of goal activation affects the central emotion, but not strongly the non-central emotion
data.mixed %>%
  group_by(kind) %>%
  summarize(r.int.E1 = cor(intensity, stateE1),
            r.int.E2 = cor(intensity, stateE2), .groups = "drop")

# # A tibble: 3 × 3
# kind  r.int.E1 r.int.E2
# <fct>    <dbl>    <dbl>
#   1 E1       0.741   -0.155
# 2 E2      -0.176    0.755
# 3 E1E2     0.682    0.665

# --> Activation of E1 correlates negatively with intensity of E2, and vice versa.

# --> B confirmed.


###c) L-shape (judged visually)
#two arm regressions (E2 ~ E1 in E1 cases, E1 ~ E2 in E2 cases) + negative overall correlation (dashed)
lshape.mixed <- test_lshape(toy.p$mixed, overall.line = FALSE)
# --- mixed: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u            p
# E1   stateE2 ~ stateE1  4500   -0.143 0.008 -0.160 -0.127 1.852084e-63
# E2   stateE1 ~ stateE2  4500   -0.147 0.008 -0.163 -0.131 1.748503e-69
# overall r(stateE1, stateE2) 10000   -0.474    NA -0.489 -0.459 0.000000e+00
# Overall correlation negative (r < 0, p < .05): TRUE  

lshape.mixed
lshape.res <- attr(lshape.mixed, "lshape")

ggsave(file.path(out.dir, "mixed_scatter_lshape.png"), plot = lshape.mixed,
       width = 7, height = 6, dpi = 300)

###additional implication: parts of a mixed emotion are weaker than the pure emotions
mean_by_kind(data.mixed, "stateE1")
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     5.39  2.16 0.0322  4500
# 2 E2     2.27  1.23 0.0183  4500
# 3 E1E2   4.89  1.97 0.0624  1000

mean_by_kind(data.mixed, "stateE2")
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.31  1.26 0.0187  4500
# 2 E2     5.46  2.16 0.0323  4500
# 3 E1E2   4.91  2.03 0.0641  1000

# --> C confirmed.

##################################################
##################################################
#Simulation for positive and negative blends
#a) emotions co-occur, especially in multi-goal situations
#d) emotion blends correlate positively
#e) correlation of emotion blends is higher for positive than negative emotions
#f) appraisal overlap (shared nodes) predicts emotion co-occurrence
##################################################
##################################################

###a) MIN is higher in E1E2 than in E1 or E2
min.positive <- mean_by_kind(data.positive, "MIN")
min.negative <- mean_by_kind(data.negative, "MIN")
min.positive
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.61  1.45 0.0216  4500
# 2 E2     2.65  1.47 0.0219  4500
# 3 E1E2   4.28  2.41 0.0762  1000

min.negative
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.65  1.33 0.0199  4500
# 2 E2     2.68  1.35 0.0201  4500
# 3 E1E2   4.24  2.17 0.0686  1000

# --> A confirmed.

###d) positive correlations, strongest in E1E2
cor.positive <- cor_by_kind(data.positive)
cor.negative <- cor_by_kind(data.negative)
cor.positive
# set     n     r  ci.l  ci.u             p
# 1 overall 10000 0.434 0.418 0.450  0.000000e+00
# 2      E1  4500 0.639 0.621 0.656  0.000000e+00
# 3      E2  4500 0.648 0.630 0.664  0.000000e+00
# 4    E1E2  1000 0.833 0.813 0.851 9.201891e-259

cor.negative
# set     n     r  ci.l  ci.u             p
# 1 overall 10000 0.236 0.217 0.254 2.530707e-126
# 2      E1  4500 0.498 0.476 0.520 8.069408e-281
# 3      E2  4500 0.475 0.452 0.497 3.687898e-252
# 4    E1E2  1000 0.758 0.730 0.783 1.731800e-187


#intensity of goal activation affects the central emotion and weakly the non-central emotion
lapply(list(positive = data.positive, negative = data.negative), function(d)
  d %>%
    group_by(kind) %>%
    summarize(r.int.E1 = cor(intensity, stateE1),
              r.int.E2 = cor(intensity, stateE2), .groups = "drop"))

# $positive
# # A tibble: 3 × 3
# kind  r.int.E1 r.int.E2
# <fct>    <dbl>    <dbl>
#   1 E1       0.761    0.442
# 2 E2       0.458    0.754
# 3 E1E2     0.814    0.824
# 
# $negative
# # A tibble: 3 × 3
# kind  r.int.E1 r.int.E2
# <fct>    <dbl>    <dbl>
#   1 E1       0.754    0.340
# 2 E2       0.331    0.762
# 3 E1E2     0.783    0.777

# --> D confirmed.

###e) r(positive blend) > r(negative blend)
#-> no significance test (any difference is significant at the simulated n);
#-> criterion: r.pos > r.neg, effect size Cohen's q (as in the empirical example)
r.pos <- cor(data.positive$stateE1, data.positive$stateE2)
r.neg <- cor(data.negative$stateE1, data.negative$stateE2)
q.pos.neg <- atanh(r.pos) - atanh(r.neg)
q.pos.neg # 0.2251428


###f) shared nodes predict MIN more strongly than ordinary (non-shared) nodes
#-> joint model per pair: MIN ~ shared + I(shared + ordinary); dB = B.shared - B.ordinary
#-> in the toy networks all shared nodes stand for appraisals
#one comparison node (node1 belongs to E1 only)
test_shared_slope(toy.p$positive, ordinary = "node1")
# #  shared ordinary     set B.shared B.ordinary    dB SE.dB         p.dB
# 1  node9    node1 overall    2.337      1.369 0.969 0.058 8.504064e-63
# 2 node10    node1 overall    2.368      1.334 1.034 0.057 1.840006e-72
# 3  node9    node1      E1    2.109      0.953 1.157 0.079 3.486232e-48
# 4 node10    node1      E1    2.115      0.962 1.153 0.076 1.059205e-50
# 5  node9    node1      E2    2.208      1.489 0.719 0.076 3.373930e-21
# 6 node10    node1      E2    2.223      1.470 0.753 0.077 6.617247e-23
# 7  node9    node1    E1E2    3.323      2.611 0.711 0.226 8.683771e-04
# 8 node10    node1    E1E2    3.368      2.289 1.079 0.242 4.675238e-06

# --> shared nodes predict co-occurrence more strongly than an ordinary node

test_shared_slope(toy.p$negative, ordinary = "node1")
# shared ordinary     set B.shared B.ordinary    dB SE.dB          p.dB
# 1  node9    node1 overall    2.084      1.096 0.988 0.054  3.776339e-73
# 2 node10    node1 overall    2.215      1.052 1.163 0.053 1.967178e-105
# 3  node9    node1      E1    1.902      0.692 1.210 0.074  3.456948e-59
# 4 node10    node1      E1    1.932      0.650 1.283 0.074  8.198434e-66
# 5  node9    node1      E2    1.908      1.302 0.607 0.073  4.164495e-17
# 6 node10    node1      E2    2.111      1.335 0.776 0.067  1.459115e-30
# 7  node9    node1    E1E2    3.101      1.881 1.220 0.225  3.607047e-08
# 8 node10    node1    E1E2    3.199      1.948 1.251 0.210  1.743960e-09

# --> shared nodes predict co-occurrence more strongly than an ordinary node

#all non-shared nodes as comparison
slope.positive <- test_shared_slope(toy.p$positive, comparison = "all")
# --- positive: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.805    1.0165  1.358
# node9      E1         16       16  0.694    0.9720  1.480
# node9      E2         16       16  0.596    1.0180  1.394
# node9    E1E2         16       16  0.526    0.8700  1.608
# node10 overall         16       16  0.880    1.0485  1.398
# node10      E1         16       16  0.618    1.0195  1.512
# node10      E2         16       16  0.607    1.0680  1.415
# node10    E1E2         16       16  0.773    1.1000  1.690

slope.negative <- test_shared_slope(toy.p$negative, comparison = "all")
# --- negative: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.901    1.1250  1.314
# node9      E1         16       16  0.556    0.9860  1.521
# node9      E2         16       16  0.522    1.0760  1.456
# node9    E1E2         16       16  0.797    1.2300  1.650
# node10 overall         16       16  1.073    1.3055  1.516
# node10      E1         16       16  0.557    1.0910  1.591
# node10      E2         16       16  0.776    1.2770  1.659
# node10    E1E2         16       16  0.870    1.2650  1.929


###additional implication: parts of a blend are stronger than the pure emotions
lapply(list(positive = data.positive, negative = data.negative), function(d)
  list(E1 = mean_by_kind(d, "stateE1"), E2 = mean_by_kind(d, "stateE2")))

# $positive
# $positive$E1
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     4.35  2.20 0.0328  4500
# 2 E2     2.79  1.46 0.0218  4500
# 3 E1E2   4.82  2.49 0.0786  1000
# 
# $positive$E2
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.75  1.45 0.0216  4500
# 2 E2     4.39  2.22 0.0331  4500
# 3 E1E2   4.89  2.49 0.0789  1000
# 
# 
# $negative
# $negative$E1
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     4.64  2.13 0.0317  4500
# 2 E2     2.85  1.36 0.0203  4500
# 3 E1E2   4.86  2.25 0.0712  1000
# 
# $negative$E2
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.80  1.34 0.0200  4500
# 2 E2     4.70  2.18 0.0325  4500
# 3 E1E2   4.86  2.28 0.0721  1000

# --> F confirmed.

##################################################
##################################################
#Summary table (cf. Table 5, toy example columns)
##################################################
##################################################

fmt <- function(m, se) sprintf("%.2f (%.2f)", m, se)

#F: per shared node (overall), median dB vs. all non-shared nodes, its range,
#   and the number of comparisons with dB > 0
fmt_f <- function(s) {
  o <- attr(s, "summary")
  o <- o[o$set == "overall", ]
  paste(sprintf("%s: median dB = %.2f [%.2f, %.2f], %d/%d > 0",
                o$shared, o$dB.median, o$dB.min, o$dB.max, o$n.dB.pos, o$n.ordinary),
        collapse = "; ")
}

table5.toy <- data.frame(
  phenomenon = c("A", "A", "A", "B", "C", "D", "E", "F"),
  criterion = c("Mean MIN score: E1", "Mean MIN score: E2", "Mean MIN score: E1E2",
                "Correlation of mixed emotion states",
                "Arm regressions: E2 ~ E1 (E1 cases), E1 ~ E2 (E2 cases)",
                "Correlation of same-valenced emotion states",
                "r positive vs. negative blend (Cohen's q)",
                "Linear regression MIN ~ each shared node"),
  mixed = c(fmt(min.mixed$mean, min.mixed$se),
            sprintf("r = %.3f", cor.mixed$r[cor.mixed$set == "overall"]),
            sprintf("E1 arm: B = %.2f (%.2f); E2 arm: B = %.2f (%.2f)",
                    lshape.res$estimate[1], lshape.res$se[1],
                    lshape.res$estimate[2], lshape.res$se[2]),
            "-", "-", "-"),
  positive = c(fmt(min.positive$mean, min.positive$se), "-", "-",
               sprintf("r = %.3f", cor.positive$r[cor.positive$set == "overall"]),
               sprintf("r = %.3f vs. %.3f, q = %.2f", r.pos, r.neg, q.pos.neg),
               fmt_f(slope.positive)),
  negative = c(fmt(min.negative$mean, min.negative$se), "-", "-",
               sprintf("r = %.3f", cor.negative$r[cor.negative$set == "overall"]),
               "-",
               fmt_f(slope.negative))
)
table5.toy
write.csv(table5.toy, file.path(out.dir, "table5_toy_example.csv"), row.names = FALSE)

##################################################
##################################################
#Figures
##################################################
##################################################

col.E1 <- "#E69F00"
col.E2 <- "#009E73"
kind.cols <- c(E1 = col.E1, E2 = col.E2, E1E2 = "#0072B2")

###Combined network and scatter plots (A-C networks, D-F scatter plots)
#as in Model/Plausibility Checks/1_0_1: networks exported to PNG, read in as
#ggplot panels, combined with the scatter plots via patchwork

#same edge-width scale in all three networks (qgraph otherwise scales per network)
w.max <- max(sapply(toy.p, function(p) max(abs(p$W))))

net.files <- setNames(file.path(out.dir, paste0("network_", names(toy.p), ".png")),
                      names(toy.p))

# net.mixed <- plot_network(toy.p$mixed, labels = seq_len(nrow(toy.p$mixed$W)),
#                           seed = seed, maximum = w.max, cut = 0, file = net.files[["mixed"]])
# net.positive <- plot_network(toy.p$positive, labels = seq_len(nrow(toy.p$positive$W)),
#                              seed = seed, maximum = w.max, cut = 0, file = net.files[["positive"]])
# #negative blend: same layout as positive blend (same nodes), as in the old figure
# net.negative <- plot_network(toy.p$negative, seed = seed, labels = seq_len(nrow(toy.p$negative$W)),
#                              layout = net.negative$layout, maximum = w.max, cut = 0,
#                              file = net.files[["negative"]])

net.mixed    <- plot_network(toy.p$mixed, labels = seq_len(nrow(toy.p$mixed$W)),
                             layout = cluster_layout(toy.p$mixed, seed = seed),
                             maximum = w.max, cut = 0, file = net.files[["mixed"]])
lay.blend    <- cluster_layout(toy.p$positive, seed = seed)
net.positive <- plot_network(toy.p$positive, labels = seq_len(nrow(toy.p$positive$W)),
                             layout = lay.blend, maximum = w.max, cut = 0,
                             file = net.files[["positive"]])
net.negative <- plot_network(toy.p$negative, labels = seq_len(nrow(toy.p$negative$W)),
                             layout = lay.blend, maximum = w.max, cut = 0,
                             file = net.files[["negative"]])

#scatter plots: black points, y-axis title only in the left column
scatter.mixed    <- scatter_states(toy.p$mixed,    colour = FALSE, title = NULL,
                                   xlab = "stateE1", ylab = "stateE2")
scatter.positive <- scatter_states(toy.p$positive, colour = FALSE, title = NULL,
                                   xlab = "stateE1", ylab = NULL)
scatter.negative <- scatter_states(toy.p$negative, colour = FALSE, title = NULL,
                                   xlab = "stateE1", ylab = NULL)

combined.net.scatter <-
  net_panel(net.files[["mixed"]]) + net_panel(net.files[["positive"]]) +
  net_panel(net.files[["negative"]]) +
  scatter.mixed + scatter.positive + scatter.negative +
  plot_layout(nrow = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 16)))

print(combined.net.scatter)
ggsave(file.path(out.dir, "combined_network_scatter_ABCDEF.png"),
       plot = combined.net.scatter, width = 15, height = 10, dpi = 300)

###MIN histograms by kind of goal activation
#coord_cartesian instead of scale limits, so that no bar is dropped at the edges
min.long <- bind_rows(lapply(names(toy), function(cell)
  data.frame(cell = cell, kind = toy[[cell]]$data$kind, MIN = toy[[cell]]$data$MIN)))
min.long$cell <- factor(min.long$cell, levels = c("mixed", "positive", "negative"),
                        labels = c("A: Mixed emotions", "B: Positive blends", "C: Negative blends"))

min.hist <- ggplot(min.long, aes(x = MIN, fill = kind)) +
  geom_histogram(position = "dodge", binwidth = 1, boundary = 0, colour = NA) +
  facet_wrap(~ cell, scales = "free_y") +
  scale_fill_manual(values = kind.cols) +
  scale_x_continuous(breaks = 0:n.nodes, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, n.nodes)) +
  labs(x = "MIN score", y = "Frequency", fill = "Activation") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())
min.hist
ggsave(file.path(out.dir, "combined_MIN_hist_ABC.png"), plot = min.hist,
       width = 14, height = 4.5, dpi = 300)
