##################################################
##################################################
#Sensitivity analyses for the toy example (preregistration, section "Sensitivity Analyses")
#
#a) number of shared nodes equal to the number of shared appraisals in Yeo and Ong (2024)
#b) relative frequency of activating situations (dist.goal.activation)
#c) cross-edge type "central" instead of "random"

#Model: model/model_1_1_3.R
##################################################
##################################################

#all paths are relative to the project root: open Analysis.Rproj before running this script

source("model/model_1_1_3.R")

library(ggplot2)
library(dplyr)
library(patchwork)
library(igraph)
library(qgraph)
library(png)
library(grid)

for (f in list.files("helpers", pattern = "\\.R$", full.names = TRUE)) source(f)

out.dir <- "supplementary analyses/output"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

#TRUE: run the simulations and save them; FALSE: load them from out.dir
run.sim <- FALSE

##################################################
##################################################
#Parameter settings
##################################################
##################################################

n.nodes <- 10

#main simulation (as in toy_example.R)
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

cell.par <- list(
  mixed    = list(n.shared.nodes = 0, n.cross.pos.edges = 0,  n.cross.neg.edges = 16),
  positive = list(n.shared.nodes = 2, n.cross.pos.edges = 16, n.cross.neg.edges = 0),
  negative = list(n.shared.nodes = 2, n.cross.pos.edges = 4,  n.cross.neg.edges = 0)
)

seed <- 4186

#sensitivity conditions: parameters that replace the main settings, per cell
#c) the model allows at most one central cross-edge per unique node
#   (mixed: 10, positive: 8).
sens.par <- list(
  a_shared = list(
    mixed    = list(n.shared.nodes = 1),
    positive = list(n.shared.nodes = 2),
    negative = list(n.shared.nodes = 2)),
  b_both20 = list(
    mixed    = list(dist.goal.activation = c(0.40, 0.40, 0.20)),
    positive = list(dist.goal.activation = c(0.40, 0.40, 0.20)),
    negative = list(dist.goal.activation = c(0.40, 0.40, 0.20))),
  b_both05 = list(
    mixed    = list(dist.goal.activation = c(0.475, 0.475, 0.05)),
    positive = list(dist.goal.activation = c(0.475, 0.475, 0.05)),
    negative = list(dist.goal.activation = c(0.475, 0.475, 0.05))),
  c_central = list(
    mixed    = list(cross.edge.type = "central", n.cross.neg.edges = 10),
    positive = list(cross.edge.type = "central", n.cross.pos.edges = 8),
    negative = list(cross.edge.type = "central")
  )
)

##################################################
##################################################
#Simulation functions
##################################################
##################################################

simulate_toy <- function(cell, condition = NULL, seed = 4186, ...) {
  args <- modifyList(c(fixed.par, cell.par[[cell]]),
                     if (is.null(condition)) list() else sens.par[[condition]][[cell]])
  args <- modifyList(args, list(...))
  set.seed(seed)
  sim <- do.call(cooccurrence.model, args)
  sim$cell <- cell
  sim$condition <- if (is.null(condition)) "main" else condition
  sim$settings <- args
  sim
}

cells      <- names(cell.par)
conditions <- names(sens.par)

if (run.sim) {
  sens <- lapply(setNames(names(sens.par), names(sens.par)), function(cond)
    lapply(setNames(cells, cells), function(cell) simulate_toy(cell, cond, seed = seed)))
  saveRDS(sens, file.path(out.dir, "toy_sensitivity_sims_1-1-3.rds"))
} else {
  sens <- readRDS(file.path(out.dir, "toy_sensitivity_sims_1-1-3.rds"))
}

#main simulation from toy_example.R as reference
sens <- c(list(main = readRDS("toy example/output/toy_sims_1-1-3.rds")), sens)
conditions <- names(sens)

sens.p <- lapply(sens, function(cond) lapply(cond, prep_pair))

###check, whether settings are correct
lapply(sens, function(cond) sapply(cond, function(s)
  unlist(s$settings[c("n.shared.nodes", "n.cross.pos.edges", "n.cross.neg.edges",
                      "weight.pos.E1.E2", "weight.neg.E1.E2", "cross.edge.type")])))
lapply(sens, function(cond) sapply(cond, function(s) table(s$data$kind)))

##################################################
##################################################
#Simulation for mixed emotions
#a) emotions co-occur, especially in multi-goal situations
#b) mixed emotions correlate negatively
#c) the negative correlation in b) is L-shaped
##################################################
##################################################

###a) MIN is higher in E1E2 than in E1 or E2
min.mixed <- bind_rows(lapply(conditions, function(cond)
  mutate(mean_by_kind(sens.p[[cond]]$mixed$data, "MIN"), condition = cond)))
min.mixed

# # A tibble: 15 × 6
# kind   mean    sd     se     n condition
# <fct> <dbl> <dbl>  <dbl> <int> <chr>    
#   1 E1     2.10  1.09 0.0163  4500 main     
# 2 E2     2.08  1.07 0.0160  4500 main     
# 3 E1E2   3.94  1.72 0.0543  1000 main     
# 4 E1     2.36  1.13 0.0168  4500 a_shared 
# 5 E2     2.38  1.14 0.0171  4500 a_shared 
# 6 E1E2   3.94  1.87 0.0592  1000 a_shared 
# 7 E1     2.11  1.09 0.0173  4000 b_both20 
# 8 E2     2.07  1.07 0.0169  4000 b_both20 
# 9 E1E2   3.96  1.78 0.0399  2000 b_both20 
# 10 E1     2.10  1.09 0.0158  4750 b_both05 
# 11 E2     2.08  1.06 0.0154  4750 b_both05 
# 12 E1E2   4.01  1.76 0.0788   500 b_both05 
# 13 E1     2.12  1.07 0.0160  4500 c_central
# 14 E2     2.13  1.09 0.0162  4500 c_central
# 15 E1E2   3.96  1.92 0.0608  1000 c_central

###b) overall negative correlation
cor.mixed <- bind_rows(lapply(conditions, function(cond)
  mutate(cor_by_kind(sens.p[[cond]]$mixed$data), condition = cond)))
cor.mixed
# set     n      r   ci.l   ci.u             p condition
# 1  overall 10000 -0.474 -0.489 -0.459  0.000000e+00      main
# 2       E1  4500 -0.247 -0.274 -0.219  1.852084e-63      main
# 3       E2  4500 -0.258 -0.285 -0.231  1.748503e-69      main
# 4     E1E2  1000  0.285  0.228  0.341  3.286506e-20      main
# 5  overall 10000 -0.275 -0.293 -0.256 2.324264e-172  a_shared
# 6       E1  4500  0.015 -0.014  0.044  3.209184e-01  a_shared
# 7       E2  4500 -0.007 -0.036  0.023  6.561963e-01  a_shared
# 8     E1E2  1000  0.448  0.397  0.496  1.833304e-50  a_shared
# 9  overall 10000 -0.381 -0.397 -0.364  0.000000e+00  b_both20
# 10      E1  4000 -0.247 -0.275 -0.217  1.879524e-56  b_both20
# 11      E2  4000 -0.248 -0.276 -0.218  6.292813e-57  b_both20
# 12    E1E2  2000  0.316  0.276  0.355  1.668039e-47  b_both20
# 13 overall 10000 -0.516 -0.530 -0.501  0.000000e+00  b_both05
# 14      E1  4750 -0.246 -0.273 -0.219  1.701840e-66  b_both05
# 15      E2  4750 -0.258 -0.284 -0.231  4.194286e-73  b_both05
# 16    E1E2   500  0.299  0.217  0.377  8.966815e-12  b_both05
# 17 overall 10000 -0.393 -0.409 -0.376  0.000000e+00 c_central
# 18      E1  4500 -0.181 -0.209 -0.153  1.445107e-34 c_central
# 19      E2  4500 -0.158 -0.187 -0.130  1.357711e-26 c_central
# 20    E1E2  1000  0.453  0.403  0.501  7.028387e-52 c_central


###c) L-shape (judged visually)
lshape.mixed <- lapply(setNames(conditions, conditions), function(cond)
  test_lshape(sens.p[[cond]]$mixed, overall.line = FALSE))
# --- mixed: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u            p
# E1   stateE2 ~ stateE1  4500   -0.143 0.008 -0.160 -0.127 1.852084e-63
# E2   stateE1 ~ stateE2  4500   -0.147 0.008 -0.163 -0.131 1.748503e-69
# overall r(stateE1, stateE2) 10000   -0.474    NA -0.489 -0.459 0.000000e+00
# Overall correlation negative (r < 0, p < .05): TRUE 
# --- mixed: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u             p
# E1   stateE2 ~ stateE1  4500    0.008 0.009 -0.008  0.025  3.209184e-01
# E2   stateE1 ~ stateE2  4500   -0.004 0.009 -0.021  0.013  6.561963e-01
# overall r(stateE1, stateE2) 10000   -0.275    NA -0.293 -0.256 2.324264e-172
# Overall correlation negative (r < 0, p < .05): TRUE 
# --- mixed: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u            p
# E1   stateE2 ~ stateE1  4000   -0.143 0.009 -0.161 -0.126 1.879524e-56
# E2   stateE1 ~ stateE2  4000   -0.139 0.009 -0.156 -0.122 6.292813e-57
# overall r(stateE1, stateE2) 10000   -0.381    NA -0.397 -0.364 0.000000e+00
# Overall correlation negative (r < 0, p < .05): TRUE 
# --- mixed: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u            p
# E1   stateE2 ~ stateE1  4750   -0.143 0.008 -0.159 -0.127 1.701840e-66
# E2   stateE1 ~ stateE2  4750   -0.146 0.008 -0.161 -0.130 4.194286e-73
# overall r(stateE1, stateE2) 10000   -0.516    NA -0.530 -0.501 0.000000e+00
# Overall correlation negative (r < 0, p < .05): TRUE 
# --- mixed: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u            p
# E1   stateE2 ~ stateE1  4500   -0.101 0.008 -0.117 -0.085 1.445107e-34
# E2   stateE1 ~ stateE2  4500   -0.091 0.008 -0.108 -0.074 1.357711e-26
# overall r(stateE1, stateE2) 10000   -0.393    NA -0.409 -0.376 0.000000e+00
# Overall correlation negative (r < 0, p < .05): TRUE 

lshape.res <- bind_rows(lapply(conditions, function(cond)
  mutate(attr(lshape.mixed[[cond]], "lshape"), condition = cond)))

lshape.res
# set               model     n estimate    se   ci.l   ci.u             p condition
# 1       E1   stateE2 ~ stateE1  4500   -0.143 0.008 -0.160 -0.127  1.852084e-63      main
# 2       E2   stateE1 ~ stateE2  4500   -0.147 0.008 -0.163 -0.131  1.748503e-69      main
# 3  overall r(stateE1, stateE2) 10000   -0.474    NA -0.489 -0.459  0.000000e+00      main
# 4       E1   stateE2 ~ stateE1  4500    0.008 0.009 -0.008  0.025  3.209184e-01  a_shared
# 5       E2   stateE1 ~ stateE2  4500   -0.004 0.009 -0.021  0.013  6.561963e-01  a_shared
# 6  overall r(stateE1, stateE2) 10000   -0.275    NA -0.293 -0.256 2.324264e-172  a_shared
# 7       E1   stateE2 ~ stateE1  4000   -0.143 0.009 -0.161 -0.126  1.879524e-56  b_both20
# 8       E2   stateE1 ~ stateE2  4000   -0.139 0.009 -0.156 -0.122  6.292813e-57  b_both20
# 9  overall r(stateE1, stateE2) 10000   -0.381    NA -0.397 -0.364  0.000000e+00  b_both20
# 10      E1   stateE2 ~ stateE1  4750   -0.143 0.008 -0.159 -0.127  1.701840e-66  b_both05
# 11      E2   stateE1 ~ stateE2  4750   -0.146 0.008 -0.161 -0.130  4.194286e-73  b_both05
# 12 overall r(stateE1, stateE2) 10000   -0.516    NA -0.530 -0.501  0.000000e+00  b_both05
# 13      E1   stateE2 ~ stateE1  4500   -0.101 0.008 -0.117 -0.085  1.445107e-34 c_central
# 14      E2   stateE1 ~ stateE2  4500   -0.091 0.008 -0.108 -0.074  1.357711e-26 c_central
# 15 overall r(stateE1, stateE2) 10000   -0.393    NA -0.409 -0.376  0.000000e+00 c_central

###f) only in condition a_shared (one shared node)
slope.mixed.a <- test_shared_slope(sens.p$a_shared$mixed, comparison = "all")
# --- mixed: shared vs. all 18 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node10 overall         18       18  0.963    1.2145  1.363
# node10      E1         18       18  0.224    0.9025  1.780
# node10      E2         18       18  0.305    1.0130  2.000
# node10    E1E2         18       18  1.287    1.8155  2.093

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
min.blends <- bind_rows(lapply(conditions, function(cond) bind_rows(
  mutate(mean_by_kind(sens.p[[cond]]$positive$data, "MIN"), cell = "positive", condition = cond),
  mutate(mean_by_kind(sens.p[[cond]]$negative$data, "MIN"), cell = "negative", condition = cond))))
print(min.blends, n = Inf)
# A tibble: 30 × 7
# kind   mean    sd     se     n cell     condition
# <fct> <dbl> <dbl>  <dbl> <int> <chr>    <chr>    
#   1 E1     2.61  1.45 0.0216  4500 positive main     
# 2 E2     2.65  1.47 0.0219  4500 positive main     
# 3 E1E2   4.28  2.41 0.0762  1000 positive main     
# 4 E1     2.65  1.33 0.0199  4500 negative main     
# 5 E2     2.68  1.35 0.0201  4500 negative main     
# 6 E1E2   4.24  2.17 0.0686  1000 negative main     
# 7 E1     2.61  1.45 0.0216  4500 positive a_shared 
# 8 E2     2.65  1.47 0.0219  4500 positive a_shared 
# 9 E1E2   4.28  2.41 0.0762  1000 positive a_shared 
# 10 E1     2.65  1.33 0.0199  4500 negative a_shared 
# 11 E2     2.68  1.35 0.0201  4500 negative a_shared 
# 12 E1E2   4.24  2.17 0.0686  1000 negative a_shared 
# 13 E1     2.60  1.44 0.0228  4000 positive b_both20 
# 14 E2     2.66  1.48 0.0233  4000 positive b_both20 
# 15 E1E2   4.30  2.37 0.0529  2000 positive b_both20 
# 16 E1     2.65  1.33 0.0210  4000 negative b_both20 
# 17 E2     2.67  1.33 0.0211  4000 negative b_both20 
# 18 E1E2   4.30  2.17 0.0485  2000 negative b_both20 
# 19 E1     2.61  1.45 0.0211  4750 positive b_both05 
# 20 E2     2.65  1.48 0.0215  4750 positive b_both05 
# 21 E1E2   4.38  2.43 0.109    500 positive b_both05 
# 22 E1     2.65  1.33 0.0193  4750 negative b_both05 
# 23 E2     2.68  1.35 0.0196  4750 negative b_both05 
# 24 E1E2   4.15  2.18 0.0975   500 negative b_both05 
# 25 E1     2.66  1.35 0.0202  4500 positive c_central
# 26 E2     2.67  1.39 0.0207  4500 positive c_central
# 27 E1E2   4.32  2.28 0.0722  1000 positive c_central
# 28 E1     2.63  1.33 0.0199  4500 negative c_central
# 29 E2     2.70  1.39 0.0207  4500 negative c_central
# 30 E1E2   4.21  2.21 0.0700  1000 negative c_central

###d) positive correlations
cor.blends <- bind_rows(lapply(conditions, function(cond) bind_rows(
  mutate(cor_by_kind(sens.p[[cond]]$positive$data), cell = "positive", condition = cond),
  mutate(cor_by_kind(sens.p[[cond]]$negative$data), cell = "negative", condition = cond))))
cor.blends
# set     n     r  ci.l  ci.u             p     cell condition
# 1  overall 10000 0.434 0.418 0.450  0.000000e+00 positive      main
# 2       E1  4500 0.639 0.621 0.656  0.000000e+00 positive      main
# 3       E2  4500 0.648 0.630 0.664  0.000000e+00 positive      main
# 4     E1E2  1000 0.833 0.813 0.851 9.201891e-259 positive      main
# 5  overall 10000 0.236 0.217 0.254 2.530707e-126 negative      main
# 6       E1  4500 0.498 0.476 0.520 8.069408e-281 negative      main
# 7       E2  4500 0.475 0.452 0.497 3.687898e-252 negative      main
# 8     E1E2  1000 0.758 0.730 0.783 1.731800e-187 negative      main
# 9  overall 10000 0.434 0.418 0.450  0.000000e+00 positive  a_shared
# 10      E1  4500 0.639 0.621 0.656  0.000000e+00 positive  a_shared
# 11      E2  4500 0.648 0.630 0.664  0.000000e+00 positive  a_shared
# 12    E1E2  1000 0.833 0.813 0.851 9.201891e-259 positive  a_shared
# 13 overall 10000 0.236 0.217 0.254 2.530707e-126 negative  a_shared
# 14      E1  4500 0.498 0.476 0.520 8.069408e-281 negative  a_shared
# 15      E2  4500 0.475 0.452 0.497 3.687898e-252 negative  a_shared
# 16    E1E2  1000 0.758 0.730 0.783 1.731800e-187 negative  a_shared
# 17 overall 10000 0.504 0.489 0.518  0.000000e+00 positive  b_both20
# 18      E1  4000 0.636 0.617 0.654  0.000000e+00 positive  b_both20
# 19      E2  4000 0.655 0.637 0.673  0.000000e+00 positive  b_both20
# 20    E1E2  2000 0.832 0.818 0.845  0.000000e+00 positive  b_both20
# 21 overall 10000 0.319 0.301 0.336 9.826391e-235 negative  b_both20
# 22      E1  4000 0.495 0.471 0.518 1.853190e-246 negative  b_both20
# 23      E2  4000 0.477 0.453 0.501 9.082602e-227 negative  b_both20
# 24    E1E2  2000 0.751 0.731 0.769  0.000000e+00 negative  b_both20
# 25 overall 10000 0.400 0.383 0.416  0.000000e+00 positive  b_both05
# 26      E1  4750 0.638 0.620 0.654  0.000000e+00 positive  b_both05
# 27      E2  4750 0.652 0.636 0.668  0.000000e+00 positive  b_both05
# 28    E1E2   500 0.838 0.809 0.862 8.186268e-133 positive  b_both05
# 29 overall 10000 0.191 0.172 0.209  1.624771e-82 negative  b_both05
# 30      E1  4750 0.501 0.479 0.522 4.949675e-300 negative  b_both05
# 31      E2  4750 0.478 0.456 0.500 6.806354e-270 negative  b_both05
# 32    E1E2   500 0.741 0.698 0.778  5.327195e-88 negative  b_both05
# 33 overall 10000 0.311 0.293 0.329 3.374129e-223 positive c_central
# 34      E1  4500 0.552 0.531 0.572  0.000000e+00 positive c_central
# 35      E2  4500 0.542 0.521 0.562  0.000000e+00 positive c_central
# 36    E1E2  1000 0.805 0.782 0.826 1.092727e-228 positive c_central
# 37 overall 10000 0.268 0.250 0.286 1.626021e-164 negative c_central
# 38      E1  4500 0.501 0.479 0.522 1.870491e-284 negative c_central
# 39      E2  4500 0.532 0.510 0.552  0.000000e+00 negative c_central
# 40    E1E2  1000 0.764 0.736 0.788 9.097700e-192 negative c_central

###e) r(positive blend) > r(negative blend), Cohen's q
e.res <- bind_rows(lapply(conditions, function(cond) {
  r.pos <- cor(sens.p[[cond]]$positive$data$stateE1, sens.p[[cond]]$positive$data$stateE2)
  r.neg <- cor(sens.p[[cond]]$negative$data$stateE1, sens.p[[cond]]$negative$data$stateE2)
  data.frame(condition = cond, r.pos = r.pos, r.neg = r.neg, q = atanh(r.pos) - atanh(r.neg))
}))
e.res
# condition     r.pos     r.neg          q
# 1      main 0.4344696 0.2357322 0.22514284
# 2  a_shared 0.4344696 0.2357322 0.22514284
# 3  b_both20 0.5038777 0.3186185 0.22438110
# 4  b_both05 0.3995856 0.1906829 0.23010998
# 5 c_central 0.3110230 0.2683552 0.04658705

###f) shared nodes predict MIN more strongly than all non-shared nodes
slope.blends <- lapply(setNames(conditions, conditions), function(cond)
  list(positive = test_shared_slope(sens.p[[cond]]$positive, comparison = "all"),
       negative = test_shared_slope(sens.p[[cond]]$negative, comparison = "all")))
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
# --- positive: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.683    0.8850  1.254
# node9      E1         16       16  0.654    0.9420  1.443
# node9      E2         16       16  0.530    1.0105  1.343
# node9    E1E2         16       16  0.579    0.8520  1.408
# node10 overall         16       16  0.760    0.9805  1.370
# node10      E1         16       16  0.651    1.0180  1.529
# node10      E2         16       16  0.585    1.0895  1.432
# node10    E1E2         16       16  0.698    0.9645  1.524
# --- negative: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.801    1.0640  1.269
# node9      E1         16       16  0.576    0.9835  1.505
# node9      E2         16       16  0.472    1.0160  1.329
# node9    E1E2         16       16  0.887    1.2020  1.643
# node10 overall         16       16  0.930    1.1855  1.423
# node10      E1         16       16  0.537    1.0505  1.522
# node10      E2         16       16  0.752    1.2790  1.662
# node10    E1E2         16       16  0.758    1.1325  1.635
# --- positive: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.867    1.0765  1.385
# node9      E1         16       16  0.677    0.9965  1.475
# node9      E2         16       16  0.613    1.0345  1.433
# node9    E1E2         16       16  0.086    0.6855  1.351
# node10 overall         16       16  0.934    1.1025  1.421
# node10      E1         16       16  0.630    1.0390  1.499
# node10      E2         16       16  0.616    1.0575  1.432
# node10    E1E2         16       16  0.641    1.1850  1.753
# --- negative: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.980    1.1725  1.358
# node9      E1         16       16  0.537    0.9735  1.505
# node9      E2         16       16  0.548    1.0975  1.461
# node9    E1E2         16       16  1.038    1.2735  2.114
# node10 overall         16       16  1.073    1.3345  1.527
# node10      E1         16       16  0.557    1.0835  1.581
# node10      E2         16       16  0.806    1.2805  1.681
# node10    E1E2         16       16  0.732    1.2390  2.339
# --- positive: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.964    1.1180  1.250
# node9      E1         16       16  0.608    0.9600  1.451
# node9      E2         16       16  0.668    1.0795  1.395
# node9    E1E2         16       16  1.071    1.4235  1.759
# node10 overall         16       16  0.988    1.1785  1.356
# node10      E1         16       16  0.487    0.9770  1.486
# node10      E2         16       16  0.777    1.1960  1.566
# node10    E1E2         16       16  0.926    1.3555  1.765
# --- negative: shared vs. all 16 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# node9 overall         16       16  0.825    1.1430  1.379
# node9      E1         16       16  0.558    1.0150  1.486
# node9      E2         16       16  0.545    1.1245  1.526
# node9    E1E2         16       16  0.589    1.0380  1.561
# node10 overall         16       16  0.995    1.3310  1.584
# node10      E1         16       16  0.514    1.1235  1.580
# node10      E2         16       16  0.749    1.2700  1.774
# node10    E1E2         16       16  1.158    1.7450  2.232

##################################################
##################################################
#Summary table (cf. Table 5, toy example columns)
##################################################
##################################################

fmt <- function(m, sd) sprintf("%.2f (%.2f)", m, sd)

fmt_f <- function(s) {
  o <- attr(s, "summary")
  o <- o[o$set == "overall", ]
  paste(sprintf("%s: median dB = %.2f [%.2f, %.2f], %d/%d > 0",
                o$shared, o$dB.median, o$dB.min, o$dB.max, o$n.dB.pos, o$n.ordinary),
        collapse = "; ")
}

fmt_min <- function(d) {
  m <- mean_by_kind(d, "MIN")
  paste(sprintf("%s: %s", m$kind, fmt(m$mean, m$sd)), collapse = "; ")
}

table.sens <- bind_rows(lapply(conditions, function(cond) {
  p  <- sens.p[[cond]]
  ls <- attr(lshape.mixed[[cond]], "lshape")
  e  <- e.res[e.res$condition == cond, ]
  data.frame(
    condition = cond,
    phenomenon = c("A", "B", "C", "D", "E", "F"),
    mixed = c(fmt_min(p$mixed$data),
              sprintf("r = %.3f", cor(p$mixed$data$stateE1, p$mixed$data$stateE2)),
              sprintf("E1 arm: B = %.2f; E2 arm: B = %.2f", ls$estimate[1], ls$estimate[2]),
              "-", "-",
              if (cond == "a_shared") fmt_f(slope.mixed.a) else "-"),
    positive = c(fmt_min(p$positive$data), "-", "-",
                 sprintf("r = %.3f", e$r.pos),
                 sprintf("r = %.3f vs. %.3f, q = %.2f", e$r.pos, e$r.neg, e$q),
                 fmt_f(slope.blends[[cond]]$positive)),
    negative = c(fmt_min(p$negative$data), "-", "-",
                 sprintf("r = %.3f", e$r.neg), "-",
                 fmt_f(slope.blends[[cond]]$negative))
  )
}))
table.sens
write.csv(table.sens, file.path(out.dir, "table_sensitivity_toy_example.csv"), row.names = FALSE)

##################################################
##################################################
#Figures
##################################################
##################################################

# A: one shared node for mixed emotions
lshape.mixed.a <- test_lshape(sens.p$a_shared$mixed, overall.line = FALSE)
lshape.mixed.a
ggsave(file.path(out.dir, "mixed_scatter_lshape_a_shared.png"), plot = lshape.mixed.a,
       width = 7, height = 6, dpi = 300)

# B: distribution of goal activations
for (cond in c("b_both20", "b_both05")) {
  ggsave(file.path(out.dir, paste0("mixed_scatter_lshape_", cond, ".png")),
         plot = lshape.mixed[[cond]], width = 7, height = 6, dpi = 300)
}

# C: cross edge type 'central'
ggsave(file.path(out.dir, "mixed_scatter_lshape_c_central.png"),
       plot = lshape.mixed$c_central, width = 7, height = 6, dpi = 300)

# networks and scatter plots
p.c <- sens.p$c_central

#same edge-width scale in all three networks
w.max <- max(sapply(p.c, function(p) max(abs(p$W))))

net.files.c <- setNames(file.path(out.dir, paste0("network_c_central_", names(p.c), ".png")),
                        names(p.c))

plot_network(p.c$mixed, labels = seq_len(nrow(p.c$mixed$W)),
             layout = cluster_layout(p.c$mixed, seed = seed),
             maximum = w.max, cut = 0, file = net.files.c[["mixed"]])
lay.blend.c <- cluster_layout(p.c$positive, seed = seed)
plot_network(p.c$positive, labels = seq_len(nrow(p.c$positive$W)),
             layout = lay.blend.c, maximum = w.max, cut = 0,
             file = net.files.c[["positive"]])
plot_network(p.c$negative, labels = seq_len(nrow(p.c$negative$W)),
             layout = lay.blend.c, maximum = w.max, cut = 0,
             file = net.files.c[["negative"]])

combined.c <-
  net_panel(net.files.c[["mixed"]]) + net_panel(net.files.c[["positive"]]) +
  net_panel(net.files.c[["negative"]]) +
  scatter_states(p.c$mixed,    colour = FALSE, title = NULL, xlab = "stateE1", ylab = "stateE2") +
  scatter_states(p.c$positive, colour = FALSE, title = NULL, xlab = "stateE1", ylab = NULL) +
  scatter_states(p.c$negative, colour = FALSE, title = NULL, xlab = "stateE1", ylab = NULL) +
  plot_layout(nrow = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 16)))

print(combined.c)
ggsave(file.path(out.dir, "combined_network_scatter_c_central_ABCDEF.png"),
       plot = combined.c, width = 15, height = 10, dpi = 300)
