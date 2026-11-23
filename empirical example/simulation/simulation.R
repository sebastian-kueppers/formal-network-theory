##################################################
##################################################
#Empirical example: simulate data from the estimated Ising networks and test
#phenomena A-F (same procedure as toy example/toy_example.R)
#
#Three emotion pairs (top co-occurring pairs in Study 2):
#  love / sadness            -> mixed emotions   (A, B, C; F for the shared appraisal)
#  enjoyment / joy           -> positive blend   (A, D, E, F)
#  anxiety / frustration     -> negative blend   (A, D, E, F)
#
#Networks: empirical example/top_emotion_networks_ising.RDS
#          (Analysis_final.R, section 1.2; run it on the anonymized data first)
#Model: model/model_1_1_3.R
##################################################
##################################################

#all paths are relative to the project root: open Analysis.Rproj before running this script
emp.dir <- "empirical example"

source("model/model_1_1_3.R") #function for simulating the model

library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(patchwork) #needed by helpers/pair_figure.R
library(qgraph)
library(png)       #needed by helpers/net_panel.R
library(grid)      #needed by helpers/net_panel.R

#shared helper functions (one file per function, also used by toy example/toy_example.R)
for (f in list.files("helpers", pattern = "\\.R$", full.names = TRUE)) source(f)

#output directories (relative to the project root)
out.dir     <- file.path(emp.dir, "simulation/results")
out.dir.fig <- file.path(emp.dir, "simulation/figures")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out.dir.fig, recursive = TRUE, showWarnings = FALSE)

#TRUE: run the three simulations (takes a while) and save them
#FALSE: load the saved simulations from out.dir
run.sim  <- TRUE
sim.file <- file.path(out.dir, "cooccurrence_sims_threshold-0-4_m1-1-3.rds")

##################################################
##################################################
#Parameter settings
##################################################
##################################################

ising_nets <- readRDS(file.path(emp.dir, "top_emotion_networks_ising.RDS"))

#node threshold (fixed for all pairs; the estimated intercepts are not used)
THRESHOLD <- -0.4

#Parameters fixed across pairs
fixed.par <- list(
  emotion.model = "existing.network",
  network.coding = "0/1",
  threshold.E1 = THRESHOLD, threshold.E2 = THRESHOLD,
  beta = 1,
  dist.goal.activation = c(0.45, 0.45, 0.1),
  n.Glauber = 10000,
  n.goal.activation = 10000
)

#one seed per pair (as in the original simulations)
seeds  <- c(love_sadness = 1, enjoyment_joy = 2, anxiety_frustration = 3)
labels <- c(love_sadness = "love / sadness", enjoyment_joy = "enjoyment / joy",
            anxiety_frustration = "anxiety / frustration")

##################################################
##################################################
#Simulation
##################################################
##################################################

#simulate one pair from its estimated network; further arguments in ... override the settings above
simulate_pair <- function(pair, seed = seeds[[pair]], ...) {
  net  <- ising_nets[[pair]]
  args <- modifyList(c(fixed.par,
                       list(existing.network = net$net.custom.target3.2,
                            names.E1 = net$names.E1, names.E2 = net$names.E2)),
                     list(...))
  set.seed(seed)
  do.call(cooccurrence.model, args)
}

#run or load
if (run.sim) {
  sims <- lapply(setNames(names(seeds), names(seeds)), simulate_pair)
  saveRDS(sims, sim.file)
} else {
  sims <- readRDS(sim.file)
}

#pair objects for the helpers (adds kind as factor and MIN)
pairs <- lapply(setNames(names(seeds), names(seeds)), function(pair)
  prep_pair(sims[[pair]], ising_nets[[pair]]$names.E1, ising_nets[[pair]]$names.E2,
            labels[[pair]]))

data.ls <- pairs$love_sadness$data
data.ej <- pairs$enjoyment_joy$data
data.af <- pairs$anxiety_frustration$data

###check, whether settings are correct
lapply(pairs, function(p) c(n.E1 = length(p$names.E1), n.E2 = length(p$names.E2),
                            n.shared = p$n.shared))
lapply(pairs, function(p) table(sign(p$W[upper.tri(p$W)])))
lapply(pairs, function(p) table(p$data$kind))
lapply(pairs, function(p) p$shared)

##################################################
##################################################
#Simulation for mixed emotions (love / sadness)
#a) emotions co-occur, especially in multi-goal situations
#b) different-valence state emotions correlate moderately negatively
#c) the moderate negative correlation is L-shaped
##################################################
##################################################

###a) MIN is higher in E1E2 than in E1 or E2
min.ls <- mean_by_kind(data.ls, "MIN")
min.ls
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     2.78  1.45 0.0216  4500
# 2 E2     2.56  1.48 0.0221  4500
# 3 E1E2   4.92  2.73 0.0865  1000

# --> A confirmed.

###b) overall negative correlation; weaker/positive within E1E2
cor.ls <- cor_by_kind(data.ls)
cor.ls
# set     n      r   ci.l   ci.u             p
# 1 overall 10000 -0.229 -0.247 -0.210 7.271180e-119
# 2      E1  4500  0.040  0.011  0.069  6.978382e-03
# 3      E2  4500  0.228  0.200  0.256  3.561642e-54
# 4    E1E2  1000  0.473  0.424  0.520  5.207630e-57

# --> B confirmed.

#intensity of goal activation affects the central emotion, but not strongly the non-central emotion
data.ls %>%
  group_by(kind) %>%
  summarize(r.int.E1 = cor(intensity, stateE1),
            r.int.E2 = cor(intensity, stateE2), .groups = "drop")
# # A tibble: 3 × 3
# kind  r.int.E1 r.int.E2
# <fct>    <dbl>    <dbl>
#   1 E1       0.804   0.0472
# 2 E2       0.175   0.749 
# 3 E1E2     0.756   0.660 


###c) L-shape (judged visually)
#two arm regressions (E2 ~ E1 in E1 cases, E1 ~ E2 in E2 cases) + negative overall correlation (dashed)
lshape.ls <- test_lshape(pairs$love_sadness, overall.line = FALSE,
                         emotions = c("love", "sadness"))
# --- love / sadness: L-shape (C) ---
#   set               model     n estimate    se   ci.l   ci.u             p
# E1   stateE2 ~ stateE1  4500    0.016 0.006  0.004  0.028  6.978382e-03
# E2   stateE1 ~ stateE2  4500    0.141 0.009  0.123  0.159  3.561642e-54
# overall r(stateE1, stateE2) 10000   -0.229    NA -0.247 -0.210 7.271180e-119
# Overall correlation negative (r < 0, p < .05): TRUE 

lshape.ls
lshape.res <- attr(lshape.ls, "lshape")

# --> C confirmed.

ggsave(file.path(out.dir.fig, "scatter_love_sadness_lshape.png"), plot = lshape.ls,
       width = 7, height = 6, dpi = 300)

###additional implication: parts of a mixed emotion are weaker than the pure emotions
mean_by_kind(data.ls, "stateE1")
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     7.35  3.80 0.0566  4500
# 2 E2     2.70  1.64 0.0244  4500
# 3 E1E2   6.82  3.83 0.121   1000

mean_by_kind(data.ls, "stateE2")
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     3.17  1.51 0.0226  4500
# 2 E2     6.60  2.65 0.0395  4500
# 3 E1E2   5.98  2.60 0.0821  1000

##################################################
##################################################
#Simulation for positive (enjoyment / joy) and negative (anxiety / frustration) blends
#a) emotions co-occur, especially in multi-goal situations
#d) emotion blends correlate highly positively
#e) correlation emotion blends is higher for positive than negative emotions
#f) appraisal overlap (shared appraisals) predicts emotion co-occurrence
##################################################
##################################################

###a) MIN is higher in E1E2 than in E1 or E2
min.ej <- mean_by_kind(data.ej, "MIN")
min.af <- mean_by_kind(data.af, "MIN")
min.ej
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     8.05  4.11 0.0613  4500
# 2 E2     7.66  3.98 0.0593  4500
# 3 E1E2   7.93  4.08 0.129   1000

min.af
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     6.79  3.67 0.0547  4500
# 2 E2     7.13  3.85 0.0574  4500
# 3 E1E2   7.39  3.85 0.122   1000

# --> A confirmed.

###d) positive correlations, strongest in E1E2
cor.ej <- cor_by_kind(data.ej)
cor.af <- cor_by_kind(data.af)
cor.ej
# set     n     r  ci.l  ci.u p
# 1 overall 10000 0.998 0.998 0.998 0
# 2      E1  4500 0.998 0.998 0.998 0
# 3      E2  4500 0.998 0.997 0.998 0
# 4    E1E2  1000 0.998 0.998 0.998 0

cor.af
# set     n     r  ci.l  ci.u p
# 1 overall 10000 0.986 0.985 0.986 0
# 2      E1  4500 0.987 0.987 0.988 0
# 3      E2  4500 0.990 0.989 0.990 0
# 4    E1E2  1000 0.989 0.988 0.991 0

# --> Near perfect correlations.

#intensity of goal activation affects the central emotion and weakly the non-central emotion
lapply(list(enjoyment_joy = data.ej, anxiety_frustration = data.af), function(d)
  d %>%
    group_by(kind) %>%
    summarize(r.int.E1 = cor(intensity, stateE1),
              r.int.E2 = cor(intensity, stateE2), .groups = "drop"))

# $enjoyment_joy
# # A tibble: 3 × 3
# kind  r.int.E1 r.int.E2
# <fct>    <dbl>    <dbl>
#   1 E1       0.854    0.852
# 2 E2       0.832    0.837
# 3 E1E2     0.856    0.854
# 
# $anxiety_frustration
# # A tibble: 3 × 3
# kind  r.int.E1 r.int.E2
# <fct>    <dbl>    <dbl>
#   1 E1       0.820    0.800
# 2 E2       0.840    0.850
# 3 E1E2     0.829    0.832

###e) r(positive blend) > r(negative blend)
r.ej <- cor(data.ej$stateE1, data.ej$stateE2)
r.af <- cor(data.af$stateE1, data.af$stateE2)

r.pos.neg <- r.test(n = nrow(data.ej), n2 = nrow(data.af), r12 = r.ej, r34 = r.af)
r.pos.neg

# Correlation tests 
# Call:r.test(n = nrow(data.ej), r12 = r.ej, r34 = r.af, n2 = nrow(data.af))
# Test of difference between two independent correlations 
# z value 64.03    with probability  0

# report Cohen's q
q.pos.neg <- atanh(r.ej) - atanh(r.af)
q.pos.neg #0.9056267

# --> Statistical difference is not informative as any non-zero difference woule become 
# significant at a given number of simulated samples.

# --> Both correlations are near perfect, so no real difference.

# --> E not confirmed.


###f) shared appraisals predict MIN more strongly than ordinary (non-shared) nodes
#-> joint model per pair: MIN ~ shared + I(shared + ordinary); dB = B.shared - B.ordinary
#-> only shared APPRAISALS are passed as `shared` 
#-> appraisal names from variable_groups.rds (saved in Analysis_final.R, section 0.4)
appraisals <- readRDS(file.path(emp.dir, "variable_groups.rds"))$appraisals

LS.app <- intersect(pairs$love_sadness$shared, appraisals)       
EJ.app <- intersect(pairs$enjoyment_joy$shared, appraisals)      
AF.app <- intersect(pairs$anxiety_frustration$shared, appraisals)
LS.app; EJ.app; AF.app
# [1] "Reality"
# [1] "Pleasantness"       "Goal.conduciveness" "Value"              "Closeness"          "Reality"           
# [1] "Difficulty" "Challenge"  "Reality"

#one comparison node per pair
# test_shared_slope(pairs$love_sadness,        shared = LS.app, ordinary = "Pleasantness")
# test_shared_slope(pairs$enjoyment_joy,       shared = EJ.app, ordinary = "Liking")
# test_shared_slope(pairs$anxiety_frustration, shared = AF.app, ordinary = "I.fell.silent.")

#all non-shared nodes as comparison
slope.ls <- test_shared_slope(pairs$love_sadness,        shared = LS.app, comparison = "all")
# --- love / sadness: shared vs. all 22 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# Reality overall         22       21 -0.395    1.2525  1.941
# Reality      E1         22       19 -0.764    0.7795  1.463
# Reality      E2         22       18 -0.364    0.7325  2.428
# Reality    E1E2         22       22  0.064    1.8315  3.840

# One ordinary/non-shared component is stronger
slope.ls[slope.ls$dB < 0, ]
# shared                                          ordinary     set B.shared B.ordinary     dB SE.dB      p.dB
# 15 Reality               I.felt.my.heartbeat.getting.faster. overall    1.444      1.839 -0.395 0.057 1.0000000
# 36 Reality                                         Challenge      E1    1.175      1.466 -0.291 0.081 0.9998267
# 37 Reality               I.felt.my.heartbeat.getting.faster.      E1    1.004      1.352 -0.348 0.073 0.9999992
# 38 Reality                                       I.frowned..      E1    1.350      2.114 -0.764 0.080 1.0000000
# 47 Reality                                         Closeness      E2    1.607      1.672 -0.064 0.072 0.8125932
# 51 Reality I.wanted.to.be.near.or.close.to.people.or.things.      E2    1.868      2.189 -0.321 0.067 0.9999992
# 52 Reality           I.wanted.to.be.tender..sweet..and.kind.      E2    1.725      2.089 -0.364 0.075 0.9999994
# 59 Reality               I.felt.my.heartbeat.getting.faster.      E2    1.563      1.705 -0.142 0.069 0.9800185

slope.ej <- test_shared_slope(pairs$enjoyment_joy,       shared = EJ.app, comparison = "all")
# --- enjoyment / joy: shared vs. all 1 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# Pleasantness overall          1        0 -4.619    -4.619 -4.619
# Pleasantness      E1          1        0 -5.289    -5.289 -5.289
# Pleasantness      E2          1        0 -4.147    -4.147 -4.147
# Pleasantness    E1E2          1        0 -5.135    -5.135 -5.135
# Goal.conduciveness overall          1        1  1.094     1.094  1.094
# Goal.conduciveness      E1          1        1  0.168     0.168  0.168
# Goal.conduciveness      E2          1        1  1.896     1.896  1.896
# Goal.conduciveness    E1E2          1        0 -0.162    -0.162 -0.162
# Value overall          1        1  2.036     2.036  2.036
# Value      E1          1        1  1.161     1.161  1.161
# Value      E2          1        1  2.746     2.746  2.746
# Value    E1E2          1        1  1.053     1.053  1.053
# Closeness overall          1        0 -6.068    -6.068 -6.068
# Closeness      E1          1        0 -6.605    -6.605 -6.605
# Closeness      E2          1        0 -5.776    -5.776 -5.776
# Closeness    E1E2          1        0 -6.368    -6.368 -6.368
# Reality overall          1        0 -0.431    -0.431 -0.431
# Reality      E1          1        0 -1.317    -1.317 -1.317
# Reality      E2          1        1  0.279     0.279  0.279
# Reality    E1E2          1        0 -1.739    -1.739 -1.739

slope.af <- test_shared_slope(pairs$anxiety_frustration, shared = AF.app, comparison = "all")
# --- anxiety / frustration: shared vs. all 4 non-shared nodes ---
#   shared     set n.ordinary n.dB.pos dB.min dB.median dB.max
# Difficulty overall          4        4  2.418    3.3390  4.299
# Difficulty      E1          4        4  3.080    3.2090  4.322
# Difficulty      E2          4        4  1.429    2.6205  5.033
# Difficulty    E1E2          4        4  2.105    3.0745  4.251
# Challenge overall          4        4  2.408    2.9945  4.404
# Challenge      E1          4        4  2.472    2.7870  4.652
# Challenge      E2          4        4  2.062    3.0475  4.180
# Challenge    E1E2          4        4  1.034    2.0290  4.105
# Reality overall          4        4  0.695    1.7280  3.071
# Reality      E1          4        4  1.369    2.0680  3.457
# Reality      E2          4        3 -0.447    1.0075  3.773
# Reality    E1E2          4        2 -0.796    0.5810  1.994

###additional implication: parts of a blend are stronger than the pure emotions
lapply(list(enjoyment_joy = data.ej, anxiety_frustration = data.af), function(d)
  list(E1 = mean_by_kind(d, "stateE1"), E2 = mean_by_kind(d, "stateE2")))

# $enjoyment_joy
# $enjoyment_joy$E1
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     8.56  4.45 0.0663  4500
# 2 E2     8.01  4.24 0.0632  4500
# 3 E1E2   8.42  4.42 0.140   1000
# 
# $enjoyment_joy$E2
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     8.05  4.11 0.0613  4500
# 2 E2     7.66  3.98 0.0593  4500
# 3 E1E2   7.93  4.08 0.129   1000
# 
# 
# $anxiety_frustration
# $anxiety_frustration$E1
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     6.87  3.64 0.0542  4500
# 2 E2     7.16  3.82 0.0570  4500
# 3 E1E2   7.43  3.82 0.121   1000
# 
# $anxiety_frustration$E2
# # A tibble: 3 × 5
# kind   mean    sd     se     n
# <fct> <dbl> <dbl>  <dbl> <int>
#   1 E1     7.34  3.96 0.0591  4500
# 2 E2     8.29  4.50 0.0671  4500
# 3 E1E2   8.42  4.40 0.139   1000

# --> F confirmed for Love/Sadness and Anxiety/Frustration.
# --> F can hardly be confirmed or refuted for Joy/Anjoyment which share all but one node.

##################################################
##################################################
#Summary table (cf. Table 5, empirical example columns)
##################################################
##################################################

#M (SD); SEs are not reported for the simulated data (they only reflect the chosen n)
fmt <- function(m, sd) sprintf("%.2f (%.2f)", m, sd)

#F: per shared appraisal (overall), median dB vs. all non-shared nodes, its range,
#   and the number of comparisons with dB > 0
fmt_f <- function(s) {
  o <- attr(s, "summary")
  o <- o[o$set == "overall", ]
  paste(sprintf("%s: median dB = %.2f [%.2f, %.2f], %d/%d > 0",
                o$shared, o$dB.median, o$dB.min, o$dB.max, o$n.dB.pos, o$n.ordinary),
        collapse = "; ")
}

table5.emp <- data.frame(
  phenomenon = c("A", "A", "A", "B", "C", "D", "E", "F"),
  criterion = c("Mean MIN score: E1", "Mean MIN score: E2", "Mean MIN score: E1E2",
                "Correlation of mixed emotion states",
                "Arm regressions: E2 ~ E1 (E1 cases), E1 ~ E2 (E2 cases)",
                "Correlation of same-valenced emotion states",
                "r positive vs. negative blend (Cohen's q)",
                "Joint regression MIN ~ shared appraisal + ordinary node: dB = B.shared - B.ordinary"),
  mixed = c(fmt(min.ls$mean, min.ls$sd),
            sprintf("r = %.3f", cor.ls$r[cor.ls$set == "overall"]),
            sprintf("E1 arm: B = %.2f; E2 arm: B = %.2f",
                    lshape.res$estimate[1], lshape.res$estimate[2]),
            "-", "-",
            fmt_f(slope.ls)),
  positive = c(fmt(min.ej$mean, min.ej$sd), "-", "-",
               sprintf("r = %.3f", cor.ej$r[cor.ej$set == "overall"]),
               sprintf("r = %.3f vs. %.3f, q = %.2f", r.ej, r.af, q.pos.neg),
               fmt_f(slope.ej)),
  negative = c(fmt(min.af$mean, min.af$sd), "-", "-",
               sprintf("r = %.3f", cor.af$r[cor.af$set == "overall"]),
               "-",
               fmt_f(slope.af))
)
table5.emp
write.csv(table5.emp, file.path(out.dir, "table5_empirical_example.csv"), row.names = FALSE)

##################################################
##################################################
#Figures
##################################################
##################################################

kind.cols <- c(E1 = "#E69F00", E2 = "#009E73", E1E2 = "#0072B2")

###node labels (L1, S1, SL1, ...)
node.labels <- list(
  love_sadness        = make_labels(pairs$love_sadness$names.E1, pairs$love_sadness$names.E2,
                                    prefix.E1 = "L", prefix.E2 = "S", prefix.shared = "SL"),
  enjoyment_joy       = make_labels(pairs$enjoyment_joy$names.E1, pairs$enjoyment_joy$names.E2,
                                    prefix.E1 = "E", prefix.E2 = "J", prefix.shared = "EJ"),
  anxiety_frustration = make_labels(pairs$anxiety_frustration$names.E1, pairs$anxiety_frustration$names.E2,
                                    prefix.E1 = "A", prefix.E2 = "F", prefix.shared = "AF")
)

###networks (one PNG per pair)
net.files <- file.path(out.dir.fig, paste0("network_", names(pairs), ".png"))
names(net.files) <- names(pairs)
for (pair in names(pairs)) {
  keep <- union(pairs[[pair]]$names.E1, pairs[[pair]]$names.E2)
  plot_network(pairs[[pair]], labels = unname(node.labels[[pair]][keep]),
               file = net.files[[pair]])
}

###Combined network and scatter plots (A-C networks, D-F scatter plots)
#col 1 love/sadness = mixed, col 2 enjoyment/joy = positive, col 3 anxiety/frustration = negative
combined.net.scatter <-
  net_panel(net.files[["love_sadness"]]) +
  net_panel(net.files[["enjoyment_joy"]]) +
  net_panel(net.files[["anxiety_frustration"]]) +
  scatter_states(pairs$love_sadness,        colour = FALSE, title = NULL,
                 xlab = "stateE1", ylab = "stateE2") +
  scatter_states(pairs$enjoyment_joy,       colour = FALSE, title = NULL,
                 xlab = "stateE1", ylab = NULL) +   # y title only on left column
  scatter_states(pairs$anxiety_frustration, colour = FALSE, title = NULL,
                 xlab = "stateE1", ylab = NULL) +
  plot_layout(nrow = 2, byrow = TRUE) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E", "F")),
    theme = theme(plot.tag = element_text(face = "bold", size = 16))
  )
combined.net.scatter
ggsave(file.path(out.dir.fig, "combined_network_scatter_ABCDEF.png"),
       plot = combined.net.scatter, width = 15, height = 10, dpi = 300)

###per-pair figures (network over scatter plot)
for (pair in names(pairs))
  pair_figure(net.files[[pair]], pairs[[pair]],
              file.path(out.dir.fig, paste0("fig_", pair, "_AB.png")))

###MIN histograms by kind of goal activation
#x-axis per panel runs to the smaller number of nodes of the pair (maximum of MIN)
min.long <- bind_rows(lapply(names(pairs), function(pair)
  data.frame(pair = pair, kind = pairs[[pair]]$data$kind, MIN = pairs[[pair]]$data$MIN)))
min.long$pair <- factor(min.long$pair, levels = names(pairs),
                        labels = c("A: Love / sadness", "B: Enjoyment / joy",
                                   "C: Anxiety / frustration"))

min.hist <- ggplot(min.long, aes(x = MIN, fill = kind)) +
  geom_histogram(position = "dodge", binwidth = 1, boundary = 0, colour = NA) +
  facet_wrap(~ pair, scales = "free") +
  scale_fill_manual(values = kind.cols) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "MIN score", y = "Frequency", fill = "Activation") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())
min.hist
ggsave(file.path(out.dir.fig, "combined_MIN_hist_ABC.png"), plot = min.hist,
       width = 14, height = 4.5, dpi = 300)
