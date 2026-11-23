## IMEC theory comparison (coding = IMEC/IMEC.docx, state 2026-11-17)
## Maier, van Dongen & Borsboom (2023). Ising Model of Explanatory Coherence.
##
## The edge table in IMEC.docx lists EFFECTIVE weights per hypothesis (i.e.,
## joint explanations are already split). The network is therefore built
## directly from that table, not via IMEC::explain() (which would split again).
## Each edge is tagged by type so that sensitivity analyses can vary it:
##   direct = follows from the hypothesis alone
##   joint  = hypothesis co-explains together with another one (weight already split)
##   aux    = needs one auxiliary assumption (rule: .5)
##   contra = contradiction (-4)

#all paths are relative to the project root: open Analysis.Rproj before running this script

library(IsingSampler)  # IsingLikelihood()
library(qgraph)
library(IMEC)

## ---- Phenomena and evidence (thresholds) ------------------------------------
phen <- c("A1", "A2", "A3", "B", "C", "D", "E", "F")
ev_default <- c(A1 = 2,  # emotion blends co-occur
                A2 = 2,  # mixed emotions co-occur
                A3 = 2,  # valence antonyms (happy-sad) co-occur (contested)
                B  = 2,  # mixed emotions correlate negatively
                C  = 1,  # L-shape (antonyms, strictly unipolar; contested)
                D  = 2,  # emotion blends correlate positively
                E  = 2,  # positive blends correlate more strongly than negative
                F  = 1)  # appraisal overlap predicts co-occurrence (contested)

## ---- Propositions ------------------------------------------------------------
theories <- list(
  NW  = c("NW1", "NW2"),   # NW1 overlapping causal networks; NW2 opposing components inhibit
  BI  = c("BI1"),          # BI1 one point in bipolar valence x arousal space
  AIM = c("AIM1", "AIM2"), # AIM1 two all-to-all self-exciting pools; AIM2 inhibition between pools
  AP  = c("AP1", "AP2")    # AP1 appraisal sets, multiple activated -> co-occurrence; AP2 conflict -> mixed
)

## ---- Edges (effective weights, as in IMEC.docx) -----------------------------
edges <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
  from  to   weight type
  NW1   A1   1      direct
  NW1   A2   1      direct
  NW1   A3   1      direct
  NW1   B    .5     joint
  NW1   C    .5     joint
  NW1   D    1      direct
  NW1   E    .5     aux
  NW1   F    .5     aux
  NW2   B    .5     joint
  NW2   C    .5     joint
  BI1   A1   .5     aux
  BI1   A2   .5     aux
  BI1   A3   -4     contra
  BI1   B    .5     aux
  BI1   C    1      direct
  BI1   D    .5     aux
  BI1   E    .5     aux
  AIM1  A1   .5     aux
  AIM1  A2   .5     joint
  AIM1  A3   .5     joint
  AIM1  C    .5     joint
  AIM1  D    .5     aux
  AIM2  A2   .5     joint
  AIM2  A3   .5     joint
  AIM2  B    1      direct
  AIM2  C    .5     joint
  AP1   A1   1      direct
  AP1   D    .5     aux
  AP1   E    .5     aux
  AP1   F    1      direct
  AP2   A2   1      direct
  AP2   A3   1      direct
  BI1   AIM1 -4     contra
")

## ---- Build the weight matrix -------------------------------------------------
## v = list of variant switches:
##   aux_factor   multiplies all 'aux' edges (e.g., .5 -> .25 or .75)
##   no_competition drops hypothesis-hypothesis contradictions
##   extra        data.frame of additional edges (from, to, weight, type)
build <- function(ev = ev_default, v = list()) {
  e <- edges
  if (!is.null(v$extra)) e <- rbind(e, v$extra)
  if (!is.null(v$aux_factor)) e$weight[e$type == "aux"] <- e$weight[e$type == "aux"] * v$aux_factor
  props <- unlist(theories, use.names = FALSE)
  if (isTRUE(v$no_competition)) e <- e[!(e$from %in% props & e$to %in% props), ]
  nodes <- c(props, phen)
  W <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (i in seq_len(nrow(e))) W[e$from[i], e$to[i]] <- W[e$to[i], e$from[i]] <- e$weight[i]
  list(W = W, ev = ev)
}

main <- build()
print(main$W)                      # check the coded network before computing anything

## ===========================================================================

## ---- 1. Figures -------------------------------------------------------------------
lab <- colnames(main$W)
grp <- factor(c(rep(c("Network theory", "Bipolar valence", "Affective Ising Model",
                      "Appraisal account of mixed emotions"), lengths(theories)),
                rep("Phenomena", length(phen))),
              levels = c("Network theory", "Bipolar valence", "Affective Ising Model",
                         "Appraisal account of mixed emotions", "Phenomena"))
lay <- rbind(cbind(c(-1, -.6), 2),                          # NW  top-left
             cbind(-.8, -2),                                # BI  bottom-left
             cbind(c(.6, 1), -2),                           # AIM bottom-right
             cbind(c(.6, 1), 2),                            # AP  top-right
             cbind(seq(-1, 1, length.out = length(phen)), 0))
png("theory comparison/imec_network.png", width = 2400, height = 1600, res = 200)
qgraph(main$W, layout = lay, labels = lab, groups = grp,
       shape = ifelse(lab %in% phen, "square", "circle"),
       color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "grey85"),
       posCol = "#0072B2", negCol = "#D55E00", vsize = 7, label.cex = .9,
       legend = TRUE, legend.cex = .5)
dev.off()


## ---- 2. Main analysis ----------------------------------------------------------
props <- unlist(theories, use.names = FALSE)
r <- computeIMEC(main$W, main$ev[phen], phen, theory1 = props)  # all propositions as "theory1"
P <- setNames(r$ExplanatoryCoherenceT1[[2]], props)
sapply(theories, function(h) round(mean(P[h]),3))
# NW    BI   AIM    AP 
# 0.934 0.007 0.992 0.988 

round(P, 3)
# NW1   NW2   BI1  AIM1  AIM2   AP1   AP2 
# 1.000 0.867 0.007 0.993 0.991 0.996 0.980 





## ---- 4. Robustness: auxiliary assumptions as explicit proposition nodes -------
## Main coding: an explanation that needs one auxiliary assumption gets .5.
## Here the auxiliary is added as its own proposition node instead (Thagard's
## Principle 2c): core and auxiliary jointly explain the phenomenon, each with .5.
## The auxiliary belongs to the core's theory and therefore ENTERS THE THEORY MEAN.
## Default: one auxiliary per 'aux' edge. To code that several edges rest on the
## SAME assumption, give them the same label (example commented out below).
aux_e <- edges[edges$type == "aux", ]
aux_e$aux <- paste0(aux_e$from, "a_", aux_e$to)          # e.g. "BI1a_A1"
# aux_e$aux[aux_e$from == "BI1" & aux_e$to %in% c("A1", "A2")] <- "BI1a_readout"
aux_ids <- unique(aux_e$aux)

## theory membership of the auxiliaries
owner <- setNames(rep(names(theories), lengths(theories)), unlist(theories))
theories_aux <- theories
for (a in aux_ids) {
  t <- owner[[aux_e$from[aux_e$aux == a][1]]]
  theories_aux[[t]] <- c(theories_aux[[t]], a)
}

## auxiliary -> phenomenon weights (core -> phenomenon edges stay at .5)
A <- matrix(0, length(aux_ids), length(phen), dimnames = list(aux_ids, phen))
for (i in seq_len(nrow(aux_e))) A[aux_e$aux[i], aux_e$to[i]] <- aux_e$weight[i]

## Exact computation without enumerating 2^(15 + k) states: auxiliary nodes have
## threshold 0 and are connected only to phenomena, so they can be summed out
## analytically. With h_a = sum_j A[a, j] * x_j:
##   P(core state) is proportional to P_main(core state) * prod_a 2*cosh(h_a)
##   P(a = +1 | core state) = plogis(2 * h_a)
props <- unlist(theories, use.names = FALSE)
thr <- c(rep(0, length(props)), main$ev[phen])
L   <- IsingLikelihood(main$W, thr, beta = 1, responses = c(-1L, 1L))
S   <- as.matrix(L[, -1]); colnames(S) <- colnames(main$W)
h   <- S[, phen, drop = FALSE] %*% t(A)                  # states x auxiliaries
w   <- L$Probability * apply(2 * cosh(h), 1, prod); w <- w / sum(w)
P_aux_model <- c(colSums((S == 1) * w), colSums(plogis(2 * h) * w))
round(P_aux_model, 3)                                    # per node
sapply(theories_aux, function(x) round(mean(P_aux_model[x]), 3))  # per theory
