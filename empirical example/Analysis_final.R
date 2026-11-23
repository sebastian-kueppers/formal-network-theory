
#### --------------------------------- ####
## CONTENTS ---------------------------- ##
#### --------------------------------- ####

# Packages

library(dplyr)
library(tidyr)
library(purrr)
library(segmented)
library(lme4)
library(doParallel)
library(doSNOW)
library(foreach)
library(corpcor)
library(Matrix)
library(bootnet)
library(IsingFit)
library(qgraph)
library(bgms)
library(glmnet)
library(qgraph)

#all paths are relative to the project root: open Analysis.Rproj before running this script
emp.dir <- "empirical example"

#### --------------------------------- ####
## 0 - PREPROCESSING ------------------- ##
#### --------------------------------- ####
# 
# data <- read.csv(file.path(emp.dir, "survey_data.csv"))
# 
# ## ------------------------------------- ##
# # 0.1 - Renaming columns ---------------- #
# ## ------------------------------------- ##
# 
# # rename col names
# colnames(data)[8] <- "codeword"
# colnames(data)[9] <- "age"
# colnames(data)[10] <- "gender"
# colnames(data)[11] <- "experience"
# 
# # emotions
# colnames(data)[12:74] <- sub(paste0(
#   "Now.think.again.about.the.emotional.situation.you.experienced..",
#   "What.emotions.did.you.experience.in.that.situation...",
#   "Click.only.on.the.ones.you.experienced...."
# ), "", colnames(data)[12:74])
# 
# colnames(data)[12:74] <- sub("\\.$", "", colnames(data)[12:74])
# 
# emotions <- tolower(gsub("[._]", "-", colnames(data)[12:74]))
# emotions[emotions == "irritated"] <- "irritation"
# colnames(data)[12:74] <- emotions
# 
# # emotions scale
# colnames(data)[75] <- "emotion.scale.instruction"
# colnames(data)[76:138] <- paste0(
#   emotions,
#   ".scale"
# )
# 
# # appraisals
# colnames(data)[139:185] <- sub(".*\\.\\.\\.\\.([A-Z][^A-Z]*)[A-Z].*", "\\1", colnames(data)[139:185])
# colnames(data)[152] <- "Accountability.Circumstances"
# colnames(data)[153] <- "Accountability.Other"
# colnames(data)[154] <- "Accountability.Self"
# colnames(data)[156] <- "Control.Circumstances"
# colnames(data)[157] <- "Control.Other"
# colnames(data)[158] <- "Control.Self"
# colnames(data)[183] <- "Self.Esteem"
# 
# appraisals <- colnames(data[139:185])
# 
# # appraisals scale
# colnames(data)[186] <- "appraisals.scale.instruction"
# colnames(data)[187:233] <- paste0(
#   appraisals,
#   '.scale'
# )
# 
# # bodily symptoms
# colnames(data)[234:250] <- sub(paste0(
#   "Considering.the.emotional.situation.you.experienced..",
#   "what.bodily.symptoms.did.you.experience.in.that.situation...",
#   "Click.only.on.the.ones.you.experienced..."
# ), "", colnames(data)[234:250])
# 
# symptoms <- colnames(data)[234:250]
# 
# # bodily symptoms scale
# colnames(data)[251] <- "bodily.symptoms.scale.instruction"
# colnames(data)[252:268] <- paste0(
#   symptoms,
#   '.scale'
# )
# 
# # expressions
# colnames(data)[269:290] <- sub(paste0(
#   "Considering.the.emotional.situation.you.experienced..",
#   "what.expressions.did.you.have.in.that.situation...",
#   "Click.only.on.the.ones.you.had..."
# ), "", colnames(data)[269:290])
# 
# expressions <- colnames(data)[269:290]
# 
# # expressions scale
# colnames(data)[291] <- "expressions.scale.instruction"
# colnames(data)[292:313] <- paste0(
#   expressions,
#   '.scale'
# )
# 
# # action tendencies
# colnames(data)[314:353] <- sub(paste0(
#   "Considering.the.emotional.situation.you.experienced..",
#   "what.did.you.want.to.do.in.that.situation...",
#   "Click.only.on.the.ones.you.wanted.to.do..."
# ), "", colnames(data)[314:353])
# 
# actions <- colnames(data)[314:353]
# 
# # action tendencies scale
# colnames(data)[354] <- "action.tendencies.scale.instruction"
# colnames(data)[355:394] <- paste0(
#   actions,
#   '.scale'
# )
# 
# # feelings
# colnames(data)[395:416] <- sub(paste0(
#   "Considering.the.emotional.situation.you.experienced..",
#   "what.did.you.feel.in.that.situation...",
#   "Click.only.on.the.ones.you.felt..."
# ), "", colnames(data)[395:416])
# 
# feelings <- colnames(data)[395:416]
# 
# # feelings scale
# colnames(data)[417] <- "feelings.scale.instruction"
# colnames(data)[418:439] <- paste0(
#   feelings,
#   '.scale'
# )
# 
# # seriousness
# colnames(data)[440] <- "seriousness"
# 
# # all components
# all_components <- c(appraisals,
#                     symptoms,
#                     expressions,
#                     actions,
#                     feelings)
# 
# 
# ## ------------------------------------- ##
# # 0.2 - Replacing values ---------------- #
# ## ------------------------------------- ##
# 
# data[data == "N/A"] <- NA
# data[data == "Yes"] <- TRUE
# data[data == "TRUE"] <- TRUE
# data[data == "No"] <- FALSE
# data[data == "FALSE"] <- FALSE
# data[data == "1 (does not apply at all)"] <- 1
# data[data == "7 (applies very much)"] <- 7
# 
# ## ------------------------------------- ##
# # 0.3 - Participant exclusion ----------- #
# ## ------------------------------------- ##
# 
# # only keep participants that have submitted
# data <- subset(data, Date.submitted != "")
# 
# # seriousness
# data$seriousness <- ifelse(
#   data$seriousness == "I have taken part seriously",
#   TRUE,
#   FALSE
# )
# 
# # omit those who have not taken part seriously
# data <- subset(data, seriousness == TRUE)
# 
# # omit those who did not described an emotional situation
# # we omit all participants who did not entered an understandable English word
# # --> None.
# 
# # make all emotion and component entries logical
# data[emotions] <- lapply(
#   data[emotions],
#   function(x) as.logical(x)
# )
# 
# data[all_components] <- lapply(
#   data[all_components],
#   function(x) as.logical(x)
# )
# 
# ## ------------------------------------- ##
# # 0.4 - Anonymize ----------------------- #
# ## ------------------------------------- ##
# 
# # remove dates
# data$Date.submitted <- NULL
# data$Date.started <- NULL
# data$Date.last.action <- NULL
# 
# # remove codeword
# data$codeword <- NULL
# 
# # remove emotional experience
# data$experience <- NULL
# 
# # check for descriptives 
# psych::describe(data$age)
# # vars   n  mean    sd median trimmed   mad min max range skew kurtosis   se
# # X1    1 617 43.12 13.53     41   42.18 13.34  19  95    76  0.6    -0.13 0.54
# 
# # make age groups
# age_breaks <- quantile(data$age, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
# age_breaks  # report these cut points in the paper
# data$age_group <- cut(
#   data$age,
#   breaks = age_breaks,
#   labels = c("young", "middle", "old"),
#   include.lowest = TRUE
# )
# 
# table(data$age_group)
# # young middle    old 
# # 224    197    196
# 
# # remove raw age
# data$age <- NULL
# 
# table(data$gender)
# # diverse       female      male no information 
# # 4            305            307              1 
# 
# # remove gender
# data$gender <- NULL
# 
# 
# # save
# saveRDS(data, file.path(emp.dir, "data_preprocessed.RDS"))
# 
# # save variable groups (column names per category) so 0.5 runs without 0.1-0.4
# saveRDS(list(emotions = emotions, appraisals = appraisals, symptoms = symptoms,
#              expressions = expressions, actions = actions, feelings = feelings),
#         file.path(emp.dir, "variable_groups.rds"))
# 

## ------------------------------------- ##
# 0.5 - Load preprocessed data ---------- #
## ------------------------------------- ##

# load
data <- readRDS(file.path(emp.dir, "data_preprocessed.RDS"))
list2env(readRDS(file.path(emp.dir, "variable_groups.rds")), envir = environment())



# `emotions` was already defined by name in 0.1 (no positional re-derivation,
# since columns were removed in 0.4)

# Keep only emotions that are investigated in one of these papers:
# - Diener & Iran-Nejad (1986)
# - Erbas et al. (2015)
# - Ketonen et al. (2018)
# - Moeller et al. (2018)
# - Tong & Jia (2017)
# - Trampe et al. (2015)
# - Vansteelandt et al. (2005)

emotions.Diener.Iran.Nejad <- c(
  "joy", "pleasure", "enjoyment", "frustration", "anger", "worry", "anxiety"
)
pairs.Diener.Iran.Nejad <- combn(emotions.Diener.Iran.Nejad, 2)

emotions.Watson.Stanton <- c(
  "guilt", "sadness", "fear"
) 
pairs.Watson.Stanton <- combn(emotions.Watson.Stanton, 2)

emotions.Moeller <- c(
  "interest", "boredom"
)
pairs.Moeller <- combn(emotions.Moeller, 2)

emotions.Ketonen <- c(
  "interest", "enthusiasm", "anxiety", "irritation"
)
pairs.Ketonen <- combn(emotions.Ketonen, 2)

emotions.Vansteelandt <- c(
  "anger", "disgust", "shame", "guilt", "fear", "sadness", "joy", "love"
)
pairs.Vansteelandt <- combn(emotions.Vansteelandt, 2)

# Vansteelandt et al. summarize anger-disgust and shame-guilt
# Remove these two pairs for which we have no data
pairs.Vansteelandt <- pairs.Vansteelandt[, -c(1, 14)]

# aggregate
pairs.in_papers <- cbind(
  pairs.Diener.Iran.Nejad,
  pairs.Watson.Stanton,
  pairs.Moeller,
  pairs.Ketonen,
  pairs.Vansteelandt
)

saveRDS(pairs.in_papers, file = file.path(emp.dir, "emotion_pairs_from_papers.rds"))

# canonical key per column, order-independent (joy/anger == anger/joy)
keys <- apply(pairs.in_papers, 2, function(col) paste(sort(col), collapse = "_"))

pairs.in_papers <- pairs.in_papers[, !duplicated(keys)]

emotions.in_papers <- sort(unique(as.vector(pairs.in_papers)))

emotions.to_drop <- setdiff(emotions, emotions.in_papers)
cols.to_drop <- c(emotions.to_drop, paste0(emotions.to_drop, ".scale"))

data <- data[, !colnames(data) %in% cols.to_drop]

# component vectors from 0.1 are name-based, so they stay valid after
# dropping columns; intersect() only guards against missing columns
appraisals <- intersect(appraisals, colnames(data))
symptoms <- intersect(symptoms, colnames(data))
expressions <- intersect(expressions, colnames(data))
actions <- intersect(actions, colnames(data))
feelings <- intersect(feelings, colnames(data))

all_components <- c(appraisals,
                    symptoms,
                    expressions,
                    actions,
                    feelings)

valence_df <- read.csv(file.path(emp.dir, "VALENCE_DF.CSV"), sep=";")




#### --------------------------------- ####
## 1 - EMOTION-COMPONENT RELATIONSHIPS - ##
#### --------------------------------- ####

# We consider that a component belongs to an emotion if that component was 
# present in more than 10% (higher than chance) of the situations in which the
# respective emotion was present.

# Get emotions with positive or negative valence
valence.non_complex <- valence_df[valence_df$valence %in% c("positive", "negative"), ]
emotions.non_complex <- valence.non_complex$emotion[valence.non_complex$emotion %in% emotions.in_papers]

# --- DF 1: emotion x component percentages ---

emo.comp.rates <- lapply(emotions.non_complex, function(emo) {
  emo_scale <- paste0(emo, ".scale")
  
  # Rows where this emotion is active
  emo_rows <- data[[emo]] == TRUE & !is.na(data[[emo]]) &
    data[[emo_scale]] > 1  & !is.na(data[[emo_scale]])
  
  n_obs_emo <- sum(emo_rows)
  
  # For each component, get % of those rows where component is active
  comp_rates <- sapply(all_components, function(comp) {
    comp_scale <- paste0(comp, ".scale")
    
    comp_active <- data[[comp]] == TRUE & !is.na(data[[comp]]) &
      data[[comp_scale]] > 1  & !is.na(data[[comp_scale]])
    
    if (n_obs_emo == 0) return(NA)
    sum(comp_active & emo_rows) / n_obs_emo
  })
  
  n_obs_comp <- sapply(all_components, function(comp) {
    comp_scale <- paste0(comp, ".scale")
    
    comp_active <- data[[comp]] == TRUE & !is.na(data[[comp]]) &
      data[[comp_scale]] > 1  & !is.na(data[[comp_scale]])
    
    sum(comp_active & emo_rows)
  })
  
  data.frame(
    emotion    = emo,
    component  = all_components,
    n_obs_emo  = n_obs_emo,
    n_obs_comp = n_obs_comp,
    pct        = comp_rates,
    row.names  = NULL
  )
})

emo.comp.rates <- do.call(rbind, emo.comp.rates)
emo.comp.rates <- merge(emo.comp.rates, valence.non_complex, by = "emotion", all.x = TRUE)

# --- DF 2: per-emotion summary ---

COMP_THRESHOLD <- 0.40

emo.num.comps <- do.call(rbind, lapply(emotions.non_complex, function(emo) {
  rows      <- emo.comp.rates[emo.comp.rates$emotion == emo, ]
  n_obs     <- rows$n_obs_emo[1]
  n_comp    <- sum(rows$pct >= COMP_THRESHOLD, na.rm = TRUE)
  
  data.frame(emotion = emo, n_obs = n_obs, n_components = n_comp)
}))

emo.num.comps <- merge(emo.num.comps, valence.non_complex, by = "emotion", all.x = TRUE)

# How many emotions have 30+ obs?
nrow(emo.num.comps) # 17
emo.num.comps.30 <- emo.num.comps[emo.num.comps$n_obs >= 30,]
nrow(emo.num.comps.30) # 16 (all but boredom) 
mean(emo.num.comps.30$n_components) # 16.6875

# Get dataframe of components belonging to an emotion
emo.comp.rates.belongs <- subset(emo.comp.rates, pct >= COMP_THRESHOLD)


## ------------------------------------- ##
# 1.1 Number of components per valence -- #
## ------------------------------------- ##

# How many components does a positive emotion (on average) have?

emo.num.comps.pos <- emo.num.comps.30[emo.num.comps.30$valence == "positive",]

comp.counts.pos <- emo.num.comps.pos %>%
  group_by(emotion) 

mean(comp.counts.pos$n_components) # 16.5
N_COMP.POS <- round(mean(comp.counts.pos$n_components)) # 16 (IEC 60559)


# How many components does a negative emotion (on average) have?

emo.num.comps.neg <- emo.num.comps.30[emo.num.comps.30$valence == "negative",]

comp.counts.neg <- emo.num.comps.neg %>%
  group_by(emotion) 

mean(comp.counts.neg$n_components) # 16.8
N_COMP.NEG <- round(mean(comp.counts.neg$n_components)) # 17

## ------------------------------------- ##
# 1.1.1 Additional Decsriptives --------- #

# Do positive and negative emotions vary in their frequency or intensity?

mean_int.df <- data %>%
  mutate(across(ends_with(".scale"), as.numeric)) %>% 
  pivot_longer(
    cols = ends_with(".scale"),
    names_to = "emotion",
    values_to = "intensity"
  ) %>%
  mutate(emotion = sub("\\.scale$", "", emotion)) %>%   
  inner_join(valence.non_complex, by = "emotion") %>%    
  filter(!is.na(intensity), intensity != 1) %>% # intensity 1 = not at all experienced                          
  group_by(valence) %>%
  summarise(
    mean_intensity = mean(intensity),
    sd_intensity = sd(intensity), 
    n = n(), 
    .groups = "drop"
  )

mean_int.df$n[mean_int.df$valence == "positive"] # 513
mean_int.df$n[mean_int.df$valence == "negative"] # 1401

mean_int.df$mean_intensity[mean_int.df$valence == "positive"] # 6.282651
mean_int.df$mean_intensity[mean_int.df$valence == "negative"] # 5.513919

mean_int.df$sd_intensity[mean_int.df$valence == "positive"] # 0.9499338
mean_int.df$sd_intensity[mean_int.df$valence == "negative"] # 1.402034

# Which emotions often co-occurred?

# emotion_pairs <- combn(emotions.non_complex, 2, simplify = FALSE)
emotion_pairs <- pairs.in_papers

# occurrence dataframe
occ <- sapply(emotions.in_papers, function(e) {
  flag  <- data[[e]]
  scale <- as.numeric(data[[paste0(e, ".scale")]])
  out <- (flag == TRUE) & (scale > 1)
  out[is.na(out)] <- FALSE            
  out
})

# co-occurrence counts
co.occ <- crossprod(occ)  

# upper tri
idx  <- which(upper.tri(co.occ), arr.ind = TRUE)
emo1 <- rownames(co.occ)[idx[, 1]]
emo2 <- colnames(co.occ)[idx[, 2]]

# valence lookup vector
val <- setNames(valence.non_complex$valence, valence.non_complex$emotion)

# co-occurrence df
val <- setNames(valence.non_complex$valence, valence.non_complex$emotion)
v1  <- val[emo1]; v2 <- val[emo2]
valence_pair <- ifelse(is.na(v1) | is.na(v2), NA_character_,
                       ifelse(v1 == v2, v1, "mixed"))

co.occ.val <- tibble(
  emo1 = emo1,
  emo2 = emo2,
  valence_pair = valence_pair,
  n = co.occ[idx]
)

# sort by n (desc)
co.occ.val <- co.occ.val %>% arrange(desc(n))

# top 5 emotion pairs
head(co.occ.val[co.occ.val$valence_pair == "negative",], n = 5)
# emo1        emo2        valence_pair     n
# <chr>       <chr>       <chr>        <dbl>
# 1 anxiety     frustration negative       143
# 2 anxiety     worry       negative       140
# 3 frustration irritation  negative       113
# 4 frustration worry       negative       107
# 5 anxiety     sadness     negative       100

head(co.occ.val[co.occ.val$valence_pair == "positive",], n = 5)
# emo1       emo2       valence_pair     n
# <chr>      <chr>      <chr>        <dbl>
# 1 enjoyment  joy        positive        71
# 2 enjoyment  enthusiasm positive        54
# 3 enjoyment  love       positive        52
# 4 joy        love       positive        50
# 5 enthusiasm joy        positive        45

head(co.occ.val[co.occ.val$valence_pair == "mixed",], n = 5)
# emo1    emo2      valence_pair     n
# <chr>   <chr>     <chr>        <dbl>
# 1 love    sadness   mixed           36
# 2 anxiety love      mixed           24
# 3 love    worry     mixed           18
# 4 anxiety enjoyment mixed           16
# 5 fear    love      mixed           14

# How often are happiness and sadness experienced together?
co.occ["joy", "sadness"] # 5
co.occ["enjoyment", "sadness"] # 6


## ------------------------------------- ##
# 1.2 Ising Networks from Top Pairs #
## ------------------------------------- ##

# Ising-network analog of make_pair_bootnet()

node_cols <- function(p) {
  v <- colnames(p$W)
  ifelse(v %in% p$shared, "#B5C400",
         ifelse(v %in% p$names.E1, "#E69F00", "#009E73"))
}

# plot network
plot_bootnet <- function(net, names.E1, names.E2, node.labels = NULL,
                         label = NULL, file = NULL, layout = "spring", ...) {
  
  W <- if (inherits(net, "bootnetResult")) net$graph
  else if (inherits(net, "qgraph"))   qgraph::getWmat(net)
  else                                 as.matrix(net)
  
  keep <- union(names.E1, names.E2)
  if (!all(keep %in% colnames(W)))
    stop("Nodes missing from the network: ",
         paste(setdiff(keep, colnames(W)), collapse = ", "))
  W <- W[keep, keep, drop = FALSE]
  
  p <- list(W = W, shared = intersect(names.E1, names.E2), names.E1 = names.E1)
  cols <- node_cols(p)                       # needs original names
  
  if (is.null(node.labels)) node.labels <- colnames(W)
  stopifnot(length(node.labels) == ncol(W))
  
  if (!is.null(file)) png(file, width = 7, height = 7, units = "in", res = 300)
  g <- qgraph(W, weighted = TRUE, theme = "colorblind", layout = layout,
              edge.width = 0.5, color = cols, labels = node.labels,
              title = label, ...)
  if (!is.null(file)) dev.off()
  invisible(g)
}


fit_ising_at_lambda <- function(node_mat, lambda) {
  X <- as.matrix(node_mat); p <- ncol(X); nm <- colnames(X)
  W  <- matrix(0, p, p, dimnames = list(nm, nm))
  b0 <- setNames(numeric(p), nm)
  for (j in seq_len(p)) {
    fit <- glmnet(X[, -j, drop = FALSE], X[, j], family = "binomial",
                  lambda = lambda, standardize = FALSE)
    b <- as.numeric(coef(fit))
    b0[j]    <- b[1]               # intercept = threshold of node j
    W[j, -j] <- b[-1]
  }
  Wsym <- (W + t(W)) / 2
  attr(Wsym, "intercepts") <- b0   # set after symmetrizing
  Wsym
}

max_abs_edge <- function(W) max(abs(W[upper.tri(W)]))

# smallest lambda (least penalization) such that max|edge| <= target
find_lambda_for_target <- function(node_mat, target,
                                   lambda_lo = 1e-4, lambda_hi = 5,
                                   tol = 1e-3, max_iter = 40) {
  W_lo <- fit_ising_at_lambda(node_mat, lambda_lo)
  if (max_abs_edge(W_lo) <= target) {
    message("Target ", target, " already met at lambda_lo = ", lambda_lo,
            " (max |edge| = ", round(max_abs_edge(W_lo), 3), "); no extra penalization needed.")
    return(list(lambda = lambda_lo, max_edge = max_abs_edge(W_lo),
                W = W_lo, already_met = TRUE))
  }
  W_hi <- fit_ising_at_lambda(node_mat, lambda_hi)
  if (max_abs_edge(W_hi) > target)
    stop("Target ", target, " not met even at lambda_hi = ", lambda_hi, "; raise lambda_hi.")
  
  for (i in seq_len(max_iter)) {
    mid <- sqrt(lambda_lo * lambda_hi)
    if (max_abs_edge(fit_ising_at_lambda(node_mat, mid)) <= target) lambda_hi <- mid
    else lambda_lo <- mid
    if (log(lambda_hi / lambda_lo) < tol) break
  }
  W <- fit_ising_at_lambda(node_mat, lambda_hi)
  list(lambda = lambda_hi, max_edge = max_abs_edge(W), W = W, already_met = FALSE)
}

# IsingFit() itself only accepts {0,1}-coded input (hard stop otherwise), so
# this transform is how the {-1,1} network is obtained.

make_pair_ising <- function(e1, e2,
                            comp_belongs = emo.comp.rates.belongs,
                            dat = data) {
  
  # 1. nodes: components belonging to either emotion (unchanged)
  names.E1 <- comp_belongs$component[comp_belongs$emotion == e1]
  names.E2 <- comp_belongs$component[comp_belongs$emotion == e2]
  nodes    <- union(names.E1, names.E2)
  
  # 2. keep only cases where BOTH emotions were experienced (flag & scale > 1)
  #    -- identical selection rule to make_pair_bootnet()
  keep <- (dat[[e1]] == TRUE) & (as.numeric(dat[[paste0(e1, ".scale")]]) > 1) &
    (dat[[e2]] == TRUE) & (as.numeric(dat[[paste0(e2, ".scale")]]) > 1)
  keep[is.na(keep)] <- FALSE
  
  # 3. node matrix of BINARY component endorsements (0/1) on that subsample
  #    -- raw TRUE/FALSE flag, no intensity threshold applied
  node_mat <- vapply(
    nodes,
    function(cmp) {
      x <- dat[[cmp]][keep]
      x[is.na(x)] <- FALSE   # treat non-response as "component not endorsed"
      as.integer(x)
    },
    integer(sum(keep))
  )
  node_mat <- as.data.frame(node_mat)
  colnames(node_mat) <- nodes
  
  # 4. drop constant components (all 0 or all 1 in this subsample):
  #    IsingFit cannot estimate edges for a node without variance
  ok      <- vapply(node_mat, function(x) length(unique(x)) > 1, logical(1))
  dropped <- nodes[!ok]
  node_mat <- node_mat[, ok, drop = FALSE]

  
  # 5. walk the lasso path directly (bypassing IsingFit's EBIC selection)
  #    to find the least-penalized lambda whose max |edge| stays under target
  # under {0;1} data, these are four times as high as under {-1;1}
  res.target3.2 <- find_lambda_for_target(node_mat, target = 3.2)
  
  list(
    net.custom.target3.2 = res.target3.2$W,
    lambda.target3.2 = res.target3.2$lambda,
    max_edge.target3.2 = res.target3.2$max_edge,
    names.E1 = intersect(names.E1, colnames(node_mat)),
    names.E2 = intersect(names.E2, colnames(node_mat)),
    n_cases  = sum(keep),
    n_nodes  = ncol(node_mat),
    dropped  = dropped
  )
}

# Same top pairs as in 1.2
pairs_top <- list(
  love_sadness        = c("love",     "sadness"),
  enjoyment_joy       = c("enjoyment","joy"),
  anxiety_frustration = c("anxiety",  "frustration")
)

ising_nets <- lapply(pairs_top, function(p) make_pair_ising(p[1], p[2]))

saveRDS(ising_nets, file.path(emp.dir, "top_emotion_networks_ising.RDS"))


## ------------------------------------- ##
# 1.3 Visualize Ising Networks ---------- #
## ------------------------------------- ##

plot_pair_ising <- function(pair_fit, e1, e2, maximum = NULL,
                            vsize = 6, label.cex = 1.4, ...) {
  
  net   <- pair_fit$net_01
  nodes <- colnames(net$graph)
  
  # group nodes: unique to e1, unique to e2, or shared between both
  group <- ifelse(
    nodes %in% pair_fit$names.E1 & nodes %in% pair_fit$names.E2, "shared",
    ifelse(nodes %in% pair_fit$names.E1, e1, e2)
  )
  group <- factor(group, levels = c("shared", e1, e2))
  
  qgraph(
    net$graph,
    labels      = nodes,
    groups      = group,
    color       = c("#B3B3B3", "#66C2A5", "#FC8D62"),
    layout      = "spring",
    maximum     = maximum,   # fix a common scale across separate plots (see note below)
    title       = paste0(e1, " \u2013 ", e2,
                         " (n = ", pair_fit$n_cases, ", nodes = ", pair_fit$n_nodes, ")"),
    negDashed   = TRUE,
    vsize       = vsize,       # smaller, fixed-size nodes ...
    label.cex   = label.cex,   # ... paired with larger, fixed-size labels
    label.scale = FALSE,       # (default TRUE shrinks labels to fit inside vsize;
    #  FALSE lets them spill outside the node instead)
    ...
  )
}

# common edge-weight scale
max_edge <- max(sapply(ising_nets, function(x) max(abs(x$net_01$graph))))

# plot each pair 
plot_pair_ising(ising_nets$love_sadness,        "love",      "sadness",     maximum = max_edge)
plot_pair_ising(ising_nets$enjoyment_joy,       "enjoyment", "joy",         maximum = max_edge)
plot_pair_ising(ising_nets$anxiety_frustration, "anxiety",   "frustration", maximum = max_edge)

