# helper: test_shared_slope()

# phenomenon F: a shared node predicts co-occurrence (MIN) more strongly than an
# ordinary (non-shared) node (p = prep_pair() output)
# For every pair of one shared and one ordinary node, one joint model per set
# (overall, E1, E2, E1E2):
#   MIN ~ shared + I(shared + ordinary)
# This is MIN ~ shared + ordinary rewritten so that the coefficient of `shared`
# is the difference B.shared - B.ordinary (dB) with its own SE; the coefficient
# of I(shared + ordinary) is B.ordinary.
#   shared     : shared node(s); default all shared nodes. In the empirical
#                pairs the caller passes only the shared APPRAISALS.
#   ordinary   : non-shared comparison node(s); needed if comparison = "given"
#   comparison : "given" = every shared node vs every node in `ordinary`
#                "all"   = every shared node vs ALL non-shared nodes
#                          (`ordinary` is ignored, `shared` is still used)
#-> p.dB: one-sided test of H0: dB <= 0 (shared node not stronger than ordinary)
#-> attr(<result>, "summary"): per shared node and set, how many comparisons
#   have dB > 0 and the range of dB (printed if comparison = "all")
test_shared_slope <- function(p, shared = p$shared, ordinary = NULL,
                              comparison = c("given", "all")) {
  comparison <- match.arg(comparison)
  d <- p$data
  non.shared <- setdiff(union(p$names.E1, p$names.E2), p$shared)
  
  if (p$n.shared == 0)
    stop(p$label, ": no shared nodes, F is not applicable.")
  if (!length(shared) || !all(shared %in% p$shared))
    stop(p$label, ": not a shared node: ",
         paste(setdiff(shared, p$shared), collapse = ", "))
  if (comparison == "all") {
    if (!is.null(ordinary)) message(p$label, ": comparison = \"all\", `ordinary` is ignored.")
    ordinary <- non.shared
  } else {
    if (!length(ordinary))
      stop(p$label, ": give `ordinary` node(s) or use comparison = \"all\".")
    if (!all(ordinary %in% non.shared))
      stop(p$label, ": not a non-shared node: ",
           paste(setdiff(ordinary, non.shared), collapse = ", "))
  }
  
  sets <- list(overall = rep(TRUE, nrow(d)),
               E1 = d$kind == "E1", E2 = d$kind == "E2", E1E2 = d$kind == "E1E2")
  pairs <- expand.grid(shared = shared, ordinary = ordinary, set = names(sets),
                       stringsAsFactors = FALSE)
  
  res <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    s  <- pairs$shared[i]
    o  <- pairs$ordinary[i]
    ds <- d[sets[[pairs$set[i]]], ]
    bs <- paste0("`", s, "`")
    bo <- paste0("`", o, "`")
    co <- coef(summary(lm(reformulate(c(bs, sprintf("I(%s + %s)", bs, bo)), "MIN"),
                          data = ds)))
    if (nrow(co) < 3) {   # a node is constant in this set -> model not estimable
      return(data.frame(shared = s, ordinary = o, set = pairs$set[i],
                        B.shared = NA, B.ordinary = NA, dB = NA, SE.dB = NA, p.dB = NA))
    }
    dB    <- co[2, ]   # B.shared - B.ordinary
    B.ord <- co[3, ]   # B.ordinary
    data.frame(shared = s, ordinary = o, set = pairs$set[i],
               B.shared = round(dB[["Estimate"]] + B.ord[["Estimate"]], 3),
               B.ordinary = round(B.ord[["Estimate"]], 3),
               dB = round(dB[["Estimate"]], 3), SE.dB = round(dB[["Std. Error"]], 3),
               p.dB = pt(dB[["t value"]], df = nrow(ds) - 3, lower.tail = FALSE))
  }))
  
  smry <- do.call(rbind, lapply(split(res, list(res$shared, res$set), drop = TRUE),
                                function(r)
    data.frame(shared = r$shared[1], set = r$set[1], n.ordinary = nrow(r),
               n.dB.pos = sum(r$dB > 0, na.rm = TRUE), dB.min = min(r$dB, na.rm = TRUE),
               dB.median = median(r$dB, na.rm = TRUE), dB.max = max(r$dB, na.rm = TRUE))))
  smry <- smry[order(match(smry$shared, shared), match(smry$set, names(sets))), ]
  rownames(smry) <- NULL
  if (comparison == "all") {
    cat("--- ", p$label, ": shared vs. all ", length(ordinary), " non-shared nodes ---\n", sep = "")
    print(smry, row.names = FALSE)
  }
  
  attr(res, "summary") <- smry
  res
}
