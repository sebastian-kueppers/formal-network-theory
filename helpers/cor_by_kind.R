# helper: cor_by_kind()

# correlation of the two emotion states, overall and by kind, with 95% CI and p
# (phenomena B and D)
cor_by_kind <- function(d) {
  sets <- list(overall = rep(TRUE, nrow(d)),
               E1 = d$kind == "E1", E2 = d$kind == "E2", E1E2 = d$kind == "E1E2")
  do.call(rbind, lapply(names(sets), function(s) {
    ds <- d[sets[[s]], ]
    ct <- cor.test(ds$stateE1, ds$stateE2)
    data.frame(set = s, n = nrow(ds), r = round(unname(ct$estimate), 3),
               ci.l = round(ct$conf.int[1], 3), ci.u = round(ct$conf.int[2], 3),
               p = ct$p.value)
  }))
}
