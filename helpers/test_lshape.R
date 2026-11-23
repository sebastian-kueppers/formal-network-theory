# helper: test_lshape()
# requires: ggplot2; scatter_states()

# phenomenon C: the negative correlation of mixed emotions is L-shaped
#   E1 arm (kind E1)   : stateE2 regressed on stateE1 (E1 predicts E2)
#   E2 arm (kind E2)   : stateE1 regressed on stateE2 (E2 predicts E1)
#   overall (all kinds): the correlation of stateE1 and stateE2 must be negative
# Prints both regressions and the overall correlation. Returns scatter_states(p)
# with the two arm regression lines (solid) and, if overall.line = TRUE, the
# overall SD line (dashed): slope sign(r) * SD(E2) / SD(E1) through the means,
# i.e. it shows the sign of r and treats both emotions symmetrically.
# The results table is attached as attr(<plot>, "lshape").
#   ... : further arguments to scatter_states(), e.g. title or
#         emotions = c("love", "sadness") for axis titles and legend
test_lshape <- function(p, overall.line = TRUE, ...) {
  d <- p$data
  
  arm <- function(k, outcome, predictor) {
    fit <- lm(reformulate(predictor, outcome), data = d[d$kind == k, ])
    co  <- coef(summary(fit))[predictor, ]
    ci  <- confint(fit)[predictor, ]
    data.frame(set = k, model = paste(outcome, "~", predictor), n = nobs(fit),
               estimate = co[["Estimate"]], se = co[["Std. Error"]],
               ci.l = ci[[1]], ci.u = ci[[2]], p = co[["Pr(>|t|)"]])
  }
  
  ct <- cor.test(d$stateE1, d$stateE2)
  res <- rbind(
    arm("E1", "stateE2", "stateE1"),
    arm("E2", "stateE1", "stateE2"),
    data.frame(set = "overall", model = "r(stateE1, stateE2)", n = nrow(d),
               estimate = unname(ct$estimate), se = NA,
               ci.l = ct$conf.int[1], ci.u = ct$conf.int[2], p = ct$p.value)
  )
  res[c("estimate", "se", "ci.l", "ci.u")] <-
    round(res[c("estimate", "se", "ci.l", "ci.u")], 3)
  
  cat("--- ", p$label, ": L-shape (C) ---\n", sep = "")
  print(res, row.names = FALSE)
  cat("Overall correlation negative (r < 0, p < .05):",
      unname(ct$estimate) < 0 && ct$p.value < .05, "\n")
  
  # E1 arm: E2 regressed on E1; E2 arm: E1 regressed on E2 (orientation = "y"),
  # so in each arm the activated emotion is the predictor
  g <- scatter_states(p, ...) +
    geom_smooth(data = d[d$kind == "E1", ],
                aes(stateE1, stateE2, colour = kind),
                method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.2) +
    geom_smooth(data = d[d$kind == "E2", ],
                aes(stateE1, stateE2, colour = kind), orientation = "y",
                method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.2)
  # caption <- "Solid: E2 ~ E1 in E1-only cases, E1 ~ E2 in E2-only cases"
  
  if (overall.line) {
    sd.slope <- sign(cor(d$stateE1, d$stateE2)) * sd(d$stateE2) / sd(d$stateE1)
    g <- g +
      geom_abline(slope = sd.slope,
                  intercept = mean(d$stateE2) - sd.slope * mean(d$stateE1),
                  colour = "grey20", linetype = "dashed", linewidth = 1.2)
    # caption <- paste0(caption, "; dashed: SD line, all cases")
  }
  
  # g <- g + labs(caption = caption)
  attr(g, "lshape") <- res
  g
}
