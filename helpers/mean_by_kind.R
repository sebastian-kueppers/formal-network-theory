# helper: mean_by_kind()
# requires: dplyr

# mean, SD and SE of a variable by kind of goal activation (phenomenon A: var = "MIN")
mean_by_kind <- function(d, var) {
  d %>%
    group_by(kind) %>%
    summarize(mean = mean(.data[[var]]),
              sd = sd(.data[[var]]),
              se = sd(.data[[var]]) / sqrt(n()),
              n = n(), .groups = "drop")
}
