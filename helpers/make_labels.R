# helper: make_labels()

# make shortened and shared emotion labels
make_labels <- function(names.E1, names.E2, prefix.E1 = "L", prefix.E2 = "S", prefix.shared = "SL") {
  shared  <- intersect(names.E1, names.E2)
  only.E1 <- setdiff(names.E1, shared)
  only.E2 <- setdiff(names.E2, shared)
  
  lab <- c(
    setNames(paste0(prefix.E1,     seq_along(only.E1)), only.E1),
    setNames(paste0(prefix.E2,     seq_along(only.E2)), only.E2),
    setNames(paste0(prefix.shared, seq_along(shared)),  shared)
  )
  lab
}
