# Conditional Permutation -------------------------------------------------

# Conditionally permutes a listw object
permute_listw <- function(listw) {
  n <- length(listw$neighbours)

  cards <- lengths(listw$neighbours)
  # Shuffle the neighbors by randomly sampling from all possible neighbors
  # except where x exists
  perm_nb <- mapply(shuffle_nbs, 1:n, n, cards, SIMPLIFY = FALSE)
  class(perm_nb) <- c("nb", "list")
  listw$neighbours <- perm_nb

  listw
}

# Internal function to shuffle neighbors
#
# Used in conditional permutation and the function `permute_listw()`.
# i is the index position of observation `i`
# n is the length of neighbor list
# card is the cardinality for observation i
shuffle_nbs <- function(i, n, card) {
  x <- 1:n
  sample(x[-i], size = card)
}
