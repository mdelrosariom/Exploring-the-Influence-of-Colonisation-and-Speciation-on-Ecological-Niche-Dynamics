library(stringr) # Make sure to load the stringr package

change_tips <- function(X_tree) {
  for (i in 1:length(X_tree$tip.label)) {
    x <- X_tree$tip.label[i]
    if (grepl("_", x)) {
      X_tree$tip.label[i] <- str_extract(x, "[^_]*_[^_]*")
    }
  }
  return(X_tree)
}



