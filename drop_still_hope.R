library(phytools)
drop_still_hope <- function (vector_still_hope,tree){ 
  representative_individual <- sample(vector_still_hope, 1)
  # the individuals that will be removed from the tree
  kill_from_tree <- vector_still_hope[vector_still_hope != representative_individual]  
  #we prune the tree 
  tree <- drop.tip(tree,  kill_from_tree)
  old_name <- representative_individual
  new_name <- str_replace_all(old_name, "_", ".")
 
  tree$tip.label[tree$tip.label==old_name] <- new_name
  return(tree)
  }