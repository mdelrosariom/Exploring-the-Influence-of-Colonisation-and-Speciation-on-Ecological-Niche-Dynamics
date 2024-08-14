new_pruned <- function(tree_analyzed, name_of_tree) {
  #
  Y_tree <- read.nexus(name_of_tree)
  #all the clades we can prune
  to_prune_ok <- tree_analyzed$to_prune
  #is a list of lists so we need to iterate over it
  for (p in 1:length(to_prune_ok)) {
   #we need to transform it into a vector because if not it will not do anything we ask for  
   prunable_vect <- as.vector(unlist(to_prune_ok[p]))
   #the one individual of the clade that will stay in the tree
   representative_individual <- sample(prunable_vect, 1)
   # the individuals that will be removed from the tree
   kill_from_tree <- prunable_vect[prunable_vect != representative_individual]  
   #we prune the tree 
   Y_tree <- drop.tip(Y_tree,  kill_from_tree)
 
  }
  
  #and save it
  tree_name <- sub('\\.tree$', '', name_of_tree) 
  file_name <- paste("C:/Users/mdrmi/OneDrive/Documents/MEGA/trees_reunion_703", tree_name, "_pruned.tree", sep = "")
  writeNexus(Y_tree, file = file_name)
  return(Y_tree)
}
