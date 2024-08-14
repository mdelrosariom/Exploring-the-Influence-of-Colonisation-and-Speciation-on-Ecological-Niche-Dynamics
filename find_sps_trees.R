setwd("C:/Users/mdrmi/Downloads/Phylogenetic trees used for global analysis of island birds/Phylogenetic trees used for global analysis of island birds/Maximum clade credibility trees/")
data_frame_names <- list.files(pattern = "*.tre")  

trees <- data.frame(FileNames = data_frame_names)

find_species_trees <- function(specie) {
  for (x in 1:dim(trees)[1]) {
    print(x)
    print(trees$FileNames[x])
    tree <- read.nexus(trees$FileNames[x])
    
    species_tree <- data.frame(species_tr = tree$tip.label)
    
    if (any(specie %in% species_tree$species_tr)) {
      print(paste("Species", specie, "found in tree", trees$FileNames[x]))
    } else {
      print(paste("Species", specie, "not found in tree", trees$FileNames[x]))
    }
  }
}




