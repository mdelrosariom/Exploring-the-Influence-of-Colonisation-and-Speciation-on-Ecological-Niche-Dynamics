library("readxl")
library("writexl")
##function to get only the specific traits of the species included on the tree
get_traits_tree <- function(tree, name_of_tree_2) { 
  #import the big ds that contains the data of all birds (more than 300)
  species_traits <- read_excel("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/traits_purified_good_version.xlsx") 
  #name of the species
  species <- tree$tip.label
  traits_tree <- list()
  #if the species is in avonet, then add their data to the data of the tree....
  for (i in 1:length(species)) {
    if (species[i] %in% species_traits$specie) {
      current_row <- species_traits[species_traits$specie == species[i], , drop = FALSE]
      traits_tree[[i]] <- current_row
    } else {#... if not add NA
      current_row <- data.frame(
        specie_sep = species[[i]],
        BLC_mean = NA, 
        BLC_var = NA, 
        BLC_sd = NA, 
        BLN_mean = NA,
        BLN_var = NA,
        BLN_sd = NA,
        BW_mean = NA, 
        BW_var = NA, 
        BW_sd = NA,
        BD_mean = NA, 
        BD_var = NA, 
        BD_sd = NA, 
        specie = species[[i]]
      )
      traits_tree[[i]] <- current_row #create dataframe
    }
  }
  
  tree_name <- sub('\\.tree$', '', name_of_tree_2)
  #export 
  file_name <- paste("C:/Users/mdrmi/OneDrive/Escritorio/tree_trait/","trait_",tree_name ,".xlsx")
  
  # Combine the list of data frames into a single data frame
  traits_tree <- do.call(rbind, traits_tree)
  write_xlsx(traits_tree, file_name)
  #return 
  return(traits_tree)
} #<3


