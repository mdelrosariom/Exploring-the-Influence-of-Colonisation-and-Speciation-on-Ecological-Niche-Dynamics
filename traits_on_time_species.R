library("phytools")
library("phylobase")

source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/tree_branch_finder.R")


traits_time_species <- function(tree, species, res_species, colonization_time,tree_name){

tree_path <- as(tree, "phylo4d")

evo_path_species <- ancestors(tree_path, species, "ALL")

evo_path_species <- as.vector(evo_path_species)
evo_path_species <- rev(evo_path_species)

nodes<-  cbind(evo_path_species[1:length(evo_path_species )-1], evo_path_species[2:length(evo_path_species )])

branches <- list()

for (n in 1:nrow(nodes)){
  
  first_node <- nodes[n,1]
  second_node <- nodes[n,2]
  
  branch <- find_branch_number(first_node, second_node, tree$edge)
  
  branches <- c(branches,branch)
  
}

branches <-as.vector(unlist(branches))



#now we make the data of the times and traits for the branches involved with the specie

#first time
times <- list()

for (i in branches){
  
  times <- c(times,unlist(res_species$times[[i]]))
}
times_species <- as.vector(unlist(times))
times_species <- round(times_species, 6)

#second traits 

traits <- list()

for (i in branches) {
  print(i)
  traits <- c(traits, list(res_species$traits[[i]]))
  print(traits)
}

# Combine the subsets into a single dataframe
species_traits <- do.call(rbind, traits)

#now we make a vector that will check the traits each certain time steps. in this case 0.1 time steps (meaning 100.000 years)
#modification 15-04, until the colonization time of the species on the island

#sampling_times <- seq(0, max(times_species), by = 0.1)
sampling_times <- seq(0, colonization_time, by = 0.1)

#as the values of times will not be exactly the same as the sampling time, we need to get the closest of times_species from 
#sampling_times
closest_indices <- sapply(sampling_times, function(x) which.min(abs(times_species - x)))

closest_times <- times_species[closest_indices]

traits_in_time_species <- list()

for (time in closest_times){
  index <- which(times_species == time)
  
  
  trait_time <- species_traits[index,]
  
  traits_in_time_species <- c(traits_in_time_species, list(trait_time))
}

combined_time_traits <- do.call(rbind, traits_in_time_species)
combined_time_traits <- as.data.frame(combined_time_traits)
combined_time_traits <- as.data.frame(combined_time_traits)
colnames(combined_time_traits) <- c("BLC_mean", "BLN_mean", "BW_mean", "BD_mean")
time_traits_species <-data.frame("BLC_mean" = combined_time_traits$BLC_mean, 
                                             "BLN_mean" = combined_time_traits$BLN_mean, 
                                             "BW_mean" = combined_time_traits$BW_mean,
                                             "BD_mean" = combined_time_traits$BD_mean,
                                             "time_tree" = closest_times,
                                             "theoretical_times" = sampling_times )


return(time_traits_species)
}

