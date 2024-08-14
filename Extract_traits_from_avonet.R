##extract traits from avonet raw (multiple individuals per species) and 
##calculate mean, standard deviation and variance

library("readxl")
library(xlsx)
avonet_raw = read_excel("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/AVONET Supplementary dataset 1 (2).xlsx", sheet="AVONET_Raw_Data") 
#create a new column to put the species name
#avonet_raw <- cbind(avonet_raw,"species_name") 
#for (i in 1:length(avonet_raw$Species1_BirdLife)){
 # bird_traits <- avonet_raw$Species1_BirdLife[[i]]
  
 # label_parts <- strsplit(bird_traits, "_") #SPLIT NAME OF TIPS BY THE UNDERSCORE (I KNOW THEY HAVE THIS STRUCTURE)
  
#  common_prefix <- sapply(label_parts, function(x) paste(x[1], sep = "_")) #WE EXTRACT THE SPS NAME (GENUS_SPNAME)
  #create a new column with the species name (wo the code that comes after)
#  avonet_raw$species_name[[i]]  <- common_prefix
  
#}
#species repeat because multiple individuals so we need to do a list to know which sps. are
unique_species<- unique(avonet_raw$Species1_BirdLife)


traits <- data.frame()  # Initialize an empty data frame

for (t in 1:length(unique_species)) {
  # Filter rows where first_column is equal to the current unique value
  #filter does not work
  #filtered_avonet_raw <- filter(avonet_raw, species_name == unique_species[2])
  #supposely %<% more tolerand and else 
  filtered_avonet_raw <- avonet_raw[avonet_raw$Species1_BirdLife %in% unique_species[t], ]

  # Calculate mean, sd, and var for beak length culmen
  mean_BLC <- mean(as.numeric(filtered_avonet_raw$Beak.Length_Culmen),na.rm=TRUE)
  sd_BLC <- sd(as.numeric(filtered_avonet_raw$Beak.Length_Culmen),na.rm=TRUE)
  var_BLC <- var(as.numeric(filtered_avonet_raw$Beak.Length_Culmen),na.rm=TRUE)
  
  #same for bl nares 
  mean_BLN <- mean(as.numeric(filtered_avonet_raw$Beak.Length_Nares),na.rm=TRUE)
  sd_BLN <- sd(as.numeric(filtered_avonet_raw$Beak.Length_Nares),na.rm=TRUE)
  var_BLN <- var(as.numeric(filtered_avonet_raw$Beak.Length_Nares),na.rm=TRUE)
  
  #same for "Beak.Width"
  
  mean_BW <- mean(as.numeric(filtered_avonet_raw$Beak.Width),na.rm=TRUE)
  sd_BW <- sd(as.numeric(filtered_avonet_raw$Beak.Width),na.rm=TRUE)
  var_BW <- var(as.numeric(filtered_avonet_raw$Beak.Width),na.rm=TRUE)
  
  #same for Beak.Depth
  
  mean_BD <- mean(as.numeric(filtered_avonet_raw$Beak.Depth),na.rm=TRUE)
  sd_BD <- sd(as.numeric(filtered_avonet_raw$Beak.Depth),na.rm=TRUE)
  var_BD <- var(as.numeric(filtered_avonet_raw$Beak.Depth),na.rm=TRUE)
  
  parts_name  <- strsplit(unique_species[t], " ") 
  unique_species[t] <- paste(parts_name[[1]][1], "_", parts_name[[1]][2], sep="")
  
  # Create a data frame for the current unique value
  current_row <- data.frame(
    specie = unique_species[t],
    #beak length culmen
    BLC_mean = mean_BLC, 
    BLC_var = var_BLC, 
    BLC_sd = sd_BLC, 
    #beak_length_narses
    BLN_mean = mean_BLN,
    BLN_var = var_BLN,
    BLN_sd = sd_BLN,
    #beak_width
    BW_mean = mean_BW, 
    BW_var = var_BW, 
    BW_sd = sd_BW,
    #beak depth 
    BD_mean = mean_BD, 
    BD_var = var_BD, 
    BD_sd = sd_BD)
  
  # Append the data frame to traits using rbind
  traits <- rbind(traits, current_row)
  
  
}



write.xlsx(
  traits, "C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/traits_purified_good_version.xlsx"
)







