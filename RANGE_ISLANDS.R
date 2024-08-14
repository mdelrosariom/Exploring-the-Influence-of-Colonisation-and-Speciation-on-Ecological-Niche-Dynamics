library("xlsx")
library("dplyr")
library("readxl")
library("spatstat")

# Define the function
RANGE_island <- function(archipelago, folder, value_to_pad, trait) {
  
  list_files <- list.files(
    path = paste("C:/Users/mdrmi/OneDrive/Escritorio/", archipelago, sep = ""), 
    pattern = ".xlsx", full.names = TRUE, recursive = FALSE
  )
  
  island_info <- island_complete_info(archipelago)
  island_info <- subset(island_info, tree != "Not_sampled_in_phylogeny")
  
  times <- seq(0, value_to_pad, by = 0.1)
  times_mnn <- data.frame(times = times, total = rep(NA, length(times)), colonization = rep(NA, length(times)), anagenetic = rep(NA, length(times)), cladogenetic = rep(NA, length(times)))
  
  for (time in times) {
    time_sps_mn_total <- data.frame()
    time_sps_mn_colo <- data.frame()
    time_sps_mn_ana <- data.frame()
    time_sps_mn_cla <- data.frame()
    
    for (file in list_files) {
      
      species_name <- unlist(strsplit(basename(file), split = '_'))
      species_name <- paste(species_name[1], species_name[2], sep = "_")
      print(species_name)
      if (grepl(" ", species_name)) {
        species_name <- gsub(" ", "", species_name)
      }
      sps_row <- island_info %>% filter(species == species_name)
      
      # Read the excel file
      data_sps <- read_excel(file)
      
      # Find the matching theoretical time and add data to the corresponding data frame
      for (t in 1:length(data_sps$theoretical_times)) {
        if (time == data_sps$theoretical_times[[t]]) {
          time_sps_mn_total <- rbind(time_sps_mn_total, data.frame(
            BLC_mean = data_sps$BLC_mean[t],
            BLN_mean = data_sps$BLN_mean[t],
            BW_mean = data_sps$BW_mean[t],
            BD_mean = data_sps$BD_mean[t],
            species = species_name
          ))
          
          if (sps_row$origin == "Non_endemic") {
            time_sps_mn_colo <- rbind(time_sps_mn_colo, data.frame(
              BLC_mean = data_sps$BLC_mean[t],
              BLN_mean = data_sps$BLN_mean[t],
              BW_mean = data_sps$BW_mean[t],
              BD_mean = data_sps$BD_mean[t],
              species = species_name
            ))
          }
          
          if (sps_row$origin == "Endemic") {
            if (grepl(",", sps_row$colonization) | sps_row$colonization == "Cladogenetic") {
              time_sps_mn_cla <- rbind(time_sps_mn_cla, data.frame(
                BLC_mean = data_sps$BLC_mean[t],
                BLN_mean = data_sps$BLN_mean[t],
                BW_mean = data_sps$BW_mean[t],
                BD_mean = data_sps$BD_mean[t],
                species = species_name
              ))
            } else {
              time_sps_mn_ana <- rbind(time_sps_mn_ana, data.frame(
                BLC_mean = data_sps$BLC_mean[t],
                BLN_mean = data_sps$BLN_mean[t],
                BW_mean = data_sps$BW_mean[t],
                BD_mean = data_sps$BD_mean[t],
                species = species_name
              ))
            }
          }
        }
      }
    }
    
   
    
    times_mnn$total[times == time] <- abs(max(time_sps_mn_total)-min(time_sps_mn_total))
    times_mnn$colonization[times == time] <- abs(max(time_sps_mn_colo)-min(time_sps_mn_colo))
    times_mnn$anagenetic[times == time] <- abs(max(time_sps_mn_ana)-min(time_sps_mn_ana))
    times_mnn$cladogenetic[times == time] <- abs(max(time_sps_mn_cla)-min(time_sps_mn_cla))
  }
  
  output_file <- paste0("C:/Users/mdrmi/OneDrive/Escritorio/", archipelago, "/","INFO/" ,archipelago,"_RANGE_",trait, "_ISLAND.xlsx")
  write.xlsx(times_mnn, output_file)
  return(times_mnn)
}
