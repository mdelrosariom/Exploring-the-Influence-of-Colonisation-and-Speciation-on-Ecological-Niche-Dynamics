library(dplyr)
library(readxl)
library(openxlsx)
source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/island_complete_info.R")
# Define the custom range function
custom_range <- function(time_sps_mn, trait) {
  if (nrow(time_sps_mn) <= 1) {
    return(0)
  } else {
    trait_range <- abs(max(time_sps_mn[[trait]], na.rm = TRUE) - min(time_sps_mn[[trait]], na.rm = TRUE))
    return(trait_range)
  }
}

# Define the main function
RANGE_island <- function(archipelago, folder, value_to_pad) {
  path <- file.path("C:/Users/mdrmi/OneDrive/Escritorio", archipelago)
  list_files <- list.files(
    path = path, 
    pattern = "\\.xlsx$", full.names = TRUE, recursive = FALSE
  )
  
  island_info <- island_complete_info(archipelago)
  island_info <- subset(island_info, tree != "Not_sampled_in_phylogeny")
  
  times <- seq(0, value_to_pad, by = 0.1)
  times_mnn <- data.frame(
    times = times, 
    total_BLC_mean = NA, 
    total_BLN_mean = NA, 
    total_BW_mean = NA, 
    total_BD_mean = NA, 
    colo_BLC_mean = NA, 
    colo_BLN_mean = NA, 
    colo_BW_mean = NA, 
    colo_BD_mean = NA, 
    ana_BLC_mean = NA, 
    ana_BLN_mean = NA, 
    ana_BW_mean = NA, 
    ana_BD_mean = NA, 
    cla_BLC_mean = NA, 
    cla_BLN_mean = NA, 
    cla_BW_mean = NA, 
    cla_BD_mean = NA, 
    speci_BLC_mean = NA, 
    speci_BLN_mean = NA, 
    speci_BW_mean = NA, 
    speci_BD_mean = NA
  )
  
  for (time in times) {
    time_sps_mn_total <- data.frame()
    time_sps_mn_colo <- data.frame()
    time_sps_mn_ana <- data.frame()
    time_sps_mn_cla <- data.frame()
    time_sps_mn_speci <- data.frame()
    
    for (file in list_files) {
      species_name <- unlist(strsplit(basename(file), split = '_'))
      species_name <- paste(species_name[1], species_name[2], sep = "_")
      if (grepl(" ", species_name)) {
        species_name <- gsub(" ", "", species_name)
      }
      
      sps_row <- island_info %>% filter(species == species_name)
      data_sps <- tryCatch(read_excel(file), error = function(e) NULL)
      if (is.null(data_sps)) next
      
      for (t in seq_along(data_sps$theoretical_times)) {
        if (time == data_sps$theoretical_times[t]) {
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
            time_sps_mn_speci <- rbind(time_sps_mn_speci, data.frame(
              BLC_mean = data_sps$BLC_mean[t],
              BLN_mean = data_sps$BLN_mean[t],
              BW_mean = data_sps$BW_mean[t],
              BD_mean = data_sps$BD_mean[t],
              species = species_name
            ))
            
            if (grepl(",", sps_row$colonization) || sps_row$colonization == "Cladogenetic") {
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
    
    index <- which(times == time)
    times_mnn$total_BLC_mean[index] <- custom_range(time_sps_mn_total, "BLC_mean")
    times_mnn$total_BLN_mean[index] <- custom_range(time_sps_mn_total, "BLN_mean")
    times_mnn$total_BW_mean[index] <- custom_range(time_sps_mn_total, "BW_mean")
    times_mnn$total_BD_mean[index] <- custom_range(time_sps_mn_total, "BD_mean")
    
    times_mnn$colo_BLC_mean[index] <- custom_range(time_sps_mn_colo, "BLC_mean")
    times_mnn$colo_BLN_mean[index] <- custom_range(time_sps_mn_colo, "BLN_mean")
    times_mnn$colo_BW_mean[index] <- custom_range(time_sps_mn_colo, "BW_mean")
    times_mnn$colo_BD_mean[index] <- custom_range(time_sps_mn_colo, "BD_mean")
    
    times_mnn$ana_BLC_mean[index] <- custom_range(time_sps_mn_ana, "BLC_mean")
    times_mnn$ana_BLN_mean[index] <- custom_range(time_sps_mn_ana, "BLN_mean")
    times_mnn$ana_BW_mean[index] <- custom_range(time_sps_mn_ana, "BW_mean")
    times_mnn$ana_BD_mean[index] <- custom_range(time_sps_mn_ana, "BD_mean")
    
    times_mnn$cla_BLC_mean[index] <- custom_range(time_sps_mn_cla, "BLC_mean")
    times_mnn$cla_BLN_mean[index] <- custom_range(time_sps_mn_cla, "BLN_mean")
    times_mnn$cla_BW_mean[index] <- custom_range(time_sps_mn_cla, "BW_mean")
    times_mnn$cla_BD_mean[index] <- custom_range(time_sps_mn_cla, "BD_mean")
    
    times_mnn$speci_BLC_mean[index] <- custom_range(time_sps_mn_speci, "BLC_mean")
    times_mnn$speci_BLN_mean[index] <- custom_range(time_sps_mn_speci, "BLN_mean")
    times_mnn$speci_BW_mean[index] <- custom_range(time_sps_mn_speci, "BW_mean")
    times_mnn$speci_BD_mean[index] <- custom_range(time_sps_mn_speci, "BD_mean")
  }
  
  output_file <- file.path("C:/Users/mdrmi/OneDrive/Escritorio", archipelago, "INFO", paste0(archipelago, "_RANGE_ISLAND.xlsx"))
  write.xlsx(times_mnn, output_file)
  return(times_mnn)
}
