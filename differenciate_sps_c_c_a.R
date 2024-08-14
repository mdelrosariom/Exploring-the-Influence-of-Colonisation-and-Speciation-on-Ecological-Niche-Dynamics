library("xlsx")
library("umx")
library("matrixStats")
library("dplyr")
library("openxlsx")
library("readxl") 

source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/island_complete_info.R")
source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/matrix_num.R")




var_range_sps <- function(archipelago, folder,value_to_pad ){
  #archipelago we want to analyze
  #foler in which we want to export
  #earlier colonization time (or maximum in colonization column, its too messy to try to code that)
  
  
  list_files <- list.files(path = paste("C:/Users/mdrmi/OneDrive/Escritorio/", archipelago, sep = ""), 
                           pattern = ".xlsx", all.files = FALSE,
                           full.names = TRUE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  island_info <- island_complete_info(archipelago)
  island_info <- subset(island_info, tree != "Not_sampled_in_phylogeny")
  
  tot_BLC <- list()
  tot_BLN <- list()
  tot_BW <- list()
  tot_BD <- list()
  
  colonization_BLC <- list()
  colonization_BLN <- list()
  colonization_BW <- list()
  colonization_BD <- list()
  
  cladogenesis_BLC <- list()
  cladogenesis_BLN <- list()
  cladogenesis_BW <- list()
  cladogenesis_BD <- list()
  
  anagenesis_BLC <- list()
  anagenesis_BLN <- list()
  anagenesis_BW <- list()
  anagenesis_BD <- list()
  
  #----- ok 
  
  for (file in list_files) {
    species_name <- unlist(strsplit(file, split = '/'))
    species_name <- tail(species_name, n = 1)
    species_name <- unlist(strsplit(species_name, split = '_'))
    species_name <- paste(species_name[1], species_name[2], sep = "_") # name of the species
    if (grepl(" ", species_name)) {
      species_name <- gsub(" ", "", species_name)
    } 
    data <- read_excel(file)  # data species
    data <- umx_pad(data, value_to_pad)
    
    # info island for species
    sps_row <- island_info %>% filter(species == species_name)
    
    # Append data to appropriate lists based on conditions
    tot_BLC <- c(tot_BLC, list(data$BLC_mean))
    tot_BLN <- c(tot_BLN, list(data$BLN_mean))
    tot_BW <- c(tot_BW, list(data$BW_mean))
    tot_BD <- c(tot_BD, list(data$BD_mean))
    
    if (sps_row$origin == "Non_endemic") {
     
      colonization_BLC <- c(colonization_BLC, list(data$BLC_mean))
      colonization_BLN <- c(colonization_BLN, list(data$BLN_mean))
      colonization_BW <- c(colonization_BW, list(data$BW_mean))
      colonization_BD <- c(colonization_BD, list(data$BD_mean))
    }
    if (sps_row$origin == "Endemic") {
      if (grepl(",", sps_row$colonization) | sps_row$colonization == "Cladogenetic") {
        cladogenesis_BLC <- c(cladogenesis_BLC, list(data$BLC_mean))
        cladogenesis_BLN <- c(cladogenesis_BLN, list(data$BLN_mean))
        cladogenesis_BW <- c(cladogenesis_BW, list(data$BW_mean))
        cladogenesis_BD <- c(cladogenesis_BD, list(data$BD_mean))
      } else {
        anagenesis_BLC <- c(anagenesis_BLC, list(data$BLC_mean))
        anagenesis_BLN <- c(anagenesis_BLN, list(data$BLN_mean))
        anagenesis_BW <- c(anagenesis_BW, list(data$BW_mean))
        anagenesis_BD <- c(anagenesis_BD, list(data$BD_mean))
      }
    }
  }
  # Convert lists to matrices and apply numeric transformation
  tot_BLC <- apply(do.call(cbind, tot_BLC), c(1, 2), matrix_to_numeric)
  tot_BLN <- apply(do.call(cbind, tot_BLN), c(1, 2), matrix_to_numeric)
  tot_BW <- apply(do.call(cbind, tot_BW), c(1, 2), matrix_to_numeric)
  tot_BD <- apply(do.call(cbind, tot_BD), c(1, 2), matrix_to_numeric)
  
  # For tot_BLC_var
  tot_BLC_var <- rowVars(as.matrix(tot_BLC), na.rm=TRUE)
  tot_BLC_var[is.na(tot_BLC_var)] <- 0
  tot_BLC_var <- rev(tot_BLC_var)
  
  # For tot_BLN_var
  tot_BLN_var <- rowVars(tot_BLN, na.rm=TRUE)
  tot_BLN_var[is.na(tot_BLN_var)] <- 0
  tot_BLN_var <- rev(tot_BLN_var)
  
  # For tot_BW_var
  tot_BW_var <- rowVars(tot_BW, na.rm=TRUE)
  tot_BW_var[is.na(tot_BW_var)] <- 0
  tot_BW_var <- rev(tot_BW_var)
  
  # For tot_BD_var
  tot_BD_var <- rowVars(tot_BD, na.rm=TRUE)
  tot_BD_var[is.na(tot_BD_var)] <- 0
  tot_BD_var <- rev(tot_BD_var)
  
  #############################----------------------------NOW WE SUM ALL VARIANCES!!!!!!!!!!!!!!!!!!!
  #############################---------------------------
  
  tot_variance_island <- cbind(tot_BLC_var,tot_BLN_var,tot_BW_var,tot_BD_var)
  tot_variance_island <- rowSums(tot_variance_island)
  tot_variance_island <- as.data.frame(tot_variance_island)
  write.xlsx(tot_variance_island, file = paste(folder, "TOT_VAR_", archipelago, ".xlsx", sep = ""))
  
  
  
  
  range_tot_BLC <- rowRanges(tot_BLC, na.rm=TRUE)
  range_tot_BLC <- range_tot_BLC[,2]-range_tot_BLC[,1]
  range_tot_BLC <- rev(range_tot_BLC) 
  
  range_tot_BLN <- rowRanges(tot_BLN, na.rm=TRUE)
  range_tot_BLN <- range_tot_BLN[,2]-range_tot_BLN[,1]
  range_tot_BLN <- rev(range_tot_BLN) 
  
  range_tot_BW <- rowRanges(tot_BW, na.rm=TRUE)
  range_tot_BW <- range_tot_BW[,2]-range_tot_BW[,1]
  range_tot_BW <- rev(range_tot_BW) 
  
  range_tot_BD <- rowRanges(tot_BD, na.rm=TRUE)
  range_tot_BD <- range_tot_BD[,2]-range_tot_BD[,1]
  range_tot_BD <- rev(range_tot_BD) 
  
  
  if (length(colonization_BD)>0){
  colonization_BLC <- apply(do.call(cbind, as.list(colonization_BLC)), c(1, 2), matrix_to_numeric)
  colonization_BLN <- apply(do.call(cbind, as.list(colonization_BLN)), c(1, 2), matrix_to_numeric)
  colonization_BW <- apply(do.call(cbind, as.list(colonization_BW)), c(1, 2), matrix_to_numeric)
  colonization_BD <- apply(do.call(cbind, as.list(colonization_BD)), c(1, 2), matrix_to_numeric)
  # For colonization variables
  col_BLC_var <- rowVars(colonization_BLC, na.rm=TRUE)
  col_BLC_var[is.na(col_BLC_var)] <- 0
  col_BLC_var <- rev(col_BLC_var)
  
  col_BLN_var <- rowVars(colonization_BLN, na.rm=TRUE)
  col_BLN_var[is.na(col_BLN_var)] <- 0
  col_BLN_var <- rev(col_BLN_var)
  
  col_BW_var <- rowVars(colonization_BW, na.rm=TRUE)
  col_BW_var[is.na(col_BW_var)] <- 0
  col_BW_var <- rev(col_BW_var)
  
  col_BD_var <- rowVars(colonization_BD, na.rm=TRUE)
  col_BD_var[is.na(col_BD_var)] <- 0
  col_BD_var <- rev(col_BD_var)
  
  #-------------------------------------------------
  #-------------------------------------------------
  col_variance_island <- cbind(col_BLC_var,col_BLN_var,col_BW_var,col_BD_var)
  col_variance_island <- rowSums(col_variance_island)
  col_variance_island <- as.data.frame(col_variance_island)
  write.xlsx(col_variance_island, file = paste(folder, "COL_VAR_", archipelago, ".xlsx", sep = ""))
  
  
  # For range of colonization variables
  range_colo_BLC <- rowRanges(colonization_BLC, na.rm=TRUE)
  range_colo_BLC <- range_colo_BLC[,2] - range_colo_BLC[,1]
  range_colo_BLC <- rev(range_colo_BLC)
  
  range_colo_BLN <- rowRanges(colonization_BLN, na.rm=TRUE)
  range_colo_BLN <- range_colo_BLN[,2] - range_colo_BLN[,1]
  range_colo_BLN <- rev(range_colo_BLN)
  
  range_colo_BW <- rowRanges(colonization_BW, na.rm=TRUE)
  range_colo_BW <- range_colo_BW[,2] - range_colo_BW[,1]
  range_colo_BW <- rev(range_colo_BW)
  
  range_colo_BD <- rowRanges(colonization_BD, na.rm=TRUE)
  range_colo_BD <- range_colo_BD[,2] - range_colo_BD[,1]
  range_colo_BD <- rev(range_colo_BD)
  
  col_BLC <- as.data.frame(cbind(col_BLC_var, range_colo_BLC ))
  col_BLN <- as.data.frame(cbind(col_BLN_var, range_colo_BLN ))
  col_BW <- as.data.frame(cbind(col_BW_var, range_colo_BW ))
  col_BD <- as.data.frame(cbind(col_BD_var, range_colo_BD ))
  
  col_BLC <- as.data.frame(cbind(col_BLC_var, range_colo_BLC))
  col_BLN <- as.data.frame(cbind(col_BLN_var, range_colo_BLN))
  col_BW <- as.data.frame(cbind(col_BW_var, range_colo_BW))
  col_BD <- as.data.frame(cbind(col_BD_var, range_colo_BD))
  
  
  write.xlsx(col_BLC, file = paste(folder, "col_BLC_", archipelago, ".xlsx", sep = ""))
  write.xlsx(col_BLN, file = paste(folder, "col_BLN_", archipelago, ".xlsx", sep = ""))
  write.xlsx(col_BW, file = paste(folder, "col_BW_", archipelago, ".xlsx", sep = ""))
  write.xlsx(col_BD, file = paste(folder, "col_BD_", archipelago, ".xlsx", sep = ""))
  
  
  
  }
  
  if (length(cladogenesis_BD)>0){
  cladogenesis_BLC <- apply(do.call(cbind, as.list(cladogenesis_BLC)), c(1, 2), matrix_to_numeric)
  cladogenesis_BLN <- apply(do.call(cbind, as.list(cladogenesis_BLN)), c(1, 2), matrix_to_numeric)
  cladogenesis_BW <- apply(do.call(cbind, as.list(cladogenesis_BW)), c(1, 2), matrix_to_numeric)
  cladogenesis_BD <- apply(do.call(cbind, as.list(cladogenesis_BD)), c(1, 2), matrix_to_numeric)
  # For cladogenesis variables
  cla_BLC_var <- rowVars(cladogenesis_BLC, na.rm=TRUE)
  cla_BLC_var[is.na(cla_BLC_var)] <- 0
  cla_BLC_var <- rev(cla_BLC_var)
  
  cla_BLN_var <- rowVars(cladogenesis_BLN, na.rm=TRUE)
  cla_BLN_var[is.na(cla_BLN_var)] <- 0
  cla_BLN_var <- rev(cla_BLN_var)
  
  cla_BW_var <- rowVars(cladogenesis_BW, na.rm=TRUE)
  cla_BW_var[is.na(cla_BW_var)] <- 0
  cla_BW_var <- rev(cla_BW_var)
  
  cla_BD_var <- rowVars(cladogenesis_BD, na.rm=TRUE)
  cla_BD_var[is.na(cla_BD_var)] <- 0
  cla_BD_var <- rev(cla_BD_var)
  
  #-----------------------------------------------------------
  #-----------------------------------------------------------
  
  cla_variance_island <- cbind(cla_BLC_var,cla_BLN_var,cla_BD_var,cla_BW_var)
  cla_variance_island <- rowSums(cla_variance_island)
  cla_variance_island <- as.data.frame(cla_variance_island)
  write.xlsx(cla_variance_island, file = paste(folder, "CLA_VAR_", archipelago, ".xlsx", sep = ""))
  
  
  
  # For range of cladogenesis variables
  range_cla_BLC <- rowRanges(cladogenesis_BLC, na.rm=TRUE)
  range_cla_BLC <- range_cla_BLC[,2] - range_cla_BLC[,1]
  range_cla_BLC <- rev(range_cla_BLC)
  
  range_cla_BLN <- rowRanges(cladogenesis_BLN, na.rm=TRUE)
  range_cla_BLN <- range_cla_BLN[,2] - range_cla_BLN[,1]
  range_cla_BLN <- rev(range_cla_BLN)
  
  range_cla_BW <- rowRanges(cladogenesis_BW, na.rm=TRUE)
  range_cla_BW <- range_cla_BW[,2] - range_cla_BW[,1]
  range_cla_BW <- rev(range_cla_BW)
  
  range_cla_BD <- rowRanges(cladogenesis_BD, na.rm=TRUE)
  range_cla_BD <- range_cla_BD[,2] - range_cla_BD[,1]
  range_cla_BD <- rev(range_cla_BD)
  
  cla_BLC <-   as.data.frame(cbind(cla_BLC_var,range_cla_BLC ))
  cla_BLN <-   as.data.frame(cbind(cla_BLN_var,range_cla_BLN ))  
  cla_BW <-   as.data.frame(cbind(cla_BW_var,range_cla_BW ))
  cla_BD <-   as.data.frame(cbind(cla_BD_var,range_cla_BD ))
  
  cla_BLC <- as.data.frame(cbind(cla_BLC_var, range_cla_BLC))
  cla_BLN <- as.data.frame(cbind(cla_BLN_var, range_cla_BLN))
  cla_BW <- as.data.frame(cbind(cla_BW_var, range_cla_BW))
  cla_BD <- as.data.frame(cbind(cla_BD_var, range_cla_BD))
  
  write.xlsx(cla_BLC, file = paste(folder, "var_range_cla_BLC_", archipelago, ".xlsx", sep = ""))
  write.xlsx(cla_BLN, file = paste(folder, "var_range_cla_BLN_", archipelago, ".xlsx", sep = ""))
  write.xlsx(cla_BW, file = paste(folder, "var_range_cla_BW_", archipelago, ".xlsx", sep = ""))
  write.xlsx(cla_BD, file = paste(folder, "var_range_cla_BD_", archipelago, ".xlsx", sep = ""))
  
  
  }
  
  if (length(anagenesis_BD)>0){
  anagenesis_BLC <- apply(do.call(cbind, as.list(anagenesis_BLC)), c(1, 2), matrix_to_numeric)
  anagenesis_BLN <- apply(do.call(cbind, as.list(anagenesis_BLN)), c(1, 2), matrix_to_numeric)
  anagenesis_BW <- apply(do.call(cbind, as.list(anagenesis_BW)), c(1, 2), matrix_to_numeric)
  anagenesis_BD <- apply(do.call(cbind, as.list(anagenesis_BD)), c(1, 2), matrix_to_numeric)
  # For anagenesis variables
  ana_BLC_var <- rowVars(anagenesis_BLC, na.rm=TRUE)
  ana_BLC_var[is.na(ana_BLC_var)] <- 0
  ana_BLC_var <- rev(ana_BLC_var)
  
  ana_BLN_var <- rowVars(anagenesis_BLN, na.rm=TRUE)
  ana_BLN_var[is.na(ana_BLN_var)] <- 0
  ana_BLN_var <- rev(ana_BLN_var)
  
  ana_BW_var <- rowVars(anagenesis_BW, na.rm=TRUE)
  ana_BW_var[is.na(ana_BW_var)] <- 0
  ana_BW_var <- rev(ana_BW_var)
  
  ana_BD_var <- rowVars(anagenesis_BD, na.rm=TRUE)
  ana_BD_var[is.na(ana_BD_var)] <- 0
  ana_BD_var <- rev(ana_BD_var)
  #-----------------------------------------------------------
  #-----------------------------------------------------------
  
  ana_variance_island <- cbind(ana_BLC_var ,ana_BLN_var ,ana_BW_var ,ana_BD_var )
  ana_variance_island <- rowSums(ana_variance_island)
  ana_variance_island <- as.data.frame(ana_variance_island)
  write.xlsx(ana_variance_island, file = paste(folder, "ANA_VAR_", archipelago, ".xlsx", sep = ""))
  
  
  
  
  
  # For range of anagenesis variables
  range_ana_BLC <- rowRanges(anagenesis_BLC, na.rm=TRUE)
  range_ana_BLC <- range_ana_BLC[,2] - range_ana_BLC[,1]
  range_ana_BLC <- rev(range_ana_BLC)
  
  range_ana_BLN <- rowRanges(anagenesis_BLN, na.rm=TRUE)
  range_ana_BLN <- range_ana_BLN[,2] - range_ana_BLN[,1]
  range_ana_BLN <- rev(range_ana_BLN)
  
  range_ana_BW <- rowRanges(anagenesis_BW, na.rm=TRUE)
  range_ana_BW <- range_ana_BW[,2] - range_ana_BW[,1]
  range_ana_BW <- rev(range_ana_BW)
  
  range_ana_BD <- rowRanges(anagenesis_BD, na.rm=TRUE)
  range_ana_BD <- range_ana_BD[,2] - range_ana_BD[,1]
  range_ana_BD <- rev(range_ana_BD)
  
  ana_BLC <- as.data.frame(cbind(ana_BLC_var,range_ana_BLC ))
  ana_BLN <- as.data.frame(cbind(ana_BLN_var,range_ana_BLN ))
  ana_BW <- as.data.frame(cbind(ana_BW_var,range_ana_BW ))
  ana_BD <- as.data.frame(cbind(ana_BD_var,range_ana_BD ))
  
  ana_BLC <- as.data.frame(cbind(ana_BLC_var, range_ana_BLC))
  ana_BLN <- as.data.frame(cbind(ana_BLN_var, range_ana_BLN))
  ana_BW <- as.data.frame(cbind(ana_BW_var, range_ana_BW))
  ana_BD <- as.data.frame(cbind(ana_BD_var, range_ana_BD))
  
  write.xlsx(ana_BLC, file = paste(folder, "var_range_ana_BLC_", archipelago, ".xlsx", sep = ""))
  write.xlsx(ana_BLN, file = paste(folder, "var_range_ana_BLN_", archipelago, ".xlsx", sep = ""))
  write.xlsx(ana_BW, file = paste(folder, "var_range_ana_BW_", archipelago, ".xlsx", sep = ""))
  write.xlsx(ana_BD, file = paste(folder, "var_range_ana_BD_", archipelago, ".xlsx", sep = ""))
  
  }


  
#---------------------combine variance and range 
 
tot_BLC <- as.data.frame(cbind(tot_BLC_var,range_tot_BLC ))
tot_BLN <- as.data.frame(cbind(tot_BLN_var,range_tot_BLN ))
tot_BW <- as.data.frame(cbind(tot_BW_var,range_tot_BW))
tot_BD <- as.data.frame(cbind(tot_BD_var,range_tot_BD ))
  


# Export to separate Excel files
write.xlsx(tot_BLC, file = paste(folder, "var_range_tot__BLC_", archipelago, ".xlsx", sep = ""))
write.xlsx(tot_BLN, file = paste(folder, "var_range_tot_BLN_", archipelago, ".xlsx", sep = ""))
write.xlsx(tot_BW, file = paste(folder, "var_range_tot_BW_", archipelago, ".xlsx", sep = ""))
write.xlsx(tot_BD, file = paste(folder, "var_range_tot_BD_", archipelago, ".xlsx", sep = ""))


}
