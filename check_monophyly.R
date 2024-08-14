check_monophyly <- function(tree_name_file) {
  
  X_tree <- read.nexus(tree_name_file) 
  
  birds <- X_tree$tip.label #WE ARE WORKING WITH THE NAME OF THE TIPS
  
  label_parts <- strsplit(birds, "_") #SPLIT NAME OF TIPS BY THE UNDERSCORE (I KNOW THEY HAVE THIS STRUCTURE)
  
  common_prefix <- sapply(label_parts, function(x) paste(x[1], x[2], sep = "_")) #WE EXTRACT THE SPS NAME (GENUS_SPNAME)
  
  common_prefix <- unique(common_prefix) #AS THEY ARE MULTIPLE INDIVIDUAL PER SPECIES. THEN WE CHECK THE UNIQUE OCCURENCES
  
  all_of_sp <- list()
  
  for (j in 1:length(common_prefix)) {
    vect <- birds[grepl(common_prefix[j], birds)]  #COMBINE ALL THE INDIVIDUALS OF THE SAME SPS IN A LIST
    all_of_sp <- c(all_of_sp, list(vect))
  }      
  
  all_of_mono <- list()
  
  for (i in 1:length(all_of_sp)) {
    if (length(all_of_sp[[i]]) > 1) {
      all_of_mono <- c(all_of_mono, list(all_of_sp[[i]]))  #IF THERE IS MORE THAN 1 INDIVIDUAL OF THE SAME SPS IN THE LIST WE NEED
      #TO CHECK FOR MONOPHYLY
    }
  }
  
  to_prune <- c() #WHICH ONES WE CAN PRUNE (IF MONO)
  to_inspect <- c() #WHICH ONES WE CANT PRUNE FOR NOW
  
  for (i in 1:length(all_of_mono)) {
    if (is.monophyletic(X_tree, all_of_mono[[i]])) { #IF THEY ARE MONOPHYLETIC, WE CAN PRUNE
      to_prune <- c(to_prune, list(all_of_mono[[i]]))
    } else {    #ELSE WE NEED TO CONTINUE THE ANALYSIS
      to_inspect <- c(to_inspect, list(all_of_mono[[i]]))
    }
  }
  
  lost_case_ask_ben <- list()   # IF THE SPECIES IS JUST DEFINED AS A SPECIES (NOT DIVIDED INTO GEOGRAPHICAL RANGE/ SUBSPECIES)
  #THEN THERE IS NOT MUCH TO DO (ASK BEN OR SEARCH IN PAPERS)
  still_hope <- list() #IF THE SPS IS DIVIDED INTO GEOGRAPHICAL RANGE/ SUBSPECIES THEN WE NEED TO CHECK FOR MONOPHYLY TAKING THAT INTO ACCOUTNT
  #BUT IS DIFFICULT TO MAKE THE CODE SO DO IT MANUALLY
  
  for (z in 1:length(to_inspect)) {
    label_parts <- strsplit(as.character(to_inspect[[z]]), "_")
    common_prefix <- sapply(label_parts, function(x) paste(x[3], sep = "_"))
    
    if (length(common_prefix) == length(unique(common_prefix))) { #THIS IS BECAUSE EACH SPECIES IS XX_YY_ZZ, IF ZZ IS DIFFERENT, THEN IS 
      #THE CODE OF THE GENENBANK AND THERE IS NOT A DIFFERENCIATION BETWEEN SUBSPECIE/GEOGRAPHICAL REGION- SO THERE IS NO MORE HOPE
      lost_case_ask_ben <- c(lost_case_ask_ben, list(to_inspect[[z]]))
    } else { #ELSE BECAUSE AT LEAST SOME ZZ WILL BE SHARED BETWEEN SPS. SO WE KNOW IS THE DESCRIPTION OF THE GEOGRAPHICAL REGION/ SUBSPECIE
      still_hope <- c(still_hope, list(to_inspect[[z]]))
    }
  }
  return(list(
    all_of_mono = all_of_mono,
    to_prune = to_prune,    
    lost_case_ask_ben = lost_case_ask_ben,
    still_hope = still_hope
  ))
}



