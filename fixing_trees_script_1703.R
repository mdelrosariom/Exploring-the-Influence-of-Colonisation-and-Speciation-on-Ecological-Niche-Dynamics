#-------------------------------------------------------------------------------

library("phytools")
library("stringr")

source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/UPDATED_check_monophyly.R")
#tp prune the prunabble (monophyletic) branches
source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/generate_pruned.R")
#to delete still hope 
source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/drop_still_hope.R")
#to change tips
source("C:/Users/mdrmi/OneDrive/Escritorio/codes_phylo/change_tips.R")

#set working directory 
setwd("C:/Users/mdrmi/OneDrive/Escritorio/trees")

#-------------------------------------------------------------------------------

#Acrocephalus_maxcred.tre

tree <- read.nexus("Acrocephalus_maxcred.tre")

#check monophyletic problems 

tree_analyzed <- check_monophyly("Acrocephalus_maxcred.tre")

tree<- new_pruned(tree_analyzed, "Acrocephalus_maxcred.tre")
#-------------------------------------------------------------------------------
#Anairetes_maxcred.tre

#Open tree 

tree <- read.nexus("Anairetes_maxcred.tre")

#check monophyletic problems 

tree_analyzed <- check_monophyly("Anairetes_maxcred.tre")

#prune prunable 
#send to trees_reunion_703
tree <- new_pruned(tree_analyzed, "Anairetes_maxcred.tre")

#check monophyly in still hope 
#is.monophyletic(tree, "Anairetes_reguloides_AF067003"   ,"Anairetes_reguloides_reguloides_JX273114" )
#TRUE
tree <- drop_still_hope(c("Anairetes_reguloides_AF067003","Anairetes_reguloides_reguloides_JX273114"), tree)


#is.monophyletic(tree, "Anairetes_reguloides_albiventris_KJ828968"   ,"Anairetes_reguloides_albiventris_KJ828967")
#TRUE

tree <- drop_still_hope(c( "Anairetes_reguloides_albiventris_KJ828968","Anairetes_reguloides_albiventris_KJ828967"), tree)

#change tips 
tree <- change_tips(tree)

#not done still need to ask Ben


write.nexus(tree, file = ".../Anairetes_maxcred_final.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#Carduelis_maxcred.tre

tree <- read.nexus("Carduelis_maxcred.tre")

tree_analyzed <- check_monophyly("Carduelis_maxcred.tre")

#prune prunable 
#send to trees_reunion_703
tree <- new_pruned(tree_analyzed, "Carduelis_maxcred.tre")

# is.monophyletic(tree,"Carduelis_carduelis_CANARIES_Va3_L380",  "Carduelis_carduelis_CANARIES_L380"    )
#[1] TRUE

#is.monophyletic(tree, "Carduelis_carduelis_IberianPeninsula_IP5_L345" ,"Carduelis_carduelis_IberianPeninsula_L345", "Carduelis_carduelis_Morocco_MO6_L374"    ,      "Carduelis_carduelis_Sweden_EU325788" , "Carduelis_carduelis_parva_Spain_L76387"      , Carduelis_carduelis_caniceps_Nepal_L76388     )
#[1] TRUE

#is.monophyletic(tree, "Carduelis_carduelis_CANARIES_Va3_L380",
+                # "Carduelis_carduelis_CANARIES_L380")
#[1] TRUE

#ALORS... STILL HOPE CAN BE SEPARATED BETWEEN CONTINENTAL SPECIES AND CANARY SPECIES

tree <- drop_still_hope(c("Carduelis_carduelis_IberianPeninsula_IP5_L345" ,"Carduelis_carduelis_IberianPeninsula_L345", "Carduelis_carduelis_Morocco_MO6_L374"    ,      "Carduelis_carduelis_Sweden_EU325788" , "Carduelis_carduelis_parva_Spain_L76387" , "Carduelis_carduelis_caniceps_Nepal_L76388"     )
, tree)

tree<- drop_still_hope(c("Carduelis_carduelis_CANARIES_Va3_L380", "Carduelis_carduelis_CANARIES_L380"), tree)

#is.monophyletic(tree, "Carduelis_chloris_IberianPeninsula_PN357_L359", "Carduelis_chloris_Morocco_GU592668" , "Carduelis_chloris_Spain_L76297"   )
#  [1] TRUE  
  
tree <- drop_still_hope(c("Carduelis_chloris_IberianPeninsula_PN357_L359", "Carduelis_chloris_Morocco_GU592668" , "Carduelis_chloris_Spain_L76297" ), tree)



 # is.monophyletic(tree,"Carduelis_chloris_CANARIES_Va5_L382"   ,        "Carduelis_chloris_CANARIES_L382"  )
#  [1] TRUE 
tree <- drop_still_hope(c("Carduelis_chloris_CANARIES_Va5_L382"   ,        "Carduelis_chloris_CANARIES_L382" ), tree)
  #ALORS... STILL HOPE CAN BE SEPARATED BETWEEN CONTINENTAL SPECIES AND CANARY SPECIES  
#change tips 

tree <- change_tips(tree)

#not done still need to ask Ben


write.nexus(tree, file = ".../Carduelis_maxcred.tre_final.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#Cisticola_maxcred.tre
tree <- read.nexus("Cisticola_maxcred.tre")

tree_analyzed <- check_monophyly("Cisticola_maxcred.tre")

tree <- new_pruned(tree_analyzed, "Cisticola_maxcred.tre")



#is.monophyletic(tree,  "Cisticola_juncidis_AJ004314","Cisticola_juncidis_HQ608851", "Cisticola_juncidis_Z73474",  "Cisticola_juncidis_juncidis_DQ008493"  )
#[1] TRUE

tree <- drop_still_hope(c("Cisticola_juncidis_AJ004314","Cisticola_juncidis_HQ608851", "Cisticola_juncidis_Z73474",  "Cisticola_juncidis_juncidis_DQ008493"    )
                        , tree)

#s.monophyletic(tree, "Cisticola_juncidis_brunniceps_AB239499" ,
#                  "Cisticola_juncidis_brunniceps_AB239500" , "Cisticola_juncidis_brunniceps_AB239501" )
#[1] TRUE

tree <- drop_still_hope(c("Cisticola_juncidis_brunniceps_AB239499",  "Cisticola_juncidis_brunniceps_AB239500" , "Cisticola_juncidis_brunniceps_AB239501"),  tree )

#is.monophyletic(tree, "Cisticola_juncidis_tinnabulans_AB239502",
#+  "Cisticola_juncidis_tinnabulans_AB239503", "Cisticola_juncidis_tinnabulans_AB239504")
#[1] TRUE

tree <- drop_still_hope(c("Cisticola_juncidis_tinnabulans_AB239502",
 "Cisticola_juncidis_tinnabulans_AB239503", "Cisticola_juncidis_tinnabulans_AB239504"), tree)

tree <- change_tips(tree)

write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Cisticola_maxcred.tre", translate = TRUE)

#-------------------------------------------------------------------------------

#Clytorhynchus_maxcred.tre

tree <- read.nexus("Clytorhynchus_maxcred.tre")

tree_analyzed <- check_monophyly("Clytorhynchus_maxcred.tre")

tree <- new_pruned(tree_analyzed, "Clytorhynchus_maxcred.tre")


#-------------------------------------------------------------------------------

#Coccyzus_maxcred.tre

tree <- read.nexus("Coccyzus_maxcred.tre")
tree_analyzed <- check_monophyly("Coccyzus_maxcred.tre")

#is.monophyletic(tree,"Coccyzus_americanus_AF204993" , "Coccyzus_americanus_AY509696", "Coccyzus_americanus_HE793186",      "Coccyzus_americanus_U09265__CAU09265_" , "Coccyzus_americanus_americanus_AY046910"     )
#[1] TRUE


#is.monophyletic(tree, "Coccyzus_americanus_Mexico_AY046908",  "Coccyzus_americanus_Mexico_AY046909"   ,              "Coccyzus_americanus_NewMexico_AY046905", "Coccyzus_americanus_americanus_Minnesota_AF249270" ,  "Coccyzus_americanus_americanus_Minnesota_AF249271"  ,
#                  "Coccyzus_americanus_occidentalis_Alaska_AF249268"  ,  "Coccyzus_americanus_occidentalis_Alaska_AF249269"   
#                , "Coccyzus_americanus_occidentalis_NewMexico_AY046906", "Coccyzus_americanus_occidentalis_NewMexico_AY04690")    

#TRUE 

#ALORS SEPARATE BETWEEN AMERICUS AMERICUS (X ) AND THE ONES IN MEXICO.... SP1 SP2 BECAUSE I DONT KNOW THE DIFFERENCE

tree <- drop_still_hope(c("Coccyzus_americanus_AF204993" , "Coccyzus_americanus_AY509696", "Coccyzus_americanus_HE793186",      "Coccyzus_americanus_U09265__CAU09265_" , "Coccyzus_americanus_americanus_AY046910"  
), tree)


tree<- drop_still_hope(c("Coccyzus_americanus_Mexico_AY046908",  "Coccyzus_americanus_Mexico_AY046909"   ,    "Coccyzus_americanus_NewMexico_AY046905", "Coccyzus_americanus_americanus_Minnesota_AF249270" ,  "Coccyzus_americanus_americanus_Minnesota_AF249271"  ,
             "Coccyzus_americanus_occidentalis_Alaska_AF249268"  ,  "Coccyzus_americanus_occidentalis_Alaska_AF249269"   
                , "Coccyzus_americanus_occidentalis_NewMexico_AY046906", "Coccyzus_americanus_occidentalis_NewMexico_AY04690"), tree)



tree <- change_tips(tree)


#-------------------------------------------------------------------------------

#Columbiformes_maxcred.tre

tree <- read.nexus("Columbiformes_maxcred.tre")
tree_analyzed <- check_monophyly("Columbiformes_maxcred.tre")


# is.monophyletic(tree, "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912609",
#+  "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912610",
#+  "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912611",
#+  "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912612",
#+  "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912613",
#+  "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912614")
#[1] TRUE

#2 problem. idk problem
# #is.monophyletic(tree,      "Treron_calvus_BIOKO_B01_L746",  
# +                  "Treron_calvus_BIOKO_B02_L747" ,    "Treron_calvus_SAOTOME_L511",   "Treron_calvus_SAOTOME_L512"  ,   "Treron_calvus_SAOTOME_L513"    , "Treron_calvus_SAOTOME_P09_L748",
# +             "Treron_calvus_SAOTOME_P10_L749", "Treron_calvus_SAOTOME_P11_L750")
# [1] TRUE
# > is.monophyletic(tree,"Treron_calvus_AY443674" ,        "Treron_calvus_Angola_L687",  "Treron_calvus_Ghana_KT023396" )
# [1] TRUE

tree <- new_pruned(tree_analyzed, "Columbiformes_maxcred.tre")

tree <- drop_still_hope(c("Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912610",
                          "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912611",
                          "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912612",
                       "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912613",
                        "Hemiphaga_novaeseelandiae_novaeseelandiae_New_Zealand_GQ912614"),tree)

tree<- drop_still_hope(c( "Treron_calvus_BIOKO_B01_L746",  
                                         "Treron_calvus_BIOKO_B02_L747" ,    "Treron_calvus_SAOTOME_L511",   "Treron_calvus_SAOTOME_L512"  ,   "Treron_calvus_SAOTOME_L513"    , "Treron_calvus_SAOTOME_P09_L748",
                                    "Treron_calvus_SAOTOME_P10_L749", "Treron_calvus_SAOTOME_P11_L750"), tree)

tree<- drop_still_hope(c("Treron_calvus_AY443674" ,        "Treron_calvus_Angola_L687",  "Treron_calvus_Ghana_KT023396"),tree)


#-------------------------------------------------------------------------------
#Copsychus_maxcred.tre
tree <- read.nexus("Copsychus_maxcred.tre")
tree_analyzed <- check_monophyly("Copsychus_maxcred.tre")
#-------------------------------------------------------------------------------

#Corvus_maxcred.tre

tree <- read.nexus("Corvus_maxcred.tre")
tree_analyzed <- check_monophyly("Corvus_maxcred.tre")
tree <- new_pruned(tree_analyzed, "Corvus_maxcred.tre")

tree <- drop_still_hope(c( "Corvus_corone_U86032","Corvus_corone_corone_Austria_HE805700"), tree)

tree <- change_tips(tree)


write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Corvus_maxcred_final.tre", translate = TRUE)


#-------------------------------------------------------------------------------
#Cuculiformes_maxcred.tre
tree <- read.nexus("Cuculiformes_maxcred.tre")
tree_analyzed <- check_monophyly("Cuculiformes_maxcred.tre")
tree <- new_pruned(tree_analyzed,"Cuculiformes_maxcred.tre")
tree <- drop_still_hope(c("Coccyzus_americanus_Mexico_AY046908"  ,"Coccyzus_americanus_Mexico_AY046909"  ,               "Coccyzus_americanus_NewMexico_AY046905" ,"Coccyzus_americanus_americanus_Minnesota_AF249270" ,  "Coccyzus_americanus_americanus_Minnesota_AF249271"  ,"Coccyzus_americanus_AF204993"  ,  "Coccyzus_americanus_AY509696","Coccyzus_americanus_HE793186", "Coccyzus_americanus_U09265__CAU09265_",               "Coccyzus_americanus_americanus_AY046910" ), tree)
tree<- drop_still_hope(c("Coccyzus_americanus_occidentalis_Alaska_AF249268",    "Coccyzus_americanus_occidentalis_Alaska_AF249269"   ,
                                         "Coccyzus_americanus_occidentalis_NewMexico_AY046906", "Coccyzus_americanus_occidentalis_NewMexico_AY046"), tree)
#-------------------------------------------------------------------------------

#Cyanistes_maxcred.tre

tree <- read.nexus("Cyanistes_maxcred.tre")
tree_analyzed <- check_monophyly("Cyanistes_maxcred.tre")


tree <- new_pruned(tree_analyzed,"Cyanistes_maxcred.tre")

# #is.monophyletic(tree,"Cyanistes_teneriffae_GranCanaria_CANARIES_DQ474055"   ,       
# +                  "Cyanistes_teneriffae_GranCanaria_CANARIES_DQ474058"        ,
# +                 "Cyanistes_teneriffae_LaGomera_CANARIES_DQ474002"             ,
# +                  "Cyanistes_teneriffae_LaGomera_CANARIES_KP759165"             ,
# +                  "Cyanistes_teneriffae_degener_CanaryIslands_CANARIES_FV1_L338",
# +                  "Cyanistes_teneriffae_hedwigii_GranCanaria_CANARIES_KP759205" ,
# +                  "Cyanistes_teneriffae_hedwigii_GranCanaria_CANARIES_KP759233" ,
# +                  "Cyanistes_teneriffae_palmensis_CANARIES_PN496_L332"          ,
# +                  "Cyanistes_teneriffae_teneriffae_LaGomeraK_CANARIES_P759166"  ,
# +                  "Cyanistes_teneriffae_teneriffae_LaGomera_CANARIES_KP759202"  ,
# +                  "Cyanistes_teneriffae_teneriffae_LaGomera_CANARIES_KP759216"  ,
# +                  "Cyanistes_teneriffae_teneriffae_Tenerife_CANARIES_KP759192"  ,
# +                  "Cyanistes_teneriffae_teneriffae_Tenerife_CANARIES_KP759210"  ,
# +                  "Cyanistes_teneriffae_teneriffae_Tenerife_CANARIES_KP759247"  ,
# +                  "Cyanistes_teneriffae_teneriffae_LaGomera_CANARIES_KC202344"  )
# [1] TRUE
# > is.monophyletic(tree, "Cyanistes_teneriffae_Algeria_KP759209_KP759209",
#                   +                 "Cyanistes_teneriffae_Algeria_KP759232"                ,
#                   +                  "Cyanistes_teneriffae_Libya_KP759204"                  ,
#                   +                  "Cyanistes_teneriffae_Libya_KP759256"                   ,
#                   +                  "Cyanistes_teneriffae_cyrenaicae_LybiaKC202348"          ,
#                   +                  "Cyanistes_teneriffae_cyrenaicae_Lybia_KC202346"          ,
#                   +                  "Cyanistes_teneriffae_cyrenaicae_Lybia_KC202347"           ,
#                   +                  "Cyanistes_teneriffae_cyrenaicae_Lybia_KC202349"            ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Algeria_KP759237"           ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Algeria_KP759254"            ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Algeria_KP759268"            ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Italy_KC202355"              ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Italy_KC202356"              ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Morocco_KC202333"            ,
#                   +                  "Cyanistes_teneriffae_teneriffae_Morocco_KC202357"            ,
#                   +                 "Cyanistes_teneriffae_teneriffae_Morocco_KC202358" )
# [1] TRUE

tree <- drop_still_hope(c("Cyanistes_teneriffae_GranCanaria_CANARIES_DQ474055"   ,       
                                            "Cyanistes_teneriffae_GranCanaria_CANARIES_DQ474058"        ,
                                           "Cyanistes_teneriffae_LaGomera_CANARIES_DQ474002"             ,
                                            "Cyanistes_teneriffae_LaGomera_CANARIES_KP759165"             ,
                                            "Cyanistes_teneriffae_degener_CanaryIslands_CANARIES_FV1_L338",
                                            "Cyanistes_teneriffae_hedwigii_GranCanaria_CANARIES_KP759205" ,
                                           "Cyanistes_teneriffae_hedwigii_GranCanaria_CANARIES_KP759233" ,
                                           "Cyanistes_teneriffae_palmensis_CANARIES_PN496_L332"          ,
                                           "Cyanistes_teneriffae_teneriffae_LaGomeraK_CANARIES_P759166"  ,
                                           "Cyanistes_teneriffae_teneriffae_LaGomera_CANARIES_KP759202"  ,
                                            "Cyanistes_teneriffae_teneriffae_LaGomera_CANARIES_KP759216"  ,
                                            "Cyanistes_teneriffae_teneriffae_Tenerife_CANARIES_KP759192"  ,
                                           "Cyanistes_teneriffae_teneriffae_Tenerife_CANARIES_KP759210"  ,
                                            "Cyanistes_teneriffae_teneriffae_Tenerife_CANARIES_KP759247"  ,
                                           "Cyanistes_teneriffae_teneriffae_LaGomera_CANARIES_KC202344"  ), tree)
tree <- drop_still_hope(c("Cyanistes_teneriffae_Algeria_KP759209_KP759209",
                                           "Cyanistes_teneriffae_Algeria_KP759232"                ,
                                           "Cyanistes_teneriffae_Libya_KP759204"                  ,
                                            "Cyanistes_teneriffae_Libya_KP759256"                   ,
                                           "Cyanistes_teneriffae_cyrenaicae_LybiaKC202348"          ,
                                           "Cyanistes_teneriffae_cyrenaicae_Lybia_KC202346"          ,
                                            "Cyanistes_teneriffae_cyrenaicae_Lybia_KC202347"           ,
                                            "Cyanistes_teneriffae_cyrenaicae_Lybia_KC202349"            ,
                                            "Cyanistes_teneriffae_teneriffae_Algeria_KP759237"           ,
                                           "Cyanistes_teneriffae_teneriffae_Algeria_KP759254"            ,
                                            "Cyanistes_teneriffae_teneriffae_Algeria_KP759268"            ,
                                           "Cyanistes_teneriffae_teneriffae_Italy_KC202355"              ,
                                            "Cyanistes_teneriffae_teneriffae_Italy_KC202356"              ,
                                            "Cyanistes_teneriffae_teneriffae_Morocco_KC202333"            ,
                                            "Cyanistes_teneriffae_teneriffae_Morocco_KC202357"            ,
                                           "Cyanistes_teneriffae_teneriffae_Morocco_KC202358" ),tree)

#-------------------------------------------------------------------------------

#Cyanolanius_maxcred.tre

tree <- read.nexus("Cyanolanius_maxcred.tre")
tree_analyzed <- check_monophyly("Cyanolanius_maxcred.tre")


tree <- new_pruned(tree_analyzed,"Cyanolanius_maxcred.tre")

tree <- drop_still_hope(c("Emberiza_citrinella_Sweden_EU325753",	"Emberiza_citrinella_Sweden_EU571277"),tree)


tree <- change_tips(tree)

write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Cyanolanius_maxcred_final.tre", translate = TRUE)


#-------------------------------------------------------------------------------
#Emberiza_maxcred.tre
tree <- read.nexus("Emberiza_maxcred.tre")
tree_analyzed <- check_monophyly("Emberiza_maxcred.tre")
tree <- new_pruned(tree_analyzed,"Emberiza_maxcred.tre")

tree <- drop_still_hope(c("Emberiza_citrinella_Sweden_EU325753",	"Emberiza_citrinella_Sweden_EU571277"),tree)

tree <- change_tips(tree)
write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Emberiza_maxcred.tre", translate = TRUE)


#-------------------------------------------------------------------------------

#Estrilda_Erythrura_maxcred.tre
tree <- read.nexus("Estrilda_Erythrura_maxcred.tre")
tree_analyzed <- check_monophyly("Estrilda_Erythrura_maxcred.tre")
tree <- new_pruned(tree_analyzed,"Estrilda_Erythrura_maxcred.tre")

tree <- drop_still_hope(c(),tree)

tree <- change_tips(tree)
write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Emberiza_maxcred.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#Finches_Galapagos_Cocos.tre

tree <- read.nexus("Finches_Galapagos_Cocos.tre")
tree_analyzed <- check_monophyly("Finches_Galapagos_Cocos.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Finches_Galapagos_Cocos.tre")


#-------------------------------------------------------------------------------
#Foudia_maxcred_Lerner_rate.tre

tree <- read.nexus("Foudia_maxcred_Lerner_rate.tre")
tree_analyzed <- check_monophyly("Foudia_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Foudia_maxcred_Lerner_rate.tre")

# > is.monophyletic(tree, "F._e._aldabrana_ALDABRA_341a" ,  "F._e._aldabrana_ALDABRA_342a","F._e._aldabrana_ALDABRA_343a")
# [1] TRUE
# > is.monophyletic(tree,"F._e._algondae_COMOROS_90c","F._e._algondae_COMOROS_94a","F._e._algondae_COMOROS_99a","F._e._anjouanensis_COMOROS_249a" , "F._e._consobrina_COMOROS_209a","F._e._consobrina_COMOROS_228" , "F._e._consobrina_COMOROS_229","F._e._eminentissima_COMOROS_135a","F._e._eminentissima_COMOROS_179a","F._e._eminentissima_COMOROS_194",  "F._e._omissa_O1711_E_Central_Ma" ,  "F._e._omissa_O1712_E_Central_Ma" )
# [1] TRUE

tree <- drop_still_hope(c("F._e._aldabrana_ALDABRA_341a" ,  "F._e._aldabrana_ALDABRA_342a","F._e._aldabrana_ALDABRA_343a"),tree)
tree <- drop_still_hope(c("F._e._algondae_COMOROS_90c","F._e._algondae_COMOROS_94a","F._e._algondae_COMOROS_99a","F._e._anjouanensis_COMOROS_249a" , "F._e._consobrina_COMOROS_209a","F._e._consobrina_COMOROS_228" , "F._e._consobrina_COMOROS_229","F._e._eminentissima_COMOROS_135a","F._e._eminentissima_COMOROS_179a","F._e._eminentissima_COMOROS_194",  "F._e._omissa_O1711_E_Central_Ma" ,  "F._e._omissa_O1712_E_Central_Ma" ),tree)

#-------------------------------------------------------------------------------
#Foudia_maxcred_rate_of_cytb.tre

tree <- read.nexus("Foudia_maxcred_rate_of_cytb.tre")
tree_analyzed <- check_monophyly("Foudia_maxcred_rate_of_cytb.tre")
tree_analyzed

#-------------------------------------------------------------------------------

#Fregilupus_maxcred_Lerner_rate.tre

tree <- read.nexus("Fregilupus_maxcred_Lerner_rate.tre")
tree_analyzed <- check_monophyly("Fregilupus_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Fregilupus_maxcred_Lerner_rate.tre")




#-------------------------------------------------------------------------------
#Fringilla_maxcred.tre

tree <- read.nexus("Fringilla_maxcred.tre")
tree_analyzed <- check_monophyly("Fringilla_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Fringilla_maxcred.tre" )
tree <- drop_still_hope(c(  "Fringilla_coelebs_moreletti_AZORES_AF002879",    
                            "Fringilla_coelebs_moreletti_AZORES_iAF002880",   
                            "Fringilla_coelebs_moreletti_AZORES_iAF002881"), tree)
tree<- drop_still_hope(c(  "Fringilla_coelebs_maderensis_MADEIRA_AF002887" , 
                           "Fringilla_coelebs_maderensis_MADEIRA_AF002888"),tree)
tree<- drop_still_hope(c(  "Fringilla_coelebs_canariensis_CANARIES_GQ330546",
                           "Fringilla_coelebs_canariensis_CANARIES_GQ330547",
                           "Fringilla_coelebs_canariensis_CANARIES_GQ330548"
                             ),tree)
tree<- drop_still_hope(c("Fringilla_coelebs_ombriosa_CANARIES_GQ330543",   
                          "Fringilla_coelebs_ombriosa_CANARIES_GQ330544" ),tree)
tree<-drop_still_hope(c("Fringilla coelebs palmae CANARIES CQ330545", "Fringilla coelebs palmae CANARIES PN493 L331"), tree)
tree <- drop_still_hope(c("Fringilla_coelebs_Madrid_L76609", "Fringilla_coelebs_Morocco_MO1_L369","Fringilla_coelebs_AF447368","Fringilla_coelebs_AY495389","Fringilla_coelebs_L76610", "Fringila_coelebs_IberianPeninsula_PN157_L352"      ),tree)


tree<- change_tips(tree)

write.nexus( tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Fringilla_maxcred.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#Honeycreepers_maxcred.tre
tree<- read.nexus("Honeycreepers_maxcred.tre")
tree_analyzed<- check_monophyly("Honeycreepers_maxcred.tre")
tree <- new_pruned(tree_analyzed,"Honeycreepers_maxcred.tre" )

#-------------------------------------------------------------------------------
#Horornis_maxcred.tre
tree<- read.nexus("Horornis_maxcred.tre")
tree_analyzed<- check_monophyly("Horornis_maxcred.tre")
tree <- new_pruned(tree_analyzed,"Horornis_maxcred.tre" )

#-------------------------------------------------------------------------------
#Humblotia_maxcred.tre

tree<- read.nexus("Humblotia_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Humblotia_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Humblotia_maxcred.tre" )

#-------------------------------------------------------------------------------

#Hypsipetes_maxcred_Lerner_rate.tre

tree<- read.nexus("Hypsipetes_maxcred_Lerner_rate.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Hypsipetes_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Hypsipetes_maxcred_Lerner_rate.tre" )
plotTree(tree)
# [1] "Hypsipetes_leucocephalus_China_GU112695"
# [2] "Hypsipetes_leucocephalus_India_KP398686"
# [3] "Hypsipetes_leucocephalus_India_KP398689"
# [4] "Hypsipetes_leucocephalus_India_KP398694"
# [5] "Hypsipetes_leucocephalus_India_KP398696"
# [6] "Hypsipetes_leucocephalus__GQ242114" 
#2,3,4,5

tree<- drop_still_hope(c("Hypsipetes_leucocephalus_India_KP398686",
           "Hypsipetes_leucocephalus_India_KP398689",
           "Hypsipetes_leucocephalus_India_KP398694",
         "Hypsipetes_leucocephalus_India_KP398696"), tree)



# [1] "Hypsipetes_madagascariensis_COMOROS_AY590717"           
# [2] "Hypsipetes_madagascariensis_COMOROS_AY590719"           
# [3] "Hypsipetes_madagascariensis_COMOROS_AY590720"           
# [4] "Hypsipetes_madagascariensis_COMOROS_AY590721"           
# [5] "Hypsipetes_madagascariensis_COMOROS_AY590725"           
# [6] "Hypsipetes_madagascariensis_COMOROS_AY590730"           
# [7] "Hypsipetes_madagascariensis_COMOROS_Mayotte_AY590723"   
# [8] "Hypsipetes_madagascariensis_COMOROS_Mayotte_AY590724"   
# [9] "Hypsipetes_madagascariensis_COMOROS__Mayotte_AY590722"  
# [10] "Hypsipetes_madagascariensis_Madagascar_AY590718"        
# [11] "Hypsipetes_madagascariensis_Madagascar_AY590726"        
# [12] "Hypsipetes_madagascariensis_Madagascar_AY590727"        
# [13] "Hypsipetes_madagascariensis_Madagascar_AY590728"        
# [14] "Hypsipetes_madagascariensis_Madagascar_AY590729"        
# [15] "Hypsipetes_madagascariensis_Madagascar_AY590731"        
# [16] "Hypsipetes_madagascariensis_Madagascar_AY590732"        
# [17] "Hypsipetes_madagascariensis_Madagascar_DQ402286"        
# [18] "Hypsipetes_madagascariensis__AY590745"                  
# [19] "Hypsipetes_madagascariensis__AY590746"                  
# [20] "Hypsipetes_madagascariensis__AY590747"                  
# [21] "Hypsipetes_madagascariensis_psaroides_India_AY590756"   
# [22] "Hypsipetes_madagascariensis_psaroides_Pakistan_AY590755"
# [23] "Hypsipetes_madagascariensis_rostratus_ALDABRA_AY590749" 
# [24] "Hypsipetes_madagascariensis_rostratus_ALDABRA_AY590750" 
# [25] "Hypsipetes_madagascariensis_rostratus_ALDABRA_AY590751" 
# 


#still hope 3 


tree <- drop_still_hope(c("Hypsipetes_parvirostris_moheliensis_COMOROS_AY590737", 
                          "Hypsipetes_parvirostris_moheliensis_COMOROS_AY590738" ,
                          "Hypsipetes_parvirostris_moheliensis_COMOROS_AY590748" ), tree)


tree <- drop_still_hope(c( "Hypsipetes_parvirostris_parvirostris_COMOROS_AY590739",
                           "Hypsipetes_parvirostris_parvirostris_COMOROS_AY590740"),tree)

#-------------------------------------------------------------------------------

#Lamprotornis_maxcred_Lerner_rate.tre

tree<- read.nexus("Lamprotornis_maxcred_Lerner_rate.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Lamprotornis_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Lamprotornis_maxcred_Lerner_rate.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Lanius_maxcred.tre

tree<- read.nexus("Lanius_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Lanius_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Lanius_maxcred.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
tree <- read.nexus("Leptosomus_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Leptosomus_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Lanius_maxcred.tre" )
plotTree(tree)

tree <- drop_still_hope(c("Lepstosomus_discolor_COMOROS_L612A", "Lepstosomus_discolor_COMOROS_L613A", 
                        "Leptosomus_discolor_Madagascar_AF407449"), tree)
tree <- change_tips(tree)

write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Leptosomus_maxcred.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#"Lonchura_maxcred.tre
tree <- read.nexus("Lonchura_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Lonchura_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Lonchura_maxcred.tre" )
plotTree(tree)


#-------------------------------------------------------------------------------

#Loxia_maxcred.tre


tree <- read.nexus("Loxia_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Loxia_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Loxia_maxcred.tre" )
plotTree(tree)


#-------------------------------------------------------------------------------
#Microeca_Eopsaltria_maxcred.tre

tree <- read.nexus("Microeca_Eopsaltria_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Microeca_Eopsaltria_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Microeca_Eopsaltria_maxcred.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Mimus_maxcred_Lerner_rate.tre

tree <- read.nexus("Mimus_maxcred_Lerner_rate.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Mimus_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Mimus_maxcred_Lerner_rate.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Monarcha_maxcred.tre

tree <- read.nexus("Monarcha_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Monarcha_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Monarcha_maxcred.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Motacilla_maxcred.tre

tree <- read.nexus("Motacilla_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Motacilla_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Motacilla_maxcred.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Myadestes_maxcred.tre

tree <- read.nexus("Myadestes_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Myadestes_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Myadestes_maxcred.tre" )
plotTree(tree)

tree <- drop_still_hope(c("Myadestes_ralloides_Bolivia_DQ463688","Myadestes_ralloides_Ecuador_DQ463687",
                          "Myadestes_ralloides_HM633343","Myadestes_ralloides_Peru_DQ463686",
                          "Myadestes_ralloides_Peru_DQ463690" ), tree)

tree <- change_tips(tree)

write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Myadestes_maxcred.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#Myiagra_maxcred_Lerner_rate.tre
tree <- read.nexus("Myiagra_maxcred_Lerner_rate.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Myiagra_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Myiagra_maxcred_Lerner_rate.tre" )
plotTree(tree)


#-------------------------------------------------------------------------------
#Myiarchus_maxcred_Lerner_rate.tre
tree <- read.nexus("Myiarchus_maxcred_rate_of_cytb.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Myiarchus_maxcred_rate_of_cytb.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Myiarchus_maxcred_rate_of_cytb.tre" )
plotTree(tree)


tree$tip.label <- str_replace_all(tree$tip.label,'M','Myiarchus')


#-------------------------------------------------------------------------------
#Petroica_maxcred.tre

tree <- read.nexus("Petroica_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Petroica_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Petroica_maxcred.tre" )
plotTree(tree)
#-------------------------------------------------------------------------------
#Pipilo_maxcred.tre

tree <- read.nexus("Pipilo_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Pipilo_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Pipilo_maxcred.tre" )
plotTree(tree)


#-------------------------------------------------------------------------------
#Phylloscopus_maxcred.tre

tree <- read.nexus("Phylloscopus_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Phylloscopus_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Phylloscopus_maxcred.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Prinia_maxcred.tre

tree <- read.nexus("Prinia_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Prinia_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Prinia_maxcred.tre" )
plotTree(tree)


#-------------------------------------------------------------------------------
#Progne_maxcred.tre
tree <- read.nexus("Progne_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Progne_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Progne_maxcred.tre" )
plotTree(tree)

tree<- drop_still_hope(c("Progne_chalybea_Panama_AY825948","Progne_chalybea_Panama_EU427738", 
                         "Progne_chalybea_CostaRica_EU427737" ),tree)

#-------------------------------------------------------------------------------
#Psittaciformes_maxcred.tre
tree <- read.nexus("Psittaciformes_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Psittaciformes_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Psittaciformes_maxcred.tre" )
plotTree(tree)



tree <- drop_still_hope(c())
#-------------------------------------------------------------------------------
#Regulus_maxcred.tre
tree <- read.nexus("Regulus_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Regulus_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Regulus_maxcred.tre" )
plotTree(tree)

tree <- drop_still_hope(c())

#-------------------------------------------------------------------------------
#Saxicola_maxcred.tre
tree <- read.nexus("Saxicola_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Saxicola_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Saxicola_maxcred.tre" )
plotTree(tree)

tree <- drop_still_hope(c("Saxicola_rubicola_Austria_CytB_EU190904_" ,         
                           "Saxicola_rubicola_IberianPeninsula_CytB_EU421082_", 
                          "Saxicola_rubicola_IberianPeninsula_CytB_EU421083_" ,
                           "Saxicola_rubicola_IberianPeninsula_CytB_EU421089_" ,
                           "Saxicola_rubicola_IberianPeninsula_CytB__EU421080_"
                              ), tree)

tree <- drop_still_hope(c( "Saxicola_torquata_Tanzania_EU421093_"   ,                    
                          "Saxicola_torquata_Tanzania__EU190921_"    ,                  
                           "Saxicola_torquata_Tanzania__EU421094"     ,                  
                           "Saxicola_torquata_Tanzania__EU421095_"), tree)



tree <- drop_still_hope(c("Saxicola_torquata_sibilla_Madagascar_EU190913"   ,           
                          "Saxicola_torquata_sibilla_Madagascar_EU190914"    ,          
                          "Saxicola_torquata_sibilla_Madagascar_EU190915"     ,         
                           "Saxicola_torquata_sibilla_Madagascar_EU190916"     ,         
                           "Saxicola_torquata_sibilla_Madagascar_L745"          ,        
                         "Saxicola_torquata_voeltzkowi_COMOROS_EU421110"     ), tree)


new_name <- str_replace_all("Saxicola_torquata_Nepal_CytB_EU421088_", "_", "")
tree$tip.label[tree$tip.label=="Saxicola_torquata_Nepal_CytB_EU421088_"] <- new_name


new_name <- str_replace_all("Saxicola_torquata_maura__AY286399_"  , "_", "")

tree$tip.label[tree$tip.label=="Saxicola_torquata_maura__AY286399_"  ] <- new_name

new_name <- str_replace_all("Saxicola_torquata_rubicola_IberianPeninsula_CytB_PN200_L354" , "_", "")

tree$tip.label[tree$tip.label=="Saxicola_torquata_rubicola_IberianPeninsula_CytB_PN200_L354"  ] <- new_name


tree<- change_tips(tree)
write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Saxicola_maxcred.tre", translate = TRUE)

#-------------------------------------------------------------------------------
#Setophaga_CytB_macred.tre
tree <- read.nexus("Setophaga_CytB_macred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Setophaga_CytB_macred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Setophaga_CytB_macred.tre" )
plotTree(tree)

#-------------------------------------------------------------------------------
#Sturnus_maxcred_Lerner_rate.tre
tree <- read.nexus("Sturnus_maxcred_Lerner_rate.tre")
plotTree(tree)
tree_analyzed <- check_monophyly("Sturnus_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed, "Sturnus_maxcred_Lerner_rate.tre")

#-------------------------------------------------------------------------------
#Sunbirds_maxcred.tre
tree <- read.nexus("Sunbirds_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Sunbirds_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Sunbirds_maxcred.tre" )
plotTree(tree)

tree <- drop_still_hope(c())


#-------------------------------------------------------------------------------
#Sylvia_maxcred.tre
tree <- read.nexus("Sylvia_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Sylvia_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Sylvia_maxcred.tre" )
plotTree(tree)

new_name <- str_replace_all("Sylvia_deserticola_JF502294"  , "_", "")
tree$tip.label[tree$tip.label=="Sylvia_deserticola_JF502294"] <- new_name


new_name <- str_replace_all("Sylvia_deserticola_Morocco_AJ534540", "_", "")
tree$tip.label[tree$tip.label=="Sylvia_deserticola_Morocco_AJ534540"] <- new_name


new_name <- str_replace_all("Sylvia_hortensis_JF502295"     , "_", "")
tree$tip.label[tree$tip.label=="Sylvia_hortensis_JF502295"   ] <- new_name

new_name <- str_replace_all("Sylvia_hortensis_Spain_AJ534534", "_", "")
tree$tip.label[tree$tip.label== "Sylvia_hortensis_Spain_AJ534534"] <- new_name

new_name <- str_replace_all("Sylvia_undata_FranceCorsica_AJ534542", "_", "")
tree$tip.label[tree$tip.label=="Sylvia_undata_FranceCorsica_AJ534542"] <- new_name

new_name <- str_replace_all("Sylvia_undata_JF502312", "_", "")
tree$tip.label[tree$tip.label=="Sylvia_undata_JF502312"] <- new_name

tree <- change_tips(tree)
plotTree(tree)

write.nexus(tree, file = "C:/Users/mdrmi/OneDrive/Escritorio/fixed_trees_703/Sylvia_maxcred.tre", translate = TRUE)


#-------------------------------------------------------------------------------
#Terpsiphone_maxcred.tre
tree <- read.nexus("Terpsiphone_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Terpsiphone_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Terpsiphone_maxcred.tre" )
plotTree(tree)

tree<- drop_still_hope(c("Terpsiphone_rufiventer_Angola_L688"  ,          
                          "Terpsiphone_rufiventer_Ghana_KP036785",         
                          "Terpsiphone_rufiventer_KP036666"       ,        
                         "Terpsiphone_rufiventer_KP036710"         ,      
                          "Terpsiphone_rufiventer_SierraLeone_KP036725",   
                          "Terpsiphone_rufiventer_mayombe_JQ885720"     ,  
                          "Terpsiphone_rufiventer_tricolor_BIOKO_JQ885721"), tree)


plotTree(tree)

tree <- drop_still_hope(c("Terpsiphone_rufocinerea_JQ885722", "Terpsiphone_rufocinerea_JQ885723"), tree)



#[1]          
# [2] "Terpsiphone_paradisi_EF081356"            
# [3]             
# [4] "Terpsiphone_paradisi_KP036641"            
# [5] "Terpsiphone_paradisi_KP036736"            
# [6]  
# [7] 
# [8] "
# [9] 
# [10]      
# [11] "Terpsiphone_paradisi_incei_JQ885716"      
# [12] "Terpsiphone_paradisi_leucogaster_JQ885717"
# [13]
# [14] "Terpsiphone_paradisi_paradisi_JQ885719"   

tree <-drop_still_hope(c("Terpsiphone_paradisi_incei_JQ885716", "Terpsiphone_paradisi_KP036641",
                         "Terpsiphone_paradisi_KP036736", "Terpsiphone_paradisi_EF081356"  ), tree)
tree<- drop_still_hope(c("Terpsiphone_paradisi_floris_JQ885715", "Terpsiphone_paradisi_nicobarica_JQ885718", "Terpsiphone_paradisi_affinis_JQ885711",  
                         "Terpsiphone_paradisi_borneensis_JQ885713" ,"Terpsiphone_paradisi_borneensis_JQ885714" ,
                         "Terpsiphone_paradisi_australis_JQ885712"  ,"Terpsiphone_paradisi_KJ456485", "Terpsiphone_paradisi_AF096466"  ), tree)

#-------------------------------------------------------------------------------

#Troglodytes_Thryomanes_maxcred_Lerner_rate.tre

tree <- read.nexus("Troglodytes_Thryomanes_maxcred_Lerner_rate.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Troglodytes_Thryomanes_maxcred_Lerner_rate.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Troglodytes_Thryomanes_maxcred_Lerner_rate.tre" )
plotTree(tree)

tree <- drop_still_hope(c("Troglodytes_aedon_Peru_KF386173",
                           "Troglodytes_aedon_Peru_KF386201","Troglodytes_aedon_AF104979"),tree)

tree <- change_tips(tree)

#-------------------------------------------------------------------------------

#Turdus_maxcred.tre

tree <- read.nexus("Turdus_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Turdus_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Turdus_maxcred.tre" )
plotTree(tree)

tree <- drop_still_hope(c("Turdus_rubrocanus_KJ456507_","Turdus_rubrocanus_gouldi_EU154670_China"), tree)

tree <- drop_still_hope(c("Turdus_olivaceus_nyikae_Tanzania_EU154649", "Turdus_olivaceus_nyikae_EU154649_Tanzania"), tree)


#-------------------------------------------------------------------------------

#Weavers_maxcred.tre

tree <- read.nexus("Weavers_maxcred.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Weavers_maxcred.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Weavers_maxcred.tre")
plotTree(tree)

tree <- drop_still_hope(c("Euplectes_orix_SouthAfrica_L737", "Euplectes_orix_SouthAfrica_L738","Euplectes_orix_South_Africa_AM710244"), tree)
tree <- drop_still_hope(c( "Ploceus_philippinus_travancoreensis__KF289832", "Ploceus_philippinus_travancoreensis__KJ442613"), tree)
tree <- drop_still_hope(c("Quelea_erythrops_SAOTOME_L494","Quelea_erythrops_SAOTOME_L495"),tree)
tree <- drop_still_hope(c("Quelea_erythrops_Gabon_L714",   "Quelea_erythrops_Gabon_L726",   "Quelea_erythrops_SAOTOME_L493"),tree)


#-------------------------------------------------------------------------------
#Zosterops_Indian_Atlantic.tre


tree <- read.nexus("Zosterops_Indian_Atlantic.tre")
plotTree(tree)
tree_analyzed<- check_monophyly("Zosterops_Indian_Atlantic.tre")
tree_analyzed
tree <- new_pruned(tree_analyzed,"Zosterops_Indian_Atlantic.tre")
plotTree(tree)
















