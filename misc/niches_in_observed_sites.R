



#==============================================================================
# check how many of each trophi niche are present within the null assemblage for each site
#==============================================================================
 species_in_null_assemblage <- data.frame(Trophic_niche = c("In","Om","Gr","Fr","Ne"))
 
 
 for(i in names(sites_tpd_diet_mapped_data_random)){
 
   data <- sites_tpd_diet_mapped_data_random[[i]]
 
   pull_species <- data.frame(Birdlife_Name = unique(unlist(str_split(data$occupying_species, pattern = "/"))))
 
   species_trophic <- pull_species %>% dplyr::left_join(Forage[,c("Birdlife_Name","Trophic_Niche")])
 
   table <- as.matrix(table(species_trophic$Trophic_Niche))
 
   trophic_df <- data.frame(Trophic_niche = rownames(table),species = table)
   colnames(trophic_df)[2] <- i
 
   species_in_null_assemblage <- species_in_null_assemblage %>% dplyr::left_join(trophic_df)
 
 }
 
 rownames(species_in_null_assemblage) <- species_in_null_assemblage$Trophic_niche
 species_in_null_assemblage <- data.frame(t(species_in_null_assemblage[,-1]))
 
 
 
 
 write_rds(species_in_null_assemblage, file = "Outputs/trophic_guilds_in_null_assemblage.rds")

#=================================
# observed species in site
#=================================
 species_in_observed_assemblage <- data.frame(Trophic_niche = c("In","Om","Gr","Fr","Ne"))
 
 
 for(i in names(sites_tpd_diet_mapped_data)){
 
   data <- sites_tpd_diet_mapped_data[[i]]
 
   pull_species <- data.frame(Birdlife_Name = unique(unlist(str_split(data$occupying_species, pattern = "/"))))
 
   species_trophic <- pull_species %>% dplyr::left_join(Forage[,c("Birdlife_Name","Trophic_Niche")])
 
   table <- as.matrix(table(species_trophic$Trophic_Niche))
 
   trophic_df <- data.frame(Trophic_niche = rownames(table),species = table)
   colnames(trophic_df)[2] <- i
 
   species_in_observed_assemblage <- species_in_observed_assemblage %>% dplyr::left_join(trophic_df)
 
 }
 
 rownames(species_in_observed_assemblage) <- species_in_observed_assemblage$Trophic_niche
 species_in_observed_assemblage <- data.frame(t(species_in_observed_assemblage[,-1]))
 
 
 write_rds(species_in_observed_assemblage, file = "Outputs/trophic_guilds_in_observed_assemblage.rds")
