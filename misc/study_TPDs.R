#===========================
# SCRIPT
#===========================


# remove environment objects and clear memory
rm(list = ls())
gc()

# load in packages 

require(tidyverse)
require(geometry)
source("Functions/tpd_computation_functions.R")
require(doParallel)
require(TPD)


# load in data 

site_comparisons <- readRDS("Outputs/land_use_centroid_shift.rds") %>%
  dplyr::filter(str_count(land_use_comp, pattern = "Primary") < 2)


# predicts TPDs

predicts_tpds <- readRDS("Outputs/PREDICTS_sites_tpds_68_unscaled.rds")

# PREDICTS assemblage information

PREDICTS <- readRDS("Outputs/refined_predicts.rds")

# species pools

species_pool <- readRDS("Outputs/predicts_sites_species_pools.rds")

# species traits

traits <- readRDS("Outputs/full_morpho_traits_list_68.rds")


Forage <- readRDS("../PREDICTS_Taxonomy/PREDICTS_imputed_BL_traits_forage.rds")

#=======================================================================


studies <- unique(site_comparisons$SS)

study_TPDs <- list()



  
  
for(study in studies){

study_species  <- unique(PREDICTS %>% dplyr::filter(SS == study) %>% pull(Birdlife_Name))
    
study_species_pool <- unique(c(Reduce(lapply(species_pool[grep(names(species_pool), pattern = study)], FUN = function(x) x[["0.9"]]), f = "c"), study_species))

  
trait_ranges <- get_species_trait_ranges(species = study_species_pool,
                                         traits = traits, buffer = 0.1)


registerDoParallel(cores = 15)

species_TPDs <- foreach(sp = study_species_pool,
                        .combine = "c",
                        .packages = c("TPD","tidyverse"), .inorder = TRUE) %dopar% {
                          
                          sp_TPD <- species_TPD(species = sp,
                                                trait_ranges = trait_ranges,
                                                traits = traits, alpha = 0.68)
                          
                          
                          
                          # which cells have a non-zero probability of occupancy
                          presence <- which(sp_TPD$TPDs[[1]] > 0)
                          # what is the probability of occupancy of these cells
                          prob <- sp_TPD$TPDs[[1]][which(sp_TPD$TPDs[[1]] > 0)]
                          
                          # store in dataframe
                          df <- data.frame(presence = presence, probability = prob)
                          
                          return(list(df))
                          
                        }

closeAllConnections()



grid_TPD <- species_TPD(species = study_species_pool[1],
                        trait_ranges = trait_ranges,
                        traits = traits, alpha = 0.68)

# extract grid
evaluation_grid <- grid_TPD$data$evaluation_grid


# which cells when considering all species are occupied 
cells <- unique(unlist(lapply(species_TPDs, FUN = function(x) x[,"presence"])))
# order cells
cells <- cells[order(cells, decreasing = FALSE)]

# filter grid to just the occupied cells
evaluation_grid <- evaluation_grid[cells,]

# create an empty matrix of the occupied cells and the species as columns
species_occupancy_matrix <- matrix(rep(0, length(cells)* length(study_species_pool)), ncol = length(study_species_pool), dimnames = list(as.character(cells),study_species_pool ))


# insert data locating cells which the species are present in and inputting the
# probability of occupancy data 
for(i in 1:length(species_TPDs)){
  
  species_occupancy_matrix[as.character(species_TPDs[[i]][["presence"]]),i] <- species_TPDs[[i]][["probability"]] 
  
}

# bind the trait values and the probability of occupancy data.
TPD_matrix <- cbind(evaluation_grid,species_occupancy_matrix)



TPD_matrix$Trophic_niche <- apply(
  TPD_matrix[, 4:ncol(TPD_matrix)],
  MARGIN = 1,
  FUN = function(x) {
    occ_sp <- names(x)[x > 0]
    
    prop_occ <- Forage %>% dplyr::filter(Birdlife_Name %in% occ_sp) %>%
      dplyr::distinct(Birdlife_Name, Trophic_Niche) %>%
      dplyr::group_by(Trophic_Niche) %>% dplyr::summarise(count = n()) %>% dplyr::ungroup () %>%
      dplyr::mutate(prop = count / sum(count)) %>%
      data.frame()
    
    
    return(ifelse(
      any(prop_occ$prop >= 0.7),
      prop_occ$Trophic_Niche[prop_occ$prop >= 0.7],
      "unclassified"
    ))
    
  }
)




study_TPDs[[study]] <- TPD_matrix

}





write_rds(study_TPDs, file = "Outputs/study_TPDs.rds")

