PREDICTSEcosystemIntegrityTPDs <- function(predicts_path,
                                       randomisations_path,
                                       species_pools_path,
                            ei_traits,
                            avonet_means_path,
                            dir_out){
  
  
  predicts <- qs2::qs_read(predicts_path)
  
  
  traits <- qs2::qs_read(ei_traits)
  
  species_pools <- qs2::qs_read(species_pools_path)
  
  all_species <- unique(unlist(species_pools))
  
  predicts_studies <- unique(as.character(predicts$SS))
  
  avonet_means <- read.csv(avonet_means_path)
  
  
  predicts_trophic <- predicts |>
    dplyr::distinct(Birdlife_Name, Trophic_Niche)
  

  rest_trophic <- avonet_means |>
    dplyr::select(Birdlife_Name = Species1, Trophic_Niche = Trophic.Niche) |>
    dplyr::filter(Birdlife_Name %in% setdiff(all_species, predicts_trophic$Birdlife_Name)) |>
    dplyr::mutate(
      Trophic_Niche = ifelse(Trophic_Niche == "Invertivore", "In", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Vertivore", "Vt", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Nectarivore", "Ne", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Granivore", "Gr", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Omnivore", "Om", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Aquatic predator", "Aq.p", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Scavenger", "Sc", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Frugivore", "Fr", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Herbivore terrestrial", "Hb.T", paste(Trophic_Niche)),
      Trophic_Niche = ifelse(Trophic_Niche == "Herbivore aquatic", "Hb.A", paste(Trophic_Niche))
      )
  
  
  species_trophic <- rbind(predicts_trophic, rest_trophic) |>
   dplyr::mutate(Trophic_Niche = ifelse(is.na(Trophic_Niche), "Om", paste(Trophic_Niche)))
  
  
  randomisations <- qs2::qs_read(randomisations_path)
  
  doParallel::registerDoParallel(cores = 10)
  

  predicts_site_tpdc <-  
    foreach::foreach(x = predicts_studies, .combine = "c", 
                     .packages = c("TPD","tidyr","purrr","dplyr")) %dopar% {
                       
                       study_sites <- predicts |>
                         dplyr::filter(SS == x) |>
                         dplyr::distinct(SSBS) |>
                         dplyr::pull() |>
                         as.character()
                       
                       study_species <- predicts |>
                         dplyr::filter(SS == x) |>
                         dplyr::distinct(Birdlife_Name) |>
                         dplyr::pull()
                       
                       random_species <- unique(unlist(species_pools[study_sites]))
                       
                       
                       study_traits <- 
                         traits$specimen_traits[traits$specimen_traits$Birdlife_Name %in% study_species,]
                       
                       study_trait_ranges <- .CalcTraitRange(traits = study_traits, buffer = 0.05)
                       
                       
                       study_bandwidths <- traits$site_bandwidths |>
                         dplyr::filter(SSBS %in% study_sites) |>
                         dplyr::select(-SSBS) |>
                         colMeans(na.rm = TRUE)
                       
                       comm_mat <- predicts |>
                         dplyr::filter(SS == x) |>
                         dplyr::group_by(SSBS, Birdlife_Name) |>
                         dplyr::summarise(abundance = RelativeStudyAbundance) |>
                         tidyr::pivot_wider(names_from = "Birdlife_Name", values_from = "abundance") |>
                         tibble::column_to_rownames("SSBS")
                       comm_mat[is.na(comm_mat)] <- 0  
                       
                       

                       random_mat <- purrr::map(study_sites, .f = function(y){
                         
                         random_mat <- randomisations[[y]] |>
                           dplyr::select(-Birdlife_Name) |>
                           tidyr::pivot_longer(names_to = "randomisation", values_to = "Birdlife_Name", cols = 2:1001) |>
                           tidyr::pivot_wider(names_from = "Birdlife_Name", values_from = "RelativeStudyAbundance") |>
                           tibble::column_to_rownames("randomisation")
                         
                         random_mat <- t(as.matrix(colMeans(random_mat, na.rm = TRUE)))
                         
                         return(as.data.frame(random_mat))
                         
                       }) |>
                         data.table::rbindlist(fill = TRUE) |>
                         as.data.frame() |>
                         magrittr::set_rownames(study_sites)
                       
                
                       random_mat[is.na(random_mat)] <- 0  
                       
                       
                       
                       random_mat[,setdiff(random_species, colnames(random_mat))] <- 0
                       
                       
                       
                       site_tpds <- eTPDs(species = random_species, 
                                          traits = traits, 
                                          trait_ranges = study_trait_ranges, 
                                          divisions = 50, 
                                          alpha = 0.95,
                                          bandwidths = study_bandwidths)
                       
                       
                      
                        
                      cell_niche <- 
                        purrr::map(1:nrow(site_tpds), .f = function(z){
                       
                        sp <- colnames(site_tpds)[-c(1:3)][which(site_tpds[z,-c(1:3)] > 0)]
                      
                        
                        out <- list(species = paste(sp, collapse = ";"))
                        
                         cell_niche <- data.frame(Birdlife_Name = sp, prob = as.numeric(site_tpds[z,sp])) |>
                           dplyr::left_join(species_trophic, by = "Birdlife_Name") |>
                           dplyr::group_by(Trophic_Niche) |>
                           dplyr::summarise(prob = sum(prob)) |>
                          dplyr::mutate(niche_prob = prob/sum(prob)) |>
                           dplyr::filter(niche_prob >= 0.7) |>
                           dplyr::pull(Trophic_Niche)
                        
                           

                       
                       if(purrr::is_empty(cell_niche)){
                         
                         out <- list(species = paste(sp, collapse = ";"),
                                     trophic_niche = "Unclassified")
                         
                         return(out)
                       }
                         
                         out <- list(species = paste(sp, collapse = ";"),
                                     trophic_niche = cell_niche)
                         
                         return(out)
                         
                       }) 
                       
                       
                      
                      
                      
                       cell_niche <- data.frame(site_tpds[,c(1:3)], 
                                                Trophic_Niche = sapply(1:length(cell_niche), FUN = function(x) cell_niche[[x]][["trophic_niche"]]),
                                                Occupying_species = sapply(1:length(cell_niche), FUN = function(x) cell_niche[[x]][["species"]]))
                       colnames(cell_niche)[1:3] <- c("locomotion","foraging","body_size")
                       
                       
                       random_site_tpdc <- eTPDc(eTPDs = site_tpds, comm = random_mat, normalise = FALSE)
                       random_site_tpdc <- random_site_tpdc |> dplyr::left_join(cell_niche)
                       
                       obvs_site_tpdc <- eTPDc(eTPDs = site_tpds, comm = comm_mat, normalise = FALSE)
                       obvs_site_tpdc <- obvs_site_tpdc |> dplyr::left_join(cell_niche)
                       
                       
                       table(obvs_site_tpdc$Trophic_Niche)
                       site_tpdc <- list(list(observed = obvs_site_tpdc,
                                         random = random_site_tpdc))
                       
                       
                       names(site_tpdc) <- x
                       
                       
                       return(site_tpdc)
                       
                     }
  
  foreach::registerDoSEQ()
  
  
  predicts_random_tpdc <- lapply(predicts_site_tpdc,FUN = function(x) x[["random"]])
  predicts_site_tpdc <- lapply(predicts_site_tpdc, FUN = function(x) x[["observed"]])
  
  
  file.out_random <- file.path(dir_out, "predicts-randomisations-tpdc.qs") 
  file.out_observed <- file.path(dir_out, "predicts-sites-tpdc.qs") 
  
  
  qs2::qs_save(predicts_random_tpdc, file.out_random)
  qs2::qs_save(predicts_site_tpdc, file.out_observed)
  
  
  return(c(file.out_observed,file.out_random))
  
}