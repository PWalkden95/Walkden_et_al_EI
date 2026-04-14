PREDICTSEIRealmTPDs <- function(predicts_path,
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
  
  
  realms <- levels(predicts$Realm)
  
  doParallel::registerDoParallel(cores = length(realms))
    
    
  predicts_realm_tpdc <-  
    foreach::foreach(x = realms, .combine = "c", 
                     .packages = c("TPD","tidyr","purrr","dplyr")) %dopar% {
                       
                       realm_sites <- predicts |>
                         dplyr::filter(Realm == x) |>
                         dplyr::distinct(SSBS) |>
                         dplyr::pull() |>
                         as.character()
                       
                       realm_species <- predicts |>
                         dplyr::filter(Realm == x) |>
                         dplyr::distinct(Birdlife_Name) |>
                         dplyr::pull()
                       
                       random_species <- unique(unlist(species_pools[realm_sites]))
                       
                       
                       realm_traits <- 
                         traits$specimen_traits[traits$specimen_traits$Birdlife_Name %in% realm_species,]
                       
                       realm_trait_ranges <- .CalcTraitRange(traits = realm_traits, buffer = 0.05)
                       
                       
                       realm_bandwidths <- traits$site_bandwidths |>
                         dplyr::filter(SSBS %in% realm_sites) |>
                         dplyr::select(-SSBS) |>
                         colMeans(na.rm = TRUE)
                       
                       comm_mat <- predicts |>
                         dplyr::filter(Realm == x) |>
                         dplyr::group_by(Predominant_habitat, Birdlife_Name) |>
                         dplyr::reframe(abundance = sum(RelativeStudyAbundance)) |>
                         tidyr::pivot_wider(names_from = "Birdlife_Name", values_from = "abundance") |>
                         tibble::column_to_rownames("Predominant_habitat")
                       comm_mat[is.na(comm_mat)] <- 0  
                       
                       
                       land_uses <- rownames(comm_mat)
                       
                       random_mat <- purrr::map_dfr(land_uses, .f = function(y){
                         
                         realm_lu_sites <- predicts |>
                           dplyr::filter(Realm == x, Predominant_habitat == y) |>
                           dplyr::distinct(SSBS) |>
                           dplyr::pull() |> as.character()
                         
                         realm_lu_mat <- purrr::map(realm_lu_sites, function(z){
                           
                           random_mat <- randomisations[[z]] |>
                             dplyr::select(-Birdlife_Name) |>
                             tidyr::pivot_longer(names_to = "randomisation", values_to = "Birdlife_Name", cols = 2:1001) |>
                             tidyr::pivot_wider(names_from = "Birdlife_Name", values_from = "RelativeStudyAbundance") |>
                             tibble::column_to_rownames("randomisation")
                           
                           random_mat <- t(as.matrix(colMeans(random_mat, na.rm = TRUE)))
                           
                           return(as.data.frame(random_mat))
                           
                         }) |>
                           data.table::rbindlist(fill = TRUE) |>
                           as.data.frame() |>
                           colSums(na.rm = TRUE)
                           
                         
                         }) |>
                         as.data.frame()
                         
                         
                       rownames(random_mat) <- land_uses
                       
                       random_mat[is.na(random_mat)] <- 0  
                       
                       
                       
                       random_mat[,setdiff(random_species, colnames(random_mat))] <- 0
                       
                       
                       
                       realm_tpds <- eTPDs(species = random_species, 
                                          traits = traits, 
                                          trait_ranges = realm_trait_ranges, 
                                          divisions = 50, 
                                          alpha = 0.95,
                                          bandwidths = study_bandwidths)
                       
                       
                       
                       
                       cell_niche <- 
                         purrr::map(1:nrow(realm_tpds), .f = function(z){
                           
                           sp <- colnames(realm_tpds)[-c(1:3)][which(realm_tpds[z,-c(1:3)] > 0)]
                           
                           
                           cell_niche <- data.frame(Birdlife_Name = sp, prob = as.numeric(realm_tpds[z,sp])) |>
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
                       
                       
               
                       cell_niche <- data.frame(realm_tpds[,c(1:3)], 
                                                Trophic_Niche = sapply(1:length(cell_niche), FUN = function(x) cell_niche[[x]][["trophic_niche"]]),
                                                Occupying_species = sapply(1:length(cell_niche), FUN = function(x) cell_niche[[x]][["species"]]))
                       colnames(cell_niche)[1:3] <- c("locomotion","foraging","body_size")
                       
                       random_realm_tpdc <- eTPDc(eTPDs = realm_tpds, comm = random_mat, normalise = FALSE)
                       random_realm_tpdc <- random_realm_tpdc |> dplyr::left_join(cell_niche)
                       
                       obvs_realm_tpdc <- eTPDc(eTPDs = realm_tpds, comm = comm_mat, normalise = FALSE)
                       obvs_realm_tpdc <- obvs_realm_tpdc |> dplyr::left_join(cell_niche)
                       
                       
                       
                       realm_tpdc <- list(list(observed = obvs_realm_tpdc,
                                              random = random_realm_tpdc))
                       
                       
                       names(realm_tpdc) <- x
                       
                       
                       return(realm_tpdc)
                       
                     }
  
  foreach::registerDoSEQ()
  
  
  predicts_random_tpdc <- lapply(predicts_realm_tpdc,FUN = function(x) x[["random"]])
  predicts_realm_tpdc <- lapply(predicts_realm_tpdc, FUN = function(x) x[["observed"]])
  
  
  file.out_random <- file.path(dir_out, "predicts-randomisations-realm-tpdc.qs") 
  file.out_observed <- file.path(dir_out, "predicts-realm-tpdc.qs") 
  
  
  qs2::qs_save(predicts_random_tpdc, file.out_random)
  qs2::qs_save(predicts_realm_tpdc, file.out_observed)
  
  
  return(c(file.out_observed,file.out_random))
  
}