SiteRandomisations <- function(predicts_path,
                           dispersal_probabilities,
                           environmental_affinity,
                           species_pools_path,
                           dir_out){
  
  
 predicts <- qs2::qs_read(predicts_path)
  
 
 species_pools <- qs2::qs_read(species_pools_path)
 
 sites <- names(species_pools)
 
  disperal_prob <- qs2::qs_read(dispersal_probabilities)

 
 environmental_prob <- qs2::qs_read(environmental_affinity)

 
 site_randomisations <-  purrr::map(sites, function(x){
   
   data <- predicts |>
     dplyr::filter(SSBS == x) |>
     dplyr::select(RelativeStudyAbundance, Birdlife_Name)
   
   
   sp_pool <- species_pools[[x]]
   
   dis_p <- disperal_prob[[x]]
   environ_p <- environmental_prob[[x]]
   
  
   species_weights <- data.frame(Birdlife_Name = sp_pool,
                                 weights = dis_p*environ_p) 
   
   
   
    randomisations <- purrr::map_dfc(1:1000, .f = function(y){
     
     
     random_df <- data.frame(randomised_spp = sample(species_weights$Birdlife_Name,
                                                         replace = FALSE,
                                                         prob = species_weights$weights, size = nrow(data)))
     colnames(random_df) <- glue::glue("randomisation_{y}")
     
     return(random_df)
   })
   
 
    randomisations <- cbind(data,randomisations)
        
    
    return(randomisations)
 })
 
 names(site_randomisations) <- sites
 
 
 file.out <- file.path(dir_out, "site-randomisations.qs")

 qs2::qs_save(site_randomisations, file.out)

 return(file.out)
   
}