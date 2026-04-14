PREDICTSEITraits <- function(predicts_path,
                             randomisations_path,
                             avonet_specimen_data,
                             avonet_mean_data,
                             dir_out
) {
  

  randomisations <- qs2::qs_read(randomisations_path)
    
  all_species <- unique(unlist(randomisations))
  
  ei_traits <- DefineTraitAxes(species = all_species,
                  predicts_avonet_path = predicts_path, 
                  avonet_specimen_data = avonet_specimen_data,
                  avonet_mean_data = avonet_mean_data)
  
  
  
  
  file.out <- file.path(dir_out, "ei-traits-prepared.qs")
  
  qs2::qs_save(ei_traits, file = file.out)
  
  return(file.out)
  
  
}