SiteCombinations <- function(predicts_path,
                             dir_out){
  
 
  
  predicts <- qs2::qs_read(predicts_path) |>
    dplyr::distinct(SS,SSBS,SSB,Predominant_habitat,Realm)
  
predicts_studies <- unique(predicts$SS)  


site_combinations <- purrr::map(predicts_studies, .f = function(x){
    
    study_sites <- predicts |>
      dplyr::filter(SS == x)
    
    primary_sites <- study_sites |>
      dplyr::filter(Predominant_habitat == "Primary vegetation")
    
    
    if(nrow(primary_sites) == 0){
      return(NA)
    }
    
    combinations <- expand.grid(primary_sites$SSBS, study_sites$SSBS) |>
      magrittr::set_colnames(c("SSBS_1","SSBS_2")) |>
      dplyr::filter(SSBS_1 != SSBS_2) |>
      dplyr::left_join(study_sites[,c("SSBS","Predominant_habitat")], by = c("SSBS_2" = "SSBS")) |>
      dplyr::rename(Predominant_habitat_SSBS_2 = Predominant_habitat) |>
      dplyr::left_join(study_sites[,c("SSBS","SS")], by = c("SSBS_2" = "SSBS")) |>
      dplyr::left_join(study_sites[,c("SSBS","SSB")], by = c("SSBS_2" = "SSBS")) |>
      dplyr::rename(SSB_SSBS_2 = SSB) |>
      dplyr::left_join(study_sites[,c("SSBS","SSB")], by = c("SSBS_1" = "SSBS")) |>
      dplyr::rename(SSB_SSBS_1 = SSB)  |>
      dplyr::left_join(study_sites[,c("SSBS","Realm")], by = c("SSBS_2" = "SSBS"))
    
    
    return(combinations)

    
  }) |>
  Reduce(f = "rbind") |>
  na.omit()


file.out <- file.path(dir_out, "site_combinations.qs")

qs2::qs_save(site_combinations, file.out)

return(file.out)
   
}
