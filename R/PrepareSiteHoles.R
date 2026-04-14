PrepareSiteHoles <- function(predicts_path,
                                         predicts_site_tpd_path,
                            predicts_randomisation_tpd_path,
                                         dir_out){
  
  predicts <- qs2::qs_read(predicts_path) |>
    dplyr::distinct(SS,SSB,SSBS, Predominant_habitat, Realm)
  
  predicts_tpds <- qs2::qs_read(predicts_site_tpd_path)
  
  predicts_random_tpds <- qs2::qs_read(predicts_randomisation_tpd_path)
  
  
  site_holes <- purrr::map2_df(.x = predicts_tpds, .y = predicts_random_tpds, .f = function(x,y){
    
    #25 min points
    #density 1.5 mini point 16
    eTPDc_holes(eTPDc = x, eTPDc_Null = y, threshold = 0.95, minimum_points = 16, density = 1.5)
    
  })
  
  
  SiteHoleMetrics <- predicts |>
    dplyr::left_join(site_holes)
  
  
  file.out <- file.path(dir_out, "site-hole-metrics.qs")
  
  
  qs2::qs_save(SiteHoleMetrics, file.out)  
  
  return(file.out)
  
  
}
