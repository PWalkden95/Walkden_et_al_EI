PrepareTrophicAlphaMetrics <- function(predicts_path,
                                         predicts_site_tpd_path,
                                         dir_out){
  
  predicts <- qs2::qs_read(predicts_path) |>
    dplyr::distinct(SS,SSB,SSBS, Predominant_habitat, Realm)
  
  predicts_tpds <- qs2::qs_read(predicts_site_tpd_path)
  
  trophic_alpha_metrics <- purrr::map_dfr(predicts_tpds, .f = function(x){
    eTPDFunctionalMetricsAlpha(eTPDc = x, trophic_niche = TRUE)
  })
  
  
  TrophicAlphaMetrics <- predicts |>
    dplyr::left_join(trophic_alpha_metrics)
  
  
  file.out <- file.path(dir_out, "trophic-alpha-metrics.qs")
  
  
  qs2::qs_save(TrophicAlphaMetrics, file.out)  
  
  return(file.out)
  
  
}
