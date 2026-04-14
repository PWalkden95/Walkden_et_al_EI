PrepareSiteLevelAlphaMetrics <- function(predicts_path,
                                         predicts_site_tpd_path,
                                         dir_out){
  
  predicts <- qs2::qs_read(predicts_path) |>
    dplyr::distinct(SS,SSB,SSBS, Predominant_habitat, Realm)
  
  predicts_tpds <- qs2::qs_read(predicts_site_tpd_path)

site_alpha_metrics <- purrr::map_dfr(predicts_tpds, .f = function(x){
  eTPDFunctionalMetricsAlpha(eTPDc = x, trophic_niche = FALSE)
})
  
  
SiteAlphaMetrics <- predicts |>
  dplyr::left_join(site_alpha_metrics)


file.out <- file.path(dir_out, "site-alpha-metrics.qs")


qs2::qs_save(SiteAlphaMetrics, file.out)  

return(file.out)

  
}
