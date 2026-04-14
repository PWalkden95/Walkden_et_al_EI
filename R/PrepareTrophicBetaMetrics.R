PrepareTrophicBetaMetrics <- function(predicts_combinations_path,
                                        predicts_site_tpd_path,
                                        dir_out){
  
  predicts_combinations <- qs2::qs_read(predicts_combinations_path)
  
  predicts_tpds <- qs2::qs_read(predicts_site_tpd_path)
  
  comb_studies <- unique(predicts_combinations$SS)
  
  

  trophic_beta_metrics <- purrr::map_df(comb_studies, .f = function(x){
    
    comb_tpd <- predicts_tpds[[x]]
    
    combs <- predicts_combinations |>
      dplyr::filter(SSBS_1 %in% colnames(comb_tpd)|
                      SSBS_2 %in% colnames(comb_tpd))
    if(nrow(combs) == 0){
      return(NULL)
    }
    
    output <- eTPDFunctionalMetricsBeta(eTPDc = comb_tpd,combinations = combs, trophic_niche = TRUE) 
    
    return(output)
  })
  
  
  
  file.out <- file.path(dir_out, "trophic-beta-metrics.qs")
  
  
  qs2::qs_save(trophic_beta_metrics, file.out)  
  
  return(file.out)
  
  
}
