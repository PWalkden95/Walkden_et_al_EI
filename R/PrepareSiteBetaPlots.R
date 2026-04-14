PrepareSiteBetaPlots <- function(site_level_analysis,
                                  beta_data,
                                  dir_out){
  
  
  
  
  if(!dir.exists(file.path(dir_out, "site_level_plots"))){
    dir.create(file.path(dir_out, "site_level_plots"))
  }
  
  
  
  analysis <- qs::qread(site_level_analysis)
  
  beta_metrics <- qs::qread(beta_data)
  
  
  fig.out <- file.path(dir_out,"site_level_plots","site-level-similarity.png") 
  
  
  .SiteLevelPlot(raw_data = beta_metrics,
                 bootstrap_data = analysis$`similarity Univariate`$Bootstraps,
                 metric = "similarity",
                 filename = fig.out, limits = c(-0.3,1.55))
  

  fig.out <- file.path(dir_out,"site_level_plots","site-level-centroid.png") 
  
  .SiteLevelPlot(raw_data = beta_metrics,
                 bootstrap_data = analysis$`centroid Univariate`$Bootstraps,
                 metric = "centroid_shift",
                 filename = fig.out)
  
    
  
  
  realms <- levels(beta_metrics$Realm)
  
  
  purrr::map(realms, .f = function(x){
    
    fig.out <- file.path(dir_out,"site_level_plots",glue::glue("site-level-similarity-{x}.png")) 
    
    realm_alpha <- beta_metrics |>
      dplyr::filter(Realm == x)
    
    
    realm_boot <- analysis$`similarity Multivariate`$Bootstraps |>
      dplyr::filter(Realm == x)
    
    .SiteLevelPlot(raw_data = realm_alpha,
                   bootstrap_data = realm_boot,
                   metric = "similarity",
                   filename = fig.out)
    
  })
  
  
  purrr::map(realms, .f = function(x){
    
    fig.out <- file.path(dir_out,"site_level_plots",glue::glue("site-level-centorid-{x}.png")) 
    
    realm_alpha <- beta_metrics |>
      dplyr::filter(Realm == x)
    
    
    realm_boot <- analysis$`centroid Multivariate`$Bootstraps |>
      dplyr::filter(Realm == x)
    
    .SiteLevelPlot(raw_data = realm_alpha,
                   bootstrap_data = realm_boot,
                   metric = "centroid_shift",
                   filename = fig.out)
    
    
  })
  
  
}