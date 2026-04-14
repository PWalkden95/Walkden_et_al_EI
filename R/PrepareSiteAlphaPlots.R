PrepareSiteAlphaPlots <- function(site_level_analysis,
                                  alpha_data,
                                  dir_out){
  
  
  
  
  if(!dir.exists(file.path(dir_out, "site_level_plots"))){
    dir.create(file.path(dir_out, "site_level_plots"))
  }
  
  
  
  analysis <- qs::qread(site_level_analysis)
  
  alpha_metrics <- qs::qread(alpha_data)
  

  fig.out <- file.path(dir_out,"site_level_plots","site-level-alpha.png") 
  
  
  .SiteLevelPlot(raw_data = alpha_metrics,
                 bootstrap_data = analysis$Univariate$Bootstraps,
                 metric = "FRich",
                 limits = c(-0.2,1.75),
                 filename = fig.out)
  
  
  
  
  realms <- levels(alpha_metrics$Realm)
  
  
    purrr::map(realms, .f = function(x){
    
    fig.out <- file.path(dir_out,"site_level_plots",glue::glue("site-level-alpha-{x}.png")) 
      
    realm_alpha <- alpha_metrics |>
      dplyr::filter(Realm == x)
    
    
    realm_boot <- analysis$Mulitvariate$Bootstraps |>
      dplyr::filter(Realm == x)
    
    .SiteLevelPlot(raw_data = realm_alpha,
                   bootstrap_data = realm_boot,
                   metric = "FRich",
                   filename = fig.out)
    
      
  })
  
      
}