SiteHolesAnalysis <- function(site_level_hole_metrics, dir_out){
  
  
  hole_metrics <- qs2::qs_read(site_level_hole_metrics) |>
    dplyr::mutate(convex_hull_volume = convex_hull_volume/max(convex_hull_volume),
                  mean_hole_volume = ifelse(is.nan(mean_hole_volume), 0, mean_hole_volume),
                  normal_mean_hole = mean_hole_volume/convex_hull_volume,
                  normal_hole_number = number_of_holes/convex_hull_volume) 
  
  
  
  
  hole_metrics$Predominant_habitat <- factor(hole_metrics$Predominant_habitat, levels = c("Primary vegetation",
                                                                                          "Secondary",
                                                                                          "Plantation forest",
                                                                                          "Pasture",
                                                                                          "Cropland",
                                                                                          "Urban"))
  
  
  
  hist(hole_metrics$normal_mean_hole)
  hist(log(hole_metrics$normal_mean_hole + 1))
  hist(sqrt(hole_metrics$normal_mean_hole))
  
  
  mod_vol <- glmmTMB::glmmTMB(normal_mean_hole ~ Predominant_habitat  +
                                (1|SS) + (1|SSB),ziformula = ~1,
                              data = hole_metrics, family = glmmTMB::ziGamma(link = "log"))
  
  
  summary(mod_vol)  
  

  #----------------------------------------------------
  # Data for model visualisations through bootstrapping 
  #------------------------------------------------------
  
  # the best way to visulaise mixed effects models is to perform bootstraping, 
  # additionally this is can assess the significance of the effect variables.
  
  # first will need blank data for all the biogeographic realm - land use 
  # combinations with sufficient data 
  
  
  boot_newdata <- hole_metrics |>
    dplyr::distinct(Predominant_habitat)
  
  
 
  
  ## run the function to bootstrap the three models

  volume_model_bootstraps <-
    .MixedModelBootstrapping(
      model = list(mod_vol),
      newdata = boot_newdata,
      transformation = c("log"),
      nboot = 250,
      n_cores = 1
    )




  volume_output <- list(Model = mod_vol,
                        Bootstraps = volume_model_bootstraps)



  #------------------------
  # hole number models
  #------------------------


  hist(hole_metrics$normal_hole_number)
  hist(sqrt(hole_metrics$normal_hole_number))
  hist(log1p(hole_metrics$normal_hole_number))




  mod_num <- glmmTMB::glmmTMB(normal_hole_number ~ Predominant_habitat  +
                                (1|SS) + (1|SSB),ziformula = ~1,
                              data = hole_metrics,family = glmmTMB::ziGamma(link = "log"))


  summary(mod_num)


  #----------------------------------------------------
  # Data for model visualisations through bootstrapping
  #------------------------------------------------------

  # the best way to visulaise mixed effects models is to perform bootstraping,
  # additionally this is can assess the significance of the effect variables.


  ## run the function to bootstrap the three models

  number_model_bootstraps <-
    .MixedModelBootstrapping(
      model = list(mod_num),
      newdata = boot_newdata,
      transformation = c("log"),
      nboot = 250,
      n_cores = 1
    )



  number_output <- list(Model = mod_num,
                        Bootstraps = number_model_bootstraps)


  output <- list(`hole volume` = volume_output,
                 `hole number` = number_output)


  file.out <- file.path(dir_out, "site-holes-analysis.qs")


  qs2::qs_save(output, file.out)


  return(file.out)
}
