RealmTPDPlots <- function(realm_tpdc,
                          dir_out){
  
  
  if(!dir.exists(file.path(dir_out,"realm_tpd_plots"))){
    dir.create(file.path(dir_out,"realm_tpd_plots"))
  }
  
  
  realm_data <- qs::qread(realm_tpdc)
  
  
  realms <- names(realm_data)
  
  
  niche_colours <- data.frame(Trophic_Niche = c("Fr","Gr","Ne","In","Vt","Aq.p","Sc","Hb.A","Hb.T","Om","Unclassified"),
                        colour = c("#ef4136",
                                   "#c9842a",
                                   "#87d1d2",
                                   "#164970",
                                   "purple1",
                                   "#aaddf9",
                                   "goldenrod4",
                                   "#007173",
                                   "#0d703b",
                                   "#b77ab4",
                                   "grey25"))
  
  
  land_uses <- c("Primary vegetation",
                 "Secondary",
                 "Plantation forest",
                 "Pasture",
                 "Cropland",
                 "Urban")
  
  
  
  files.out <- purrr::map(realms, function(x){
    
    
    if(!dir.exists(file.path(dir_out,"realm_tpd_plots",x))){
      dir.create(file.path(dir_out,"realm_tpd_plots",x))
    }
    
    
    data <- realm_data[[x]]
    
    plot_limits =    list(xmin = min(data[,1]) - 0.25,
                          xmax = max(data[,1]) + 0.25,
                          ymin = min(data[,2]) - 0.25,
                          ymax = max(data[,2]) + 0.25,
                          zmin = min(data[,3]) - 0.25,
                          zmax = max(data[,3]) + 0.25)
    
    
    
    
    all_probs <- as.matrix(data[,c(land_uses)])
    all_probs <- all_probs[all_probs > 0]
    
    cap <- quantile(all_probs, 0.95)
    
    
    size_values <- scale(all_probs)
    
    data[,land_uses] <- sapply(land_uses, FUN = function(z){
      out <- ifelse(data[,z] > cap, cap, data[,z])
      
      
      out <- (out - attr(size_values, "scaled:center"))/attr(size_values,"scaled:scale")
      
      out <- (out + (abs(min(out))) + 0.5) /1.4
      
      return(out)
      
    })
    
 
    min_val <- min(as.matrix(data[,land_uses]))
    
    
    out <- purrr::map_vec(land_uses[land_uses %in% colnames(data)], function(y){
      
    
      lu_data <- data |>
        dplyr::select(1:3, dplyr::contains(y), "Trophic_Niche") |>
        dplyr::left_join(niche_colours, by = "Trophic_Niche") |>
        dplyr::filter(!!rlang::sym(y) > min_val) |>
        dplyr::filter(Trophic_Niche %in% c("Gr","Fr","Ne","In","Om"))
      
      save.out <- file.path(dir_out,"realm_tpd_plots",x,paste0(y,".png"))
      
      
   
      
      TPD3DPlot(eTPDc = lu_data,
                save = TRUE,
                filename = save.out,
                limits = plot_limits,
                scale_size = FALSE)
      
      
      
      return(save.out)
    })
    
    return(out)
    
    
    
  }) |> unlist()
  
  
  return(files.out)
  
}