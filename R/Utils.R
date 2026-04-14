
########################## create blank matrix row
.BlankMatrix <- function(rows, columns) {
  mat <-
    matrix(
      rep(NA, length(rows) * length(columns)),
      nrow = length(rows),
      ncol = length(columns),
      dimnames = list(as.character(rows), as.character(columns))
    )
  
  return(mat)
}


########## combine spatial polygons ready for analysis

.SpatialCombinePolygons <-
  function(geometry) {
    if (any(
      class(geometry$Shape)[1] == "sfc_MULTISURFACE",
      class(geometry$Shape)[1] == "sfc_GEOMETRY"
    )) {
      for (k in 1:NROW(geometry)) {
        geometry$Shape[[k]] <- st_cast(geometry$Shape[[k]], "MULTIPOLYGON")
      }
    }
    
    
    shape <-
      st_combine(geometry$Shape)
    
    
    
    st_crs(shape) <-
      "+proj=longlat +datum=WGS84 +no_defs"
    
    return(shape)
  }


####################### filter polygons to be the right presence, origin and seasonal categories

.FilterPolygons <-
  function(BL_data,
           presence_codes,
           seasonal_codes,
           origin_codes) {
    require(sf)
    data <-
      BL_data |> dplyr::filter(presence %in% presence_codes,
                                origin %in% origin_codes,
                                seasonal %in% seasonal_codes)
    return(data)
  }




## calculate similarity of species between cells
.SpeciesCellSimilarity <- function(species,cell, species_in_all_cells){
  
  if(purrr::is_empty(species_in_all_cells[[cell]])) {
    return(NA)
  } else {
    
    other_spp <- species_in_all_cells[[cell]]
    
    ### if the cell is terrestrial but has no species in it then the other species is none
    
    if (purrr::is_empty(other_spp)) {
      other_spp <- "none"
    }
    
    ### how many species are shared between the two cells
    
    sim <- length(which(species %in% other_spp))
    
    # if none similarity is 0
    
    
    return(ifelse(sim != 0, sim / length(species), 0))
  }
}



.SiteFind <- function(species,species_pools){
  
  
  sp_sites <- names(species_pools)[purrr::map_vec(1:length(species_pools), .f = function(x){
    species %in% species_pools[[x]]
  })]
  
  return(sp_sites)
}


#=============================================================================
# FUNCTION: Bootstraps mixed effects models 
#
#INPUT: Model to be bootstrapped; newdata to be bootstrapped; transformation
# of the response variable if needed; how many bootstraps
# OUTPUTS: The bootstrapped model for each of the 
#=============================================================================

.MixedModelBootstrapping <-
  function(model = list(),
           newdata,
           transformation = c("none", "log", "sqrt","logit","log+1"),
           nboot,
           n_cores) {
    
    
    # the newdata to be bootstrapped
    
    data <- newdata
    

    ## if multiple models given to the function iterate through them 
  
      full_boot_models <- purrr::map_dfr(1:length(model), .f = function(i){
      
      # extract the response variable
      
      FD_met <-
        ifelse(class(model[[i]]) == "glmmTMB",
               as.character(formula(model[[i]]))[2],
               unlist(stringr::str_split(as.character(model[[i]]@call)[2], pattern = " ~ "))[1]) |>
        stringr::str_squish()
      
      data$estimate <- predict(model[[i]], newdata = data, re.form = NA)
      
if(class(model[[i]]) == "glmmTMB"){
  mod_frame <- model[[i]]$frame 
  
  # perform the bootstraps 
  bootstraps <-  purrr::map_dfr(1:nboot, function(y){
    
    boot_dat <- mod_frame[sample(nrow(mod_frame), size=nrow(mod_frame),replace=TRUE),]  
    
    
    boot_mod <- update(model[[i]], data = boot_dat)   
    
    
    boot <- as.data.frame(matrix(predict(boot_mod, newdata = data, re.form = NA), nrow = 1))
    
    return(boot)
    
  })
  
} else {
  mod_frame <- model[[i]]@frame
  
  ## function to do the prediction - this is required for teh bootstrapping
  
  predict_function <- function(x) {
    as.numeric(predict(x, newdata = data, re.form = NA,allow.new.levels=TRUE))
  }
  # perform the bootstraps 
  booted_mod <-
    lme4::bootMer(
      model[[i]],
      FUN = function(x)
        predict_function(x),
      nsim = nboot, seed = 1234, parallel = "multicore", ncpus = n_cores
    )
  
  bootstraps <- booted_mod$t
  
}

      
      
      ## depending on the transformation back transform the output
      
        if (transformation[i] == "log") {
          data$estimate <- exp(data$estimate)
          bootstraps <- exp(bootstraps)
        } else {
          if (transformation[i] == "log+1") {
            data$estimate <- exp(data$estimate) -1
            bootstraps <- exp(bootstraps) -1
          } else {
            if (transformation[i] == "sqrt") {
              data$estimate <- data$estimate ^ 2
              bootstraps <- bootstraps ^ 2
            } else {
              if(transformation[i] == "logit"){
                data$estimate <- arm::invlogit(data$estimate)
                bootstraps <- arm::invlogit(bootstraps)
              }
            }
            
          }
        }
      
      
      ## we don't want all the bootstraps so we are just going to take 95% 
      ## of the data so removing the top and bottom 2.5% of data
      
      
      ninety_five_boot <- apply(bootstraps, MARGIN = 2, FUN = function(x){
        
        
        range <- as.numeric(quantile(x,c(0.025,0.975)))
        

        return(x[dplyr::between(x,range[1],range[2])])
           
      })
      

      

      
      present_columns <- grep(
        pattern = "Predominant_habitat|realm|Trophic_niche|Predominant_habitat_SSBS_2",
        colnames(data),
        ignore.case = TRUE, value = TRUE
      )
      
      data_column_assignments <- which(colnames(data) %in% present_columns)
      
      data <- data.frame(data)
      if(length(data_column_assignments) == 1){
        
        
        
        boot_colnames <- as.character(data[, data_column_assignments])
      } else {
        boot_colnames <- paste(data[, data_column_assignments[1]],
                               data[, data_column_assignments[2]], sep = ";")
        
      }
      
      
      
      ## assign the column names indicating the realm and land use combo
      ninety_five_boot <-
        data.frame(ninety_five_boot) |> 
        magrittr::set_colnames(boot_colnames) 
      
      
      
      # pivot long to have a single column for the bootstrapped values enabling easier plotting
      boot_plot <-
        ninety_five_boot |>
        tidyr::pivot_longer(cols = colnames(ninety_five_boot),
                                          names_to = "levels") |>
        data.frame()
      
      # split name to isolate realm
      
    
      
      if(length(present_columns) > 1){
        
        if("Realm" %in% present_columns|
           "realm" %in% present_columns){
          
          boot_plot[,grep(present_columns, pattern = "Realm|realm", value = TRUE)] <-
            unlist(stringr::str_split(boot_plot[, "levels"], pattern = ";"))[seq(grep(present_columns, pattern = "realm|Realm"), nrow(boot_plot) *
                                                                          2, 2)]
        } else {
          boot_plot$Trophic_niche <-
            unlist(stringr::str_split(boot_plot[, "levels"], pattern = ";"))[seq(grep(present_columns,pattern = "Trophic_niche"), nrow(boot_plot) *
                                                                          2, 2)]
        }
        
        
        boot_plot[,grep(present_columns, pattern = "land_use|Predominant_habitat|Predominant_habitat_SSBS_2", value = TRUE)] <-
          unlist(stringr::str_split(boot_plot[, "levels"], pattern = ";"))[seq(grep(present_columns, pattern = "land_use|Predominant_habitat|Predominant_habitat_SSBS_2"), nrow(boot_plot) *
                                                                        2, 2)]
        
      } else {
        colnames(boot_plot)[1] <- present_columns
      }
      
      
      
      
      
      boot_plot$model <- FD_met
      
      return(boot_plot)
      
    })
    
    return(full_boot_models)
  }

## site level plot




.SiteLevelPlot <- function(raw_data,
                           bootstrap_data,
                           metric,
                           limits = c(),
                           filename){
  
  
  land_uses <-
    c(
      "Primary vegetation",
      "Secondary",
      "Plantation forest",
      "Pasture",
      "Cropland",
      "Urban"
    )
  
  land_use_colours <-
    data.frame(
      land_use = land_uses,
      colours = rev(c(
        "#6d6e71",
        "#d1cc1d",
        "#fbb040",
        "#8dc63f",
        "#39b54a",
        "#225816"
      )
      ))
  rownames(land_use_colours) <- land_uses
  
  
  colnames(raw_data)[grep(colnames(raw_data), pattern = "Predominant_habitat")] <-
    "Predominant_habitat"
  
  colnames(raw_data)[grep(colnames(raw_data), pattern = "SSB_SSBS_2")] <-
    "SSB"
  
  background_points <-  raw_data |>
    dplyr::group_by(SS,SSB) |>
    dplyr::mutate(!!rlang::sym(metric) := as.numeric(!!rlang::sym(metric)),
                  med = median(as.numeric(!!rlang::sym(metric))[Predominant_habitat == "Primary vegetation"])) |>
    dplyr::ungroup() |>
    dplyr::mutate(relative_metric = as.numeric(!!rlang::sym(metric))/ med,
                  Predominant_habitat = factor(Predominant_habitat,
                                               levels = land_uses)) |>
    na.omit() |> 
    dplyr::group_by(Predominant_habitat) |>
    dplyr::slice_sample(n = 500)
  
  
  ninety_five <- as.numeric(quantile(background_points$relative_metric, 0.95, na.rm = TRUE))
  
  background_points <- background_points |>
    dplyr::filter(relative_metric < ninety_five)
  
  
  colnames(bootstrap_data)[grep(colnames(bootstrap_data), pattern = "Predominant_habitat")] <-
    "Predominant_habitat"
  
  bootstrap_data <- bootstrap_data |>
    dplyr::mutate(relative_value = value/median(value[Predominant_habitat == "Primary vegetation"]))
  
  
  bootstrap_data$Predominant_habitat <- factor(bootstrap_data$Predominant_habitat,
                                               levels = land_uses)
  
  
  metric_plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
    ggplot2::geom_point(
      data = background_points,
      ggplot2::aes(
        x = Predominant_habitat,
        y = relative_metric,
        group = Predominant_habitat,
        fill = Predominant_habitat
      ),
      size =  3,
      alpha = 0.8,
      position = ggplot2::position_jitter(0.3),
      show.legend = FALSE,
      pch = 21,
      colour = "grey50"
    ) +
    ggplot2::geom_boxplot(
      data = bootstrap_data,
      width = 0.4,
      ggplot2::aes(fill = Predominant_habitat, x = Predominant_habitat, y = relative_value),
      size = 0.9,
      show.legend = FALSE,
      outlier.shape = NA,
      colour = "black"
    ) +
    ggplot2::scale_colour_manual(values = land_use_colours[land_uses, "colours"]) +
    ggplot2::scale_fill_manual(values = land_use_colours[land_uses, "colours"]) +
    # ggplot2::scale_y_continuous(limits = c(-0.5, 3.5),
    #                    breaks = c(0.5, 1, 2,3)) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = 'white', color = 'white'),
      panel.grid.major = ggplot2::element_line(color = 'white'),
      panel.grid.minor = ggplot2::element_line(color = 'white'),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linetype = "solid"),
      #axis.ticks.x = element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
  
  
  if(!purrr::is_null(limits)){
    metric_plot <- metric_plot +
      ggplot2::scale_y_continuous(limits = limits)
  }
  
  
  ggplot2::ggsave(
    metric_plot,
    filename =
      filename,
    width = 200,
    height = 150,
    units = "mm",
    dpi = 1000
  )
  
}

