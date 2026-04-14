PrepareTrophicBetaPlots <- function(trophic_beta_analysis,
                                    dir_out){
  

  if(!dir.exists(file.path(dir_out, "trophic_plots"))){
    dir.create(file.path(dir_out, "trophic_plots"))
  }
  

  trophic_beta_data <- qs::qread(trophic_beta_analysis)
  
  
  effect_sizes <- trophic_beta_data$`Marginal effects`
  

  
  # store land uses as vector
  
  land_uses <-
    c(
      "Primary vegetation",
      "Secondary",
      "Plantation forest",
      "Pasture",
      "Cropland",
      "Urban"
    )
  
  # store realms as vector
  
  trophic_niches <-
    c("Om","In","Ne","Fr","Gr")
  
  # store land use colours
  
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
  
  
  ### significance
  
  contrasts <- trophic_beta_data$Contrasts
  
  
  signif <- purrr::map_dfr(1:nrow(contrasts), function(x){
    
    lu <- stringr::str_split(contrasts[x,"contrast"], pattern = " - ", simplify = TRUE) |>
  stringr::str_remove(pattern = paste(paste0(" ",trophic_niches), collapse = "|"))  
    
    lu <- lu[lu != "Primary vegetation"]
    
    
    niche <- stringr::str_split(contrasts[x,"contrast"], pattern = " - ", simplify = TRUE) |>
      stringr::str_remove(pattern = paste(paste0(land_uses, " "), collapse = "|")) |>
      unique()
   
    p.value <- contrasts[x,"p.value"]
    
  
    out <- data.frame(Predominant_habitat_SSBS_2 = lu,
                      Trophic_Niche = niche,
                      p_value = p.value)
       
    return(out)
    
  })
  
  
  signif$alpha <- ifelse(signif$p_value <= 0.05, 1, 0.8)
  
  #convert land use to factor so order variables on plot
  
  effect_sizes <- effect_sizes |>
    dplyr::left_join(signif)
  effect_sizes[is.na(effect_sizes)] <- 1
  
  effect_sizes$Predominant_habitat_SSBS_2 <- factor(effect_sizes$Predominant_habitat_SSBS_2,
                                                    levels = land_uses)
  
  
  effect_sizes$Trophic_Niche <- factor(effect_sizes$Trophic_Niche,
                                                    levels = trophic_niches)
  
  
  
  
  trophic_similarity_plot <- ggplot2::ggplot(data = effect_sizes, ggplot2::aes(x = Predominant_habitat_SSBS_2,
                                                                               y = rel_similarity, alpha = alpha)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
    ggplot2::geom_point(
      ggplot2::aes(colour = Predominant_habitat_SSBS_2),
      size = 5,
      position = ggplot2::position_dodge(1),
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower,
        ymax = upper,
        colour = Predominant_habitat_SSBS_2
      ),
      width = 0.5,
      linewidth = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = land_use_colours[land_uses, "colours"]) +
    ggplot2::scale_y_continuous(limits = c(0, 1.35),
                       breaks = c(0, 0.5, 1, 1.25)) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::facet_grid(. ~ Trophic_Niche) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black", linetype = "solid"),
      axis.title.y = ggplot2::element_blank(),
      #axis.ticks.y = element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      # axis.ticks.x = element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = 'white', color = 'white'),
      panel.grid.major = ggplot2::element_line(color = 'white'),
      panel.grid.minor = ggplot2::element_line(color = 'white'),
      panel.spacing = ggplot2::unit(2, "lines")
    )
  
    
  
  fig.out <- file.path(file.path(dir_out, "trophic_plots"), "trophic-similarity.png") 
  
  
  
  ggplot2::ggsave(plot = trophic_similarity_plot,
         filename = 
           fig.out,
         width = 350,
         height = 75,
         units = "mm",
         dpi = 500
  )
  

    
    
    
    }