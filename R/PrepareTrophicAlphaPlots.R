PrepareTrophicAlphaPlots <- function(trophic_alpha_analysis, 
                                    dir_out){
  
  
  if(!dir.exists(file.path(dir_out, "trophic_plots"))){
    dir.create(file.path(dir_out, "trophic_plots"))
  }
  
  trophic_alpha_data <- qs::qread(trophic_alpha_analysis)
  
  effect_sizes <- trophic_alpha_data$`Marginal effects`
  
  trophic_niches <-
    c("Om","In","Ne","Fr","Gr")
  
  niche_colours <-
    data.frame(
      land_use = trophic_niches,
      colours = c(
        "#b77ab4",
        "#164970",
        "#87d1d2",
        "#ef4136",
        "#c9842a"
      ))
  rownames(niche_colours) <- trophic_niches
  
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
  
  contrasts <- trophic_alpha_data$Contrasts
  
  
  signif <- purrr::map_dfr(1:nrow(contrasts), function(x){
    
    lu <- stringr::str_split(contrasts[x,"contrast"], pattern = " - ", simplify = TRUE) |>
      stringr::str_remove(pattern = paste(paste0(" ",trophic_niches), collapse = "|"))  
    
    lu <- lu[lu != "Primary vegetation"]
    
    
    niche <- stringr::str_split(contrasts[x,"contrast"], pattern = " - ", simplify = TRUE) |>
      stringr::str_remove(pattern = paste(paste0(land_uses, " "), collapse = "|")) |>
      unique()
    
    p.value <- contrasts[x,"p.value"]
    
    
    out <- data.frame(Predominant_habitat = lu,
                      Trophic_Niche = niche,
                      p_value = p.value)
    
    return(out)
    
  })
  
  
  signif$alpha <- ifelse(signif$p_value <= 0.05, 1, 0.8)
  
  
  # Absolute plot 
  
  effect_sizes <- effect_sizes |>
    dplyr::left_join(signif)
  effect_sizes[is.na(effect_sizes)] <- 1
  
  effect_sizes$Predominant_habitat <- factor(effect_sizes$Predominant_habitat,
                                             levels = land_uses)
  
  
  effect_sizes$Trophic_Niche <- factor(effect_sizes$Trophic_Niche,
                                       levels = trophic_niches)
  
  
  site_rao_plot <- ggplot2::ggplot(data = effect_sizes,
                                   ggplot2::aes(x = Predominant_habitat,
                                                y = RaoQ, group = Trophic_Niche,
                                                alpha = alpha)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(colour = Trophic_Niche), show.legend = FALSE, linewidth = 0.5, position = ggplot2::position_dodge(0.5)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper, colour = Trophic_Niche, alpha = alpha),width = 0.4,size =1, show.legend = FALSE, position = ggplot2::position_dodge(0.5)) +
    ggplot2::geom_point(ggplot2::aes(fill = Trophic_Niche, alpha = alpha),colour = "black",pch = 21, size = 4,  show.legend = FALSE, position = ggplot2::position_dodge(0.5)) +
    ggplot2::scale_colour_manual(values = niche_colours[trophic_niches, "colours"]) +
    ggplot2::scale_fill_manual(values = niche_colours[trophic_niches, "colours"]) +
    ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,-0.25, -0.5)) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    #facet_grid(.~Trophic_niche) +
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
      panel.grid.minor = ggplot2::element_line(color = 'white')
    )
  
  site_rao_plot
  
    fig.out <- file.path(file.path(dir_out, "trophic_plots"), "trophic-alpha-plot-absolute.png")
  
  
  ggplot2::ggsave(plot = site_rao_plot,
         filename = 
           fig.out,
         width = 120,
         height = 30,
         units = "mm",
         dpi = 500
  )
  
  
  # relative plot
  
  effect_sizes$RaoQ <- effect_sizes$RaoQ - effect_sizes$pri
  effect_sizes$upper <- effect_sizes$upper - effect_sizes$pri
  effect_sizes$lower <- effect_sizes$lower - effect_sizes$pri


  site_rao_plot <- ggplot2::ggplot(data = effect_sizes,
                                   ggplot2::aes(x = Predominant_habitat,
                                                y = RaoQ, alpha = alpha)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_point(
      ggplot2::aes(colour = Predominant_habitat),
      size = 5,
      position = ggplot2::position_dodge(1),
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower,
        ymax = upper,
        colour = Predominant_habitat
      ),
      width = 0.5,
      linewidth = 1,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = land_use_colours[land_uses, "colours"]) +
    ggplot2::scale_y_continuous(breaks = c(0, -0.5,-0.25, 0.25, 0.5)) +
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
  
  
  site_rao_plot
  
    fig.out <- file.path(file.path(dir_out, "trophic_plots"), "trophic-alpha-plot-relative.png")
  
  
  
    ggplot2::ggsave(plot = site_rao_plot,
                    filename = 
                      fig.out,
                    width = 350,
                    height = 75,
                    units = "mm",
                    dpi = 500
    )
    

  
  
    
}


