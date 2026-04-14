TrophicRichnessAnalysis <- function(trophic_alpha_metrics,
                                    guilds_in_block,
                                    dir_out){
  
  t_alpha_metrics <- qs2::qs_read(trophic_alpha_metrics)
  
  t_alpha_metrics$Predominant_habitat <- factor(t_alpha_metrics$Predominant_habitat, levels = c("Primary vegetation",
                                                                                                "Secondary",
                                                                                                "Plantation forest",
                                                                                                "Pasture",
                                                                                                "Cropland",
                                                                                                "Urban"))
  
  guild_block <- qs2::qs_read(guilds_in_block)  
  
  t_alpha_metrics <- t_alpha_metrics |>
    dplyr::left_join(guild_block) |>
    dplyr::filter(!is.na(count))
  
  t_alpha_metrics <- t_alpha_metrics |>
    dplyr::filter(Trophic_Niche %in% c("Fr","Gr","In","Ne","Om"))
  
  
  t_alpha_metrics <- t_alpha_metrics |>
    dplyr::group_by(SS,SSB) |>
    dplyr::mutate(RaoQ = RaoQ/max(RaoQ, na.rm = TRUE))
  
  
  
  
  hist(t_alpha_metrics$FRich)
  hist(log1p(t_alpha_metrics$FRich))
  hist(sqrt(t_alpha_metrics$FRich))
  
  hist(t_alpha_metrics$RaoQ)
  
  
  mod_rich <- glmmTMB::glmmTMB(RaoQ ~ Predominant_habitat * Trophic_Niche  + (1|SS) +(1|SSB),
                               data = t_alpha_metrics,
                               ziformula = ~1,
                               # weights = count,
                               family = glmmTMB::ziGamma(link = "log"),
                               control = glmmTMB::glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
  
  
  
  
  summary(mod_rich)
  
  BIC(mod_rich) # 4559.911
  
  
  # resid panel
  plot(performance::simulate_residuals(mod_rich))
  # check the significance of the interaction
  car::Anova(mod_rich)
  # get the R squared
  performance::r2(mod_rich)
  # get the pairwise contrasts
  rao_contrasts <-
    summary(emmeans::emmeans(mod_rich, revpairwise ~ Predominant_habitat * Trophic_Niche))
  
  
  emmeans <- rao_contrasts$emmeans
  
  
  emmeans$RaoQ <- exp(emmeans$emmean)
  emmeans$lower <- exp(emmeans$asymp.LCL)
  emmeans$upper <- exp(emmeans$asymp.UCL)
  
  
  emmeans <- emmeans |>
    dplyr::group_by(Trophic_Niche) |>
    dplyr::mutate(pri = RaoQ[Predominant_habitat == "Primary vegetation"],
                  rel_RaoQ = RaoQ / pri,
                  rel_upper = upper - pri,
                  rel_lower = lower - pri)
  
  
  
  contrasts <- rao_contrasts$contrasts 
  
  
  # contrast to get pairwise p-values with tukey HSD
  contrasts <- purrr::map_dfr(unique(t_alpha_metrics$Trophic_Niche), function(x){
    cons <- contrasts |>
      dplyr::filter(grepl(contrast, pattern = paste("Primary vegetation", x)),
                    stringr::str_count(contrast, pattern = x) == 2)
    return(cons)
  })
  
  
  
  output <- list(Model = mod_rich,
                 `Marginal effects` = emmeans,
                 `Contrasts` = contrasts)
  
  
  file.out <- file.path(dir_out,"trophic-alpha-analysis.qs")
  

  qs2::qs_save(output,file.out)
  
  
  return(file.out)  
}

