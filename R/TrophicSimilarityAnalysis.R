TrophicSimilarityAnalysis <- function(trophic_beta_metrics,
                                    guilds_in_site,
                                    guilds_in_block,
                                    dir_out){
  
  t_beta_metrics <- qs2::qs_read(trophic_beta_metrics)
  
  t_beta_metrics$Predominant_habitat_SSBS_2 <- factor(t_beta_metrics$Predominant_habitat_SSBS_2, levels = c("Primary vegetation",
                                                                                                            "Secondary",
                                                                                                            "Plantation forest",
                                                                                                            "Pasture",
                                                                                                            "Cropland",
                                                                                                            "Urban"))
  
  guild_site <- qs2::qs_read(guilds_in_site)  
  
  t_beta_metrics <- t_beta_metrics |>
    dplyr::left_join(guild_site, by = c("SSBS_2" = "SSBS", "Trophic_Niche")) |>
    dplyr::rename(count_SSBS_2 = count) |>
    dplyr::left_join(guild_site, by = c("SSBS_1" = "SSBS", "Trophic_Niche")) |>
    dplyr::rename(count_SSBS_1 = count) 
  
  t_beta_metrics <- t_beta_metrics |>
    dplyr::filter(Trophic_Niche %in% c("Fr","Gr","In","Ne","Om"),
                  SSB_SSBS_1 == SSB_SSBS_2) |>
    dplyr::filter(!is.na(count_SSBS_2)|
                    !is.na(count_SSBS_1))
  
  
  guild_block <- qs2::qs_read(guilds_in_block)
  
  t_beta_metrics <- t_beta_metrics |>
    dplyr::left_join(guild_block, by = c("SSB_SSBS_2" = "SSB", "Trophic_Niche"))
  
  
  t_beta_metrics <- t_beta_metrics |>
    dplyr::group_by(Trophic_Niche) |>
    dplyr::mutate(count = count/max(count))
  
  
  t_beta_metrics$logitSimilarity <- car::logit(t_beta_metrics$Similarity, adjust = 0.01)
  
  t_beta_metrics[t_beta_metrics$Similarity == 1, "Similarity"] <- 
    0.99
  

  trophic_overlap <- glmmTMB::glmmTMB(Similarity ~ Predominant_habitat_SSBS_2 * Trophic_Niche + (Trophic_Niche|SS) + (1|SSB_SSBS_2),
                                      data = t_beta_metrics,ziformula = ~1,weights = count, family = glmmTMB::beta_family(link = "logit"),
                                      control = glmmTMB::glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))


  
  # trophic_overlap <- lmerTest::lmer(logitSimilarity ~ Predominant_habitat_SSBS_2 * Trophic_Niche + (Trophic_Niche|SS) + (1|SSB_SSBS_2),
  #                                   data = t_beta_metrics,weights = count)
  
  # , control=lme4::lmerControl(optimizer="bobyqa",
  #                             optCtrl=list(maxfun=100000))
  
  summary(trophic_overlap)
  
  BIC(trophic_overlap) # 4559.911
  
  
  # resid panel
  #plot(performance::simulate_residuals(trophic_overlap))
  # check the significance of the interaction
  car::Anova(trophic_overlap)
  # get the R squared
  performance::r2(trophic_overlap)
  # get the pairwise contrasts
  similarity_contrasts <-
    summary(emmeans::emmeans(trophic_overlap, revpairwise ~ Predominant_habitat_SSBS_2 * Trophic_Niche))
  
  
  emmeans <- similarity_contrasts$emmeans
  
  
  emmeans$similarity <- arm::invlogit(emmeans$emmean)
  emmeans$lower <- arm::invlogit(emmeans$asymp.LCL)
  emmeans$upper <- arm::invlogit(emmeans$asymp.UCL)
  
  
  emmeans <- emmeans |>
    dplyr::group_by(Trophic_Niche) |>
    dplyr::mutate(pri = similarity[Predominant_habitat_SSBS_2 == "Primary vegetation"],
                  rel_similarity = similarity/pri,
                  upper = upper/pri,
                  lower = lower/pri)
  
  
  
  contrasts <- similarity_contrasts$contrasts 
  
  
  # contrast to get pairwise p-values with tukey HSD
  contrasts <- purrr::map_dfr(unique(t_beta_metrics$Trophic_Niche), function(x){
    cons <- contrasts |>
      dplyr::filter(grepl(contrast, pattern = paste("Primary vegetation", x)),
                    stringr::str_count(contrast, pattern = x) == 2)
    return(cons)
  })
  
  
  
  output <- list(Model = trophic_overlap,
                 `Marginal effects` = emmeans,
                 `Contrasts` = contrasts)
  
  
  file.out <- file.path(dir_out,"trophic-beta-analysis.qs")
  
  
  qs2::qs_save(output,file.out)
  
  
  return(file.out)  
  
}
