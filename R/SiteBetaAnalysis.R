SiteBetaAnalysis <- function(site_beta_metrics, dir_out) {
  beta_metrics <- qs2::qs_read(site_beta_metrics)
  
  beta_metrics$Predominant_habitat_SSBS_2 <-
    factor(
      beta_metrics$Predominant_habitat_SSBS_2,
      levels = c(
        "Primary vegetation",
        "Secondary",
        "Plantation forest",
        "Pasture",
        "Cropland",
        "Urban"
      )
    )
  
  
  hist(beta_metrics$similarity)
  hist(car::logit(beta_metrics$similarity, adjust = 0.01))
  
  beta_metrics$logitSim <- car::logit(beta_metrics$similarity, adjust = 0.01)
  
  
  beta_metrics$similarity[beta_metrics$similarity > 0.99] <- 0.99
  beta_metrics$similarity[beta_metrics$similarity < 0.01] <- 0.01
  
  
  # Multivare GLMMs including realm and it's interaction with land use
  # dissim_model <- glmmTMB::glmmTMB(similarity ~ Predominant_habitat_SSBS_2 + (1|SS) + (1|SSB_SSBS_2),ziformula = ~1,
  #                                  data = beta_metrics, family = glmmTMB::beta_family(link = "logit"))
  #
  similarity_model <- lmerTest::lmer(
    logitSim ~ Predominant_habitat_SSBS_2 + (1 | SS) + (1 | SSB_SSBS_2),
    data = beta_metrics,
    control = lme4::lmerControl(optimizer = "bobyqa")
  )
  
  summary(similarity_model)
  performance::r2(similarity_model)
  

  #----------------------------------------------------
  # Data for model visualisations through bootstrapping
  #------------------------------------------------------
  
  # the best way to visulaise mixed effects models is to perform bootstraping,
  # additionally this is can assess the significance of the effect variables.
  
  # first will need blank data for all the biogeographic realm - land use
  # combinations with sufficient data
  
  
  boot_newdata <- beta_metrics |>
    dplyr::distinct(Predominant_habitat_SSBS_2)
  
  
  ## run the function to bootstrap the three models
  
  similarity_model_bootstraps <-
    .MixedModelBootstrapping(
      model = list(similarity_model),
      newdata = boot_newdata,
      transformation = c("logit"),
      nboot = 250,
      n_cores = 1
    )
  
  
  similarity_univariate_output <- list(`Model` = similarity_model,
                                       `Bootstraps` = similarity_model_bootstraps)
  
  
  ## -----------------------------------------------------------------------
  # Multivariate modelling looking at the interaction between land use and
  # biogerographic realm
  #-------------------------------------------------------------------------
  
  
  
  similarity_model_realm <-
    lmerTest::lmer(
      logitSim ~ Predominant_habitat_SSBS_2 * Realm + (1 |
                                                         SS) + (1 | SSB_SSBS_2),
      data = beta_metrics,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )
  
  summary(similarity_model_realm)
  #ggResidpanel::resid_panel(mod_rich_realm)
  
  car::Anova(similarity_model_realm)
  
  performance::r2(similarity_model_realm)
  
  
  
  # get the pairwise contrasts
  similarity_model_realm_contrasts <-
    summary(
      emmeans::emmeans(
        similarity_model_realm,
        revpairwise ~ Predominant_habitat_SSBS_2 * Realm,
        adjust = "tukey"
      )
    )
  
  
  
  similarity_contrasts <- similarity_model_realm_contrasts$contrasts
  
  
  # contrast to get pairwise p-values with tukey HSD
  similarity_contrasts <- purrr::map_dfr(levels(beta_metrics$Realm), function(x) {
    cons <- similarity_contrasts |>
      dplyr::filter(grepl(contrast, pattern = paste("Primary vegetation", x)),
                    stringr::str_count(contrast, pattern = x) == 2)
    return(cons)
  })
  
  # pairwise estimates
  similarity_emmeans <- similarity_model_realm_contrasts$emmeans
  
  
  #----------------------------------------------------
  # Data for model visualisations through bootstrapping
  #------------------------------------------------------
  
  # the best way to visulaise mixed effects models is to perform bootstraping,
  # additionally this is can assess the significance of the effect variables.
  
  # first will need blank data for all the biogeographic realm - land use
  # combinations with sufficient data
  
  
  boot_newdata <- beta_metrics |>
    dplyr::distinct(Realm, Predominant_habitat_SSBS_2)
  
  
  
  ## run the function to bootstrap the three models
  
  similarity_model_bootstraps_realm <-
    .MixedModelBootstrapping(
      model = list(similarity_model_realm),
      newdata = boot_newdata,
      transformation = c("logit"),
      nboot = 250,
      n_cores = 1
    )
  
  
  
  similarity_multivariate_output <- list(
    `Model` = similarity_model_realm,
    `Marginal effects` = similarity_emmeans,
    `Contrasts` = similarity_contrasts,
    `Bootstraps` = similarity_model_bootstraps_realm
  )
  
  
  
  #======================
  #-----------
  ## centroid model
  #-----------
  #=======================
  
  hist(beta_metrics$centroid_shift)
  hist(sqrt(beta_metrics$centroid_shift))
  
  
  centroid_model <- 
    lmerTest::lmer(sqrt(centroid_shift) ~ Predominant_habitat_SSBS_2 + 
                     (1|SS) + (1|SSB_SSBS_2), data = beta_metrics)
  
  
  summary(centroid_model)
  performance::r2(centroid_model)
  
  
  
  #----------------------------------------------------
  # Data for model visualisations through bootstrapping
  #------------------------------------------------------
  
  # the best way to visulaise mixed effects models is to perform bootstraping,
  # additionally this is can assess the significance of the effect variables.
  
  ## run the function to bootstrap the three models
  
  centroid_model_bootstraps <-
    .MixedModelBootstrapping(
      model = list(centroid_model),
      newdata = boot_newdata,
      transformation = c("sqrt"),
      nboot = 250,
      n_cores = 1
    )
  
  
  
  centroid_univariate_output <- list(`Model` = centroid_model,
                                       `Bootstraps` = centroid_model_bootstraps)
  
  ## -----------------------------------------------------------------------
  # Multivariate modelling looking at the interaction between land use and
  # biogerographic realm
  #-------------------------------------------------------------------------
  
  
  
  centroid_model_realm <-
    lmerTest::lmer(
      sqrt(centroid_shift) ~ Predominant_habitat_SSBS_2 * Realm + (1 |
                                                         SS) + (1 | SSB_SSBS_2),
      data = beta_metrics,
      control = lme4::lmerControl(optimizer = "bobyqa")
    )
  
  summary(centroid_model_realm)
  #ggResidpanel::resid_panel(mod_rich_realm)
  
  car::Anova(centroid_model_realm)
  
  performance::r2(centroid_model_realm)
  
  
  
  # get the pairwise contrasts
  centroid_model_realm_contrasts <-
    summary(
      emmeans::emmeans(
        centroid_model_realm,
        revpairwise ~ Predominant_habitat_SSBS_2 * Realm,
        adjust = "tukey"
      )
    )
  
  
  
  centroid_contrasts <- centroid_model_realm_contrasts$contrasts
  
  
  # contrast to get pairwise p-values with tukey HSD
  centroid_contrasts <- purrr::map_dfr(levels(beta_metrics$Realm), function(x) {
    cons <- centroid_contrasts |>
      dplyr::filter(grepl(contrast, pattern = paste("Primary vegetation", x)),
                    stringr::str_count(contrast, pattern = x) == 2)
    return(cons)
  })
  
  # pairwise estimates
  centroid_emmeans <- centroid_model_realm_contrasts$emmeans
  
  
  #----------------------------------------------------
  # Data for model visualisations through bootstrapping
  #------------------------------------------------------
  
  # the best way to visulaise mixed effects models is to perform bootstraping,
  # additionally this is can assess the significance of the effect variables.
  
  # first will need blank data for all the biogeographic realm - land use
  # combinations with sufficient data
  
  
  boot_newdata <- beta_metrics |>
    dplyr::distinct(Realm, Predominant_habitat_SSBS_2)
  
  
  
  ## run the function to bootstrap the three models
  
  centroid_model_bootstraps_realm <-
    .MixedModelBootstrapping(
      model = list(centroid_model_realm),
      newdata = boot_newdata,
      transformation = c("sqrt"),
      nboot = 250,
      n_cores = 1
    )
  
  
  
  centroid_multivariate_output <- list(
    `Model` = centroid_model_realm,
    `Marginal effects` = centroid_emmeans,
    `Contrasts` = centroid_contrasts,
    `Bootstraps` = centroid_model_bootstraps_realm
  )
  
  
  
  
  output <- list(
    `similarity Univariate` = similarity_univariate_output,
    `similarity Multivariate` = similarity_multivariate_output,
    `centroid Univariate` = centroid_univariate_output,
    `centroid Multivariate` = centroid_multivariate_output
  )
  
  
  file.out <- file.path(dir_out, "site-beta-analysis.qs")
  
  
  qs2::qs_save(output, file.out)
  
  
  return(file.out)
  
  
  
}
