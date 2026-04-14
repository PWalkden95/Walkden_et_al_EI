SiteAlphaAnalysis <- function(site_level_alpha_metrics,
                                         dir_out){
 
  
  logger::log_info(
    "Read in site level alpha metrics"
  )
  
  alpha_metrics <- qs2::qs_read(site_level_alpha_metrics)
  
  
alpha_metrics$Predominant_habitat <- factor(alpha_metrics$Predominant_habitat, levels = c("Primary vegetation",
                                                                            "Secondary",
                                                                            "Plantation forest",
                                                                            "Pasture",
                                                                            "Cropland",
                                                                            "Urban"))
  

# hist(alpha_metrics$FRich)
# hist(log1p(alpha_metrics$FRich))
# hist(log(alpha_metrics$FRich))
# hist(sqrt(alpha_metrics$FRich))
#   # 
  # hist(alpha_metrics$RaoQ)
  # hist(log(alpha_metrics$RaoQ))
  

logger::log_info(
  "Square root transformation of funcitonal richness"
)

  alpha_metrics$sqrt_frich <- sqrt(alpha_metrics$FRich)
  
  
  logger::log_info(
    "Model functional richness and a function of land use using GLMMs"
  )
  
  mod_rich <- lmerTest::lmer(sqrt_frich ~ Predominant_habitat + (1|SS) + (1|SSB), data = alpha_metrics)

  
 # summary(mod_rich)


#ggResidpanel::resid_panel(mod_rich)

#----------------------------------------------------
# Data for model visualisations through bootstrapping 
#------------------------------------------------------

# the best way to visulaise mixed effects models is to perform bootstraping, 
# additionally this is can assess the significance of the effect variables.

# first will need blank data for all the biogeographic realm - land use 
# combinations with sufficient data 


boot_newdata <- alpha_metrics |>
  dplyr::distinct(Predominant_habitat)


## run the function to bootstrap the three models

 model_bootstraps <-
  .MixedModelBootstrapping(
    model = list(mod_rich),
    newdata = boot_newdata,
    transformation = c("sqrt"),
    nboot = 250,
    n_cores = 1
  )


# save data for visualisations


univariate_output <- list(`Model` = mod_rich,
               `Bootstraps` = model_bootstraps)


## -----------------------------------------------------------------------
# Multivariate modelling looking at the interaction between land use and 
# biogerographic realm
#-------------------------------------------------------------------------



mod_rich_realm <-
  lmerTest::lmer(sqrt_frich ~ Predominant_habitat * Realm + (1|SS) + (1|SSB)  , data = alpha_metrics)

summary(mod_rich_realm)
#ggResidpanel::resid_panel(mod_rich_realm)

car::Anova(mod_rich_realm)

performance::r2(mod_rich_realm)



# get the pairwise contrasts
mod_rich_realm_contrasts <-
  summary(emmeans::emmeans(mod_rich_realm, revpairwise ~ Predominant_habitat * Realm, adjust = "tukey"))



contrasts <- mod_rich_realm_contrasts$contrasts 


# contrast to get pairwise p-values with tukey HSD
contrasts <- purrr::map_dfr(levels(alpha_metrics$Realm), function(x){
  cons <- contrasts |>
    dplyr::filter(grepl(contrast, pattern = paste("Primary vegetation", x)),
                  stringr::str_count(contrast, pattern = x) == 2)
 return(cons)
})

# pairwise estimates
emmeans <- mod_rich_realm_contrasts$emmeans


#----------------------------------------------------
# Data for model visualisations through bootstrapping 
#------------------------------------------------------

# the best way to visulaise mixed effects models is to perform bootstraping, 
# additionally this is can assess the significance of the effect variables.

# first will need blank data for all the biogeographic realm - land use 
# combinations with sufficient data 


boot_newdata <- alpha_metrics |>
  dplyr::distinct(Realm,Predominant_habitat)


## run the function to bootstrap the three models

model_bootstraps_realm <-
  .MixedModelBootstrapping(
    model = list(mod_rich_realm),
    newdata = boot_newdata,
    transformation = c("sqrt"),
    nboot = 250,
    n_cores = 1
  )



multivariate_output <- list(`Model` = mod_rich_realm,
                            `Marginal effects` = emmeans,
                            `Contrasts` = contrasts,
                            `Bootstraps` = model_bootstraps_realm)


file.out <- file.path(dir_out, "site-alpha-analysis.qs")

output <- list(Univariate = univariate_output,
               Mulitvariate = multivariate_output)


qs2::qs_save(output, file.out)

return(file.out)
         
}



