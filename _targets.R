# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

data_path <- config::get()$data_path
project <- config::get()$project


# Set target options:
packages <- c(
  "terra",
  "stringr",
  "tidyr",
  "purrr",
  "glue",
  "data.table",
  "foreach",
  "sp",
  "dplyr",
  "glmmTMB"
) # Packages that your targets need for their tasks.
# format = "qs", # Optionally set the default storage format. qs is fast.
#
# Pipelines that take a long time to run may benefit from
# optional distributed computing. To use this capability
# in tar_make(), supply a {crew} controller
# as discussed at https://books.ropensci.org/targets/crew.html.
# Choose a controller that suits your needs. For example, the following
# sets a controller that scales up to a maximum of two workers
# which run as local R processes. Each worker launches when there is work
# to do and exits if 60 seconds pass with no tasks to run.
#
controller = crew::crew_controller_local(
  workers = parallel::detectCores() - 1,
  seconds_idle = 60
)

targets::tar_source(file.path(data_path, "predicts.tpds/R"))
targets::tar_source("R")
#
# Alternatively, if you want workers to run on a high-performance computing
# cluster, select a controller from the {crew.cluster} package.
# For the cloud, see plugin packages like {crew.aws.batch}.
# The following example is a controller for Sun Grid Engine (SGE).
#
#   controller = crew.cluster::crew_controller_sge(
#     # Number of workers that the pipeline can scale up to:
#     workers = 10,
#     # It is recommended to set an idle time so workers can shut themselves
#     # down if they are not running tasks.
#     seconds_idle = 120,
#     # Many clusters install R as an environment module, and you can load it
#     # with the script_lines argument. To select a specific verison of R,
#     # you may need to include a version string, e.g. "module load R/4.3.2".
#     # Check with your system administrator if you are unsure.
#     script_lines = "module load R"
#   )
#
# Set other options as needed.)

# Run the R scripts in the R/ folder with your custom functions:

# tar_source("other_functions.R") # Source other scripts as needed.

# Set target options:
targets::tar_option_set(
  # packages that your targets need to run
  packages = packages,
  # default storage format
  # quick to check targets
  trust_timestamps = TRUE,
  # remove objects from memory as soon as possible
  memory = "transient",
  garbage_collection = TRUE,
  # use crew to run the pipeline as a multisession
  controller = controller,
  # make sure the data can be retrieved by each worker
  storage = "worker",
  retrieval = "worker",
  error = "continue"
)

list(
  tar_target(
    name = EI_predicts_data,
    command = PreparePredictsEIData(
      predicts_avonet_path = file.path(
        data_path,
        "predicts-avonet-harmonisation/predicts.avonet-out/predicts-avonet-prepared-data.qs"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = predicts_site_overlap,
    command = PREDICTSSiteOverlap(
      predicts_path = EI_predicts_data,
      species_maps_dir = file.path(
        data_path,
        "datasets/Birdlife_maps/PREDICTS_BL/"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = inventory_matrix,
    command = SpeciesInventoryMatrix(
      species_maps_dir = file.path(
        data_path,
        "datasets/Birdlife_maps/PREDICTS_BL/"
      ),
      template_map = file.path(data_path, project, "ei-in/blank_map_1deg.tif"),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = adf_matrix,
    command = FocalSimilarity(
      predicts_path = EI_predicts_data,
      predicts_species_overlap = predicts_site_overlap,
      inventory_matrix = inventory_matrix,
      template_map = file.path(data_path, project, "ei-in/blank_map_1deg.tif"),
      island_polys = file.path(
        data_path,
        project,
        "ei-in/assembly_islands.rds"
      ),
      island_spp = file.path(
        data_path,
        project,
        "ei-in/assembly_island_spp.rds"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = regional_species_pool,
    command = RegionalSpeciesPool(
      predicts_path = EI_predicts_data,
      predicts_site_cells = file.path(
        data_path,
        project,
        "ei-out/predicts-sites-cells.qs"
      ),
      adf_matrix = adf_matrix,
      inventory_matrix = inventory_matrix,
      island_spp = file.path(
        data_path,
        project,
        "ei-in/assembly_island_spp.rds"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = dispersal_probability,
    command = DispersalProbabilities(
      predicts_path = EI_predicts_data,
      species_maps_dir = file.path(
        data_path,
        "datasets/Birdlife_maps/PREDICTS_BL/"
      ),
      species_pools = regional_species_pool,
      avonet_means = file.path(
        data_path,
        "datasets/AVONET/avonet-birdlife-mean-traits.csv"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = environmental_affinity,
    command = EnvironmentalAffinity(
      predicts_path = EI_predicts_data,
      SDM_probabilities_dir = file.path(
        data_path,
        "datasets/Bird_SDM_Projections_Probability/"
      ),
      SDM_miss_sp = file.path(data_path, project, "ei-in/Miss_sp.csv"),
      species_maps_dir = file.path(
        data_path,
        "datasets/Birdlife_maps/PREDICTS_BL/"
      ),
      crosswalk_path = file.path(
        data_path,
        "datasets/AVONET/BL_Jetz crosswalk v3.csv"
      ),
      species_pools_path = regional_species_pool,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_randomisations,
    command = SiteRandomisations(
      predicts_path = EI_predicts_data,
      dispersal_probabilities = dispersal_probability,
      environmental_affinity = environmental_affinity,
      species_pools_path = regional_species_pool,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = ei_traits,
    command = PREDICTSEITraits(
      predicts_path = EI_predicts_data,
      randomisations_path = site_randomisations,
      avonet_specimen_data = file.path(
        data_path,
        "datasets/AVONET/GBD_BiometricsRaw_combined_15_Sept_2021_MASTER.csv"
      ),
      avonet_mean_data = file.path(
        data_path,
        "datasets/AVONET/avonet-birdlife-mean-traits.csv"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),

  tar_target(
    name = trait_probability_densities,
    command = PREDICTSEcosystemIntegrityTPDs(
      predicts_path = EI_predicts_data,
      randomisations_path = site_randomisations,
      species_pools_path = regional_species_pool,
      ei_traits = ei_traits,
      avonet_means_path = file.path(
        data_path,
        "datasets/AVONET/avonet-birdlife-mean-traits.csv"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = realm_tpds,
    command = PREDICTSEIRealmTPDs(
      predicts_path = EI_predicts_data,
      randomisations_path = site_randomisations,
      species_pools_path = regional_species_pool,
      ei_traits = ei_traits,
      avonet_means_path = file.path(
        data_path,
        "datasets/AVONET/avonet-birdlife-mean-traits.csv"
      ),
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_combinations,
    command = SiteCombinations(
      predicts_path = EI_predicts_data,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),

  tar_target(
    name = site_level_alpha,
    command = PrepareSiteLevelAlphaMetrics(
      predicts_path = EI_predicts_data,
      predicts_site_tpd_path = trait_probability_densities[1],
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_level_beta,
    command = PrepareSiteLevelBetaMetrics(
      predicts_combinations_path = site_combinations,
      predicts_site_tpd_path = trait_probability_densities[1],
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),

  tar_target(
    name = site_holes,
    command = PrepareSiteHoles(
      predicts_path = EI_predicts_data,
      predicts_site_tpd_path = trait_probability_densities[1],
      predicts_randomisation_tpd_path = trait_probability_densities[2],
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = trophic_level_alpha,
    command = PrepareTrophicAlphaMetrics(
      predicts_path = EI_predicts_data,
      predicts_site_tpd_path = trait_probability_densities[1],
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = trophic_level_beta,
    command = PrepareTrophicBetaMetrics(
      predicts_combinations_path = site_combinations,
      predicts_site_tpd_path = trait_probability_densities[1],
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_alpha_analysis,
    command = SiteAlphaAnalysis(
      site_level_alpha_metrics = site_level_alpha,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_beta_analysis,
    command = SiteBetaAnalysis(
      site_beta_metrics = site_level_beta,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_hole_analysis,
    command = SiteHolesAnalysis(
      site_level_hole_metrics = site_holes,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = guilds_in_site,
    command = TrophicNicheInSite(
      predicts_path = EI_predicts_data,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = guilds_in_block,
    command = TrophicNicheInBlock(
      predicts_path = EI_predicts_data,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = trophic_alpha_analysis,
    command = TrophicRichnessAnalysis(
      trophic_alpha_metrics = trophic_level_alpha,
      guilds_in_block = guilds_in_block,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = trophic_beta_analysis,
    command = TrophicSimilarityAnalysis(
      trophic_beta_metrics = trophic_level_beta,
      guilds_in_site = guilds_in_site,
      guilds_in_block = guilds_in_block,
      dir_out = file.path(data_path, project, "ei-out")
    )
  ),
  tar_target(
    name = site_alpha_figure,
    command = PrepareSiteAlphaPlots(
      site_level_analysis = site_alpha_analysis,
      alpha_data = site_level_alpha,
      dir_out = file.path(data_path, project, "figures")
    )
  ),
  tar_target(
    name = site_beta_figure,
    command = PrepareSiteBetaPlots(
      site_level_analysis = site_beta_analysis,
      beta_data = site_level_beta,
      dir_out = file.path(data_path, project, "figures")
    )
  ),
  tar_target(
    name = site_hole_figure,
    command = PrepareSiteHolePlots(
      site_level_analysis = site_hole_analysis,
      hole_data = site_holes,
      dir_out = file.path(data_path, project, "figures")
    )
  )
)
