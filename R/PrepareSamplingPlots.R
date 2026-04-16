PrepareSamplingPlots <- function(predicts_path, dir_out) {
  predicts <- qs2::qs_read(predicts_path)

  site_info <-
    predicts |>
    dplyr::distinct(SSBS, SS, Realm, Predominant_habitat, Longitude, Latitude)

  # view distribution of land uses in biogeographic realms
  site_tbl <- table(site_info$Predominant_habitat, site_info$Realm)

  colSums(as.matrix(site_tbl))

  site_tbl

  # store realms for looping
  realms <- as.character(unique(site_info$Realm))

  # store land uses
  land_uses <-
    c(
      "Urban",
      "Cropland",
      "Pasture",
      "Plantation forest",
      "Secondary",
      "Primary vegetation"
    )

  # colours for the land use bars
  land_use_colours <-
    c(
      Urban = "#6d6e71",
      Cropland = "#d1cc1d",
      Pasture = "#fbb040",
      'Plantation forest' = "#8dc63f",
      'Secondary' = "#39b54a",
      'Primary vegetation' = "#225816"
    )

  # land use categorey data

  land_use_plot_data <-
    #group by land use and realm
    site_info |>
    dplyr::group_by(Predominant_habitat, Realm) |>
    # count the number in each category and realm
    dplyr::summarise(count = dplyr::n()) |>
    # pivot table wider to get correct format
    tidyr::pivot_wider(names_from = Realm, values_from = count) |>
    #convert land use to factor to ease plotting
    dplyr::mutate(
      Predominant_habitat = factor(Predominant_habitat, levels = land_uses)
    ) |>
    #convert to data frame
    as.data.frame()

  land_use_plot_data[is.na(land_use_plot_data)] <- 0

  rownames(land_use_plot_data) <- land_use_plot_data$Predominant_habitat
  land_use_plot_data <- land_use_plot_data[land_uses, ]

  # we want figures to be in the same scale so we will cap the limit to the most
  # numerous land use realm combination (Neotropic Primary vegetation)

  limits <- max(
    land_use_plot_data[, c(2:ncol(land_use_plot_data))],
    na.rm = TRUE
  ) +
    20

  if (!dir.exists(file.path(dir_out, "sampling_plots"))) {
    dir.create(file.path(dir_out, "sampling_plots"))
  }

  #----------------------------
  # loop to create plot and save
  #-----------------------------

  # for each realm
  purrr::map(realms, function(x) {
    realm_land_use_data <- land_use_plot_data |>
      dplyr::select(Predominant_habitat, dplyr::all_of(x))

    colnames(realm_land_use_data)[2] <- "realm"
    # create the bar graph indicating the number of sites designated as each land
    # use category
    sample_plot <-
      ggplot2::ggplot(
        data = realm_land_use_data,
        ggplot2::aes(x = Predominant_habitat, y = realm, label = realm)
      ) +
      ggplot2::ylim(0, limits) +
      ggplot2::geom_bar(
        stat = "identity",
        show.legend = FALSE,
        fill = land_use_colours
      ) +
      ggplot2::geom_text(
        colour = "black",
        ggplot2::aes(y = (realm + 2)),
        size = 20,
        fontface = "bold",
        hjust = 0
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()
      )

    fig.out <- file.path(
      dir_out,
      "sampling_plots",
      glue::glue("LU-sample-plot-{x}.png")
    )
    # save the plot

    ggplot2::ggsave(
      sample_plot,
      filename = fig.out,
      width = 250,
      height = 200,
      units = "mm",
      dpi = 500
    )
  })

  #==============================================================================
  # Dietary guild plots
  #
  # Plots to show teh distribution of species belonging to different dietary guild
  # in the study and across biogeographic realms
  #==============================================================================

  # change species of other dietary guilds to other as they are not relevant for the
  # study

  diet_data <- predicts |>
    dplyr::mutate(
      Trophic_Niche = ifelse(
        Trophic_Niche %in%
          c("Aq.p", "Hb.T", "Hb.A", "Vt"),
        "Other",
        Trophic_Niche
      )
    )

  # store guilds and the colours used to denote guilds in figure
  guilds <- c("Om", "In", "Ne", "Fr", "Gr")
  guild_colours <-
    c("#b77ab4", "#164970", "#87d1d2", "#ef4136", "#c9842a")

  #-------------------------------------------------------------------------------
  # get data for table of species recorded
  #
  # data will be for supplementary table 1 to show the number of species recorded
  # in each land use category in each biogeographic realm and worldwide. table
  # will also partion species into their dietary guilds.
  #-------------------------------------------------------------------------------

  # get data for number of species in each dietary guild in each land use with
  # each biogeographic realm
  diet_table_data <- diet_data |>
    # distinct species in each land use and realm
    dplyr::distinct(
      Realm,
      Birdlife_Name,
      Predominant_habitat,
      .keep_all = TRUE
    ) |>
    # group by land use, realm and dietary niche
    dplyr::group_by(Realm, Predominant_habitat, Trophic_Niche) |>
    # how many species of each dietary guild are present in each group
    dplyr::summarise(count = dplyr::n()) |>
    # convert to data frame
    data.frame()

  # get the total number of species in each land use class in each realm
  habitat_totals <- diet_table_data |>
    # group by realm and land use
    dplyr::group_by(Realm, Predominant_habitat) |>
    # sum the number of species
    dplyr::summarise(value = sum(count, na.rm = TRUE)) |>
    # convert to data frame
    data.frame()

  #----------------------------------------------------------------------------
  # How many unique species belonging to each dietary guild are present in each
  # biogeogrpahic realm
  #----------------------------------------------------------------------------

  diet_totals <- diet_data |>
    # distinct species in each realm
    dplyr::distinct(Realm, Birdlife_Name, .keep_all = TRUE) |>
    # group by realm and dietary niche
    dplyr::group_by(Realm, Trophic_Niche) |>
    # sum teh number of species of each dietary guild in each realm
    dplyr::summarise(value = dplyr::n()) |>
    # convert to data frame
    data.frame()

  # total number of species in each realm

  species_totals <- diet_totals |>
    dplyr::group_by(Realm) |>
    dplyr::summarise(count = sum(value, na.rm = TRUE)) |>
    # convert to data frame
    data.frame()

  #----------------------------------------
  # worldiwde species
  #----------------------------------------

  # total number of species belonging to each dietary guild that inhabit each
  #land use category worldwide
  worldwide_total <- diet_data |>
    # distinct species in each land use category
    dplyr::distinct(Birdlife_Name, Predominant_habitat, .keep_all = TRUE) |>
    # group by land use and dietary guild
    dplyr::group_by(Predominant_habitat, Trophic_Niche) |>
    # count the number of species of each dietary guild in each land use
    dplyr::summarise(count = dplyr::n()) |>
    # convert to data frame
    data.frame()

  # total number of species in each land use category
  worldwide_overalls <- worldwide_total |>
    # group by land use
    dplyr::group_by(Predominant_habitat) |>
    # sum teh number of species
    dplyr::summarise(value = sum(count, na.rm = TRUE)) |>
    # convert to data frame
    data.frame()

  # what is the overall number of species taht belong to each dietary category
  overall_total <- diet_data |>
    # distinct bordlife name in the analysis data
    dplyr::distinct(Birdlife_Name, .keep_all = TRUE) |>
    # group by dietary guild
    dplyr::group_by(Trophic_Niche) |>
    # count the number of species in each dietary guild
    dplyr::summarise(count = dplyr::n()) |>
    # convert to data frame
    data.frame()

  total_number_of_species <- sum(overall_total$count, na.rm = TRUE)

  # create limits so all plots are in the same scale

  diet_limits <- max(diet_totals$value) + 100

  #---------------------------------------------
  # loop to create and save dietary guild plots
  #--------------------------------------------

  diet_plot_data <- diet_totals |>
    tidyr::pivot_wider(values_from = value, names_from = Realm) |>
    dplyr::filter(Trophic_Niche %in% guilds) |>
    dplyr::mutate(Trophic_Niche = factor(Trophic_Niche, levels = guilds)) |>
    as.data.frame()

  rownames(diet_plot_data) <- diet_plot_data$Trophic_Niche

  diet_plot_data <- diet_plot_data[guilds, ]

  diet_plot_data[is.na(diet_plot_data)] <- 0

  # for each realm
  purrr::map(realms, function(x) {
    # create the bar graph indicating the number of species classified in each
    # dietary guild
    diet_plot <-
      ggplot2::ggplot(
        data = diet_plot_data[, c("Trophic_Niche", x)],
        ggplot2::aes(
          x = c(0.3, 0.65, 1, 1.35, 1.70),
          y = as.numeric(diet_plot_data[, x]),
          label = as.numeric(diet_plot_data[, x])
        )
      ) +
      ggplot2::geom_bar(
        stat = "identity",
        show.legend = FALSE,
        colour = c("black", "black", "black", "black", "black"),
        fill = guild_colours,
        width = 0.3
      ) +
      ggplot2::geom_text(
        colour = "black",
        ggplot2::aes(y = (as.numeric(diet_plot_data[, x]) + 10)),
        size = 15,
        fontface = "bold",
        angle = 90,
        hjust = 0
      ) +
      #geom_vline(xintercept = 1.95, lwd = 2) +
      # xlim(0, 2) +
      ggplot2::ylim(0, diet_limits) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()
      )

    fig.out <- file.path(
      dir_out,
      "sampling_plots",
      glue::glue("diet-sample-plot-{x}.png")
    )
    # save the plot
    ggplot2::ggsave(
      plot = diet_plot,
      filename = file.path(
        dir_out,
        "sampling_plots",
        glue::glue("diet-sample-plot-{x}.png")
      ),
      width = 100,
      height = 200,
      units = "mm",
      dpi = 500
    )
  })

  worldwide_diet_plot <- overall_total |>
    dplyr::filter(Trophic_Niche %in% guilds)

  worldwide_diet_plot$Trophic_Niche <- factor(
    worldwide_diet_plot$Trophic_Niche,
    levels = c(
      "Om",
      "In",
      "Ne",
      "Fr",
      "Gr"
    )
  )

  rownames(worldwide_diet_plot) <- worldwide_diet_plot$Trophic_Niche
  worldwide_diet_plot <- worldwide_diet_plot[
    c(
      "Om",
      "In",
      "Ne",
      "Fr",
      "Gr"
    ),
  ]

  # create the bar graph indicating the number of species classified in each
  # dietary guild
  plot <-
    ggplot2::ggplot(
      data = worldwide_diet_plot,
      ggplot2::aes(
        x = c(0.3, 0.65, 1, 1.35, 1.70),
        y = as.numeric(worldwide_diet_plot[, 2]),
        label = as.numeric(worldwide_diet_plot[, 2])
      )
    ) +
    ggplot2::geom_bar(
      stat = "identity",
      show.legend = FALSE,
      colour = c("black", "black", "black", "black", "black"),
      fill = guild_colours,
      width = 0.3
    ) +
    ggplot2::geom_text(
      colour = "black",
      ggplot2::aes(y = (as.numeric(worldwide_diet_plot[, 2]) + 20)),
      size = 15,
      fontface = "bold",
      angle = 90,
      hjust = 0
    ) +
    #geom_vline(xintercept = 1.95, lwd = 2) +
    # xlim(0, 2) +
    ggplot2::ylim(0, (max(worldwide_diet_plot$count) + 300)) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )

  fig.out <- file.path(
    dir_out,
    "sampling_plots",
    glue::glue("diet-sample-plot-worldwide.png")
  )
  # save the plot
  ggplot2::ggsave(
    plot = plot,
    filename = fig.out,
    width = 100,
    height = 200,
    units = "mm",
    dpi = 500
  )
}
