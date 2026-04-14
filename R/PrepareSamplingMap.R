PrepareSamplingMap <- function(predicts_path,
                               realm_path,
                               dir_out){
  
  
  
  
  #================================
  # load in relevant data
  #================================
  
  
  # PREDICTS data
  predicts <- qs::qread(predicts_path)
  
  # extract info pertaining to sites geographic location
  site_locations <-
    predicts  |>
    dplyr::distinct(SSBS, SS, Longitude, Latitude)
  
  # load in polygons that will be able to visualise the delimination of 
  # biogeographic realms
  biome_polygons <-
    sf::st_read(realm_path)
  
  
  # filter just for polygons in the relevant biogeographic realms
  biome_polygons <-
    biome_polygons |>
    dplyr::filter(REALM %in% c("NT", "PA", "NA", "IM", "AT", "AA"))

  # load in map of the world
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  # combine all polygons
  world_polygon <- sf::st_as_sf(sf::st_combine(world$geometry)) |>
    ggplot2::fortify()
  
  #=========================
  # make sampling map
  #=========================
  
  # store the colours for each biogeographic realm in the map
  
  realm_fills <-
    c("white", "grey35", "grey50", "grey65", "grey80", "grey95")
  
  
  
  # create base map plot
  raster_plot <- ggplot2::ggplot() +
    ggplot2::coord_fixed() +
    ggplot2::geom_sf(data = world_polygon,
            fill = "grey40",
            linewidth = NA)
  
  
  
  # turn off spherical geometry which sometimes causes problems
  sf::sf_use_s2(FALSE)
  i <- 1
  
  
  # iteratively add teh biogeographic realm polygons to the map with the assign colours
  
  for (realm in unique(biome_polygons$REALM)) {
    
    # filter for the realm 
    realm_poly <- biome_polygons |>
      dplyr::filter(REALM == realm)
    # join the polygon together  
    realm_poly <- sf::st_as_sf(sf::st_union(realm_poly$geometry))
    
    # add polygon to the map with the assign colour
    raster_plot <- raster_plot  +
      ggplot2::geom_sf(data = realm_poly,
              colour = "black",
              fill = realm_fills[i])
    
    
    i <- i + 1
  }
  
      ## create a dataframe that combines the site points for each study and assigns 
  ## a size depending on the number of sites each study contains
  
  site_points <-
    # group site info dataframe by study 
    site_locations |>
    dplyr::group_by(SS) |> 
    # summarise each study to get mean long lat and the number of encompassing sites
    dplyr::reframe(
      latitude = mean(Latitude),
      longitude = mean(Longitude),
      assemblages = dplyr::n()
    ) |>
    # create column to make a size variable for the point on the map
    dplyr::mutate(
      point_size = ifelse(assemblages > 0, 10 , 0),
      point_size = ifelse(assemblages > 25, 15 , point_size),
      point_size = ifelse(assemblages > 100, 30, point_size)
    )
  
  
  # add points to map and change theme
  
  raster_plot <- raster_plot +
    ggplot2::geom_point(
      data = site_points,
      ggplot2::aes(x = longitude, y = latitude, size = point_size),
      alpha = 0.7,
      fill = "#ffee80",
      colour = "black",
      pch = 21,
      show.legend = FALSE
    ) +
    ggplot2::scale_size_continuous(range = c(3, 12)) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    ggplot2::xlim(-155, 171) +
    ggplot2::ylim(-50, 90)
  
  
  plot(raster_plot)
  #save
  
  fig.out <- file.path(dir_out,"sampling_plots","sampling_map.png")
  
  ggplot2::ggsave(
    filename = fig.out,
    plot = raster_plot,
    device = "png",
    dpi = 1200,
    width = 8.73,
    height = 5.46
  )
  
  
}