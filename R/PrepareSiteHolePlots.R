PrepareSiteHolePlots <- function(site_level_analysis, hole_data, dir_out) {
  if (!dir.exists(file.path(dir_out, "site_level_plots"))) {
    dir.create(file.path(dir_out, "site_level_plots"))
  }

  analysis <- qs2::qs_read(site_level_analysis)

  hole_metrics <- qs2::qs_read(hole_data) |>
    dplyr::mutate(
      convex_hull_volume = convex_hull_volume / max(convex_hull_volume),
      mean_hole_volume = ifelse(is.nan(mean_hole_volume), 0, mean_hole_volume),
      normal_mean_hole = mean_hole_volume / convex_hull_volume,
      normal_hole_number = number_of_holes / convex_hull_volume
    )

  fig.out <- file.path(
    dir_out,
    "site_level_plots",
    "site-level-hole-volume.png"
  )

  .SiteLevelPlot(
    raw_data = hole_metrics,
    bootstrap_data = analysis$`hole volume`$Bootstraps,
    metric = "normal_mean_hole",
    limits = c(-0.9, 5),
    filename = fig.out
  )

  fig.out <- file.path(
    dir_out,
    "site_level_plots",
    "site-level-hole-number.png"
  )

  .SiteLevelPlot(
    raw_data = hole_metrics,
    bootstrap_data = analysis$`hole number`$Bootstraps,
    metric = "normal_hole_number",
    limits = c(-0.2, 2.5),
    filename = fig.out
  )
}
