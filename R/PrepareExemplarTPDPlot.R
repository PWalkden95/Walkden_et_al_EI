PrepareExemplarTPDPlot <- function(predicts_site_tpd_path,
                                   predicts_randomisation_tpd_path,
                                   dir_out) {
  
  if (!dir.exists(file.path(dir_out, "example_tpd"))) {
    dir.create(file.path(dir_out, "example_tpd"))
  }
  
  predicts_tpds <- qs::qread(predicts_site_tpd_path)
  
  predicts_random_tpds <- qs::qread(predicts_randomisations_tpd_path)
  
  #==========================
  # study plot limits so porjecting into same trait space
  #=======================================
  
  random_tpd <- predicts_random_tpds$`2020_Lees_RAS`[, c("locomotion",
                                                         "foraging",
                                                         "body_size",
                                                         "2020_Lees_RAS RAS_357 4")]
  
  xyz <- random_tpd[, 1:3]
  
  limits <- list(
    xmin = min(xyz[, 1]) - 0.25,
    xmax = max(xyz[, 1]) + 0.25,
    ymin = min(xyz[, 2]) - 0.25,
    ymax = max(xyz[, 2]) + 0.25,
    zmin = min(xyz[, 3]) - 0.25,
    zmax = max(xyz[, 3]) + 0.25
  )
  
  random_tpd <- random_tpd[random_tpd$`2020_Lees_RAS RAS_357 4` > 0, ]
  
  #======================================
  
  prob <- random_tpd$`2020_Lees_RAS RAS_357 4`
  
  prob <- ifelse(prob > quantile(prob, 0.99), quantile(prob, 0.99), prob)
  
  
  x <- xyz[prob > 0, 1] # locomotion
  y <- xyz[prob > 0, 2] # foraging
  z <- xyz[prob > 0, 3] # boy size
  
  
  size_values <- scale(prob)[, 1]
  
  size_values <- (size_values + (abs(min(size_values))) + 3) / 7
  

  rgl::clear3d()
  rgl::view3d(theta = 33,
              phi = 10,
              zoom = 0.7)
  
  rgl::plot3d(
    x = y,
    y = z,
    z = x,
    box = FALSE,
    xlab = "",
    ylab = "",
    zlab = "",
    type = "s",
    size = size_values,
    #  radius = scale,
    alpha = 1,
    zlim = c(limits[["xmin"]], limits[["xmax"]]),
    xlim = c(limits[["ymin"]], limits[["ymax"]]),
    ylim = c(limits[["zmin"]], limits[["zmax"]]),
    col = "#859D80",
    lwd = 0.1,
    axes = FALSE,
    ann = FALSE
  )
  bb <- rgl::par3d("bbox")
  x <- bb[1:2]
  y <- bb[3:4]
  z <- bb[5:6]
  
  # Draw all 12 box edges manually
  rgl::segments3d(rbind(
    c(x[1], y[1], z[1]),
    c(x[2], y[1], z[1]),
    c(x[1], y[2], z[1]),
    c(x[2], y[2], z[1]),
    c(x[1], y[1], z[2]),
    c(x[2], y[1], z[2]),
    # c(x[1], y[2], z[2]), c(x[2], y[2], z[2]),
    c(x[1], y[1], z[1]),
    c(x[1], y[2], z[1]),
    c(x[2], y[1], z[1]),
    c(x[2], y[2], z[1]),
    c(x[1], y[1], z[2]),
    c(x[1], y[2], z[2]),
    #c(x[2], y[1], z[2]), c(x[2], y[2], z[2]),
    c(x[1], y[1], z[1]),
    c(x[1], y[1], z[2]),
    c(x[2], y[1], z[1]),
    c(x[2], y[1], z[2]),
    c(x[1], y[2], z[1]),
    c(x[1], y[2], z[2])
    #c(x[2], y[2], z[1]), c(x[2], y[2], z[2])
  ),
  col = "black")
  
  
  
  fig.out <- file.path(dir_out, "example_tpd", "null_site_tpd.png")
  
  widget <- rgl::rglwidget(width = 1000, height = 1000)
  htmlwidgets::saveWidget(widget, "temp.html", selfcontained = TRUE)
  
  # Save high-res PNG (specify zoom)
  webshot2::webshot(
    "temp.html",
    fig.out,
    vwidth = 1000,
    vheight = 1000,
    zoom = 4
  )
  file.remove("temp.html")
  
  #===========================
  #==========================
  
  site_tpd <-
    predicts_tpds$`2020_Lees_RAS`[, c("locomotion",
                                      "foraging",
                                      "body_size",
                                      "2020_Lees_RAS RAS_357 4")]
  
  site_tpd <- site_tpd[site_tpd$`2020_Lees_RAS RAS_357 4` > 0, ]
  
  
  combine_tpd <- random_tpd |>
    dplyr::rename(random = `2020_Lees_RAS RAS_357 4`) |>
    dplyr::left_join(site_tpd) |>
    dplyr::rename(obvs = `2020_Lees_RAS RAS_357 4`)
  
  
  convex_hull <- t(geometry::convhulln(site_tpd[, 1:3]))
  
  
  combine_tpd <-
    combine_tpd[geometry::inhulln(ch = convex_hull, p = as.matrix(combine_tpd[, 1:3])), ]
  
  
  combine_tpd <- combine_tpd |>
    dplyr::mutate(
      random = ifelse(!is.na(obvs), NA, random),
      colour = ifelse(is.na(random), "#859D80", "grey")
    )
  
  
  
  
  prob <- as.numeric(t(as.matrix(combine_tpd[, c("random", "obvs")]))) |>
    na.omit()
  
  prob <- ifelse(prob > quantile(prob, 0.99), quantile(prob, 0.99), prob)
  
  
  
  x <- combine_tpd[, 1] # locomotion
  y <- combine_tpd[, 2] # foraging
  z <- combine_tpd[, 3] # boy size
  
  
  size_values <- scale(prob)[, 1]
  
  size_values <- (size_values + (abs(min(size_values))) + 3) / 7
  
  
  
  rgl::clear3d()
  rgl::view3d(theta = 33,
              phi = 10,
              zoom = 0.7)
  
  rgl::plot3d(
    x = y,
    y = z,
    z = x,
    box = FALSE,
    xlab = "",
    ylab = "",
    zlab = "",
    type = "s",
    size = size_values,
    #  radius = scale,
    alpha = 1,
    zlim = c(limits[["xmin"]], limits[["xmax"]]),
    xlim = c(limits[["ymin"]], limits[["ymax"]]),
    ylim = c(limits[["zmin"]], limits[["zmax"]]),
    col = combine_tpd$colour,
    lwd = 0.1,
    axes = FALSE,
    ann = FALSE
  )
  
  rgl::triangles3d(site_tpd[, 1:3][convex_hull, 2],
                   site_tpd[, 1:3][convex_hull, 3],
                   site_tpd[, 1:3][convex_hull, 1],
                   col = "transparent",
                   alpha = 0.3)
  
  bb <- rgl::par3d("bbox")
  x <- bb[1:2]
  y <- bb[3:4]
  z <- bb[5:6]
  
  # Draw all 12 box edges manually
  rgl::segments3d(rbind(
    c(x[1], y[1], z[1]),
    c(x[2], y[1], z[1]),
    c(x[1], y[2], z[1]),
    c(x[2], y[2], z[1]),
    c(x[1], y[1], z[2]),
    c(x[2], y[1], z[2]),
    # c(x[1], y[2], z[2]), c(x[2], y[2], z[2]),
    c(x[1], y[1], z[1]),
    c(x[1], y[2], z[1]),
    c(x[2], y[1], z[1]),
    c(x[2], y[2], z[1]),
    c(x[1], y[1], z[2]),
    c(x[1], y[2], z[2]),
    #c(x[2], y[1], z[2]), c(x[2], y[2], z[2]),
    c(x[1], y[1], z[1]),
    c(x[1], y[1], z[2]),
    c(x[2], y[1], z[1]),
    c(x[2], y[1], z[2]),
    c(x[1], y[2], z[1]),
    c(x[1], y[2], z[2])
    #c(x[2], y[2], z[1]), c(x[2], y[2], z[2])
  ),
  col = "black")
  #rgl::axes3d(edges = c("x--","x-+","x++","y++","y+-","y-+","z--","z++","z+-"), tick = FALSE, labels = FALSE)
  
  
  
  fig.out <- file.path(dir_out, "example_tpd", "metric_site_tpd.png")
  

    widget <- rgl::rglwidget(width = 1000, height = 1000)
    htmlwidgets::saveWidget(widget, "temp.html", selfcontained = TRUE)
    
    # Save high-res PNG (specify zoom)
    webshot2::webshot(
      "temp.html",
      fig.out,
      vwidth = 1000,
      vheight = 1000,
      zoom = 4
    )
    file.remove("temp.html")
  
  
  
}
