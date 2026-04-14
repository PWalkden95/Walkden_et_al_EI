PREDICTSSiteOverlap <- function(predicts_path,
                                species_maps_dir,
                                   dir_out){
  
  
  require(foreach)
  
  
  
  species_files <-
    list.files(species_maps_dir,
               full.names = TRUE)
  
  
  
  
  predicts <- qs2::qs_read(predicts_path)
  
  
  sites <-
    predicts |>
    dplyr::distinct(SS, SSBS, Longitude, Latitude) |>
    data.frame()
  
  
  site_coordinates <- 
    terra::vect(as.matrix(sites[,c("Longitude","Latitude")]),
                crs = "+proj=longlat +datum=WGS84 +no_defs")
  
  ################### function for generating which species overlap with each of the predicts sites
  ###################
  
  doParallel::registerDoParallel(cores = 10)
    predicts_site_overlap <- foreach(range = species_files,
                              .combine = "cbind",
                              .packages = c(
                                "dplyr",
                                "terra",
                                "sf"),
                              .inorder = FALSE) %dopar% {
                                
                                
                                data <- readRDS(range) # load in species range
                                
                                sp_name <- data$binomial[1]
                                
                                #### filter data so that distirbution polygons are those that reprsent extant, probably extnat and possibly extant presence,
                                #### native, reintroduced and introduced origin and those that are resident, breeding or non-breeding rnages
                                
                                data <-
                                  .FilterPolygons(
                                    BL_data = data,
                                    presence_codes = c(1, 2, 3),
                                    origin_codes = c(1, 2, 3),
                                    seasonal_codes = c(1, 2, 3)
                                  )
                                
                                if (nrow(data) == 0) {
                                  mat <- .BlankMatrix(columns = sp_name, rows = 1:nrow(sites))
                                  
                                  return(mat)
                                } else {
                                  
                                  
                                  shape <- terra::vect(.SpatialCombinePolygons(data))
                                  
                                  
                                  
                                  ### convert polygon into a presence raster on the blank world map
                                  
                                  #
                                  # try(ras <-  raster::rasterize(x = shape, y = blank, getCover = TRUE,background = NA),silent = TRUE)
                                  # ras@data@values <- ifelse(ras@data@values >= 0.1, 1,0)
                                  
                                  does.overlap <- as.numeric(terra::relate(shape, site_coordinates, "contains"))
                                  
                                  #### create matrix where the row is the species and each column is a cell in the raster - presence indicated by a 1
                                  
                                  mat <- matrix(does.overlap, ncol = 1)
                                  colnames(mat) <- sp_name
                                  
                                  return(mat)
                                }
                                
                              }
  
  foreach::registerDoSEQ()
  
  
  
  
  
  predicts_site_overlap <- predicts_site_overlap[, colSums(predicts_site_overlap, na.rm = TRUE)  > 0]
  
  
  rownames(predicts_site_overlap) <- sites$SSBS
  ## save output map data that can then be put given to the blank raster
  
  file.out <- file.path(dir_out, "predicts_site_species_overlap.qs")
  
  qs2::qs_save(predicts_site_overlap, file.out)
  
  return(file.out)
  
  
}

