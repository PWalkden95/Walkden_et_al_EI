SpeciesInventoryMatrix <- function(species_maps_dir,
                                template_map,
                                dir_out){
  
  
require(foreach)
  
 blank <- terra::rast(template_map)
   

 species_files <-
   list.files(species_maps_dir,
              full.names = TRUE)
  
 
 
 
 ################### function for generating which species overlap with each of the predicts sites
 ###################
 
doParallel::registerDoParallel(cores = 10)
 
    inventory_matrix <- foreach(range = species_files,
                                  .combine = "rbind",
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
       mat <- .BlankMatrix(columns = 1:terra::ncell(blank), rows = sp_name)
       
       return(mat)
     } else {
       
       
       shape <- terra::vect(.SpatialCombinePolygons(data))
       
       
       
       ### convert polygon into a presence raster on the blank world map
       
       #
       # try(ras <-  raster::rasterize(x = shape, y = blank, getCover = TRUE,background = NA),silent = TRUE)
       # ras@data@values <- ifelse(ras@data@values >= 0.1, 1,0)
       
       ras <- terra::rasterize(x = shape,y = blank, cover = TRUE) |>
         terra::classify(cbind(0.1,1)) |>
         terra::classify(rcl = matrix(c(0, 0.1, 0,
                                        0.1, 1, 1), ncol = 3, byrow = TRUE)) 
       
       #### create matrix where the row is the species and each column is a cell in the raster - presence indicated by a 1
       
       ras <- terra::values(ras)
        ras[is.na(ras)] <- 0
       
       mat <-
         matrix(
           ras,
           nrow = 1,
           ncol = length(ras),
           dimnames = list(sp_name, c(1:length(ras)))
         )
       
       return(mat)
     }
     
   }
 
    foreach::registerDoSEQ()
 
 
 
  inventory_matrix <- inventory_matrix[, colSums(inventory_matrix, na.rm = TRUE)  > 0]
    
    
 ## save output map data that can then be put given to the blank raster
 
 file.out <- file.path(dir_out, "species-inventory-matrix.qs")
 
 qs2::qs_save(inventory_matrix, file.out)
 
 return(file.out)
 
 
}

