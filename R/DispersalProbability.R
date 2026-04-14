DispersalProbabilities <- function(predicts_path,
                                   species_maps_dir,
                                   species_pools,
                                   avonet_means,
                                   dir_out){
  
  
  require(foreach)
  
  predicts <- qs2::qs_read(predicts_path)
  
  predicts_sites <- predicts |>
    dplyr::distinct(SS,SSBS,Latitude,Longitude)
  
 bio_species_pool <-  qs2::qs_read(species_pools)

 all_species <- grep(unique(unlist(bio_species_pool)),
                     pattern = "_",
                     invert = TRUE,
                     value = TRUE)
   
 
 species_maps <- list.files(species_maps_dir, full.names = TRUE)

 avonet_means <- read.csv(avonet_means)

 

 site_coordinates <- terra::vect(as.matrix(predicts_sites[,c("Longitude","Latitude")]),
                                 crs = "+proj=longlat +datum=WGS84 +no_defs", )
 
   
 
 doParallel::registerDoParallel(cores = 10)

distance_to_range_edge <- foreach::foreach(range = species_maps, 
                  .combine = "rbind",
                  .packages = c("terra","dplyr","sf")) %dopar% {
                    
                    data <- readRDS(range) # load in species range
                    
                    sp_name <- data$binomial[1]
                    
                    
                    if(!(sp_name %in% all_species)){
                      mat <- .BlankMatrix(columns = 1:nrow(predicts_sites), rows = sp_name)
                      
                      return(mat)
                    }
                    
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
                      mat <- .BlankMatrix(columns = 1:nrow(predicts_sites), rows = sp_name)
                      
                      return(mat)
                    } else {
                      
                      
                      
                      sp_sites <- which(predicts_sites$SSBS %in% .SiteFind(species = sp_name, species_pools = bio_species_pool))
                      
                      shape <- terra::vect(.SpatialCombinePolygons(data))
                      
                    
                      does.overlap <- as.numeric(terra::relate(shape, site_coordinates[sp_sites], "contains"))
                      
                      distance <-  terra::distance(x = shape, y = site_coordinates[sp_sites])/1000
                      
              
                      distance <- ifelse(does.overlap == 1, 0, distance)
                        
                      mat <- .BlankMatrix(columns = 1:nrow(predicts_sites), rows = sp_name)
                      mat[,sp_sites] <- round(distance, digits = 0)
                     
                     return(mat)
                  }
 
                  }



foreach::registerDoSEQ()


pseudosp <- grep(unique(unlist(bio_species_pool)),
     pattern = "_",
     value = TRUE)


pseudo_mat <- matrix(rep(0,length(pseudosp)*ncol(distance_to_range_edge)),
       ncol = ncol(distance_to_range_edge), dimnames = list(pseudosp, 1:ncol(distance_to_range_edge)))

distance_to_range_edge <- rbind(distance_to_range_edge,pseudo_mat)






species_pools_dispersal <- purrr::map(1:length(bio_species_pool), .f = function(x){
  
sp_pool <- bio_species_pool[[x]]  
  



sp_dists <- distance_to_range_edge[sp_pool,x]
scale_dist <- sp_dists/max(sp_dists, na.rm = TRUE)

scale_HWI <- avonet_means |>
  dplyr::select(Birdlife_Name = Species1,
                Hand.Wing.Index) |>
  dplyr::filter(Birdlife_Name %in% sp_pool) |>
  dplyr::mutate(scale.HWI = Hand.Wing.Index/max(Hand.Wing.Index, na.rm = TRUE)) |>
  dplyr::pull(scale.HWI)


dispersal_probs <- 1 - pgamma(q = scale_dist, shape = scale_HWI)

return(dispersal_probs)

})


names(species_pools_dispersal) <- predicts_sites$SSBS


file.out <- file.path(dir_out, "species-pool-dispersal-probabilities.qs")


qs2::qs_save(file.out, object = species_pools_dispersal)


return(file.out)

}