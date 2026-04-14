EnvironmentalAffinity <- function(predicts_path,
                                   SDM_probabilities_dir,
                                  SDM_miss_sp,
                                  species_maps_dir,
                                  crosswalk_path,
                                   species_pools_path,
                                   dir_out){
  
  
  
  species_maps <- list.files(species_maps_dir, full.names = TRUE)
  
  predicts_sites <- qs2::qs_read(predicts_path) |>
    dplyr::distinct(SS,SSBS,Latitude,Longitude)
  

  crosswalk <-
    read.csv(crosswalk_path)
  
  # get the list of files that have the model outputs - in the methodology species that occupied fewer that 10 grid cells were excluded from
  # the analysis so I'm going to need to come up with something to do about those.
  
  SDM_sp <-
    list.files(SDM_probabilities_dir, full.names = TRUE)
  
  species_pools <- qs2::qs_read(species_pools_path)
  
  
  ## get just species names
  sp <-
    gsub(
      x = gsub(
        x = basename(SDM_sp),
        replacement = "",
        pattern = "\\.rds"
      ),
      pattern = "_",
      replacement = " "
    )
  
  sdm_species <- data.frame(Birdlife_Name = sp)
  

  resolve_species <-
    read.csv(SDM_miss_sp)
  
  
  ##### I'm just going to be looking at the environmental affinities for those species in the 0.9 species pools
  ##### this is because this is the species pool I am initially going on with and then I had to do less manual resolutions of mismatched species
  
  assembly_species <- unique(unlist(species_pools))
  
  
  ######
  ###### So given that the SDM files does not fully cover the assembly species with the update taxonomy the three years since the inital analysis was performed there are a few different pathways to get the environmental affinty of a species to a site. First its easy when the species is present just need to transform site coordinates and extract the probability of presense given the environmental variables.
  ###### Second if teh taxonomy has changed since the stduy was perfromed I was able to detct changes and point the species to a potential surrogate and then just load that SDM and extract the probability
  ###### Third if a species has split into multiple since I've taken the mean probability of presence of the three sister taxon
  ###### Fourth if there is no map, often due to too few cells being occupied to perform the analysis, probability of presence doesn't really come into it, if the site overlaps with the study site the probability is one and else zero.
  
  ## for species in teh assembly pool
  
  SDM_data_from_species <- function(species) {
    string <- grep(
     SDM_sp,
      pattern = paste(gsub(
        species,
        pattern = " ",
        replacement = "_"
      ),
      ".rds",
      sep = ""),
      value = TRUE,
      ignore.case = FALSE
    )
    
    
    return(string)
  }

  #####################
  
  
  
  too_few_cell_species <- function(sp, species_pools,sites) {
    
    data <- readRDS(grep(species_maps, pattern = paste("/",sp,"Birdlife",sep = ""), value = TRUE))
    
    sp_poly <-
      .SpatialCombinePolygons(
        .FilterPolygons(
          BL_data =  data,
          presence_codes = c(1, 2, 3),
          seasonal_codes = c(1, 2, 3),
          origin_codes = c(1, 2, 3)
        )
      ) |>
      terra::vect()
    
    ## filter ranges
    
    
    sp_sites <- .SiteFind(species = sp,
                          species_pools = species_pools)
    
    # get coordinates
    
    sp_sites <-
      terra::vect(
        sites |> 
          dplyr::filter(SSBS %in% sp_sites) |> 
          dplyr::distinct(Longitude, Latitude) |>
          as.matrix(), crs = "+proj=longlat +datum=WGS84 +no_defs"
      )
    
    
    
    return(as.numeric(terra::relate(sp_poly, sp_sites, "contains")))
    
  }
  
  ##########################
  ###########################
  ###########################
  
  
  mean_species_environmental_affinity <- function(sp, crosswalk, species_pools, sites) {
    
    split_sp <-
      crosswalk |> dplyr::filter(BirdLife.name == sp) |> dplyr::pull(Jetz.name)
    split_sp <-
      crosswalk |> dplyr::filter(Jetz.name == split_sp) |> dplyr::pull(BirdLife.name)
    split_sp <-
      split_sp[which(split_sp %in% sdm_species$Birdlife_Name)]
    
    ## find the sites that the original species occurs in the species pools
    
    sp_sites <- .SiteFind(species = sp, species_pools = species_pools)
    
    ## for each of those sites
    
    
    
    #get coords
    
    s_coords <- sites |> 
      dplyr::filter(SSBS %in% sp_sites) |> 
      dplyr::select(Longitude,Latitude) |> 
      as.matrix() |>
      terra::vect(crs = "+proj=longlat +datum=WGS84 +no_defs") |>
      terra::project("+proj=cea +lat_ts=30")
    
    # empty vector for probability of presence vals
    
    mean_vals <- c()
    
    # for each species

    for (i in 1:length(split_sp)) {
      # load the map
      
    
      sp_file <- SDM_data_from_species(split_sp[i])
      
      
      data <- readRDS(sp_file)
      
      # make raster
      sp_ras <- terra::rast(data[, c("lon", "lat", "current")],crs = "+proj=cea +lat_ts=30" )
      

      
      # extract probability of presence
      
      val <- matrix(terra::extract(x = sp_ras, y = s_coords)[,2], ncol = 1)
      
      # join mean to vector
      
      mean_vals <- cbind(mean_vals, val)
    }
    
    return(rowMeans(mean_vals, na.rm = TRUE))
    
    
  }
  
  
  
  #####################################
  
  
  environmental_probability_matrix <- function(species, species_pools, sites, crosswalk){
    
    mat <-
      matrix(
        rep(NA, nrow(sites)),
        nrow = nrow(sites),
        ncol = 1,
        dimnames = list(sites$SSBS, species)
      )
    
    
    if(grepl(species, pattern = "_")){
      sp_sites <- .SiteFind(species, species_pools = species_pools)
      
      mat[sp_sites,species] <- 1
      
      return(mat)
      
    }
    
    
    if (species %in% resolve_species$Birdlife_Name) {
      ## get that species data
      
      res_dat <-
        resolve_species |> dplyr::filter(Birdlife_Name == species)
      
      ## if the species is set to be dropped due to too few cells often or for two species quite mysterious reasons
      
      ###########################
      ################ SCENAIO 4
      ###########################
      
      if (res_dat$DROP) {
        ### load in the species poly
        
        sp_sites <- .SiteFind(species, species_pools = species_pools)
        
        
        mat[sp_sites, species] <-
          purrr::map_vec(sp_sites, .f = function(x){
            as.numeric(species %in% species_pools[[x]])
          })
        
        
        return(mat)
      }
      
      
      
      
      ###########################
      ################ SCENAIO 2
      ###########################
      
      # if the species wasn't do be dropped and there was a single species it could translate to then just point to that species file in the
      # folder
      
      if (!res_dat$DROP & res_dat$potential_species != "Mean") {
        # extract the file name
        
        sp_file <- SDM_data_from_species(res_dat$potential_species)
        
      }
      
      
      
      ###########################
      ################ SCENAIO 3
      ###########################
      
      ###### if it is the mean species then get the sister species from the crosswalk and those that SDM data is available.
      
      else {
        sp_sites <- .SiteFind(species, species_pools = species_pools)
        
        mat[sp_sites, species] <-
          mean_species_environmental_affinity(sp = species,crosswalk = crosswalk,
                                              species_pools = species_pools, sites = sites)
        
        
        
        return(mat)
      }
    } else {
      sp_file <- SDM_data_from_species(species)
    }
    
    
    data <- readRDS(sp_file)
    
    
    sp_ras <-
      terra::rast(data[, c("lon", "lat", "current")], crs = "+proj=cea +lat_ts=30")
    
    
    sp_sites <- .SiteFind(species, species_pools = species_pools)
    
    s_coords <- sites[,c("SSBS","Longitude","Latitude")]
    
    rownames(s_coords) <- s_coords$SSBS
    
    s_coords <-
      s_coords[sp_sites,c("Longitude","Latitude")] |>
      as.matrix() |>
      terra::vect(crs = "+proj=longlat +datum=WGS84 +no_defs") |>
      terra::project("+proj=cea +lat_ts=30")
    
    
 
    
    environ_p <- terra::extract(x = sp_ras, y = s_coords)[,2]
    
  
    
    mat[sp_sites, species] <- 
      environ_p
    
    ######################
    ######################
    ######################
    
    
    if (any(environ_p == 0)) {
      zero_site <- sp_sites[which(environ_p == 0)]
      
      
      ## if the site overlaps with the species range code as 1 else 0
      
      
     
      
      
      mat[zero_site, species] <-
        purrr::map_vec(zero_site, .f = function(x){
          
          as.numeric(species %in% species_pools[[x]])
          
        })
      
      
  
    }
    
    return(mat)
    
  }
  
  
  
  
  
  doParallel::registerDoParallel(cores = 10)
  
  environmental_probabilities <- foreach(species = assembly_species,
                                         .combine = "cbind",
                                         .inorder = FALSE,
                                         .packages = c("tidyverse", "terra", "sf")) %dopar% {
    
    environmental_prob <- environmental_probability_matrix(species, species_pools = species_pools, sites = predicts_sites, crosswalk = crosswalk)
    return(environmental_prob)
    
  }
  
  registerDoSEQ()

  
  
  sp_pool_environmental_affinities <- purrr::map(1:length(species_pools), .f = function(x){
    
    site <- names(species_pools)[x]
    sp_pool <- species_pools[[x]]
    
    return(environmental_probabilities[site,sp_pool])
    
    
  })
  
  
  names(sp_pool_environmental_affinities) <- predicts_sites$SSBS
  
file.out <- file.path(dir_out, "species-pool-environmental-affinties.qs")
  
  qs2::qs_save(sp_pool_environmental_affinities, file.out)
  
  return(file.out)
}