FocalSimilarity <- function(predicts_path,
                            predicts_species_overlap,
                            inventory_matrix,
                            template_map,
                            island_polys,
                            island_spp,
                            dir_out) {
  
  options(future.globals.maxSize = 891289600000000000)
  
  predicts <- qs2::qs_read(predicts_path)
  blank <- terra::rast(template_map)
  
  range_matrix <- qs2::qs_read(inventory_matrix)
  
  predicts_species_overlap_matrix <- qs2::qs_read(predicts_species_overlap)
  
  
  
  island_polys <- readRDS(island_polys)
  island_spp <- readRDS(island_spp)
  
  
  sites <-
    predicts |>
    dplyr::distinct(SS, SSBS, Longitude, Latitude) |>
    data.frame()
  
  sites[, "cell"] <- terra::cellFromXY(blank, sites[, c("Longitude", "Latitude")])
  
  
  #### Which cell is each site located in
  .GetSpeciesFromCell <- function(inventory_matrix, cell) {
    species <-  names(which(inventory_matrix[, as.numeric(cell)] > 0))
    return(species)
  }
  
  ## This loop is going to create a matrix with rows being the cell in which a site is located and columns being all other raster cells
  ## the value of each cell will the proportion of species shared between each cell
  
  site_cells <- unique(sites$cell)
  
  
  species_in_all_cells <- purrr::map(
    1:ncol(range_matrix),
    .f = function(x) {
      .GetSpeciesFromCell(inventory_matrix = range_matrix, cell = x)
    }
  )
  names(species_in_all_cells) <- colnames(range_matrix)
  
  

  
  ADFSimilarity <- function(inventory_matrix,
                            cell,
                            island = FALSE,
                            island_polys,
                            island_spp) {
    mat <-
      matrix(
        rep(NA, ncol(inventory_matrix)),
        nrow = 1,
        ncol = ncol(inventory_matrix),
        dimnames = list(cell, colnames(inventory_matrix))
      )
    
    #####
    
    
    if (island) {
      coords <-
        sites[as.character(sites$new_cell) == as.character(cell),] |>
        dplyr::distinct(Longitude, Latitude) |>
        slice(1) 
        
      coords <- terra::vect(as.matrix(coords),crs = "+proj=longlat +datum=WGS84 +no_defs")
      
      
        ## extract species that are present within a specified cell for the island form the focal cell species will be all species overlaping with that island
      
      for (isle_name in names(island_polys)) {
        
        
       isle_poly <-  terra::vect(island_polys[[isle_name]])
       terra::crs(isle_poly) <-  "+proj=longlat +datum=WGS84 +no_defs"
       
       
      
       
        if(as.numeric(terra::distance(x = coords, y = isle_poly,unit = "m")/1000) < 100) {
          focal_spp <- island_spp[[island]]
          break()
        }
      }
    } else {
      ## extract species that are present within a specified cell
      focal_spp <- species_in_all_cells[[as.character(cell)]]
    }
    ### get all other cells minus the focal cell
    
    other_cells <-
      as.character(setdiff(names(species_in_all_cells), as.character(cell)))
    
    
    ## cell shares all species with itself
    
    mat[, as.character(cell)] <- 1
    
    
    ### for all other cells extract species and calculate proportion of pseices shared with focal cell
    
    mat[, other_cells] <- purrr::map_vec(
      other_cells,
      .f = function(x) {
        .SpeciesCellSimilarity(species = focal_spp, cell = x, species_in_all_cells = species_in_all_cells)
      }
    )
    
    
    return(mat)
  }
  
  
  
  doParallel::registerDoParallel(cores = 10)
  
  adf_similarity_matrix <- foreach::foreach(
    cell = site_cells,
    .combine = "rbind",
    .packages = c("purrr", "dplyr")
  ) %dopar% {
    adf_sim <- ADFSimilarity(inventory_matrix = range_matrix,
                             island = FALSE,
                             cell = cell)
  }
  
  
  foreach::registerDoSEQ()
  
  
  future::plan(future::multisession(workers = 10))
  
  
  adf_richness <- furrr::future_map(
    .x = site_cells,
    .f = function(x) {
      sequence <- seq(0, 1, by = 0.1)
      
      
      adf_nspp <- purrr::map_vec(
        sequence,
        .f = function(y) {
          prop_cells <- which(adf_similarity_matrix[as.character(x), ] > y)
          
          rich <- length(unique(unlist(species_in_all_cells[prop_cells])))
          
        }
      )
      
      rich_df <- matrix(adf_nspp,
                        nrow = 1,
                        dimnames = list(x, sequence))
      
    }
  ) |>
    Reduce(f = "rbind")
  
  
  future::plan(future::sequential())
  
  
  problem_cells <- c(31241, 31868, 28739, 36584, 30779, 30778, 37311, 32227)
  
  
  sites$new_cell <- sites$cell
  
  for (problem_cell in problem_cells) {
    problem_site <- sites %>%
      dplyr::filter(cell == problem_cell) %>%
      dplyr::distinct(SSBS) %>%
      dplyr::pull() %>%
      as.character()
    problem_site <- problem_site[1]
    
    ## what are the species overlapping the site
    
    prob_spp <- names(which(predicts_species_overlap_matrix[problem_site, ] == 1))
    
    ### get the surround cells
    
    surround <- terra::adjacent(blank, problem_cell, directions = 8)[1, ]
    
    
    ### for each of the surrounding cells
    
    sound_prob <- c()
    
    for (s_cell in surround) {
      ### what species are within
      
      sound_spp <-  species_in_all_cells[[as.character(s_cell)]]
      
      ### calculate similarity between surrounding cell and site
      
      sound_sim <- length(which(prob_spp %in% sound_spp)) / length(prob_spp)
      
      # store
      
      sound_prob <- c(sound_prob, sound_sim)
      
    }
    
    # which of the surrounding cells has teh greatest similarity
    
    solve_cell <- surround[which.max(sound_prob)]
    
    ## reassign
    
    
    sites <- sites |>
      dplyr::mutate(new_cell = ifelse(cell == problem_cell, as.numeric(solve_cell), new_cell))
  }
  
  
  
  new_cells <- setdiff(sites$new_cell, sites$cell)
  
  
  future::plan(future::multisession(workers = 10))
  
  
  new_adf <- furrr::future_map(
    new_cells,
    .options = furrr::furrr_options(seed = NULL),
    .f = function(x) {
      mat <- ADFSimilarity(inventory_matrix = range_matrix,
                           cell = x,
                           island = FALSE)
    }
  ) |>
    Reduce(f = "rbind")
  
  
  future::plan(future::sequential())
  
  
  adf_similarity_matrix <- rbind(adf_similarity_matrix, new_adf)
  
  for (i in rownames(new_adf)) {
    adf <- blank
    
    terra::values(adf)[as.numeric(colnames(new_adf))] <-
      new_adf["36224", ]
    
    terra::plot(adf)
    title(i)
    
  }
  
  
  persistent_problems <- c(36224,32226,32229)
  
  
  for(prob in persistent_problems) {
    m <- ADFSimilarity(inventory_matrix = range_matrix,
                        cell = prob, 
                       island = TRUE,
                       island_polys = island_polys,
                       island_spp = island_spp)
    adf_similarity_matrix[as.character(prob), ] <- m
  }
  
  
  
  qs2::qs_save(sites, file.path(dir_out, "predicts-sites-cells.qs"))
  
  
  
  file.out <- file.path(dir_out, "adf-similarity-matrix.qs")
  
  
  qs2::qs_save(adf_similarity_matrix, file.out)
  
  return(file.out)
  
}
