RegionalSpeciesPool <- function(predicts_path,
                                predicts_site_cells,
                                adf_matrix, 
                                inventory_matrix,
                                island_spp,
                                dir_out){
  
  
  island_species <- readRDS(island_spp)
  
  range_matrix <- qs2::qs_read(inventory_matrix)
  
  predicts <- qs2::qs_read(predicts_path)

  adf <- qs2::qs_read(adf_matrix)
  
  site_cells <- qs2::qs_read(predicts_site_cells)
  

  species <- rownames(range_matrix)
  
sites <- as.character(site_cells$SSBS)  
 
 biogeographic_species_pool <-  purrr::map(sites, .f = function(x){
 
    
    block <- unique(as.character(predicts[predicts$SSBS == x, "SSB"]))
    
    block_species <- predicts |>
      dplyr::filter(SSB == block) |>
      dplyr::distinct(Birdlife_Name) |>
      dplyr::pull()
  
    cell <- site_cells[site_cells$SSBS == x, "new_cell"]
    
   regional_cells <- colnames(adf)[which(adf[as.character(cell),]  > 0.9)]
    
   
   regional_species_pool <- unique(c(block_species,
                              species[which(rowSums(as.matrix(range_matrix[,regional_cells])) > 0)]))
    

   
   if(cell == 36224){
     regional_species_pool <- unique(c(regional_species_pool, island_species$Comoros))
   }

if(cell == 32226){
  regional_species_pool <- unique(c(regional_species_pool, island_species$Sao_Tome))
}
 
if(cell == 32229){
  regional_species_pool <- unique(c(regional_species_pool, island_species$Principe))
}
  
   
  return(regional_species_pool) 
   
  })
  
 names(biogeographic_species_pool) <- site_cells$SSBS
 
 
 file.out <- file.path(dir_out,"biogeographic_regional_species_pools.qs")
  

 qs2::qs_save(biogeographic_species_pool, file.out)
 
 return(file.out)
   
}