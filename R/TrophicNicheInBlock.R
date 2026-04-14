TrophicNicheInBlock <- function(predicts_path,
                               dir_out){
  
  predicts <- qs2::qs_read(predicts_path)
  
  
  
  data <- predicts |>
    dplyr::distinct(SSB,Birdlife_Name,Trophic_Niche) |>
    dplyr::group_by(SSB,Trophic_Niche) |>
    dplyr::summarise(count = dplyr::n())
  
  
  file.out <- file.path(dir_out, "guilds-in-observed-block.qs")  
  
  
  qs2::qs_save(data,file.out)
  
  return(file.out)
  
}