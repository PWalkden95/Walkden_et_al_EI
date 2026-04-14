TrophicNicheInSite <- function(predicts_path,
                                dir_out){
  
  predicts <- qs2::qs_read(predicts_path)
  
  
  study <- unique(predicts$SS)
  
  
  data <- predicts |>
    dplyr::distinct(SSBS,Birdlife_Name,Trophic_Niche) |>
    dplyr::group_by(SSBS,Trophic_Niche) |>
    dplyr::summarise(count = dplyr::n())
  

file.out <- file.path(dir_out, "guilds-in-observed-sites.qs")  


qs2::qs_save(data,file.out)

return(file.out)
  
}