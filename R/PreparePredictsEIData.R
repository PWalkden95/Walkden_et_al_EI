PreparePredictsEIData <- function(predicts_avonet_path, dir_out) {
  # Read in the PREDICTS-AVONET dataset from a .qs file
  predicts <- qs2::qs_read(predicts_avonet_path)

  # List of study names to exclude from the dataset
  drop_studies <- c(
    "Gagne",
    "CarvajalCastro",
    "Tekeuchi",
    "Deguchi",
    "Lorenzon",
    "Hastedt",
    "Humphrey",
    "Pena",
    "Silvetti"
  )

  # Filter and clean the dataset
  ei_predicts <- predicts |>
    # Remove rows where SS matches any study in drop_studies
    dplyr::filter(!grepl(SS, pattern = paste(drop_studies, collapse = "|"))) |>
    # Standardise 'Predominant_habitat' values for 'secondary' habitats
    dplyr::mutate(
      Predominant_habitat = ifelse(
        grepl(Predominant_habitat,
          pattern = "secondary", ignore.case = TRUE
        ),
        "Secondary",
        paste(Predominant_habitat)
      ),
      # Standardise 'Predominant_habitat' values for 'Primary forest'
      Predominant_habitat = ifelse(
        grepl(Predominant_habitat,
          pattern = "Primary forest", ignore.case = TRUE
        ),
        "Primary vegetation",
        paste(Predominant_habitat)
      )
    ) |>
    # Remove rows with unwanted habitat types
    dplyr::filter(!(Predominant_habitat %in% c(
      "Primary non-forest",
      "Cannot decide"
    )))

  # Define output file path
  file.out <- file.path(dir_out, "predicts-ei-data-prepared.qs")

  # Save the cleaned dataset as a .qs file
  qs2::qs_save(ei_predicts, file.out)
  
  return(file.out)
}
