fit_model2d_bayesian_woody_seperated <- function(data_model) {

  data_woody <- data_model %>% filter(as.character(woody) == "woody")
  data_non_woody <- data_model %>% filter(as.character(woody) != "woody")

  fit_group <- function(dat) {
    fit_MAP_stan <- brm(
      log_value ~ std_log_annual_precipitation_mn + std_log_annual_precipitation_c +
        (1 | taxon_name) + (0 + std_log_annual_precipitation_c | taxon_name),
      data = dat,
      chains = 4,
      iter = 10000,
      control = list(adapt_delta = 0.95),
      cores = 4
    )

    fit_MAT_stan <- brm(
      log_value ~ std_average_temperature_mn + std_average_temperature_c +
        (1 | taxon_name) + (0 + std_average_temperature_c | taxon_name),
      data = dat,
      chains = 4,
      iter = 10000,
      control = list(adapt_delta = 0.95),
      cores = 4
    )

    list(fit_MAP_stan = fit_MAP_stan, fit_MAT_stan = fit_MAT_stan)
  }

  list_to_save <- list(
    woody = if (nrow(data_woody) > 0) fit_group(data_woody) else NULL,
    non_woody = if (nrow(data_non_woody) > 0) fit_group(data_non_woody) else NULL,
    data_model = data_model
  )

  output_dir <- "outputs/model_fits/fit2d_bayesian_woody_seperated"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  saveRDS(list_to_save, paste0(output_dir, "/", data_model$trait_name[1], ".RDS"))
}