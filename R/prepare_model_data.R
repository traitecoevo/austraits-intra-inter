# huber_value and wood_density are only recorded for woody species, so they have
# no non-woody subset and cannot support a woody/non-woody interaction term.
WOODY_ONLY_TRAITS <- c("huber_value", "wood_density")

# For a given trait/growth-form subset, builds the model-ready data frame used
# by find_ITV() and the fit_model2d_* fitting functions.
prepare_model_data <- function(data, trait, woody_filter = NULL) {
  data_model <- data |>
    filter(trait_name == trait)

  if (!is.null(woody_filter)) {
    data_model <- data_model |>
      filter(woody == woody_filter)
  }

  data_model |>
    select(trait_name, log_value, std_average_temperature, std_log_annual_precipitation,
           std_log_annual_precipitation_mn, std_log_annual_precipitation_c,
           std_average_temperature_mn, std_average_temperature_c,
           log_annual_precipitation_mn, log_annual_precipitation_c,
           average_temperature_mn, average_temperature_c, taxon_name, woody) |>
    drop_na() |>
    mutate(std_log_value = scale(log_value))
}
