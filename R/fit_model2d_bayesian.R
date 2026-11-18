# iter is per chain (half warmup). Leaf area gets more iterations (see
# BAYESIAN_ITER_OVERRIDES): with ~1,260 taxa its global intercept mixes slowly against the
# random intercepts: bulk ESS was 360 at 7500 and 396 (MAT) at 10000.
fit_model2d_bayesian <- function(data_model, iter = 7500) {
    fit_MAP_stan <- brm(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c) + (1 | taxon_name) + (0 + std_log_annual_precipitation_c | taxon_name), data = data_model, chains = 4, iter = iter, control = list(adapt_delta = 0.98), cores = 4)
    fit_MAT_stan <- brm(log_value ~ (std_average_temperature_mn + std_average_temperature_c) + (1 | taxon_name) + (0 + std_average_temperature_c | taxon_name), data = data_model, chains = 4, iter = iter, control = list(adapt_delta = 0.98), cores = 4)
    print(summary(fit_MAP_stan))
    print(summary(fit_MAT_stan))
    list_to_save <- list(fit_MAP_stan, fit_MAT_stan, data_model)
  saveRDS(list_to_save,paste0("outputs/model_fits/fit2d_bayesian/",data_model$trait_name[1],".RDS"))
}

source("R/prepare_model_data.R")

# Per-trait iterations, overriding the default 7500
BAYESIAN_ITER_OVERRIDES <- c(leaf_area = 15000)

# Fits fit_model2d_bayesian() for every trait present in data (or just `traits`).
fit_model2d_bayesian_all_traits <- function(data, traits = unique(data$trait_name)) {
  for (trait in traits) {
    cat(trait, ":  ")
    iter <- if (trait %in% names(BAYESIAN_ITER_OVERRIDES)) BAYESIAN_ITER_OVERRIDES[[trait]] else 7500
    fit_model2d_bayesian(prepare_model_data(data, trait), iter = iter)
  }
}
