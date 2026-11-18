source("R/prepare_model_data.R")

source("R/compare_random_effect_structures.R")  # re_compare_formula()

# Frequentist (glmmTMB) counterpart of fit_model2d_bayesian_woody(): the combined
# model with a *woody interaction term, one model per climate driver.
#
# By default (by_group = TRUE) the random-intercept SD and the random within-species
# slope SD are both estimated separately for woody and non-woody taxa, matching
# fit_model2d_bayesian_woody(): diag(0 + woody | taxon_name) +
# diag(0 + std_<var>_c:woody | taxon_name), independent between growth forms (each
# taxon only has data in its own group's columns). This structure was best supported
# by AIC (see "Random-effect structure comparison - woody" in Analysis.qmd).
# by_group = FALSE fits the original model with both SDs pooled across growth forms.
# Formulas come from re_compare_formula() so they match the comparison exactly.
#
# Fixed effects are named as before. The woody random slope is named
# std_<var>_c:woodywoody, the same as the fixed interaction, so per-species slopes
# must be taken from fixef() + ranef() (species_within_slopes()), not coef().
fit_model2d_frequentist_woody <- function(data_model, by_group = TRUE) {
  structure <- if (by_group) "both_grp" else "pooled"
  model_formula <- function(var) re_compare_formula(var, structure, "glmmTMB")

  fit_MAP <- glmmTMB(model_formula("log_annual_precipitation"), data = data_model)
  fit_MAT <- glmmTMB(model_formula("average_temperature"), data = data_model)
  print(summary(fit_MAP))
  print(summary(fit_MAT))

  list_to_save <- list(fit_MAP, fit_MAT, data_model)
  saveRDS(list_to_save, paste0("outputs/model_fits/fit2d_frequentist_woody/", data_model$trait_name[1], ".RDS"))
}

# Fits fit_model2d_frequentist_woody() for every trait present in data,
# excluding huber_value/wood_density (woody-only traits with no non-woody
# data to support the *woody interaction). by_group is passed through (default TRUE).
fit_model2d_frequentist_woody_all_traits <- function(data, by_group = TRUE) {
  traits <- setdiff(unique(data$trait_name), WOODY_ONLY_TRAITS)
  for (trait in traits) {
    cat(trait, ":  ")
    fit_model2d_frequentist_woody(prepare_model_data(data, trait), by_group = by_group)
  }
}
