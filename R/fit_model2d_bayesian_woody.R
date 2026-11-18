source("R/prepare_model_data.R")

source("R/compare_random_effect_structures.R")  # re_compare_formula()

# Same model as fit_model2d_bayesian() (chains, iter, adapt_delta), with a
# *woody interaction term added to compare woody vs non-woody slopes.
#
# By default (by_group = TRUE) the random-intercept SD and the random within-species
# slope SD are both estimated separately for woody and non-woody taxa, via
# (1 | gr(taxon_name, by = woody)) + (0 + std_<var>_c | gr(taxon_slope, by = woody)).
# This structure was best supported by AIC (see "Random-effect structure comparison -
# woody" in Analysis.qmd); pooling lets the more numerous woody taxa set a single SD for
# both groups. by_group = FALSE fits the original model with both SDs pooled.
# taxon_slope is a copy of taxon_name so the intercept and slope are separate,
# uncorrelated grouping terms. gr(..., by = woody) requires each taxon to belong to a
# single growth form. Formulas come from re_compare_formula() so they match the
# glmmTMB models and the comparison exactly.
#
# Fixed-effect names and order are the same under both options, so everything
# downstream that reads outputs/model_fits/fit2d_bayesian_woody/ (B_r, S_w) is
# unaffected; only the random-effect SD names change, to
# sd_taxon_name__Intercept:woody<group> and sd_taxon_slope__std_<var>_c:woody<group>.
fit_model2d_bayesian_woody <- function(data_model, by_group = TRUE) {
  data_model <- data_model |> mutate(taxon_slope = taxon_name)

  structure <- if (by_group) "both_grp" else "pooled"
  formula_MAP <- re_compare_formula("log_annual_precipitation", structure, "brms")
  formula_MAT <- re_compare_formula("average_temperature", structure, "brms")

  fit_MAP_stan <- brm(formula_MAP, data = data_model, chains = 4, iter = 7500, control = list(adapt_delta = 0.98), cores = 4)
  fit_MAT_stan <- brm(formula_MAT, data = data_model, chains = 4, iter = 7500, control = list(adapt_delta = 0.98), cores = 4)
  print(summary(fit_MAP_stan))
  print(summary(fit_MAT_stan))
  list_to_save <- list(fit_MAP_stan, fit_MAT_stan, data_model)
  saveRDS(list_to_save, paste0("outputs/model_fits/fit2d_bayesian_woody/", data_model$trait_name[1], ".RDS"))
}

# Fits fit_model2d_bayesian_woody() for every trait present in data, excluding
# huber_value/wood_density (woody-only traits with no non-woody data to
# support the *woody interaction). by_group is passed through (default TRUE).
fit_model2d_bayesian_woody_all_traits <- function(data, by_group = TRUE) {
  traits <- setdiff(unique(data$trait_name), WOODY_ONLY_TRAITS)
  for (trait in traits) {
    cat(trait, ":  ")
    fit_model2d_bayesian_woody(prepare_model_data(data, trait), by_group = by_group)
  }
}
