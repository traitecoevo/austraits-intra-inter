source("R/prepare_model_data.R")

# Compares random-effect structures for the woody-split models behind Figure 4,
# for each woody-split trait and climate variable (fixed effects identical in all):
#
#   pooled     (1 | taxon) + (0 + Env_w | taxon)
#              one random-intercept SD and one random-slope SD shared by both growth forms
#   slope_grp  (1 | taxon) + slope SD per growth form
#              the previous default (random-slope SD split by growth form only)
#   both_grp   intercept SD and slope SD each estimated separately per growth form
#
# AIC: glmmTMB fits of all three structures (fast).
# LOO: brms fits of pooled vs both_grp (PSIS-LOO; slow, so results are cached).
#
# Env_m = std_<var>_mn, Env_w = std_<var>_c. Per-group SDs are independent between
# growth forms (each taxon belongs to one growth form, so no correlation is estimable):
# glmmTMB diag(0 + ... : woody | taxon), brms gr(taxon, by = woody).

RE_COMPARE_ENV_VARS <- c(MAP = "log_annual_precipitation", MAT = "average_temperature")
RE_COMPARE_STRUCTURES <- c("pooled", "slope_grp", "both_grp")

re_compare_formula <- function(var, structure, engine = c("glmmTMB", "brms")) {
  engine <- match.arg(engine)
  mn <- paste0("std_", var, "_mn")
  cc <- paste0("std_", var, "_c")
  fixed <- paste0("log_value ~ (", mn, " + ", cc, ") * woody")

  random <- switch(paste(engine, structure),
    "glmmTMB pooled"    = paste0("(1 | taxon_name) + (0 + ", cc, " | taxon_name)"),
    "glmmTMB slope_grp" = paste0("(1 | taxon_name) + diag(0 + ", cc, ":woody | taxon_name)"),
    "glmmTMB both_grp"  = paste0("diag(0 + woody | taxon_name) + diag(0 + ", cc, ":woody | taxon_name)"),
    # brms: taxon_slope is a copy of taxon_name so the intercept and slope terms are
    # separate (uncorrelated) grouping terms, as in fit_model2d_bayesian_woody()
    "brms pooled"       = paste0("(1 | taxon_name) + (0 + ", cc, " | taxon_slope)"),
    "brms slope_grp"    = paste0("(1 | taxon_name) + (0 + ", cc, " | gr(taxon_slope, by = woody))"),
    "brms both_grp"     = paste0("(1 | gr(taxon_name, by = woody)) + (0 + ", cc, " | gr(taxon_slope, by = woody))")
  )
  as.formula(paste(fixed, "+", random))
}

# AIC table for all traits x climate variables x structures (glmmTMB).
compare_re_aic <- function(data, traits) {
  tidyr::expand_grid(trait = traits, env = names(RE_COMPARE_ENV_VARS)) |>
    purrr::pmap_dfr(function(trait, env) {
      d <- prepare_model_data(data, trait)
      purrr::map_dfr(RE_COMPARE_STRUCTURES, function(structure) {
        fit <- glmmTMB::glmmTMB(re_compare_formula(RE_COMPARE_ENV_VARS[[env]], structure, "glmmTMB"), data = d)
        tibble(trait = trait, env = env, structure = structure,
               AIC = AIC(fit), df = attr(logLik(fit), "df"),
               converged = isTRUE(fit$sdr$pdHess))
      })
    }) |>
    group_by(trait, env) |>
    mutate(dAIC = AIC - min(AIC)) |>
    ungroup()
}

# LOO for one trait x climate variable x structure (brms), cached so the slow fits
# only run once. Only the loo object and a few diagnostics are kept, not the fit.
re_compare_loo_one <- function(data, trait, env, structure, iter = 4000,
                               cache_dir = "outputs/model_fits/re_structure_loo") {
  cache_path <- file.path(cache_dir, paste0(trait, "_", env, "_", structure, ".RDS"))
  if (file.exists(cache_path)) return(readRDS(cache_path))

  d <- prepare_model_data(data, trait) |> mutate(taxon_slope = taxon_name)
  fit <- brms::brm(re_compare_formula(RE_COMPARE_ENV_VARS[[env]], structure, "brms"),
                   data = d, chains = 4, iter = iter, cores = 4,
                   control = list(adapt_delta = 0.98), seed = 1, refresh = 0)
  np <- brms::nuts_params(fit)
  out <- list(
    loo = brms::loo(fit),
    diagnostics = tibble(trait = trait, env = env, structure = structure,
                         max_rhat = max(brms::rhat(fit), na.rm = TRUE),
                         divergences = sum(np$Value[np$Parameter == "divergent__"]))
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache_path)
  rm(fit); gc(verbose = FALSE)
  out
}

# LOO comparison of pooled vs both_grp for each trait x climate variable: elpd of each,
# the difference (both_grp - pooled, positive = by-group better) with its SE, Pareto-k
# diagnostics and sampler diagnostics.
compare_re_loo <- function(data, traits, structures = c("pooled", "both_grp"), iter = 4000) {
  tidyr::expand_grid(trait = traits, env = names(RE_COMPARE_ENV_VARS)) |>
    purrr::pmap_dfr(function(trait, env) {
      res <- purrr::map(structures, \(s) re_compare_loo_one(data, trait, env, s, iter = iter))
      names(res) <- structures
      comp <- brms::loo_compare(res[[1]]$loo, res[[2]]$loo)
      # loo_compare rows are ordered best first; elpd_diff is relative to the best model
      diff <- comp[, "elpd_diff"]; se <- comp[, "se_diff"]
      second <- rownames(comp)[2]
      sign <- if (second == "model1") 1 else -1   # +: structures[2] better than structures[1]
      tibble(
        trait = trait, env = env,
        elpd_pooled   = res[[1]]$loo$estimates["elpd_loo", "Estimate"],
        elpd_both_grp = res[[2]]$loo$estimates["elpd_loo", "Estimate"],
        elpd_diff     = sign * abs(diff[2]),
        se_diff       = se[2],
        high_k_pooled   = sum(res[[1]]$loo$diagnostics$pareto_k > 0.7),
        high_k_both_grp = sum(res[[2]]$loo$diagnostics$pareto_k > 0.7),
        max_rhat      = max(res[[1]]$diagnostics$max_rhat, res[[2]]$diagnostics$max_rhat),
        divergences   = res[[1]]$diagnostics$divergences + res[[2]]$diagnostics$divergences
      )
    })
}

# Word table (flextable, same layout as the other tables) of compare_re_aic() output:
# one block per trait (bold header row), one row per random-effect structure, and for
# MAP and MAT the AIC, delta AIC (relative to the best structure for that trait and
# climate variable) and number of parameters. The best structure is in bold.
export_re_aic_table <- function(re_aic, path = "outputs/tables/random_effect_structure_AIC_table.docx") {
  source("R/export_figure4_table.R")  # FIGURE_TRAIT_ORDER, TRAIT_TABLE_NAMES

  traits <- intersect(FIGURE_TRAIT_ORDER, unique(re_aic$trait))
  structure_labels <- c(pooled    = "Both SDs pooled",
                        slope_grp = "Slope SD by growth form",
                        both_grp  = "Intercept and slope SDs by growth form")

  table_data <- re_aic |>
    transmute(trait, env, structure,
              AIC  = sprintf("%.1f", AIC),
              dAIC = sprintf("%.1f", dAIC),
              df   = as.character(df)) |>
    pivot_wider(names_from = env, values_from = c(AIC, dAIC, df), names_glue = "{env}_{.value}") |>
    arrange(match(trait, traits), match(structure, names(structure_labels))) |>
    transmute(Trait = factor(TRAIT_TABLE_NAMES[trait], levels = TRAIT_TABLE_NAMES[traits]),
              Structure = structure_labels[structure],
              MAP_AIC, MAP_dAIC, MAP_df, MAT_AIC, MAT_dAIC, MAT_df)

  ft <- flextable::as_grouped_data(table_data, groups = "Trait") |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::set_header_labels(
      Structure = "Random-effect variances",
      MAP_AIC = "AIC", MAP_dAIC = "ΔAIC", MAP_df = "df",
      MAT_AIC = "AIC", MAT_dAIC = "ΔAIC", MAT_df = "df"
    ) |>
    flextable::add_header_row(values = c("", "MAP", "MAT"), colwidths = c(1, 3, 3)) |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = 2:7, part = "body", align = "center") |>
    flextable::bold(part = "header") |>
    flextable::bold(i = ~ !is.na(Trait), part = "body") |>
    # best structure (delta AIC = 0) in bold, for each climate variable
    flextable::bold(i = ~ MAP_dAIC == "0.0", j = c("MAP_AIC", "MAP_dAIC", "MAP_df"), part = "body") |>
    flextable::bold(i = ~ MAT_dAIC == "0.0", j = c("MAT_AIC", "MAT_dAIC", "MAT_df"), part = "body") |>
    flextable::add_footer_lines(paste(
      "AIC of the woody-split glmmTMB models (fixed effects identical) with the random-intercept",
      "and random within-species slope SDs pooled across growth forms, the slope SD estimated",
      "separately for woody and non-woody taxa, or both estimated separately. ΔAIC is relative",
      "to the best structure (bold) for each trait and climate variable; df = number of parameters.")) |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::autofit()

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}
