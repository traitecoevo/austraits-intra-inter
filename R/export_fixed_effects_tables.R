source("R/prepare_model_data.R")

# Extract fixed effects from a frequentist glmmTMB model.
# Returns a tidy tibble with z-statistic and p-value columns.
extract_freq_fixed <- function(fit, climate_var, trait) {
  tidy(fit, effects = "fixed") %>%
    select(term, estimate, std.error, statistic, p.value) %>%
    mutate(
      estimate         = round(estimate,  5),
      std.error        = round(std.error, 5),
      statistic        = round(statistic, 5),
      p.value          = round(p.value,   5),
      climate_variable = climate_var,   # label which climate driver this model used
      trait            = trait
    )
}

# Extract fixed effects from a Bayesian brms model.
# brms tidy() returns posterior mean/SD + credible interval; no z/p columns.
extract_bayes_fixed <- function(fit, climate_var, trait) {
  tidy(fit, effects = "fixed") %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    mutate(
      estimate         = round(estimate,  5),
      std.error        = round(std.error, 5),
      conf.low         = round(conf.low,  5),
      conf.high        = round(conf.high, 5),
      climate_variable = climate_var,
      trait            = trait
    )
}

# Writes a single Word table (flextable) of the fixed effects of the combined
# (non-split) models behind Figure 3, styled like export_figure4_table(): one block
# per trait present in data (in Figure 3's order) with the trait name as a header
# row, one row per variable, and estimate / std. error / p value for the
# precipitation (MAP) and temperature (MAT) models side by side.
#
# huber_value: frequentist model did not converge reliably; use the Bayesian
# (brms) model instead — estimate = posterior mean, std. error = posterior SD,
# p value = two-sided posterior tail probability, 2 * min(P(b <= 0), P(b >= 0)),
# as used for the Figure 3 significance keys (noted in a table footnote).
# All other traits: use the frequentist glmmTMB model.
# Both are saved as list(fit_MAP, fit_MAT, data_model).
export_fixed_effects_tables_all_traits <- function(data, path = "outputs/tables/Figure3_fixed_effects_table.docx") {
  source("R/export_figure4_table.R")  # build_trait_flextable(), format_env_columns()

  traits <- intersect(FIGURE_TRAIT_ORDER, unique(data$trait_name))

  variable_labels <- c(
    "(Intercept)" = "Intercept",
    "Intercept"   = "Intercept",
    "mn"          = "Across-species slope",
    "c"           = "Within-species slope"
  )

  extract_model <- function(fit, climate_var, env) {
    if (inherits(fit, "brmsfit")) {
      draws <- posterior::as_draws_df(fit)
      fixed <- brms::fixef(fit)
      tbl <- tibble(
        term      = rownames(fixed),
        estimate  = fixed[, "Estimate"],
        std.error = fixed[, "Est.Error"],
        p.value   = map_dbl(term, \(t) {
          d <- draws[[paste0("b_", t)]]
          2 * min(mean(d <= 0), mean(d >= 0))
        })
      )
    } else {
      tbl <- tidy(fit, effects = "fixed")
    }
    tbl |>
      mutate(Variable = variable_labels[str_remove(term, paste0("std_", climate_var, "_"))]) |>
      format_env_columns(env)
  }

  ft <- build_trait_flextable(
    traits,
    function(trait) {
      model_dir <- if (trait == "huber_value") "fit2d_bayesian" else "fit2d_frequentist"
      fits <- readRDS(paste0("outputs/model_fits/", model_dir, "/", trait, ".RDS"))
      stopifnot(length(fits) == 3)
      left_join(
        extract_model(fits[[1]], "log_annual_precipitation", "MAP"),
        extract_model(fits[[2]], "average_temperature", "MAT"),
        by = "Variable"
      ) |>
        arrange(match(Variable, unique(variable_labels)))
    },
    footer = if ("huber_value" %in% traits)
      "SA:LA (Huber value) is from the Bayesian model: estimate = posterior mean, std. error = posterior SD, p value = two-sided posterior tail probability."
  )

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}

# Writes a "<trait>_fixed_effects_table_woody.docx" for every trait present in
# data (excluding huber_value/wood_density, which have no *woody fit), pulling
# fixed effects out of the saved fit2d_bayesian_woody/fit2d_frequentist_woody RDS.
#
# huber_value: use the Bayesian woody model for the same reason as above
# (frequentist convergence issues); all other traits use the frequentist woody model.
export_fixed_effects_tables_woody_all_traits <- function(data) {
  traits <- setdiff(unique(data$trait_name), WOODY_ONLY_TRAITS)

  for (trait in traits) {

    if (trait == "huber_value") {

      model <- readRDS(paste0("outputs/model_fits/fit2d_bayesian_woody/", trait, ".RDS"))

      # Bayesian woody models are saved as:
      #   length 3: list(fit_MAP_stan, fit_MAT_stan, data_model)
      #             model[[1]] = precipitation (MAP), model[[2]] = temperature (MAT)
      #   length 2: list(fit_stan, data_model) — single combined model
      if (length(model) == 3) {
        fixed_tbl <- bind_rows(
          extract_bayes_fixed(model[[1]], "precipitation", trait),
          extract_bayes_fixed(model[[2]], "temperature",   trait)
        )
      } else if (length(model) == 2) {
        fixed_tbl <- extract_bayes_fixed(model[[1]], "combined", trait)
      }

    } else {

      model <- readRDS(paste0("outputs/model_fits/fit2d_frequentist_woody/", trait, ".RDS"))

      # Frequentist woody models are saved as:
      #   length 3: list(fit_MAP, fit_MAT, data_model)
      #             model[[1]] = precipitation (MAP), model[[2]] = temperature (MAT)
      #   length 2: list(fit_lme, data_model) — single combined model
      if (length(model) == 3) {
        fixed_tbl <- bind_rows(
          extract_freq_fixed(model[[1]], "precipitation", trait),
          extract_freq_fixed(model[[2]], "temperature",   trait)
        )
      } else if (length(model) == 2) {
        fixed_tbl <- extract_freq_fixed(model[[1]], "combined", trait)
      }

    }

    doc <- read_docx()
    doc <- body_add_table(doc, value = fixed_tbl, style = "table_template")
    print(doc, target = paste0(trait, "_fixed_effects_table_woody.docx"))
  }
}
