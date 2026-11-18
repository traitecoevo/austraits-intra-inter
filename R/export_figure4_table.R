source("R/prepare_model_data.R")
source("R/build_compressed_trait_figure.R")
source("R/test_group_slope_zero.R")

# Plain-text trait names for the table's trait header rows
TRAIT_TABLE_NAMES <- c(
  leaf_delta13C      = "Ci:Ca",
  leaf_N_per_area    = "Leaf N per area",
  leaf_mass_per_area = "Leaf mass per area (LMA)",
  leaf_area          = "Leaf area",
  huber_value        = "SA:LA (Huber value)",
  wood_density       = "Wood density",
  seed_dry_mass      = "Seed mass"
)

format_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# Builds one table with a block of rows per trait (in the order given).
# trait_rows(trait) must return one row per variable with columns Variable,
# MAP_estimate, MAP_std_error, MAP_p, MAT_estimate, MAT_std_error, MAT_p.
# Each trait gets a bold merged header row with its name above its variables,
# with MAP / MAT headers spanning their estimate, std. error and p-value columns.
# footer: optional footnote lines added under the table.
build_trait_flextable <- function(traits, trait_rows, footer = NULL) {
  table_data <- map_dfr(traits, function(trait) {
    trait_rows(trait) |>
      mutate(Trait = TRAIT_TABLE_NAMES[[trait]], .before = 1)
  }) |>
    mutate(Trait = factor(Trait, levels = TRAIT_TABLE_NAMES[traits]))

  # Trait names become merged group rows above each trait's variables
  ft <- flextable::as_grouped_data(table_data, groups = "Trait") |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::set_header_labels(
      Variable      = "Variable",
      MAP_estimate  = "Estimate", MAP_std_error = "Std. error", MAP_p = "p value",
      MAT_estimate  = "Estimate", MAT_std_error = "Std. error", MAT_p = "p value"
    ) |>
    flextable::add_header_row(values = c("", "MAP", "MAT"), colwidths = c(1, 3, 3)) |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = 2:7, part = "body", align = "center") |>
    flextable::bold(part = "header") |>
    flextable::bold(i = ~ !is.na(Trait), part = "body")

  if (!is.null(footer)) ft <- flextable::add_footer_lines(ft, footer)

  ft |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::autofit()
}

# Woody-split traits (Figure 4's traits, in its order), reading the saved
# frequentist woody fits, list(fit_MAP, fit_MAT, data_model); trait_rows(fits).
build_woody_trait_flextable <- function(trait_rows) {
  traits <- setdiff(FIGURE_TRAIT_ORDER, WOODY_ONLY_TRAITS)

  build_trait_flextable(traits, function(trait) {
    fits <- readRDS(paste0("outputs/model_fits/fit2d_frequentist_woody/", trait, ".RDS"))
    stopifnot(length(fits) == 3)
    trait_rows(fits)
  })
}

# Formats estimates / SEs / p-values for one climate model into the MAP_* or MAT_* columns.
format_env_columns <- function(tbl, env) {
  tbl |>
    transmute(
      Variable,
      !!paste0(env, "_estimate")  := sprintf("%.3f", estimate),
      !!paste0(env, "_std_error") := sprintf("%.3f", std.error),
      !!paste0(env, "_p")         := format_p(p.value)
    )
}

# Writes a single Word table (flextable) of the fixed effects behind Figure 4: the
# frequentist woody-split glmmTMB models, one row per model term.
export_figure4_table <- function(path = "outputs/tables/Figure4_fixed_effects_table.docx") {
  # Terms are named after the climate variable; strip it so MAP and MAT rows line up.
  variable_labels <- c(
    "(Intercept)"   = "Intercept",
    "mn"            = "Across-species slope",
    "c"             = "Within-species slope",
    "woodywoody"    = "Woody",
    "mn:woodywoody" = "Across-species slope × woody",
    "c:woodywoody"  = "Within-species slope × woody"
  )

  extract_model <- function(fit, climate_var, env) {
    tidy(fit, effects = "fixed") |>
      mutate(Variable = variable_labels[str_remove(term, paste0("std_", climate_var, "_"))]) |>
      format_env_columns(env)
  }

  ft <- build_woody_trait_flextable(function(fits) {
    left_join(
      extract_model(fits[[1]], "log_annual_precipitation", "MAP"),
      extract_model(fits[[2]], "average_temperature", "MAT"),
      by = "Variable"
    ) |>
      arrange(match(Variable, variable_labels))
  })

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}

# Writes a single Word table (flextable) testing whether each growth form's slope
# differs from zero (NOT whether woody differs from non-woody), for the across-
# (β_mn) and within-species (β_c) slopes of the Figure 4 woody-split models.
#   non-woody slope = β_base
#   woody slope     = β_base + β_interaction, with
#                     se = sqrt(Var(base) + Var(interaction) + 2 Cov(base, interaction))
# (see test_group_slope_zero()).
export_group_slope_table <- function(path = "outputs/tables/Figure4_group_slope_vs_zero_table.docx") {
  slope_specs <- tribble(
    ~component, ~group,      ~Variable,
    "mn",       "non-woody", "Across-species slope, non-woody",
    "mn",       "woody",     "Across-species slope, woody",
    "c",        "non-woody", "Within-species slope, non-woody",
    "c",        "woody",     "Within-species slope, woody"
  )

  test_model <- function(fit, climate_var, env) {
    slope_specs |>
      mutate(test = map2(component, group, \(component, group)
        test_group_slope_zero(fit, paste0("std_", climate_var, "_", component), group) |>
          select(estimate, std.error, p.value))) |>
      unnest(test) |>
      format_env_columns(env)
  }

  ft <- build_woody_trait_flextable(function(fits) {
    left_join(
      test_model(fits[[1]], "log_annual_precipitation", "MAP"),
      test_model(fits[[2]], "average_temperature", "MAT"),
      by = "Variable"
    )
  })

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}
