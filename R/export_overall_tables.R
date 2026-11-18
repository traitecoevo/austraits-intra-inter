source("R/export_figure4_table.R")  # build_trait_flextable(), format_env_columns(), TRAIT_TABLE_NAMES

# Overall (OLS) responses: the simple linear models drawn as the grey "overall
# response" lines in Figures 3 and 4, ignoring species identity. Fitted to the same
# data as the mixed models (prepare_model_data()), one model per climate driver:
#   Figure 3: log_value ~ std_X
#   Figure 4: log_value ~ std_X * woody
# Each table has the same layout as the fixed-effects tables: a bold header row per
# trait, one row per term, and estimate / std. error / p value for MAP and MAT.

overall_variable_labels <- c(
  "(Intercept)" = "Intercept",
  "slope"       = "Slope",
  "woodywoody"  = "Woody",
  "slope:woodywoody" = "Slope × woody"
)

# Tidies one OLS fit, renaming the climate term (std_X) to "slope" so the MAP
# and MAT rows line up.
tidy_overall <- function(fit, climate_var, env) {
  broom::tidy(fit) |>
    mutate(Variable = overall_variable_labels[str_replace(term, paste0("std_", climate_var), "slope")]) |>
    format_env_columns(env)
}

overall_rows <- function(data, trait, woody) {
  data_model <- prepare_model_data(data, trait)
  rhs <- if (woody) " * woody" else ""
  fit_MAP <- lm(as.formula(paste0("log_value ~ std_log_annual_precipitation", rhs)), data = data_model)
  fit_MAT <- lm(as.formula(paste0("log_value ~ std_average_temperature", rhs)), data = data_model)

  left_join(
    tidy_overall(fit_MAP, "log_annual_precipitation", "MAP"),
    tidy_overall(fit_MAT, "average_temperature", "MAT"),
    by = "Variable"
  ) |>
    arrange(match(Variable, overall_variable_labels))
}

# Figure 3: overall response for every trait (all growth forms combined).
export_overall_table_figure3 <- function(data, path = "outputs/tables/Figure3_overall_response_table.docx") {
  traits <- intersect(FIGURE_TRAIT_ORDER, unique(data$trait_name))
  ft <- build_trait_flextable(traits, \(trait) overall_rows(data, trait, woody = FALSE))
  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}

# Figure 4: overall response with a woody interaction, for the woody-split traits.
# Woody is the shift from the non-woody intercept; Slope × woody is the difference
# between the woody and non-woody slopes (Slope is the non-woody slope).
export_overall_table_figure4 <- function(data, path = "outputs/tables/Figure4_overall_response_table.docx") {
  traits <- setdiff(intersect(FIGURE_TRAIT_ORDER, unique(data$trait_name)), WOODY_ONLY_TRAITS)
  ft <- build_trait_flextable(traits, \(trait) overall_rows(data, trait, woody = TRUE))
  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}
