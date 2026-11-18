# Trait ranges for the observations plotted as grey dots in Figure 3: the per-cell species
# means in data_traits_climate (Analysis.qmd), which are exactly the rows the models are
# fitted to. The range is the 0.1% and 99.9% quantiles of the trait values (robust to a
# few implausible extreme records), and the fold variation is their ratio. Values are on
# the plotted scale: exp(log_value) (the geometric cell mean), except leaf_delta13C, whose
# log_value is already ci/ca.
build_trait_range_table <- function(data) {
  data |>
    ungroup() |>
    mutate(value = if_else(trait_name == "leaf_delta13C", log_value, exp(log_value))) |>
    group_by(trait_name) |>
    summarise(n = n(),
              q_lo = quantile(value, 0.001),
              q_hi = quantile(value, 0.999),
              fold = q_hi / q_lo,
              .groups = "drop")
}

# Word table of the trait ranges from build_trait_range_table().
export_trait_range_table <- function(trait_range_table, path = "outputs/tables/trait_range_table.docx") {
  source("R/export_figure4_table.R")  # FIGURE_TRAIT_ORDER, TRAIT_TABLE_NAMES

  # Units after the transformations in Analysis.qmd (huber_value x 1e2, delta13C -> ci/ca)
  trait_units <- c(
    leaf_delta13C      = "unitless",
    leaf_N_per_area    = "g m⁻²",
    leaf_mass_per_area = "g m⁻²",
    leaf_area          = "mm²",
    huber_value        = "mm² cm⁻²",
    wood_density       = "g cm⁻³",
    seed_dry_mass      = "mg"
  )
  fmt3 <- function(x) vapply(x, \(v) format(signif(v, 3), big.mark = ",", scientific = FALSE,
                                             drop0trailing = TRUE), character(1))

  traits <- intersect(FIGURE_TRAIT_ORDER, trait_range_table$trait_name)

  table_data <- trait_range_table |>
    filter(trait_name %in% traits) |>
    arrange(match(trait_name, traits)) |>
    transmute(
      Trait = TRAIT_TABLE_NAMES[trait_name],
      Unit  = trait_units[trait_name],
      q_lo  = fmt3(q_lo),
      q_hi  = fmt3(q_hi),
      fold  = fmt3(fold)
    )

  ft <- flextable::flextable(table_data) |>
    flextable::set_header_labels(Trait = "Trait", Unit = "Unit",
                                 q_lo = "0.1% quantile", q_hi = "99.9% quantile",
                                 fold = "Fold variation") |>
    flextable::add_footer_lines(
      "Values are the 0.1% and 99.9% quantiles of the species mean trait values per climate cell plotted in Figure 3; fold variation is the 99.9% quantile divided by the 0.1% quantile."
    ) |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = 2:5, part = "body", align = "center") |>
    flextable::bold(part = "header") |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::autofit()

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}
