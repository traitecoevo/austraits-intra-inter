source("R/find_ITV.r")
source("R/prepare_model_data.R")

# Fraction of intra-specific variation (ITV_frac) per trait, for combined,
# woody-only and non-woody-only subsets. Replaces three separate per-group loops
# over find_ITV() with one table of results: one row per trait x group (the
# non-woody row is omitted for WOODY_ONLY_TRAITS).
#
# Also reports, for each trait x group, how the climate space was sampled:
#   n_taxa                          — taxa with this trait
#   median_cells_per_taxon          — median number of climate cells per taxon
#   Median range in MAP (mm) / MAT (°C) — median, across taxa, of each taxon's
#                                     climate range (max - min) over its sites
#   centroid_range_MAP/_MAT         — range (max - min) spanned by the taxa's
#                                     climate centroids
# MAP is in mm (centroids back-transformed from the log scale, i.e. geometric
# means) and MAT in degrees C. Within-taxon ranges use the sites where this trait
# was measured; centroids are each taxon's mean climate across all traits, as
# used in the models.
summarise_ITV_frac <- function(data) {
  traits <- unique(data$trait_name)
  groups <- list(combined = NULL, woody = "woody", `non-woody` = "non-woody")

  tidyr::expand_grid(trait_name = traits, group = names(groups)) |>
    filter(!(group == "non-woody" & trait_name %in% WOODY_ONLY_TRAITS)) |>
    purrr::pmap_dfr(function(trait_name, group) {
      trait <- trait_name
      woody_filter <- groups[[group]]

      data_model <- prepare_model_data(data, trait, woody_filter)
      capture.output(result <- find_ITV(data_model))

      data_group <- data |>
        ungroup() |>
        filter(trait_name == trait)
      if (!is.null(woody_filter)) data_group <- filter(data_group, woody == woody_filter)

      taxon_summary <- data_group |>
        group_by(taxon_name) |>
        summarise(
          n_cells       = n_distinct(cell),
          range_MAP     = diff(range(annual_precipitation, na.rm = TRUE)),
          range_MAT     = diff(range(average_temperature, na.rm = TRUE)),
          centroid_MAP  = exp(first(log_annual_precipitation_mn)),
          centroid_MAT  = first(average_temperature_mn),
          .groups = "drop"
        )

      taxon_summary |>
        summarise(
          trait_name                    = trait,
          group                         = group,
          ITV_frac                      = result$ITV_frac,
          n_taxa                        = n(),
          median_cells_per_taxon        = median(n_cells),
          `Median range in MAP (mm)`    = median(range_MAP),
          `Median range in MAT (°C)`    = median(range_MAT),
          centroid_range_MAP            = diff(range(centroid_MAP, na.rm = TRUE)),
          centroid_range_MAT            = diff(range(centroid_MAT, na.rm = TRUE))
        )
    })
}

# Writes the summarise_ITV_frac() output as a single Word table (flextable),
# styled like the Figure 3 / Figure 4 fixed-effects tables: one block per trait
# (in Figure 3's trait order) with the trait name as a bold merged header row,
# one row per growth-form group (Combined / Woody / Non-woody), and spanning
# headers over the within-taxon and centroid climate ranges.
export_ITV_frac_table <- function(itv_summary, path = "outputs/tables/ITV_frac_table.docx") {
  source("R/export_figure4_table.R")  # FIGURE_TRAIT_ORDER, TRAIT_TABLE_NAMES

  fmt <- function(x, digits) formatC(x, format = "f", digits = digits, big.mark = ",")
  group_labels <- c(combined = "Combined", woody = "Woody", `non-woody` = "Non-woody")
  traits <- intersect(FIGURE_TRAIT_ORDER, itv_summary$trait_name)

  table_data <- itv_summary |>
    filter(trait_name %in% traits) |>
    arrange(match(trait_name, traits), match(group, names(group_labels))) |>
    transmute(
      Trait        = factor(TRAIT_TABLE_NAMES[trait_name], levels = TRAIT_TABLE_NAMES[traits]),
      Group        = group_labels[group],
      itv_frac     = fmt(ITV_frac, 2),
      n_taxa       = fmt(n_taxa, 0),
      cells        = fmt(median_cells_per_taxon, 1),
      within_MAP   = fmt(`Median range in MAP (mm)`, 0),
      within_MAT   = fmt(`Median range in MAT (°C)`, 1),
      centroid_MAP = fmt(centroid_range_MAP, 0),
      centroid_MAT = fmt(centroid_range_MAT, 1)
    )

  # Trait names become merged group rows above each trait's growth-form rows
  ft <- flextable::as_grouped_data(table_data, groups = "Trait") |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::set_header_labels(
      Group = "Group", itv_frac = "ITV fraction",
      n_taxa = "Taxa", cells = "Median cells per taxon",
      within_MAP = "MAP (mm)", within_MAT = "MAT (°C)",
      centroid_MAP = "MAP (mm)", centroid_MAT = "MAT (°C)"
    ) |>
    flextable::add_header_row(
      values = c("", "Median range in", "Centroid range"),
      colwidths = c(4, 2, 2)
    ) |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = 2:8, part = "body", align = "center") |>
    flextable::bold(part = "header") |>
    flextable::bold(i = ~ !is.na(Trait), part = "body") |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::autofit()

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}
