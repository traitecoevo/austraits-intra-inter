source("R/prepare_model_data.R")

# Trait display order (top to bottom) shared by Figure 3 and the woody vs
# non-woody comparison figure, so both stay in the same order and can't
# silently diverge. seed_dry_mass is last because it's the only trait with a
# populated x-axis label, i.e. the intended bottom row.
FIGURE_TRAIT_ORDER <- c(
  "leaf_delta13C", "leaf_N_per_area", "leaf_mass_per_area",
  "leaf_area", "huber_value", "wood_density", "seed_dry_mass"
)

# Leaf traits fill the left half (MAP | MAT columns) of the compressed figures; the
# remaining traits fill the right half (MAP | MAT) beneath the legend.
FIGURE_LEAF_TRAITS <- c("leaf_delta13C", "leaf_N_per_area", "leaf_mass_per_area", "leaf_area")

FIGURE_ENV_LABELS <- list(
  log_annual_precipitation = "MAP (mm)",
  average_temperature      = expression(paste(MAT, " (", degree * C, ")"))
)

# Scatter panels for one climate driver, named by trait, in FIGURE_TRAIT_ORDER.
# woody = TRUE plots the woody vs non-woody comparison
# (fit_model2d_frequentist_4panel_woody_compressed()) and drops
# WOODY_ONLY_TRAITS (huber_value, wood_density) entirely, since they have no
# non-woody data to compare against. woody = FALSE plots every trait with the
# combined (non-split) fit. remove_ns = TRUE drops the OLS overall,
# interspecific (β_mn), and average intraspecific (β_c) lines wherever they
# are not (marginally) significant at P < 0.10 — see fit_model2d_frequentist_
# 4panel_compressed()/_woody_compressed() for exactly what that filters.
build_compressed_trait_panels <- function(data, ind_variable, woody = FALSE, remove_ns = FALSE) {
  traits <- FIGURE_TRAIT_ORDER
  if (woody) traits <- setdiff(traits, WOODY_ONLY_TRAITS)

  panels <- list()
  for (trait in traits) {
    cat(trait, ":  ")
    data_model <- data |> filter(trait_name == trait)

    if (woody) {
      out <- fit_model2d_frequentist_4panel_woody_compressed(data_model, ind_variable = ind_variable, remove_ns = remove_ns)
    } else {
      out <- fit_model2d_frequentist_4panel_compressed(data_model, ind_variable = ind_variable, remove_ns = remove_ns)
    }
    # out[[2]] is scatter | density, where the combined figure's scatter carries its
    # significance key stacked on top (the woody figure has none); [[1]] keeps the scatter
    panels[[trait]] <- out[[2]][[1]]
  }
  panels
}

# Line legend for the top-right cell of the compressed figures, matching the
# fitted-line styling in the panels: a black outline under a coloured line, solid for
# the within-species response and dashed for the among-species response. Drawn
# manually (rather than via a ggplot guide) so the black outline stays solid under the
# dashed line, exactly as in the panels. Keys are stacked one per row.
#   woody = FALSE: blue solid (within), orange dashed (among), grey solid (overall
#                  OLS response) (fit_model2d_frequentist_4panel_compressed()).
#   woody = TRUE:  colour is growth form (non-woody green, woody vermillion, as in
#                  fit_model2d_frequentist_4panel_woody_compressed()), so the
#                  within/among keys are black and differ by linetype only, followed
#                  by two rows naming the colours.
build_line_legend <- function(woody = FALSE) {
  labels <- c("Within-species response to the env.",
              "Among-species response to the env.")

  if (!woody) {
    keys <- tibble(colour = c("#0072B2", "#E69F00", "grey60"),
                   type   = c("solid", "dashed", "solid"),
                   label  = c(labels, "Overall response to the env."))
  } else {
    keys <- tibble(colour = c("black", "black", "#009E73", "#D55E00"),
                   type   = c("solid", "dashed", "solid", "solid"),
                   label  = c(labels, "Non-woody", "Woody"))
  }

  key_x <- 0.12
  key_length <- 0.12
  keys <- keys |>
    mutate(y = rev(seq_len(n())), x = key_x, xend = key_x + key_length)

  ggplot(keys) +
    # Black outline only under coloured keys (it would fill the gaps of a black dashed key)
    geom_segment(data = \(d) filter(d, colour != "black"),
                 aes(x = x, xend = xend, y = y, yend = y), colour = "black", linewidth = 1.7) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = y, colour = colour, linetype = type), linewidth = 1.1) +
    geom_text(aes(x = xend + 0.03, y = y, label = label), hjust = 0, size = 17.5 / .pt) +
    scale_colour_identity() +
    scale_linetype_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.25, nrow(keys) + 0.75), expand = FALSE, clip = "off") +
    theme_void()
}

# Compressed trait figure on a four-column grid: MAP | MAT | MAP | MAT. The leaf traits
# (FIGURE_LEAF_TRAITS) run down the left two columns; the line legend sits in the top
# right, with the remaining traits below it in the right two columns (Figure 3:
# SA:LA, wood density, seed mass; Figure 4: seed mass only, as huber_value and
# wood_density are woody-only). MAT panels drop the trait y-axis title (it is on the
# MAP panel alongside), and the climate x-axis title is shown only on the bottom
# panel of each half. `left_traits` sets which traits go in the left half (the rest go
# in the right half under the legend, in FIGURE_TRAIT_ORDER), e.g. Figure 4 moves
# leaf_area to the right. Saved to `path`; see build_compressed_trait_panels() for
# `woody` and `remove_ns`.
build_compressed_trait_figure <- function(data, path, woody = FALSE, remove_ns = FALSE,
                                          left_traits = FIGURE_LEAF_TRAITS,
                                          height = 11, width = 13.33 * (4 / 3)) {
  env_vars <- names(FIGURE_ENV_LABELS)
  panels <- lapply(set_names(env_vars), \(v) build_compressed_trait_panels(data, v, woody = woody, remove_ns = remove_ns))

  traits <- names(panels[[1]])
  left   <- intersect(left_traits, traits)
  right  <- setdiff(traits, left_traits)
  n_rows <- max(length(left), length(right) + 1)

  # One panel, with its x title only when it is the bottom panel of its half.
  # `&` reaches the scatter inside the combined figure's key + scatter stack (the key
  # is theme_void, so the labels don't show there).
  panel <- function(trait, v, bottom) {
    p <- panels[[v]][[trait]] & labs(x = if (bottom) FIGURE_ENV_LABELS[[v]] else NULL)
    if (v == "average_temperature") p <- p & labs(y = NULL)
    p
  }
  pair <- function(trait, bottom) lapply(env_vars, \(v) panel(trait, v, bottom))

  # Cells in reading order; "#" marks empty cells in the design
  plots <- list()
  design <- character(n_rows)
  next_id <- 1
  add <- function(p) {
    plots[[length(plots) + 1]] <<- p
    id <- LETTERS[next_id]
    next_id <<- next_id + 1
    id
  }
  for (i in seq_len(n_rows)) {
    row <- character(0)
    # Left half: leaf traits
    if (i <= length(left)) {
      ids <- vapply(pair(left[i], i == length(left)), add, character(1))
      row <- c(row, ids)
    } else row <- c(row, "#", "#")
    # Right half: legend on the first row, then the remaining traits
    if (i == 1) {
      id <- add(build_line_legend(woody))
      row <- c(row, id, id)
    } else if (i - 1 <= length(right)) {
      ids <- vapply(pair(right[i - 1], i - 1 == length(right)), add, character(1))
      row <- c(row, ids)
    } else row <- c(row, "#", "#")
    design[i] <- paste(row, collapse = "")
  }

  combined <- patchwork::wrap_plots(plots, design = paste(design, collapse = "\n"))
  ggsave(path, plot = combined, height = height, width = width)
  combined
}
