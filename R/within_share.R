source("R/cache_is_stale.R")

# Within-species share of climate-associated trait variation, S_w (Appendix A):
#
#   S_w = V_w / (V_w + V_m) = k * r / (k * r + 1)
#   r   = (B_w / B_m)^2                 squared slope ratio, per posterior draw
#   k   = var(Env_w) / var(Env_m)       sampling term, at its observed value (k_obs)
#
# with Env_m = std_<var>_mn (taxon climate centroid) and Env_w = std_<var>_c
# (deviation from it). k is computed from the rows each model was fitted to (the
# data saved with the fit), separately for each growth form. S_w is bounded in
# [0, 1] (0.5 = equal contributions) and, unlike the unbounded variance-weighted
# ratio it replaces, stays stable when B_m is near zero.
#
# Groups follow Figure 5: for the woody-split traits, non-woody and woody slopes come
# from the woody interaction model (main effect, and main effect + :woodywoody) and
# the combined value from the all-species model; huber_value and wood_density are
# woody-only and use the all-species model.

WITHIN_SHARE_ENV_VARS <- c(MAP = "log_annual_precipitation", MAT = "average_temperature")
WITHIN_SHARE_WOODY_SPLIT_TRAITS <- c("leaf_mass_per_area", "leaf_N_per_area", "leaf_delta13C",
                                     "leaf_area", "seed_dry_mass")

# Trait order (top to bottom) for the S_w figures and table: the woody-split traits
# first (seed mass included), then the woody-only traits. Unlike FIGURE_TRAIT_ORDER
# (Figure 3), seed mass doesn't need to be last here, as no panel carries its own x-axis.
WITHIN_SHARE_TRAIT_ORDER <- c("leaf_delta13C", "leaf_N_per_area", "leaf_mass_per_area",
                              "leaf_area", "seed_dry_mass", "huber_value", "wood_density")

share_within <- function(r, k) k * r / (k * r + 1)

# S_w draws and diagnostics for one fit / climate variable / group.
#   group       — "combined", "woody" or "non-woody"
#   woody_shift — add the :woodywoody interaction terms (woody slopes from the woody model)
#   data_group  — the rows of the fitted data belonging to this group (for k_obs)
within_share_one <- function(fit, env, group, data_group, woody_shift = FALSE) {
  var  <- WITHIN_SHARE_ENV_VARS[[env]]
  mn   <- paste0("std_", var, "_mn")
  cc   <- paste0("std_", var, "_c")
  dr   <- posterior::as_draws_df(fit) |> as_tibble()

  B_m <- dr[[paste0("b_", mn)]]
  B_w <- dr[[paste0("b_", cc)]]
  m_hyp <- mn
  w_hyp <- cc
  if (woody_shift) {
    B_m <- B_m + dr[[paste0("b_", mn, ":woodywoody")]]
    B_w <- B_w + dr[[paste0("b_", cc, ":woodywoody")]]
    m_hyp <- paste0(mn, " + ", mn, ":woodywoody")
    w_hyp <- paste0(cc, " + ", cc, ":woodywoody")
  }

  k_obs <- var(data_group[[cc]], na.rm = TRUE) / var(data_group[[mn]], na.rm = TRUE)

  # As in Figure 5 (B_r): unfilled when BOTH slopes' 90% intervals include zero
  ns <- function(h) {
    h <- brms::hypothesis(fit, paste0(h, " = 0"), alpha = 0.1)$hypothesis
    h$CI.Lower <= 0 && h$CI.Upper >= 0
  }

  list(
    # r is kept so S_w can be recomputed at any k (within_share_at_k())
    draws = tibble(env = env, group = group, r = (B_w / B_m)^2,
                   value = share_within((B_w / B_m)^2, k_obs)),
    diagnostics = tibble(
      env = env, group = group,
      k_obs       = k_obs,
      cor_m_w     = cor(data_group[[mn]], data_group[[cc]], use = "complete.obs"),
      prop_single = mean(table(data_group$taxon_name) == 1),
      no_fill     = ns(m_hyp) && ns(w_hyp)
    )
  )
}

# All S_w draws and diagnostics for one trait (both climate variables, all groups),
# cached in outputs/model_fits/within_share/ and recomputed whenever a source fit is
# newer than the cache.
within_share_trait <- function(trait) {
  combined_path <- paste0("outputs/model_fits/fit2d_bayesian/", trait, ".RDS")
  woody_path    <- paste0("outputs/model_fits/fit2d_bayesian_woody/", trait, ".RDS")
  is_split      <- trait %in% WITHIN_SHARE_WOODY_SPLIT_TRAITS
  cache_path    <- paste0("outputs/model_fits/within_share/", trait, ".RDS")

  stale <- cache_is_stale(cache_path, combined_path) ||
    (is_split && cache_is_stale(cache_path, woody_path)) ||
    !("r" %in% names(readRDS(cache_path)$draws))  # caches from before r was stored
  if (!stale) return(readRDS(cache_path))

  # Saved as list(fit_MAP, fit_MAT, data_model)
  fits_combined <- readRDS(combined_path)
  parts <- list()
  for (i in 1:2) {
    env <- names(WITHIN_SHARE_ENV_VARS)[i]
    if (is_split) {
      fits_woody <- if (i == 1) readRDS(woody_path) else fits_woody
      d_w <- fits_woody[[3]]
      parts <- c(parts, list(
        within_share_one(fits_woody[[i]],    env, "non-woody", filter(d_w, woody == "non-woody")),
        within_share_one(fits_woody[[i]],    env, "woody",     filter(d_w, woody == "woody"), woody_shift = TRUE),
        within_share_one(fits_combined[[i]], env, "combined",  fits_combined[[3]])
      ))
    } else {
      # Woody-only traits: the all-species model is the woody estimate
      parts <- c(parts, list(within_share_one(fits_combined[[i]], env, "woody", fits_combined[[3]])))
    }
  }

  result <- list(
    draws       = map_dfr(parts, "draws") |> mutate(trait_name = trait, .before = 1),
    diagnostics = map_dfr(parts, "diagnostics") |> mutate(trait_name = trait, .before = 1)
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, cache_path)
  rm(fits_combined); if (is_split) rm(fits_woody); gc(verbose = FALSE)
  result
}

# S_w draws for the given traits at a chosen k: "observed" (each group's k_obs) or a
# single number applied to every group (e.g. 0.5, or 1 for the ceiling). Adds factor
# columns `trait` and `group_label` for plotting, and each group's no_fill flag
# (TRUE when neither slope differs from zero — see within_share_one()).
within_share_at_k <- function(traits, k = "observed") {
  source("R/export_figure4_table.R")  # TRAIT_TABLE_NAMES
  traits <- intersect(WITHIN_SHARE_TRAIT_ORDER, traits)

  map_dfr(traits, function(t) {
    ws <- within_share_trait(t)
    ws$draws |>
      left_join(select(ws$diagnostics, env, group, k_obs, no_fill), by = c("env", "group")) |>
      mutate(k = if (identical(k, "observed")) k_obs else k,
             value = share_within(r, k))
  }) |>
    mutate(trait = factor(TRAIT_TABLE_NAMES[trait_name], levels = rev(TRAIT_TABLE_NAMES[traits])),
           group_label = factor(group, levels = c("non-woody", "woody", "combined"),
                                labels = c("Non-woody", "Woody", "Combined")))
}

# Logit axis shared by the S_w figures, from p_min (smallest lower 95% bound at the
# observed k — S_w only rises with k — rounded down to a power of ten) to 0.99999,
# which holds every 95% interval up to k = 1 (the largest upper bound is ~0.99995).
# Breaks at every decade, with every break above 0.5 mirrored below it (0.9/0.1,
# 0.99/0.01, ..., 0.99999/1e-5). Values from 0.001 up are written as decimals; smaller
# ones as powers of ten. Anything beyond the limits runs to the edge (squished).
# p_min (a power of ten) and p_max can be set directly to zoom either end, e.g. 1e-4
# for the curve plot; they are then also labelled, and anything beyond runs to the edge.
within_share_logit_scale <- function(traits, scale_fn = scale_x_continuous, p_min = NULL,
                                     p_max = 0.99999) {
  lo_obs <- within_share_at_k(traits, "observed") |>
    group_by(trait_name, env, group) |>
    summarise(lo = quantile(value, 0.025), .groups = "drop")
  fixed_min <- !is.null(p_min)
  if (!fixed_min) p_min <- 10^floor(log10(min(lo_obs$lo)))
  decade <- seq(log10(p_min) + if (fixed_min) 0 else 1, -1)
  upper  <- c(0.9, 0.99, 0.999, 0.9999, 0.99999)
  upper  <- upper[upper <= p_max + 1e-12]
  breaks <- c(10^decade, 0.5, upper)
  labels <- as.expression(c(
    lapply(decade, \(e) if (e >= -3) format(10^e, scientific = FALSE) else bquote(10^.(e))),
    list("0.5"), as.list(format(upper, scientific = FALSE, drop0trailing = TRUE))))

  logit_tr <- scales::trans_new("logit_clamped",
                                \(x) qlogis(pmin(pmax(x, p_min), p_max)), plogis)
  scale_fn(transform = logit_tr, breaks = breaks, labels = labels,
           limits = c(p_min, p_max), oob = scales::oob_squish,
           # small padding (in logit units) so the end labels aren't clipped
           expand = expansion(add = c(0.3, 0.6)))
}

WITHIN_SHARE_COLOURS <- c("Combined" = "#0072B2", "Woody" = "#D55E00", "Non-woody" = "#009E73")

within_share_theme <- function() {
  theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = 14),
          text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(1.2, "lines"))
}

# Figure 5 (alternative): forest plot of S_w for all traits — one panel per climate
# variable, traits as rows, growth forms dodged within each row (colours as in
# Figure 5). Points are posterior medians, thick bars 50% and thin lines 95%
# (2.5-97.5%) intervals. As in Figure 5 (B_r), points are filled white when neither
# slope (B_m nor B_w) differs from zero (both 90% intervals include zero).
# k is "observed" (each group's k_obs, the floor), or a single
# value for every group (e.g. 0.5 intermediate, 1 the ceiling). All versions share the
# same logit axis (within_share_logit_scale()) so they can be compared directly.
# Returns a patchwork: the forest plot with a direction arrow strip under each panel.
#
# metric = "B_r" draws Figure 5 (the manuscript figure) in the same format instead:
# B_r = |B_w| / (|B_w| + |B_m|) = sqrt(r) / (1 + sqrt(r)) from the same posterior draws,
# on a linear 0-1 axis (k does not enter B_r).
plot_within_share_forest <- function(traits, k = "observed", metric = c("S_w", "B_r")) {
  metric <- match.arg(metric)
  draws <- within_share_at_k(traits, k)
  if (metric == "B_r") draws <- mutate(draws, value = sqrt(r) / (1 + sqrt(r)))

  summary <- draws |>
    group_by(trait, env, group_label, no_fill) |>
    summarise(median = median(value),
              lo   = quantile(value, 0.025), hi   = quantile(value, 0.975),
              lo50 = quantile(value, 0.25),  hi50 = quantile(value, 0.75), .groups = "drop") |>
    mutate(point_fill = if_else(no_fill, "white", WITHIN_SHARE_COLOURS[as.character(group_label)]))

  k_text <- if (identical(k, "observed")) "observed k" else paste("k =", k)
  dodge  <- position_dodge(width = 0.7)
  if (metric == "S_w") {
    # axis from 1e-5 to 0.9999: the few 95% intervals reaching beyond run to the panel edge
    x_scale <- within_share_logit_scale(traits, p_min = 1e-5, p_max = 0.9999)
    lims    <- c(1e-5, 0.9999)
    x_label <- bquote(S[w] ~ "(" * .(k_text) * ")")
    x_text  <- element_text(size = 9.5, angle = 45, hjust = 1)
  } else {
    x_scale <- scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                                  expand = expansion(add = 0.02))
    lims    <- c(0, 1)
    x_label <- bquote(B[r])
    x_text  <- element_text(size = 10)
  }

  forest <- ggplot(summary, aes(y = trait, x = median, colour = group_label)) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6, position = dodge) +
    geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 2, position = dodge) +
    # group = growth form so points dodge with their own intervals (not by fill colour)
    geom_point(aes(fill = point_fill, group = group_label), size = 2.2, shape = 21,
               colour = "black", position = dodge) +
    facet_grid(~env) +
    x_scale +
    scale_colour_manual(values = WITHIN_SHARE_COLOURS, breaks = c("Combined", "Woody", "Non-woody"), name = NULL) +
    scale_fill_identity() +
    within_share_theme() +
    # S_w: angled so the closely spaced decade labels at both ends don't overlap
    theme(axis.text.x = x_text) +
    # the x-axis title goes under the direction arrow instead (below)
    labs(y = NULL, x = NULL)

  # Direction arrow under each panel: towards 0 = among-species (turnover) response
  # stronger, towards 1 = within-species response stronger. Drawn as its
  # own strip with the same facets and x scale, so patchwork aligns each arrow with its
  # panel; a short tick marks S_w = 0.5 (equal contributions).
  arrow_data <- tibble(env = unique(summary$env))
  arrow_strip <- ggplot(arrow_data) +
    annotate("segment", x = lims[1], xend = lims[2], y = 0, yend = 0,
             arrow = arrow(ends = "both", length = unit(0.1, "inches"), type = "closed"),
             linewidth = 0.5) +
    annotate("segment", x = 0.5, xend = 0.5, y = -0.35, yend = 0.35, linewidth = 0.5) +
    annotate("text", x = lims[1], y = 0.9, label = "Among-species stronger", hjust = 0, size = 3.6) +
    annotate("text", x = lims[2], y = 0.9, label = "Within-species stronger", hjust = 1, size = 3.6) +
    facet_grid(~env) +
    x_scale +
    scale_y_continuous(limits = c(-0.5, 1.4)) +
    theme_void() +
    theme(strip.text = element_blank(),
          panel.spacing.x = unit(1.2, "lines"),
          axis.title.x = element_text(size = 12, margin = margin(t = 6))) +
    labs(x = x_label)

  patchwork::wrap_plots(forest, arrow_strip, ncol = 1, heights = c(1, 0.09)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# Curve plot (Appendix A): S_w as k rises from each group's observed value (the floor)
# to 1 (the ceiling), on a dense grid of n_k values of k. Lines are posterior medians,
# ribbons 95% intervals; one panel per trait (rows) and climate variable (columns),
# with the same logit S_w axis as the forest plots. Ribbons are 95% credible intervals
# of S_w at each fixed k (posterior uncertainty in the slopes; k itself is the value
# being varied). Points mark n_points evenly spaced k values from the observed k
# (circles) to k = 1 (squares), with the intermediate ones as small dots, plus k_mid
# (triangles), the intermediate value used in the forest-plot version.
plot_within_share_curve <- function(traits, n_k = 100, k_mid = 0.5, n_points = 8) {
  # Facet rows top-to-bottom in Figure 5 order (the forest plots list traits bottom-up)
  draws <- within_share_at_k(traits, "observed") |>
    mutate(trait = factor(trait, levels = rev(levels(trait))))

  curve <- draws |>
    group_by(trait, env, group_label) |>
    group_modify(function(d, key) {
      k_grid <- seq(d$k_obs[1], 1, length.out = n_k)
      s <- outer(d$r, k_grid, share_within)  # draws x k
      tibble(k = k_grid,
             median = apply(s, 2, median),
             lo = apply(s, 2, quantile, 0.025),
             hi = apply(s, 2, quantile, 0.975))
    }) |>
    ungroup()

  point_types <- c("Observed k", "Evenly spaced k", paste("k =", k_mid), "k = 1")
  marks <- draws |>
    group_by(trait, env, group_label) |>
    group_modify(function(d, key) {
      ks <- seq(d$k_obs[1], 1, length.out = n_points)
      which <- c("Observed k", rep("Evenly spaced k", n_points - 2), "k = 1")
      tibble(k = c(ks, k_mid), which = c(which, paste("k =", k_mid)),
             median = sapply(c(ks, k_mid), \(kk) median(share_within(d$r, kk))))
    }) |>
    ungroup() |>
    mutate(which = factor(which, levels = point_types))

  ggplot(curve, aes(x = k, y = median, colour = group_label, fill = group_label)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_vline(xintercept = k_mid, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(data = marks, aes(shape = which, size = which), colour = "black") +
    facet_grid(trait ~ env, switch = "y") +
    # lower end at 1e-4: the curves themselves all sit above it; only the lower edges
    # of some ribbons run below and are cut at the panel edge
    within_share_logit_scale(traits, scale_fn = scale_y_continuous, p_min = 1e-4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = WITHIN_SHARE_COLOURS, breaks = c("Combined", "Woody", "Non-woody"), name = NULL) +
    scale_fill_manual(values = WITHIN_SHARE_COLOURS, breaks = c("Combined", "Woody", "Non-woody"), name = NULL) +
    # Values keyed by name; "Evenly spaced k" dots are left out of the legend (breaks)
    scale_shape_manual(values = setNames(c(21, 21, 24, 22), point_types), name = NULL,
                       breaks = point_types[-2]) +
    scale_size_manual(values = setNames(c(2, 1.2, 2, 2), point_types), name = NULL,
                      breaks = point_types[-2]) +
    within_share_theme() +
    theme(strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1),
          strip.placement = "outside",
          legend.box = "vertical",  # colour and point-type legends on separate rows
          axis.text.y = element_text(size = 8)) +
    labs(x = "k = var(Env_w) / var(Env_m)", y = bquote(S[w]))
}

# ---- Supplementary tables -------------------------------------------------------
# All in the same layout as the other tables (flextable): one block per trait (bold
# merged header row) with Combined / Woody / Non-woody rows and MAP / MAT spanning
# headers. Values are posterior medians [2.5%, 97.5%] to 4 significant figures, from the
# cached draws in within_share_trait() (the same draws, groups and slopes as the figures).

WITHIN_SHARE_GROUP_LABELS <- c(combined = "Combined", woody = "Woody", `non-woody` = "Non-woody")

# 4 significant figures, trailing zeros kept, no scientific notation
sig4 <- function(x) formatC(signif(x, 4), digits = 4, format = "fg", flag = "#")
sig4_interval <- function(x) {
  sprintf("%s [%s, %s]", sig4(median(x)), sig4(quantile(x, 0.025)), sig4(quantile(x, 0.975)))
}

# rows: one row per trait_name x env x group, with the value columns in `cols`
# (named vector: column -> header label). Pivots MAP / MAT side by side and builds the
# flextable; footer is added under the table.
build_within_share_flextable <- function(rows, traits, cols, footer) {
  source("R/export_figure4_table.R")  # TRAIT_TABLE_NAMES

  table_data <- rows |>
    pivot_wider(names_from = env, values_from = all_of(names(cols)),
                names_glue = "{env}__{.value}") |>
    arrange(match(trait_name, traits), match(group, names(WITHIN_SHARE_GROUP_LABELS))) |>
    mutate(Trait = factor(TRAIT_TABLE_NAMES[trait_name], levels = TRAIT_TABLE_NAMES[traits]),
           Group = WITHIN_SHARE_GROUP_LABELS[group]) |>
    select(Trait, Group, all_of(c(paste0("MAP__", names(cols)), paste0("MAT__", names(cols)))))

  n <- length(cols)
  header_labels <- as.list(c(Group = "Group",
                             setNames(rep(cols, 2), c(paste0("MAP__", names(cols)), paste0("MAT__", names(cols))))))

  flextable::as_grouped_data(table_data, groups = "Trait") |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::set_header_labels(values = header_labels) |>
    flextable::add_header_row(values = c("", "MAP", "MAT"), colwidths = c(1, n, n)) |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = 2:(1 + 2 * n), part = "body", align = "center") |>
    flextable::bold(part = "header") |>
    flextable::bold(i = ~ !is.na(Trait), part = "body") |>
    flextable::add_footer_lines(footer) |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::autofit()
}

# B_r = |B_w| / (|B_w| + |B_m|) (Figure 5, the manuscript figure). Computed from the
# cached slope-ratio draws: with r = (B_w / B_m)^2, B_r = sqrt(r) / (1 + sqrt(r)).
export_br_table <- function(traits, path = "outputs/tables/Br_table.docx") {
  traits <- intersect(WITHIN_SHARE_TRAIT_ORDER, traits)
  rows <- map_dfr(traits, function(trait) {
    ws <- within_share_trait(trait)
    ws$draws |>
      mutate(B_r = sqrt(r) / (1 + sqrt(r))) |>
      group_by(env, group) |>
      summarise(B_r = sig4_interval(B_r), .groups = "drop") |>
      mutate(trait_name = trait, .before = 1)
  })
  ft <- build_within_share_flextable(rows, traits, c(B_r = "B_r"), paste(
    "B_r = |B_w| / (|B_w| + |B_m|): within-species slope relative to the sum of the",
    "within- and across-species slopes (0.5 = equal); median [2.5%, 97.5%] of the posterior.",
    "Woody slopes include the woody interaction terms; SA:LA and wood density (woody-only)",
    "use the all-species model."))
  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}

# S_w at a given k: "observed" (each group's k_obs; also reports k and the Env_m-Env_w
# correlation, the orthogonality check) or a single value for every group (e.g. 0.5, 1).
export_within_share_table <- function(traits, k = "observed",
                                      path = "outputs/tables/within_share_table.docx") {
  traits <- intersect(WITHIN_SHARE_TRAIT_ORDER, traits)
  observed <- identical(k, "observed")

  rows <- map_dfr(traits, function(trait) {
    ws <- within_share_trait(trait)
    ws$draws |>
      left_join(ws$diagnostics, by = c("env", "group")) |>
      mutate(S_w = share_within(r, if (observed) k_obs else k)) |>
      group_by(env, group) |>
      summarise(k = sig4(first(k_obs)), cor = sig4(first(cor_m_w)),
                S_w = sig4_interval(S_w), .groups = "drop") |>
      mutate(trait_name = trait, .before = 1)
  })

  if (observed) {
    cols <- c(k = "k (observed)", cor = "cor(Env_m, Env_w)", S_w = "S_w")
    k_text <- "the observed k = var(Env_w) / var(Env_m)"
  } else {
    rows <- select(rows, -k, -cor)
    cols <- c(S_w = "S_w")
    k_text <- paste0("k = ", k, " for every group")
  }
  ft <- build_within_share_flextable(rows, traits, cols, paste(
    "S_w: within-species share of climate-associated trait variation, median [2.5%, 97.5%]",
    "of the posterior at", paste0(k_text, ".")))
  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}

# Table for the S_w-vs-k curve (sensitivity) figure: S_w at the n_points evenly spaced
# values of k from each group's observed k (point 1) to k = 1 (point n_points), as marked
# on the curves. One row per group x climate variable x point, grouped by trait.
export_within_share_curve_table <- function(traits, n_points = 8,
                                            path = "outputs/tables/within_share_curve_table.docx") {
  source("R/export_figure4_table.R")  # TRAIT_TABLE_NAMES
  traits <- intersect(WITHIN_SHARE_TRAIT_ORDER, traits)

  table_data <- map_dfr(traits, function(trait) {
    ws <- within_share_trait(trait)
    ws$draws |>
      left_join(select(ws$diagnostics, env, group, k_obs), by = c("env", "group")) |>
      group_by(env, group) |>
      group_modify(function(d, key) {
        ks <- seq(d$k_obs[1], 1, length.out = n_points)
        tibble(point = seq_along(ks), k = sig4(ks),
               S_w = sapply(ks, \(kk) sig4_interval(share_within(d$r, kk))))
      }) |>
      ungroup() |>
      mutate(trait_name = trait, .before = 1)
  }) |>
    arrange(match(trait_name, traits), match(env, c("MAP", "MAT")),
            match(group, names(WITHIN_SHARE_GROUP_LABELS)), point) |>
    transmute(Trait = factor(TRAIT_TABLE_NAMES[trait_name], levels = TRAIT_TABLE_NAMES[traits]),
              Climate = env, Group = WITHIN_SHARE_GROUP_LABELS[group],
              Point = if_else(point == 1, "1 (observed k)",
                              if_else(point == n_points, paste0(n_points, " (k = 1)"), as.character(point))),
              k, S_w)

  ft <- flextable::as_grouped_data(table_data, groups = "Trait") |>
    flextable::as_flextable(hide_grouplabel = TRUE) |>
    flextable::set_header_labels(values = list(Climate = "Climate", Group = "Group",
                                               Point = "Point", k = "k", S_w = "S_w")) |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = 2:5, part = "body", align = "center") |>
    flextable::bold(part = "header") |>
    flextable::bold(i = ~ !is.na(Trait), part = "body") |>
    flextable::add_footer_lines(paste(
      "S_w (median [2.5%, 97.5%] of the posterior) at", n_points, "evenly spaced values of",
      "k = var(Env_w) / var(Env_m), from each group's observed k (point 1) to k = 1",
      paste0("(point ", n_points, "); the points marked on the S_w-vs-k curves."))) |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::autofit()

  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}

# ---- Orthogonality check (Appendix A) -----------------------------------------------
# S_w treats the turnover and within-species components as additive, dropping the
# cross term 2 * B_m * B_w * cov(Env_m, Env_w) from
#   var(fitted) = B_m^2 var(Env_m) + B_w^2 var(Env_w) + 2 B_m B_w cov(Env_m, Env_w).
# Env_w is a deviation from each taxon's centroid, but centroids were computed from all
# traits pooled (one record per grid cell), so within one trait's data Env_w need not
# average to zero within each taxon and the covariance need not be exactly zero.
#
# For each trait x climate variable x growth form (rows the models were fitted to:
# prepare_model_data(), split by woody for the growth-form groups) this reports:
#   cor          cor(Env_m, Env_w) across observations (the "|r| < X" in the Methods)
#   cross_share  |cross term| / (V_m + V_w), per posterior draw (median [2.5%, 97.5%]).
#                Dividing through by B_m^2: 2 sqrt(r) |cov| / (var_m + r var_w), with
#                r = (B_w / B_m)^2 from the cached draws; < ~0.05 means negligible.
#   taxon_offset mean over taxa of |taxon mean of Env_w|, relative to sd(Env_w): how far
#                each taxon's deviations are from averaging to zero (0 = perfectly centred)
within_share_orthogonality <- function(data, traits) {
  source("R/prepare_model_data.R")  # prepare_model_data(), the rows the models were fitted to
  traits <- intersect(WITHIN_SHARE_TRAIT_ORDER, traits)

  map_dfr(traits, function(trait) {
    ws <- within_share_trait(trait)
    d_all <- prepare_model_data(data, trait)
    groups <- if (trait %in% WITHIN_SHARE_WOODY_SPLIT_TRAITS) c("combined", "woody", "non-woody") else "woody"

    tidyr::expand_grid(env = names(WITHIN_SHARE_ENV_VARS), group = groups) |>
      purrr::pmap_dfr(function(env, group) {
        var <- WITHIN_SHARE_ENV_VARS[[env]]
        d <- if (group %in% c("woody", "non-woody") && trait %in% WITHIN_SHARE_WOODY_SPLIT_TRAITS)
          filter(d_all, woody == group) else d_all
        m <- d[[paste0("std_", var, "_mn")]]
        w <- d[[paste0("std_", var, "_c")]]
        r <- ws$draws$r[ws$draws$env == env & ws$draws$group == group]
        cross_share <- 2 * sqrt(r) * abs(cov(m, w)) / (var(m) + r * var(w))
        taxon_offset <- tibble(taxon = d$taxon_name, w = w) |>
          group_by(taxon) |> summarise(mean_w = mean(w), .groups = "drop")
        tibble(trait_name = trait, env = env, group = group, n = nrow(d),
               cor = cor(m, w),
               cross_share_median = median(cross_share),
               cross_share_lo = quantile(cross_share, 0.025),
               cross_share_hi = quantile(cross_share, 0.975),
               taxon_offset = mean(abs(taxon_offset$mean_w)) / sd(w))
      })
  })
}

# Word table of within_share_orthogonality() in the same layout as the other tables.
export_orthogonality_table <- function(ortho, path = "outputs/tables/orthogonality_table.docx") {
  traits <- intersect(WITHIN_SHARE_TRAIT_ORDER, unique(ortho$trait_name))
  rows <- ortho |>
    transmute(trait_name, env, group,
              cor = sig4(cor),
              cross = sprintf("%s [%s, %s]", sig4(cross_share_median), sig4(cross_share_lo), sig4(cross_share_hi)),
              offset = sig4(taxon_offset))
  ft <- build_within_share_flextable(
    rows, traits,
    c(cor = "cor(Env_m, Env_w)", cross = "|cross term| / (V_m + V_w)", offset = "Taxon offset"),
    paste("cor(Env_m, Env_w): correlation between taxon climatic centroids and within-taxon",
          "deviations across observations. |cross term| / (V_m + V_w): size of the omitted",
          "covariance term 2 B_m B_w cov(Env_m, Env_w) relative to the turnover and",
          "within-species components, median [2.5%, 97.5%] of the posterior. Taxon offset:",
          "mean absolute taxon mean of Env_w relative to sd(Env_w) (0 = deviations average",
          "to zero within every taxon)."))
  flextable::save_as_docx(ft, path = path)
  invisible(ft)
}
