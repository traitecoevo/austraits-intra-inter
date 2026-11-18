fit_model2d_bayesian_3panel_variance_metric <- function(data_model, y_limits = c(1e-6, 1e4), simple = FALSE) {

  # Plots the variance-weighted slope metric (as in the Appendix A chunks of Analysis.qmd):
  #   metric = (intra slope / inter slope)^2 * (var(deviation from centroid) / var(centroids))
  # i.e. the variance in the trait explained by the intraspecific component relative to that
  # explained by the interspecific component. metric = 1 means equal contributions.
  # Violins and 95% intervals come from the metric computed on every posterior draw;
  # the point and label show the metric computed from the posterior median slopes,
  # which is the value tabulated in Analysis.qmd.

  trait_name <- data_model$trait_name[1]

  # Axis label for the trait (used as x-label on the flipped plot)
  trait_label <- switch(trait_name,
    leaf_mass_per_area = expression(paste(LMA,             " (g ",  m^{-2}, ")")),
    huber_value        = expression(paste(SA:LA,            " (", mm^{2}, cm^{-2}, ")")),
    leaf_area          = expression(paste(Leaf~area,        " (", mm^{-2}, ")")),
    leaf_delta13C      = expression(paste(C[i], ":", C[a])),
    leaf_N_per_area    = expression(paste(Leaf~N~per~area,  " (g ", m^{-2}, ")")),
    wood_density       = expression(paste(Wood~density,     " (g ", cm^{-3}, ")")),
    leaf_P_per_area    = expression(paste(P[area],          " (g ", m^{-2}, ")")),
    seed_dry_mass = expression(paste(Seed~mass,        " (mg)"))
  )

  env_vars <- c(MAP = "std_log_annual_precipitation", MAT = "std_average_temperature")

  # Returns TRUE when the 90% credible interval spans zero.
  is_not_different_from_zero <- function(fit, hypothesis_string) {
    h <- brms::hypothesis(fit, hypothesis_string, alpha = 0.1)$hypothesis
    h$CI.Lower <= 0 && h$CI.Upper >= 0
  }

  # Per-draw metric for one model / environment / group.
  # 'woody_shift' adds the :woodywoody interaction terms to get the woody slopes.
  # 'data_group' is the data used for the variance terms.
  metric_draws <- function(fit, env, name, data_group, woody_shift = FALSE) {
    var <- env_vars[[env]]
    draws <- posterior::as_draws_df(fit) %>% as_tibble()

    inter <- draws[[paste0("b_", var, "_mn")]]
    intra <- draws[[paste0("b_", var, "_c")]]
    inter_hyp <- paste0(var, "_mn")
    intra_hyp <- paste0(var, "_c")
    if (woody_shift) {
      inter <- inter + draws[[paste0("b_", var, "_mn:woodywoody")]]
      intra <- intra + draws[[paste0("b_", var, "_c:woodywoody")]]
      inter_hyp <- paste0(inter_hyp, " + ", var, "_mn:woodywoody")
      intra_hyp <- paste0(intra_hyp, " + ", var, "_c:woodywoody")
    }

    var_ratio <- var(data_group[[paste0(var, "_c")]], na.rm = TRUE) /
      var(data_group[[paste0(var, "_mn")]], na.rm = TRUE)

    # Group is marked non-significant only when BOTH inter and intra slopes include zero.
    no_fill <- is_not_different_from_zero(fit, paste0(inter_hyp, " = 0")) &&
      is_not_different_from_zero(fit, paste0(intra_hyp, " = 0"))
if(simple == TRUE){
    list(
      draws = tibble(env = env, name = name, value = (intra / inter)^2),
      summary = tibble(
        env = env,
        name = name,
        point_value = (median(intra) / median(inter))^2,
        no_fill = no_fill
      )
    )
    } else{
         list(
      draws = tibble(env = env, name = name, value = (intra / inter)^2 * var_ratio),
      summary = tibble(
        env = env,
        name = name,
        point_value = (median(intra) / median(inter))^2 * var_ratio,
        no_fill = no_fill
      )
    ) 
    }
  }

  # Traits with a woody / non-woody split use the growth-form interaction model.
  # The remaining traits (huber_value, wood_density, leaf_P_per_area) are
  # woody-only so a single all-species model is used without a growth-form split.
  woody_split_traits <- c("leaf_mass_per_area", "leaf_N_per_area", "leaf_delta13C", "leaf_area","seed_dry_mass")
  is_woody_split <- trait_name %in% woody_split_traits

  # Saved as list(fit_MAP_stan, fit_MAT_stan, data_model) [length 3].
  fits_combined <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))

  if (is_woody_split) {
    fits <- readRDS(paste0("outputs/model_fits/fit2d_bayesian_woody/", trait_name, ".RDS"))
    stopifnot(length(fits) == 3)

    data_nonwoody <- data_model %>% filter(woody == "non-woody")
    data_woody    <- data_model %>% filter(woody == "woody")

    results <- list(
      metric_draws(fits[[1]],          "MAP", "diff",          data_nonwoody),
      metric_draws(fits[[1]],          "MAP", "diff_woody",    data_woody, woody_shift = TRUE),
      metric_draws(fits_combined[[1]], "MAP", "diff_combined", data_model),
      metric_draws(fits[[2]],          "MAT", "diff",          data_nonwoody),
      metric_draws(fits[[2]],          "MAT", "diff_woody",    data_woody, woody_shift = TRUE),
      metric_draws(fits_combined[[2]], "MAT", "diff_combined", data_model)
    )
  } else {
    results <- list(
      metric_draws(fits_combined[[1]], "MAP", "diff", data_model),
      metric_draws(fits_combined[[2]], "MAT", "diff", data_model)
    )
  }

  cat("extracting fitted vals -> ")

  # x-positions for each quantity on the coord_flip axis
  spacing <- c(diff = 1.6, diff_woody = 1.7, diff_combined = 1.8)

  group_code <- function(name) {
    case_when(
      str_detect(name, "combined") ~ "c",
      str_detect(name, "woody")    ~ "a",
      TRUE                         ~ "b"
    )
  }

  data_for_plot_group <- map_dfr(results, "summary") %>%
    left_join(
      map_dfr(results, "draws") %>%
        group_by(env, name) %>%
        summarise(
          upp_quantile = quantile(value, 0.975),
          median = quantile(value, 0.5),
          low_quantile = quantile(value, 0.025),
          .groups = "drop"
        ),
      by = c("env", "name")
    ) %>%
    mutate(xpos = spacing[name], woody = group_code(name))

  data_for_plot <- map_dfr(results, "draws") %>%
    left_join(data_for_plot_group %>% select(env, name, no_fill), by = c("env", "name")) %>%
    mutate(xpos = spacing[name],
           woody = group_code(name),
           fill_status = if_else(no_fill, "no_fill", "fill"))

  y_axis_label <- expression(paste("(", beta[intra], "/", beta[inter], ")"^2,
                                   " × ", sigma[within]^2, "/", sigma[between]^2))
  ref_line <- 1

  # Woody-split traits show three groups; woody-only traits show a single (woody) group.
  if (is_woody_split) {
    violin_nudge <- 0
    label_nudge  <- 0.08
    colour_order <- c(6, 3, 5)
    colour_breaks <- c("a", "b", "c")
    violin_layer <- geom_violinhalf(aes(alpha = fill_status), flip = FALSE)
  } else {
    violin_nudge <- 0.08
    label_nudge  <- 0.17
    colour_order <- c(6)
    colour_breaks <- c("a")
    violin_layer <- geom_violinhalf(aes(alpha = fill_status), flip = FALSE, width = 0.09,
                                    position = position_nudge(x = violin_nudge))
    data_for_plot       <- data_for_plot       %>% mutate(woody = "a")
    data_for_plot_group <- data_for_plot_group %>% mutate(woody = "a")
  }

  a <- data_for_plot %>%
    ggplot(aes(y = value, x = xpos,
               group = interaction(name, woody), colour = woody, fill = woody)) +
    geom_hline(yintercept = ref_line, linetype = "dashed") +
    violin_layer +
    geom_segment(data = data_for_plot_group,
                 aes(y = low_quantile, yend = upp_quantile, x = xpos, xend = xpos),
                 colour = "black", size = 1.2,
                 position = position_nudge(x = violin_nudge)) +
    geom_point(data = data_for_plot_group, aes(y = median),
               colour = "black", position = position_nudge(x = violin_nudge)) +
    # Metric label placed above the violin half (positive x nudge in coord_flip = upward)
    geom_text(data = data_for_plot_group,
              aes(y = median, label = signif(median, 2)),
              colour = "black", size = 3, position = position_nudge(x = label_nudge)) +
    facet_grid(~env) +
    # The metric is unbounded and right-skewed, so it is shown on a log10 scale.
    # coord limits (rather than scale limits) zoom without dropping draws.
    coord_flip(ylim = y_limits) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 12),
          axis.text = element_text(size = 10.5),
          axis.text.y = element_blank(),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = 16)) +
    labs(y = y_axis_label) +
    xlab(trait_label) +
    scale_y_log10(breaks = 10^seq(-10, 10, by = 2),
                  labels = scales::label_math(10^.x, format = log10)) +
    scale_colour_okabeito(order = colour_order, breaks = colour_breaks) +
    scale_fill_okabeito(order = colour_order,   breaks = colour_breaks) +
    scale_alpha_manual(values = c("fill" = 0.5, "no_fill" = 0))

  if (!is_woody_split) {
    a <- a + scale_x_continuous(limits = c(1.55, 1.85))
  }

  rm(list = intersect(c("fits", "fits_combined"), ls()))
  gc(verbose = FALSE)
  a
}
