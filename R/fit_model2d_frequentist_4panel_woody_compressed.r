source("R/species_within_slopes.R")
source("R/cache_is_stale.R")
source("R/test_group_slope_zero.R")
source("R/pick_interval_col.R")

# Produces a 2-panel figure summarising trait responses to a single climate driver,
# split by woody / non-woody growth form (compressed version of
# fit_model2d_frequentist_4panel_woody.R).
#
# Returns list(a, c + d) where:
#   a — Raw observations coloured by woody/non-woody + per-group OLS regression lines
#   c — Per-species mixed-model fitted lines (intraspecific component, solid) +
#         average intraspecific slope per group (solid) + species-mean points +
#         per-group interspecific fitted line (dashed) — all coloured by woody/non-woody
#   d — Density of per-species intraspecific slopes (coloured by woody); vertical lines:
#         OLS per-group slopes (dashed) and fixed interspecific β_mn per group (solid)
#
# Colour is reserved for woody vs non-woody throughout panel c; the intraspecific
# (solid) and interspecific (dashed) fitted lines are distinguished by linetype only.
#
# Bᵣ = |β_c| / (|β_mn| + |β_c|) computed separately for each growth-form group.
# β_mn = interspecific (between-species) slope; β_c = intraspecific (within-species) slope.
# For traits with a woody-split Bayesian model, β are extracted from brms fixed effects:
#   non-woody: rows 2–3; woody: rows 2–3 plus woody interaction rows 5–6.
# For huber_value and wood_density (fallback to all-species model), a single Bᵣ is used.
#
# Arguments:
#   data_model   — data frame for one trait; must contain trait_name and woody columns
#   ind_variable — climate driver: "log_annual_precipitation" (default) or
#                  "average_temperature"
#   remove_ns    — if TRUE, drop the OLS overall line, the interspecific (β_mn) line,
#                  and the average intraspecific (β_c) line from panels a/c for
#                  whichever growth-form group (woody / non-woody) is not
#                  (marginally) significant at P < 0.10. Each group is tested
#                  independently, so e.g. a significant woody slope can be kept
#                  while a non-significant non-woody slope is dropped. Per-species
#                  intraspecific lines and panel d are always kept regardless of
#                  this setting.
fit_model2d_frequentist_4panel_woody_compressed <- function(data_model, ind_variable = "log_annual_precipitation", remove_ns = FALSE) {
  # --- Derive helper variables from ind_variable ---
  # f:             back-transform function (exp for log-variables, identity for temperature)
  # standard_name: untransformed variable name used for axis labels and column naming
  # x_transform:   ggplot scale transform for the x-axis
  # off_variable:  the climate variable held fixed (opposite to ind_variable)
  if (str_detect(ind_variable, "log")) { f <- exp } else { f <- I }
  if (str_detect(ind_variable, "log")) { standard_name <- str_remove(ind_variable, "log_") } else { standard_name <- ind_variable }
  if (str_detect(ind_variable, "log")) { x_transform <- "log10" } else { x_transform <- "identity" }
  off_variable <- if (ind_variable == "log_annual_precipitation") "average_temperature" else "log_annual_precipitation"
  # Axis labels: only populated for leaf_area (bottom row in combined multi-panel figures)
  environment_label <- ""
  b1_label          <- ""
  if (data_model$trait_name[1] == "seed_dry_mass") {
    environment_label <- if (ind_variable == "log_annual_precipitation") "MAP (mm)" else expression(paste(MAT, " (", degree * C, ")"))
    b1_label          <- if (ind_variable == "log_annual_precipitation") expression(paste("Intraspecific response to MAP")) else expression(paste("Intraspecific response to MAT"))
  }

  # Trait y-axis labels
  trait_name <- data_model$trait_name[1]
  if (trait_name == "leaf_mass_per_area") trait_label <- expression(paste(LMA,                   " (g ",  m^{-2},   ")"))
  if (trait_name == "huber_value")        trait_label <- expression(paste(SA:LA,                 " (", mm^{2}, cm^{-2},  ")"))
  if (trait_name == "leaf_area")          trait_label <- expression(paste(Leaf ~ area,            " (",    mm^{-2},  ")"))
  if (trait_name == "leaf_delta13C")      trait_label <- expression(paste(C[i],                  ":",     C[a]))
  if (trait_name == "leaf_N_per_area")    trait_label <- expression(paste(Leaf ~ N ~ per ~ area, " (g ",  m^{-2},   ")"))
  if (trait_name == "wood_density")       trait_label <- expression(paste(Wood ~ density,         " (g ",  cm^{-3},  ")"))
  if (trait_name == "leaf_P_per_area")    trait_label <- expression(paste(P[area],               " (g ",  m^{-2},   ")"))
  if (trait_name == "seed_dry_mass")    trait_label <- expression(paste(Seed~mass,        " (mg)"))
  # Temperature plots share the y-axis with the precipitation plots alongside them, so
  # suppress the trait y-axis title to avoid repeating it.
  if (ind_variable == "average_temperature") trait_label <- ""

  # --- Load frequentist model fit ---
  # Length-2 RDS: list(fitmod_MAP_or_combined, data_model)
  # Length-3 RDS: list(fitmod_MAP, fitmod_MAT, data_model)
  fits <- readRDS(paste0("outputs/model_fits/fit2d_frequentist_woody/", trait_name, ".RDS"))
  bayes_ratio_path <- paste0("outputs/model_fits/fit2d_bayesian_woody/bayes_ratios/", trait_name, "_", ind_variable, ".RDS")

  # --- Compute or load Bayesian ratio (Bᵣ) ---
  # For huber_value and wood_density the woody Bayesian model did not converge;
  # fall back to the all-species Bayesian model for a single shared Bᵣ.
  # For other traits, Bᵣ is computed separately for woody and non-woody from the
  # woody Bayesian brms model. Results are cached to disk and recomputed
  # whenever the source fit is newer than the cache (see cache_is_stale()).
  if (trait_name %in% c("huber_value", "wood_density")) {
    if (cache_is_stale(bayes_ratio_path, paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))) {
      fits_bayes  <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))
      bayes_coefs <- if (ind_variable == "log_annual_precipitation") fixef(fits_bayes[[1]])[2:3, 1] else fixef(fits_bayes[[2]])[2:3, 1]
      bayes_denom <- abs(bayes_coefs[1]) + abs(bayes_coefs[2])
      bayes_ratio <- if (bayes_denom == 0) NA_real_ else abs(bayes_coefs[2]) / bayes_denom
      bayes_ratio <- round(bayes_ratio, 2)
      if (sign(bayes_coefs[1]) != sign(bayes_coefs[2])) bayes_ratio <- paste0(bayes_ratio, "*")
      dir.create(dirname(bayes_ratio_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(bayes_ratio, bayes_ratio_path)
    } else {
      bayes_ratio <- readRDS(bayes_ratio_path)
    }
  } else {
    if (cache_is_stale(bayes_ratio_path, paste0("outputs/model_fits/fit2d_bayesian_woody/", trait_name, ".RDS"))) {
      fits_bayes  <- readRDS(paste0("outputs/model_fits/fit2d_bayesian_woody/", trait_name, ".RDS"))
      # Rows 2–3: β_mn and β_c for non-woody; rows 5–6: woody interaction terms
      bayes_coefs <- if (ind_variable == "log_annual_precipitation") fixef(fits_bayes[[1]])[c(2:3, 5:6), 1] else fixef(fits_bayes[[2]])[c(2:3, 5:6), 1]
      # Non-woody: β_mn = bayes_coefs[1], β_c = bayes_coefs[2]
      bayes_denom_nonwoody <- abs(bayes_coefs[1]) + abs(bayes_coefs[2])
      # Woody: β_mn = bayes_coefs[1] + bayes_coefs[3], β_c = bayes_coefs[2] + bayes_coefs[4]
      bayes_denom_woody    <- abs(bayes_coefs[1] + bayes_coefs[3]) + abs(bayes_coefs[2] + bayes_coefs[4])
      bayes_nonwoody_ratio <- if (bayes_denom_nonwoody == 0) NA_real_ else abs(bayes_coefs[2]) / bayes_denom_nonwoody
      bayes_woody_ratio    <- if (bayes_denom_woody    == 0) NA_real_ else abs(bayes_coefs[2] + bayes_coefs[4]) / bayes_denom_woody
      bayes_nonwoody_ratio <- round(bayes_nonwoody_ratio, 2)
      bayes_woody_ratio    <- round(bayes_woody_ratio,    2)
      if (sign(bayes_coefs[1])                  != sign(bayes_coefs[2]))                  bayes_nonwoody_ratio <- paste0(bayes_nonwoody_ratio, "*")
      if (sign(bayes_coefs[1] + bayes_coefs[3]) != sign(bayes_coefs[2] + bayes_coefs[4])) bayes_woody_ratio    <- paste0(bayes_woody_ratio, "*")
      bayes_ratio <- list(bayes_nonwoody_ratio, bayes_woody_ratio)
      dir.create(dirname(bayes_ratio_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(bayes_ratio, bayes_ratio_path)
    } else {
      bayes_ratio <- readRDS(bayes_ratio_path)
    }
  }

  # Select the correct fitted model and recover the data used for fitting
  if (length(fits) == 2) {
    fitmod     <- fits[[1]]
    data_model <- fits[[2]]
  }

  if (length(fits) == 3) {
    fitmod     <- if (ind_variable == "log_annual_precipitation") fits[[1]] else fits[[2]]
    data_model <- fits[[3]]
  }

  # --- Global scaling parameters (defined in the calling environment) ---
  if (ind_variable == "log_annual_precipitation") {
    scale_mn <- log_annual_precipitation_scaler_scale
    mean_mn  <- log_annual_precipitation_scaler_center
    scale_c  <- log_annual_precipitation_scaler_scale
    mean_c   <- log_annual_precipitation_scaler_center
  } else {
    scale_mn <- average_temperature_scaler_scale
    mean_mn  <- average_temperature_scaler_center
    scale_c  <- average_temperature_scaler_scale
    mean_c   <- average_temperature_scaler_center
  }

  cat("extracting fitted vals -> ")

  # Fix the off-variable at each species' mean to isolate the ind_variable response.
  # The woody column is retained so estimates are split by growth form.
  data_for_fit <- data_model %>%
    select(std_log_annual_precipitation_mn, std_log_annual_precipitation_c,
           std_average_temperature_mn,      std_average_temperature_c,
           taxon_name, woody) %>%
    group_by(taxon_name) %>%
    mutate(
      !!paste0("std_", off_variable, "_c")  := mean(get(paste0("std_", off_variable, "_c")),  na.rm = TRUE),
      !!paste0("std_", off_variable, "_mn") := mean(get(paste0("std_", off_variable, "_mn")), na.rm = TRUE)
    ) %>%
    drop_na()

  # Per-species fitted values (with random effects) over the observed climate range
  data_fit2_re <- estimate_expectation(fitmod, data = data_for_fit) |> as_tibble()

  # Interspecific prediction: population-level fitted values varying only the
  # species-mean climate (_mn), separately for each woody group, with _c fixed at 0.
  data_fit2_global_raw <- estimate_means(fitmod, by = c("woody", paste0("std_", ind_variable, "_mn"))) |>
    as_tibble()

  data_fit2_global <- data_fit2_global_raw |>
    mutate(
      log_value     = Mean,
      log_value_low = pick_interval_col(data_fit2_global_raw, c("CI_low", "CI", "Conf.Low", "conf.low")),
      log_value_hi  = pick_interval_col(data_fit2_global_raw, c("CI_high", "CI", "Conf.High", "conf.high"))
    ) |>
    select(!!(paste0("std_", ind_variable)) := paste0("std_", ind_variable, "_mn"), log_value, log_value_low, log_value_hi, woody) |>
    mutate(
      !!standard_name := as.numeric(f(get(paste0("std_", ind_variable)) * scale_mn + mean_mn)),
      value      = exp(log_value),
      value_low  = exp(log_value_low),
      value_hi   = exp(log_value_hi),
      taxon_name = "global average"
    )

  # Species-mean predictions: set both _c deviations to zero to obtain the predicted
  # trait value at each species' mean climate position (no within-species deviation).
  data_fit2_mn <- estimate_expectation(
    fitmod,
    data = data_fit2_re |>
      mutate(std_log_annual_precipitation_c = 0, std_average_temperature_c = 0) |>
      select(-Predicted) |>
      distinct()
  ) |>
    as_tibble() |>
    distinct()

  # Build a 4-row prediction dataset (2 woody groups × 2 endpoints) for the average
  # intraspecific slope line of each growth form: the group's within-species response
  # (B_w, or B_w + ΔB_w for woody) of a taxon centred at the group's median taxon
  # centroid, drawn across the group's observed climate range x (min and max of
  # _mn + _c). Prediction uses _mn = group median centroid and _c = x - that median, so
  # each line crosses its group's turnover line at the median centroid; plotted at x.
  if (length(fits) > 2) {
    data_c_line <- data_fit2_re %>%
      group_by(woody) %>%
      summarise(
        !!paste0("std_max_05_", ind_variable, "_mn") := max(get(paste0("std_", ind_variable, "_mn"))),
        !!paste0("std_min_05_", ind_variable, "_mn") := min(get(paste0("std_", ind_variable, "_mn"))),
        !!paste0("std_med_",    ind_variable, "_mn") := median(get(paste0("std_", ind_variable, "_mn"))),
        !!paste0("std_min_05_", ind_variable)        := min(get(paste0("std_", ind_variable, "_mn")) + get(paste0("std_", ind_variable, "_c"))),
        !!paste0("std_max_05_", ind_variable)        := max(get(paste0("std_", ind_variable, "_mn")) + get(paste0("std_", ind_variable, "_c")))
      ) %>%
      pivot_longer(
        cols      = c(paste0("std_min_05_", ind_variable), paste0("std_max_05_", ind_variable)),
        values_to = paste0("std_", ind_variable, "_c_plot")
      ) %>%
      mutate(
        # Deviation from the median centroid, so the line passes through the turnover
        # line at the median taxon's centroid (prediction uses _mn = median, _c = x - median)
        !!paste0("std_", ind_variable, "_c")  := get(paste0("std_", ind_variable, "_c_plot")) - get(paste0("std_med_", ind_variable, "_mn")),
        !!paste0("std_", ind_variable, "_mn") := get(paste0("std_", "med_", ind_variable, "_mn")),
        !!paste0(ind_variable, "_mn")         := get(paste0("std_", ind_variable, "_mn")) * scale_c + mean_c,
        # x position for plotting: the absolute climate value (not the deviation)
        !!paste0(ind_variable, "_c")          := get(paste0("std_", ind_variable, "_c_plot")) * scale_c + mean_c,
        taxon_name = "global"
      )
  } else {
    data_c_line <- data_fit2_re %>%
      group_by(woody) %>%
      summarise(
        !!paste0("std_max_05_", ind_variable, "_mn") := max(get(paste0("std_", ind_variable, "_mn"))),
        !!paste0("std_min_05_", ind_variable, "_mn") := min(get(paste0("std_", ind_variable, "_mn"))),
        !!paste0("std_med_",    ind_variable, "_mn") := median(get(paste0("std_", ind_variable, "_mn"))),
        !!paste0("std_min_05_", ind_variable)        := min(get(paste0("std_", ind_variable, "_mn")) + get(paste0("std_", ind_variable, "_c"))),
        !!paste0("std_max_05_", ind_variable)        := max(get(paste0("std_", ind_variable, "_mn")) + get(paste0("std_", ind_variable, "_c"))),
        !!paste0("std_", off_variable, "_c")  := 0,
        !!paste0("std_", off_variable, "_mn") := mean(get(paste0("std_", off_variable, "_mn")))
      ) %>%
      pivot_longer(
        cols      = c(paste0("std_min_05_", ind_variable), paste0("std_max_05_", ind_variable)),
        values_to = paste0("std_", ind_variable, "_c_plot")
      ) %>%
      mutate(
        # Deviation from the median centroid, so the line passes through the turnover
        # line at the median taxon's centroid (prediction uses _mn = median, _c = x - median)
        !!paste0("std_", ind_variable, "_c")  := get(paste0("std_", ind_variable, "_c_plot")) - get(paste0("std_med_", ind_variable, "_mn")),
        !!paste0("std_", ind_variable, "_mn") := get(paste0("std_", "med_", ind_variable, "_mn")),
        !!paste0(ind_variable, "_mn")         := get(paste0("std_", ind_variable, "_mn")) * scale_c + mean_c,
        # x position for plotting: the absolute climate value (not the deviation)
        !!paste0(ind_variable, "_c")          := get(paste0("std_", ind_variable, "_c_plot")) * scale_c + mean_c,
        taxon_name = "global"
      )
  }

  # Predict trait values along the intraspecific slope line and back-transform
  data_c_line_est <- estimate_expectation(fitmod, data_c_line) |> as_tibble()
  data_c_line$est_Exp     <- data_c_line_est$Predicted
  data_c_line$est_Exp_low <- pick_interval_col(data_c_line_est, c("CI_low", "CI", "Conf.Low", "conf.low"))
  data_c_line$est_Exp_hi  <- pick_interval_col(data_c_line_est, c("CI_high", "CI", "Conf.High", "conf.high"))
  if (trait_name != "leaf_delta13C") {
    data_c_line$est_Exp     <- exp(data_c_line$est_Exp)
    data_c_line$est_Exp_low <- exp(data_c_line$est_Exp_low)
    data_c_line$est_Exp_hi  <- exp(data_c_line$est_Exp_hi)
  }

  # Back-transform x-axis: standard_name_c = f(ind_variable_c) gives actual environmental units
  data_c_line <- data_c_line |>
    mutate(!!paste0(standard_name, "_c") := f(get(paste0(ind_variable, "_c"))))

  # Combine per-species fitted values with their species-mean predictions.
  # Adds back-transformed climate columns and woody column for plotting.
  data_fit2 <- data_fit2_re |>
    rename(log_value = Predicted) |>
    left_join(
      data_fit2_mn |> select(paste0("std_", ind_variable, "_mn"), taxon_name, log_value_mn = Predicted, woody)
    ) |>
    mutate(
      !!ind_variable                 := (get(paste0("std_", ind_variable, "_mn")) + get(paste0("std_", ind_variable, "_c"))) * scale_c + mean_c,
      !!standard_name                := as.numeric(f(get(paste0(ind_variable)))),
      !!paste0(standard_name, "_mn") := as.numeric(f(get(paste0("std_", ind_variable, "_mn")) * scale_c + mean_c))
    )

  if (trait_name == "leaf_delta13C") {
    data_fit2_global <- data_fit2_global |>
      mutate(
        value    = log_value,
        value_low = log_value_low,
        value_hi  = log_value_hi
      )
  }

  cat("cooking up a sweet plot -> ")

  # Back-transform data_model climate columns to original scale for plotting
  data_model <- data_model %>%
    mutate(
      log_annual_precipitation = std_log_annual_precipitation * log_annual_precipitation_scaler_scale + log_annual_precipitation_scaler_center,
      annual_precipitation     = exp(std_log_annual_precipitation * log_annual_precipitation_scaler_scale + log_annual_precipitation_scaler_center),
      average_temperature      = std_average_temperature * average_temperature_scaler_scale + average_temperature_scaler_center
    )

  # --- OLS overall regression, split by woody (naive, ignoring species identity) ---
  # Used only for panel a; includes interaction with woody.
  if (length(fits) > 2) {
    if (ind_variable == "log_annual_precipitation") {
      fit_overall <- lm(log_value ~ std_log_annual_precipitation * woody, data = data_model, na.action = na.omit)
      # doc <- read_docx()
      # doc <- body_add_table(doc, value = as.data.frame(summary(fit_overall)$coef), style = "table_template")
      # print(doc, target = paste0(trait_name, "_MAP_woody_overall.docx"))
    }
    if (ind_variable == "average_temperature") {
      fit_overall <- lm(log_value ~ std_average_temperature * woody, data = data_model, na.action = na.omit)
      # doc <- read_docx()
      # doc <- body_add_table(doc, value = as.data.frame(summary(fit_overall)$coef), style = "table_template")
      # print(doc, target = paste0(trait_name, "_MAT_woody_overall.docx"))
    }
  } else {
    fit_overall <- lm(log_value ~ (std_log_annual_precipitation + std_average_temperature) * woody, data = data_model, na.action = na.omit)
  }
  print(summary(fit_overall))

  # OLS predictions for plotting (marginalised over the focal climate variable, per woody group)
  data_fit2_re_overall <- estimate_expectation(fit_overall, by = c("woody", paste0("std_", ind_variable))) |>
    as_tibble() |>
    mutate(!!ind_variable := f(get(paste0("std_", ind_variable)) * scale_c + mean_c))

  # --- Optional significance filtering (per growth-form group) ---
  # When remove_ns = TRUE, drop a group's OLS overall / interspecific (β_mn) /
  # intraspecific (β_c) line from panels a/c if that group's slope is not
  # (marginally) significant at ns_alpha. Groups are tested independently.
  ns_alpha     <- 0.10
  overall_term <- paste0("std_", ind_variable)
  mn_term      <- paste0("std_", ind_variable, "_mn")
  c_term       <- paste0("std_", ind_variable, "_c")
  groups       <- c("non-woody", "woody")

  sig_overall <- setNames(rep(TRUE, 2), groups)
  sig_mn      <- setNames(rep(TRUE, 2), groups)
  sig_c       <- setNames(rep(TRUE, 2), groups)

  if (remove_ns) {
    for (grp in groups) {
      sig_overall[grp] <- test_group_slope_zero(fit_overall, overall_term, grp)$p.value < ns_alpha
      sig_mn[grp]      <- test_group_slope_zero(fitmod,      mn_term,      grp)$p.value < ns_alpha
      sig_c[grp]       <- test_group_slope_zero(fitmod,      c_term,       grp)$p.value < ns_alpha
    }
  }

  data_fit2_re_overall <- data_fit2_re_overall |> filter(woody %in% names(sig_overall)[sig_overall])
  data_fit2_global     <- data_fit2_global     |> filter(woody %in% names(sig_mn)[sig_mn])
  data_c_line          <- data_c_line          |> filter(woody %in% names(sig_c)[sig_c])

  # Per-species intraspecific slopes (glmmTMB): non-woody fixed β_c + each species'
  # random slope, plus the fixed woody interaction (added for woody taxa below).
  # Built from fixef() + ranef() — see species_within_slopes().
  re_strength <- species_within_slopes(fitmod, ind_variable)

  # ---- Panel a: raw data + per-group OLS lines ----
  if (trait_name != "leaf_delta13C") {
    a <- data_model |>
      mutate(!!ind_variable := as.numeric(f(get(ind_variable)))) |>
      ggplot(aes(x = get(ind_variable), y = exp(log_value), colour = woody, group = woody)) +
      geom_point(alpha = 0.75) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 1.5 * 1.5, colour = "black") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 1.5) +
      scale_x_continuous(transform = x_transform) +
      scale_y_log10() +
      theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none",
            text         = element_text(size = 15),
            axis.text    = element_text(size = 15),
            axis.title   = element_text(size = 16.25),
            legend.text  = element_text(size = 15),
            legend.title = element_text(size = 16.25)) +
      scale_colour_okabeito(order = c(3, 6))
  } else {
    a <- data_model |>
      mutate(!!ind_variable := as.numeric(f(get(ind_variable)))) |>
      ggplot(aes(x = get(ind_variable), y = log_value, colour = woody, group = woody)) +
      geom_point(alpha = 0.75) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 1.5 * 1.5, colour = "black") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 1.5) +
      scale_x_continuous(transform = x_transform) +
      theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none",
            text         = element_text(size = 15),
            axis.text    = element_text(size = 15),
            axis.title   = element_text(size = 16.25),
            legend.text  = element_text(size = 15),
            legend.title = element_text(size = 16.25)) +
      scale_colour_okabeito(order = c(3, 6))
  }

  # ---- Panel c: per-species intraspecific lines (solid) + average intraspecific
  #               slope per group (solid) + species-mean points + per-group
  #               interspecific fitted line (dashed) ----
  # Colour is reserved for woody vs non-woody; linetype (solid vs dashed) is the
  # only thing distinguishing the intraspecific line from the interspecific line.
  #
  # fill is mapped only by the two geom_ribbon() layers below (data_c_line /
  # data_fit2_global), which remove_ns can filter down to a single group. The
  # positional scale_*_okabeito(order = ...) assigns colour by position within
  # whatever levels are actually present in the fill-mapped layers, so a lone
  # surviving group would silently get the *other* group's colour. Look the
  # colours up by name instead so they stay tied to "woody"/"non-woody" identity
  # regardless of which group's ribbon/line survives the significance filter.
  woody_levels  <- c("non-woody", "woody")
  # palette_okabeito()'s own default (palette = "full_amber") isn't a valid name in
  # okabeito_palettes and falls back to "full" with a warning; pass it explicitly to
  # match scale_colour_okabeito()/scale_fill_okabeito()'s default and avoid the warning.
  woody_colours <- setNames(palette_okabeito(palette = "full", order = c(3, 6))(2), woody_levels)
  if (trait_name != "leaf_delta13C") {
    c <- data_fit2 |>
      ggplot(aes(get(standard_name), exp(log_value), group = taxon_name)) +
      geom_point(aes(colour = woody), alpha = 0.5, shape = 1, data = data_model) +
      # One point per taxon at its climatic centroid (x = taxon mean climate, y = predicted
      # trait value there, i.e. with no within-taxon deviation; data_fit2_mn)
      geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
               aes(x = get(paste0(standard_name, "_mn")), y = exp(log_value_mn), colour = woody, group = woody), alpha = 0.75) +

      geom_line(aes(group = interaction(woody, taxon_name), colour = woody, alpha = 0.75), linewidth = 0.75) +
      geom_ribbon(data = data_c_line,
          aes(x = as.numeric(get(paste0(standard_name, "_c"))), ymin = est_Exp_low, ymax = est_Exp_hi, group = woody, fill = woody),
                  alpha = 0.15, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody), col = "black", linewidth = 1.5 * 1.5) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody, colour = woody), linewidth = 1.5) +
      geom_ribbon(data = data_fit2_global,
          aes(x = as.numeric(get(standard_name)), ymin = value_low, ymax = value_hi, group = woody, fill = woody),
                  alpha = 0.15, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_fit2_global, aes(group = woody), linewidth = 1.5 * 1.5, colour = "black", linetype = "dashed") +
      geom_line(data = data_fit2_global, aes(group = woody, colour = woody), linewidth = 1.5, linetype = "dashed") +
      scale_x_continuous(transform = x_transform) +
      scale_y_log10() +
      theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none",
            plot.title   = element_text(hjust = 0.5, size = 13.75),
            text         = element_text(size = 15),
            axis.text    = element_text(size = 15),
            axis.title   = element_text(size = 16.25),
            legend.text  = element_text(size = 15),
            legend.title = element_text(size = 16.25)) +
      scale_colour_manual(values = woody_colours, breaks = woody_levels) +
      scale_fill_manual(values = woody_colours, breaks = woody_levels)
  } else {
    c <- data_fit2 |>
      ggplot(aes(get(standard_name), log_value, group = taxon_name)) +
      geom_point(aes(colour = woody), alpha = 0.5, shape = 1, data = data_model) +
      # One point per taxon at its climatic centroid (x = taxon mean climate, y = predicted
      # trait value there, i.e. with no within-taxon deviation; data_fit2_mn)
      geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
               aes(x = get(paste0(standard_name, "_mn")), y = log_value_mn, colour = woody, group = woody), alpha = 0.75) +

      geom_line(aes(group = interaction(woody, taxon_name), colour = woody, alpha = 0.75), linewidth = 0.75) +
      geom_ribbon(data = data_c_line,
          aes(x = as.numeric(get(paste0(standard_name, "_c"))), ymin = est_Exp_low, ymax = est_Exp_hi, group = woody, fill = woody),
                  alpha = 0.15, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody), col = "black", linewidth = 1.5 * 1.5) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody, colour = woody), linewidth = 1.5) +
        geom_ribbon(data = data_fit2_global,
          aes(x = as.numeric(get(standard_name)), ymin = value_low, ymax = value_hi, group = woody, fill = woody),
                      alpha = 0.15, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_fit2_global, aes(group = woody), linewidth = 1.5 * 1.5, colour = "black", linetype = "dashed") +
      geom_line(data = data_fit2_global, aes(group = woody, colour = woody), linewidth = 1.5, linetype = "dashed") +
      scale_x_continuous(transform = x_transform) +
      theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none",
            plot.title   = element_text(hjust = 0.5, size = 13.75),
            text         = element_text(size = 15),
            axis.text    = element_text(size = 15),
            axis.title   = element_text(size = 16.25),
            legend.text  = element_text(size = 15),
            legend.title = element_text(size = 16.25)) +
      scale_colour_manual(values = woody_colours, breaks = woody_levels) +
      scale_fill_manual(values = woody_colours, breaks = woody_levels)
  }

  # ---- Panel d: density of per-species intraspecific slopes, split by woody ----
  # For woody species the total intraspecific slope = re_strength + woody_effect.
  # Vertical lines (twodash): OLS per-group slopes.
  # Vertical lines (solid): fixed interspecific β_mn per group from the mixed model.
  # The β_mn reference lines are intentional — they anchor the inter- vs intraspecific
  # comparison alongside the density of intraspecific slopes.
  data_model %>%
    filter(taxon_name %in% re_strength$taxon_name) %>%
    distinct(taxon_name, .keep_all = TRUE) %>%
    select(taxon_name, woody) -> woody_key

  d <- re_strength |>
    left_join(woody_key) |>
    mutate(re_strength = ifelse(woody == "woody", re_strength + woody_effect, re_strength)) %>%
    ggplot(aes(x = re_strength, colour = woody, group = woody, fill = woody)) +
    geom_vline(xintercept = coef(fit_overall)[paste0("std_", ind_variable)],
               col = "#009E73", linewidth = 1, linetype = "twodash") +
    geom_vline(xintercept = coef(fit_overall)[paste0("std_", ind_variable)] + coef(fit_overall)[paste0("std_", ind_variable, ":woodywoody")],
               col = "#D55E00", linewidth = 1, linetype = "twodash") +
    geom_vline(xintercept = fixef(fitmod)$cond[paste0("std_", ind_variable, "_mn")],
               col = "#009E73", linewidth = 1.25) +
    geom_vline(xintercept = fixef(fitmod)$cond[paste0("std_", ind_variable, "_mn")] + fixef(fitmod)$cond[paste0("std_", ind_variable, "_mn:woodywoody")],
               col = "#D55E00", linewidth = 1.25) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_density(linewidth = 1.5, alpha = 0.1) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = c(-max(abs(re_strength$re_strength)) * 1.1, max(abs(re_strength$re_strength)) * 1.1)) +
    theme_classic() +
    labs(x = b1_label, y = "Density") +
    theme(legend.position = "none",
          text         = element_text(size = 15),
          axis.text    = element_text(size = 15),
          axis.title   = element_text(size = 16.25),
          legend.text  = element_text(size = 15),
          legend.title = element_text(size = 16.25)) +
    scale_colour_okabeito(order = c(3, 6)) +
    scale_fill_okabeito(order = c(3, 6))

  return(list(a, c + d + plot_layout(nrow = 1)))
}
