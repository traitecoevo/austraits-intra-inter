# Produces a 4-panel figure summarising trait responses to a single climate driver,
# split by woody / non-woody growth form.
#
# Returns list(a, c + b + d) where:
#   a — Raw observations coloured by woody/non-woody + per-group OLS regression lines
#   b — Per-species mixed-model fitted lines (intraspecific component) +
#         average intraspecific slope per group + Bᵣ annotation for each group
#   c — Species-mean points (coloured by woody) + per-group interspecific fitted lines +
#         per-group OLS lines for visual comparison
#   d — Density of per-species intraspecific slopes (coloured by woody); vertical lines:
#         OLS per-group slopes (dashed) and fixed interspecific β_mn per group (solid)
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
fit_model2d_frequentist_4panel_woody <- function(data_model, ind_variable = "log_annual_precipitation") {
  # --- Derive helper variables from ind_variable ---
  # f:             back-transform function (exp for log-variables, identity for temperature)
  # standard_name: untransformed variable name used for axis labels and column naming
  # x_transform:   ggplot scale transform for the x-axis
  # off_variable:  the climate variable held fixed (opposite to ind_variable)
  if (str_detect(ind_variable, "log")) { f <- exp } else { f <- I }
  if (str_detect(ind_variable, "log")) { standard_name <- str_remove(ind_variable, "log_") } else { standard_name <- ind_variable }
  if (str_detect(ind_variable, "log")) { x_transform <- "log10" } else { x_transform <- "identity" }
  off_variable <- if (ind_variable == "log_annual_precipitation") "average_temperature" else "log_annual_precipitation"

  # Axis labels: only populated for leaf_N_per_area (bottom row in combined multi-panel figures)
  environment_label <- ""
  b1_label          <- ""
  if (data_model$trait_name[1] == "leaf_N_per_area") {
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

  # --- Load frequentist model fit ---
  # Length-2 RDS: list(fitmod_MAP_or_combined, data_model)
  # Length-3 RDS: list(fitmod_MAP, fitmod_MAT, data_model)
  fits <- readRDS(paste0("outputs/model_fits/fit2d_frequentist_woody/", trait_name, ".RDS"))
  bayes_ratio_path <- paste0("outputs/model_fits/fit2d_bayesian_woody/bayes_ratios/", trait_name, "_", ind_variable, ".RDS")

  # --- Compute or load Bayesian ratio (Bᵣ) ---
  # For huber_value and wood_density the woody Bayesian model did not converge;
  # fall back to the all-species Bayesian model for a single shared Bᵣ.
  # For other traits, Bᵣ is computed separately for woody and non-woody from the
  # woody Bayesian brms model. Results are cached to disk.
  if (trait_name %in% c("huber_value", "wood_density")) {
    if (!file.exists(bayes_ratio_path)) {
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
    if (!file.exists(bayes_ratio_path)) {
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

  # Helper: normalise interval column names returned by different prediction methods.
  # Returns NA when the prediction object does not include intervals.
  pick_interval_col <- function(df, candidates) {
    col_name <- candidates[candidates %in% names(df)][1]
    if (is.na(col_name)) return(rep(NA_real_, nrow(df)))
    as.numeric(df[[col_name]])
  }

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
  # intraspecific slope line. Fix _mn at each group's median climate; set _c to span
  # the full observed range (min and max of _mn + _c combined) within each group.
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
        !!paste0("std_", ind_variable, "_c")  := get(paste0("std_", ind_variable, "_c_plot")),
        !!paste0("std_", ind_variable, "_mn") := get(paste0("std_", "med_", ind_variable, "_mn")),
        !!paste0(ind_variable, "_mn")         := get(paste0("std_", ind_variable, "_mn")) * scale_c + mean_c,
        !!paste0(ind_variable, "_c")          := get(paste0("std_", ind_variable, "_c"))  * scale_c + mean_c,
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
        !!paste0("std_", ind_variable, "_c")  := get(paste0("std_", ind_variable, "_c_plot")),
        !!paste0("std_", ind_variable, "_mn") := get(paste0("std_", "med_", ind_variable, "_mn")),
        !!paste0(ind_variable, "_mn")         := get(paste0("std_", ind_variable, "_mn")) * scale_c + mean_c,
        !!paste0(ind_variable, "_c")          := get(paste0("std_", ind_variable, "_c"))  * scale_c + mean_c,
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
  # Used as a per-group reference line; includes interaction with woody.
  if (length(fits) > 2) {
    if (ind_variable == "log_annual_precipitation") {
      fit_overall <- lm(log_value ~ std_log_annual_precipitation * woody, data = data_model, na.action = na.omit)
      doc <- read_docx()
      doc <- body_add_table(doc, value = as.data.frame(summary(fit_overall)$coef), style = "table_template")
      print(doc, target = paste0(trait_name, "_MAP_woody_overall.docx"))
    }
    if (ind_variable == "average_temperature") {
      fit_overall <- lm(log_value ~ std_average_temperature * woody, data = data_model, na.action = na.omit)
      doc <- read_docx()
      doc <- body_add_table(doc, value = as.data.frame(summary(fit_overall)$coef), style = "table_template")
      print(doc, target = paste0(trait_name, "_MAT_woody_overall.docx"))
    }
  } else {
    fit_overall <- lm(log_value ~ (std_log_annual_precipitation + std_average_temperature) * woody, data = data_model, na.action = na.omit)
  }
  print(summary(fit_overall))

  # OLS predictions for plotting (marginalised over the focal climate variable, per woody group)
  data_fit2_re_overall <- estimate_expectation(fit_overall, by = c("woody", paste0("std_", ind_variable))) |>
    as_tibble() |>
    mutate(!!ind_variable := f(get(paste0("std_", ind_variable)) * scale_c + mean_c))

  # Per-species intraspecific slopes (glmmTMB): each species' β_c random effect, plus
  # the species-specific woody interaction effect (std_X_c:woodywoody random slope).
  re_strength <- tibble(
    taxon_name   = rownames(coef(fitmod)$cond$taxon_name),
    re_strength  = coef(fitmod)$cond$taxon_name[[paste0("std_", ind_variable, "_c")]],
    woody_effect = coef(fitmod)$cond$taxon_name[[paste0("std_", ind_variable, "_c:woodywoody")]]
  )

  # Format Bᵣ annotation: woody (orange) above non-woody (green)
  bayes_ratio_pasted <- paste0("Bᵣ = ", unlist(bayes_ratio))
  bayes_ratio_pasted <- paste0(
    "<span style='color:#D55E00;'>", bayes_ratio_pasted[2], "</span>",
    "<br>",
    "<span style='color:#009E73;'>", bayes_ratio_pasted[1], "</span>"
  )

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
            text         = element_text(size = 12),
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 13),
            legend.text  = element_text(size = 12),
            legend.title = element_text(size = 13)) +
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
            text         = element_text(size = 12),
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 13),
            legend.text  = element_text(size = 12),
            legend.title = element_text(size = 13)) +
      scale_colour_okabeito(order = c(3, 6))
  }

  # ---- Panel b: per-species fitted lines (intraspecific) + average intraspecific slope + Bᵣ ----
  if (trait_name != "leaf_delta13C") {
    b <- data_fit2 |>
      left_join(re_strength) %>%
      ggplot(aes(get(standard_name), exp(log_value))) +
      geom_point(aes(colour = woody, group = woody), alpha = 0.75, col = "grey", data = data_model) +
      geom_line(aes(group = interaction(woody, taxon_name), colour = woody, alpha = 0.75), linewidth = 0.75) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted), group = woody), alpha = 0.8, colour = "grey20", linewidth = 1.5) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted), group = woody, colour = woody), linewidth = 0.9, linetype = "twodash") +
      geom_ribbon(data = data_c_line,
          aes(x = as.numeric(get(paste0(standard_name, "_c"))), ymin = est_Exp_low, ymax = est_Exp_hi, group = woody, fill = woody),
                  alpha = 0.3, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody, colour = woody), col = "black", linewidth = 1.5 * 1.5) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody, colour = woody), linewidth = 1.5) +
      # annotate("richtext", x = Inf, y = Inf, label = bayes_ratio_pasted,
      #          hjust = 1.1, vjust = 1, size = 5, label.size = NA, fill = NA) +
      scale_x_continuous(transform = x_transform) +
      scale_y_log10() +
      theme_classic() +
      labs(y = "", x = environment_label) +
      theme(legend.position = "none",
            text         = element_text(size = 12),
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 13),
            legend.text  = element_text(size = 12),
            legend.title = element_text(size = 13)) +
          scale_colour_okabeito(order = c(3, 6)) +
          scale_fill_okabeito(order = c(3, 6))
  } else {
    b <- data_fit2 |>
      left_join(re_strength) %>%
      ggplot(aes(get(standard_name), log_value)) +
      geom_point(aes(colour = woody, group = woody), alpha = 0.75, col = "grey", data = data_model) +
      geom_line(aes(group = interaction(woody, taxon_name), colour = woody, alpha = 0.75), linewidth = 0.75) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted, group = woody), alpha = 0.8, colour = "grey20", linewidth = 1.5) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted, group = woody, colour = woody), linewidth = 0.9, linetype = "twodash") +
          geom_ribbon(data = data_c_line,
              aes(x = as.numeric(get(paste0(standard_name, "_c"))), ymin = est_Exp_low, ymax = est_Exp_hi, group = woody, fill = woody),
                      alpha = 0.3, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody, colour = woody), col = "black", linewidth = 1.5 * 1.5) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, group = woody, colour = woody), linewidth = 1.5) +
      # annotate("richtext", x = Inf, y = Inf, label = bayes_ratio_pasted,
      #          hjust = 1.1, vjust = 1, size = 5, label.size = NA, fill = NA) +
      scale_x_continuous(transform = x_transform) +
      theme_classic() +
      labs(y = "", x = environment_label) +
      theme(legend.position = "none",
            text         = element_text(size = 12),
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 13),
            legend.text  = element_text(size = 12),
            legend.title = element_text(size = 13)) +
          scale_colour_okabeito(order = c(3, 6)) +
          scale_fill_okabeito(order = c(3, 6))
  }

  # ---- Panel c: species-mean points + per-group interspecific fitted line + OLS lines ----
  if (trait_name != "leaf_delta13C") {
    c <- data_fit2 |>
      ggplot(aes(get(standard_name), exp(log_value), group = taxon_name)) +
      geom_point(alpha = 0.5, col = "grey", data = data_model) +
      geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
                 aes(x = get(standard_name), y = exp(log_value), colour = woody, group = woody), alpha = 0.75) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted), group = woody), alpha = 0.8, colour = "grey20", linewidth = 1.5) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted), group = woody, colour = woody), linewidth = 0.9, linetype = "twodash") +
      geom_ribbon(data = data_fit2_global,
          aes(x = as.numeric(get(standard_name)), ymin = value_low, ymax = value_hi, group = woody, fill = woody),
                  alpha = 0.3, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_fit2_global, aes(group = woody), linewidth = 1.5 * 1.5, colour = "black") +
      geom_line(data = data_fit2_global, aes(group = woody, colour = woody), linewidth = 1.5) +
      scale_x_continuous(transform = x_transform) +
      scale_y_log10() +
      theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none",
            text         = element_text(size = 12),
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 13),
            legend.text  = element_text(size = 12),
            legend.title = element_text(size = 13)) +
          scale_colour_okabeito(order = c(3, 6)) +
          scale_fill_okabeito(order = c(3, 6))
  } else {
    c <- data_fit2 |>
      ggplot(aes(get(standard_name), log_value, group = taxon_name)) +
      geom_point(alpha = 0.5, col = "grey", data = data_model) +
      geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
                 aes(x = get(standard_name), y = log_value, colour = woody, group = woody), alpha = 0.75) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted, group = woody), alpha = 0.8, colour = "grey20", linewidth = 1.5) +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted, group = woody, colour = woody), linewidth = 0.9, linetype = "twodash") +
          geom_ribbon(data = data_fit2_global,
              aes(x = as.numeric(get(standard_name)), ymin = value_low, ymax = value_hi, group = woody, fill = woody),
                      alpha = 0.3, colour = NA, inherit.aes = FALSE) +
      geom_line(data = data_fit2_global, aes(group = woody), linewidth = 1.5 * 1.5, colour = "black") +
      geom_line(data = data_fit2_global, aes(group = woody, colour = woody), linewidth = 1.5) +
      scale_x_continuous(transform = x_transform) +
      theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none",
            text         = element_text(size = 12),
            axis.text    = element_text(size = 12),
            axis.title   = element_text(size = 13),
            legend.text  = element_text(size = 12),
            legend.title = element_text(size = 13)) +
          scale_colour_okabeito(order = c(3, 6)) +
          scale_fill_okabeito(order = c(3, 6))
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
          text         = element_text(size = 12),
          axis.text    = element_text(size = 12),
          axis.title   = element_text(size = 13),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 13)) +
    scale_colour_okabeito(order = c(3, 6)) +
    scale_fill_okabeito(order = c(3, 6))

  return(list(a, c + b + d + plot_layout(nrow = 1)))
}

