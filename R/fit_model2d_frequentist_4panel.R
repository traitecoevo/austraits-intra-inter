# Produces a 2-panel figure summarising trait responses to a single climate driver
# across all species (no woody/non-woody split).
#
# Returns list(a, c + b + d) where:
#   a — Raw observations + OLS regression line
#   b — Per-species mixed-model fitted lines (intraspecific component) +
#         average intraspecific slope line (data_c_line) + Bᵣ annotation
#   c — Species-mean points + interspecific fitted line (data_fit2_global) +
#         OLS line for visual comparison
#   d — Density of per-species intraspecific slopes; vertical lines at
#         OLS overall slope (grey) and fixed intraspecific effect β_c (yellow)
#
# Bᵣ (Bayesian ratio) = |β_c| / (|β_mn| + |β_c|) from the Bayesian brms model.
# β_mn = interspecific (between-species) slope; β_c = intraspecific (within-species) slope.
# Asterisk (*) appended when β_mn and β_c have opposite signs.
#
# Arguments:
#   data_model     — data frame for one trait; must contain trait_name column
#   ind_variable   — climate driver: "log_annual_precipitation" (default) or
#                    "average_temperature"
#   woody_colouring — if TRUE, use trait-specific okabeito colour (for wood_density /
#                     huber_value rows in combined woody figures)
fit_model2d_frequentist_4panel <- function(data_model, ind_variable = "log_annual_precipitation", woody_colouring = FALSE) {
  # --- Derive helper variables from ind_variable ---
  # f:             back-transform function (exp for log-variables, identity for temperature)
  # standard_name: untransformed variable name used for axis labels and column naming
  # x_transform:   ggplot scale transform for the x-axis
  # off_variable:  the climate variable held fixed (opposite to ind_variable)
  if (str_detect(ind_variable, "log")) { f <- exp } else { f <- I }
  if (str_detect(ind_variable, "log")) { standard_name <- str_remove(ind_variable, "log_") } else { standard_name <- ind_variable }
  if (str_detect(ind_variable, "log")) { x_transform <- "log10" } else { x_transform <- "identity" }
  off_variable <- if (ind_variable == "log_annual_precipitation") "average_temperature" else "log_annual_precipitation"

  # Axis labels: only populated for wood_density (bottom row in combined multi-panel figures)
  environment_label <- ""
  b1_label          <- ""
  if (data_model$trait_name[1] == "wood_density") {
    environment_label <- if (ind_variable == "log_annual_precipitation") "MAP (mm)" else expression(paste(MAT, " (", degree * C, ")"))
    b1_label          <- if (ind_variable == "log_annual_precipitation") expression(paste("Intraspecific response to MAP")) else expression(paste("Intraspecific response to MAT"))
  }

  # Trait y-axis labels
  trait_name <- data_model$trait_name[1]
  if (trait_name == "leaf_mass_per_area") trait_label <- expression(paste(LMA,            " (g ",  m^{-2},   ")"))
  if (trait_name == "huber_value")        trait_label <- expression(paste(SA:LA,          " (", mm^{2}, cm^{-2},  ")"))
  if (trait_name == "leaf_area")          trait_label <- expression(paste(Leaf ~ area,    " (",    mm^{-2},  ")"))
  if (trait_name == "leaf_delta13C")      trait_label <- expression(paste(C[i],           ":",     C[a]))
  if (trait_name == "leaf_N_per_area")    trait_label <- expression(paste(Leaf ~ N ~ per ~ area, " (g ", m^{-2}, ")"))
  if (trait_name == "wood_density")       trait_label <- expression(paste(Wood ~ density, " (g ",  cm^{-3}, ")"))
  if (trait_name == "seed_dry_mass")    trait_label <- expression(paste(Seed~mass,        " (mg)"))

  # --- Load frequentist model fit ---
  # Length-2 RDS: list(fitmod_MAP_or_combined, data_model)
  # Length-3 RDS: list(fitmod_MAP, fitmod_MAT, data_model)


  fits           <- readRDS(paste0("outputs/model_fits/fit2d_frequentist/", trait_name, ".RDS"))

  # --- Compute or load Bayesian ratio (Bᵣ) ---
  # Bᵣ = |β_c| / (|β_mn| + |β_c|): relative contribution of intraspecific variation.
  # Extracted from the brms model fixed effects (rows 2 = β_mn, 3 = β_c).
  # Result is cached to disk after the first computation.
  bayes_ratio_path <- paste0("outputs/model_fits/fit2d_bayesian/bayes_ratios/", trait_name, "_", ind_variable, ".RDS")
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

  # Select the correct fitted model and recover the data used for fitting
  if (length(fits) == 2) {
    fitmod     <- fits[[1]]
    data_model <- fits[[2]]
  }
  if (length(fits) == 3) {
    fitmod     <- if (ind_variable == "log_annual_precipitation") fits[[1]] else fits[[2]]
    data_model <- fits[[3]]
  }

  # For huber_value the frequentist model has convergence issues; substitute the brms
  # model for all predictions used in panels b, c, and d. Panel a (OLS) is unaffected.
  use_bayes_fitmod <- trait_name == "huber_value"
  # use_bayes_fitmod <- trait_name == "asd"

  if (use_bayes_fitmod) {
    fits_bayes <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))
    fitmod     <- if (ind_variable == "log_annual_precipitation") fits_bayes[[1]] else fits_bayes[[2]]
    data_model <- fits_bayes[[3]]
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

  # Fix the off-variable at each species' mean to isolate the ind_variable response
  data_for_fit <- data_model %>%
    select(std_log_annual_precipitation_mn, std_log_annual_precipitation_c,
           std_average_temperature_mn,      std_average_temperature_c,
           taxon_name) %>%
    group_by(taxon_name) %>%
    mutate(
      !!paste0("std_", off_variable, "_c")  := mean(get(paste0("std_", off_variable, "_c")),  na.rm = TRUE),
      !!paste0("std_", off_variable, "_mn") := mean(get(paste0("std_", off_variable, "_mn")), na.rm = TRUE)
    ) %>%
    drop_na()

  # Per-species fitted values (with random effects) over the observed climate range
  data_fit2_re <- estimate_expectation(fitmod, data = data_for_fit) |> as_tibble()


  # Interspecific prediction: population-level fitted values varying only the
  # species-mean climate (_mn), with _c fixed at 0.
  # For brms: estimate_means only uses a few grid points for continuous predictors,
  # giving a near-flat line. Instead, build a 100-point grid and use fitted() with
  # re_formula = NA to get true population-level (no random effects) predictions.
  if (use_bayes_fitmod) {
    brms_grid <- tibble(
      !!paste0("std_", ind_variable, "_mn") := seq(
        min(data_for_fit[[paste0("std_", ind_variable, "_mn")]]),
        max(data_for_fit[[paste0("std_", ind_variable, "_mn")]]),
        length.out = 100
      ),
      !!paste0("std_", ind_variable, "_c")  := 0,
      !!paste0("std_", off_variable, "_mn") := mean(data_for_fit[[paste0("std_", off_variable, "_mn")]], na.rm = TRUE),
      !!paste0("std_", off_variable, "_c")  := 0
    )
    data_fit2_global <- brms_grid |>
      mutate(
        log_value  = fitted(fitmod, newdata = brms_grid, re_formula = NA)[, "Estimate"],
        !!standard_name := as.numeric(f(get(paste0("std_", ind_variable, "_mn")) * scale_mn + mean_mn)),
        value      = exp(log_value),
        taxon_name = "global average"
      )
  } else {
    data_fit2_global <- estimate_means(fitmod, by = paste0("std_", ind_variable, "_mn")) |>
      as_tibble() |>
      mutate(log_value = Mean) |>
      select(!!(paste0("std_", ind_variable)) := paste0("std_", ind_variable, "_mn"), log_value) |>
      mutate(
        !!standard_name := as.numeric(f(get(paste0("std_", ind_variable)) * scale_mn + mean_mn)),
        value      = exp(log_value),
        taxon_name = "global average"
      )
  }

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



  # Build a 2-row prediction dataset for the average intraspecific slope line.
  # Fix _mn at the median species' climate position; set _c to span the full
  # observed range (min and max of _mn + _c combined).
  # This visually represents how a typical species' trait changes with within-species
  # climate variation, extrapolated across the full data range.
  if (length(fits) > 2) {
    data_c_line <- data_fit2_re %>%
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
  data_c_line$est_Exp <- estimate_expectation(fitmod, data_c_line)$Predicted
  if (trait_name != "leaf_delta13C") data_c_line$est_Exp <- exp(data_c_line$est_Exp)

  # Back-transform x-axis: standard_name_c = f(ind_variable_c) gives actual environmental units
  data_c_line <- data_c_line |>
    mutate(!!paste0(standard_name, "_c") := f(get(paste0(ind_variable, "_c"))))

  # Combine per-species fitted values with their species-mean predictions.
  # Adds back-transformed climate columns for plotting on the original scale.
  data_fit2 <- data_fit2_re |>
    rename(log_value = Predicted) |>
    left_join(
      data_fit2_mn |> select(paste0("std_", ind_variable, "_mn"), taxon_name, log_value_mn = Predicted)
    ) |>
    mutate(
      !!ind_variable                 := (get(paste0("std_", ind_variable, "_mn")) + get(paste0("std_", ind_variable, "_c"))) * scale_c + mean_c,
      !!standard_name                := as.numeric(f(get(paste0(ind_variable)))),
      !!paste0(standard_name, "_mn") := as.numeric(f(get(paste0("std_", ind_variable, "_mn")) * scale_c + mean_c))
    )

  cat("cooking up a sweet plot -> ")

  # Back-transform data_model climate columns to original scale for plotting
  data_model <- data_model %>%
    mutate(
      log_annual_precipitation = std_log_annual_precipitation * log_annual_precipitation_scaler_scale + log_annual_precipitation_scaler_center,
      annual_precipitation     = exp(std_log_annual_precipitation * log_annual_precipitation_scaler_scale + log_annual_precipitation_scaler_center),
      average_temperature      = std_average_temperature * average_temperature_scaler_scale + average_temperature_scaler_center
    )

  # --- OLS overall regression (naive, ignoring species identity) ---
  # Used as a reference line: represents the combined inter+intraspecific slope.
  if (length(fits) > 2) {
    if (ind_variable == "log_annual_precipitation") {
      fit_overall <- lm(log_value ~ std_log_annual_precipitation, data = data_model, na.action = na.omit)
      print(summary(fit_overall))
      doc <- read_docx()
      doc <- body_add_table(doc, value = as.data.frame(summary(fit_overall)$coef), style = "table_template")
      print(doc, target = paste0(trait_name, "_MAP_overall.docx"))
    }
    if (ind_variable == "average_temperature") {
      fit_overall <- lm(log_value ~ std_average_temperature, data = data_model, na.action = na.omit)
      print(summary(fit_overall))
      doc <- read_docx()
      doc <- body_add_table(doc, value = as.data.frame(summary(fit_overall)$coef), style = "table_template")
      print(doc, target = paste0(trait_name, "_MAT_overall.docx"))
    }
  } else {
    fit_overall <- lm(log_value ~ std_log_annual_precipitation + std_average_temperature, data = data_model, na.action = na.omit)
    print(summary(fit_overall))
  }

  # OLS predictions for plotting (marginalised over the focal climate variable)
  data_fit2_re_overall <- estimate_expectation(fit_overall, by = paste0("std_", ind_variable)) |>
    as_tibble() |>
    mutate(!!ind_variable := f(get(paste0("std_", ind_variable)) * scale_c + mean_c))

  # Per-species intraspecific slopes: species-specific slope for std_X_c
  # (fixed effect + random deviation), used for the density plot (panel d)
  # brms stores random effects as a 3D array [species, statistic, coefficient];
  # glmmTMB stores them under $cond$taxon_name.
  if (use_bayes_fitmod) {
    re_strength <- tibble(
      taxon_name  = dimnames(coef(fitmod)$taxon_name)[[1]],
      re_strength = coef(fitmod)$taxon_name[, "Estimate", paste0("std_", ind_variable, "_c")]
    )
  } else {
    re_strength <- tibble(
      taxon_name  = rownames(coef(fitmod)$cond$taxon_name),
      re_strength = coef(fitmod)$cond$taxon_name[[paste0("std_", ind_variable, "_c")]]
    )
  }

  # ---- Panel a: raw data + OLS line ----
  if (trait_name != "leaf_delta13C") {
    if (woody_colouring) {
      a <- data_model |>
        mutate(!!ind_variable := as.numeric(f(get(ind_variable)))) |>
        ggplot(aes(x = get(ind_variable), y = exp(log_value))) +
        geom_point(aes(colour = trait_name), alpha = 0.75) +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 1.5 * 1.5, colour = "black") +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted), colour = trait_name), linewidth = 1.5) +
        scale_x_continuous(transform = x_transform) + scale_y_log10() + theme_classic() +
        labs(y = trait_label, x = environment_label) +
        theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
        scale_colour_okabeito(order = c(6))
    } else {
      a <- data_model |>
        mutate(!!ind_variable := as.numeric(f(get(ind_variable)))) |>
        ggplot(aes(x = get(ind_variable), y = exp(log_value))) +
        geom_point(alpha = 0.75, col = "grey") +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 1.5) +
        scale_x_continuous(transform = x_transform) + scale_y_log10() + theme_classic() +
        labs(y = trait_label, x = environment_label) +
        theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    }
  } else {
    a <- data_model |>
      mutate(!!ind_variable := as.numeric(f(get(ind_variable)))) |>
      ggplot(aes(x = get(ind_variable), y = log_value)) +
      geom_point(alpha = 0.75, col = "grey") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 1.5) +
      scale_x_continuous(transform = x_transform) + theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
  }

  # ---- Panel b: per-species fitted lines (intraspecific component) + average intraspecific slope + Bᵣ ----
  if (trait_name != "leaf_delta13C") {
    if (woody_colouring) {
      bayes_ratio_pasted <- paste0("Bᵣ = ", unlist(bayes_ratio))
      bayes_ratio_pasted <- paste0("<span style='color:#D55E00;'>", bayes_ratio_pasted[1], "</span>")
      b <- data_fit2 |>
        left_join(re_strength) %>%
        ggplot(aes(get(standard_name), exp(log_value), group = taxon_name)) +
        geom_point(alpha = 0.75, col = "grey", data = data_model) +
        geom_line(aes(group = taxon_name, alpha = 0.75, colour = trait_name), linewidth = 0.75) +
        geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp), col = "black", linewidth = 1.5 * 1.5) +
        geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp, colour = trait_name), linewidth = 1.5) +
        # annotate("richtext", x = Inf, y = Inf, label = bayes_ratio_pasted, hjust = 1.1, vjust = 1, size = 5, label.size = NA, fill = NA) +
        scale_x_continuous(transform = x_transform) + scale_y_log10() +
        scale_colour_okabeito(order = c(6)) + theme_classic() +
        labs(y = "", x = environment_label) +
        theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    } else {
      b <- data_fit2 |>
        left_join(re_strength) %>%
        ggplot(aes(get(standard_name), exp(log_value))) +
        geom_point(alpha = 0.75, col = "grey", data = data_model) +
        geom_line(aes(group = taxon_name, alpha = 0.75), linewidth = 0.75, col = "#0072B2") +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 1.5,  colour = "black") +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 0.9, colour = "grey60") +
        geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp), colour = "black",   linewidth = 1.5) +
        geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp), colour = "#0072B2", linewidth = 0.9) +
        # geom_text(data = data.frame(x = Inf, y = Inf), aes(x = x, y = y, label = label),
                  # label = paste0("Bᵣ = ", bayes_ratio), hjust = 1.1, vjust = 1, size = 5, parse = FALSE, inherit.aes = FALSE) +
        scale_x_continuous(transform = x_transform) + scale_y_log10() + theme_classic() +
        labs(y = "", x = environment_label) +
        theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    }
  } else {
    b <- data_fit2 |>
      left_join(re_strength) %>%
      ggplot(aes(get(standard_name), log_value)) +
      geom_point(alpha = 0.75, col = "grey", data = data_model) +
      geom_line(aes(group = taxon_name, alpha = 0.75), linewidth = 0.75, col = "#0072B2") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 1.5,  colour = "black") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 0.9, colour = "grey60") +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp), colour = "black",   linewidth = 1.5) +
      geom_line(data = data_c_line, aes(x = as.numeric(get(paste0(standard_name, "_c"))), y = est_Exp), colour = "#0072B2", linewidth = 0.9) +
      # geom_text(data = data.frame(x = Inf, y = Inf), aes(x = x, y = y, label = label),
                # label = paste0("Bᵣ = ", bayes_ratio), hjust = 1.1, vjust = 1, size = 5, parse = FALSE, inherit.aes = FALSE) +
      scale_x_continuous(transform = x_transform) + theme_classic() +
      labs(y = "", x = environment_label) +
      theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
  }

  # ---- Panel c: species-mean points + interspecific fitted line + OLS line ----
  if (trait_name != "leaf_delta13C") {
    if (woody_colouring) {
      c <- data_fit2 |>
        ggplot(aes(get(standard_name), exp(log_value))) +
        geom_point(alpha = 0.5, col = "grey", data = data_model) +
        geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
                   aes(x = get(standard_name), y = exp(log_value), colour = trait_name), alpha = 0.75) +
        geom_line(data = data_fit2_global, col = "black",   linewidth = 1.5 * 1.5) +
        geom_line(data = data_fit2_global, col = "#D55E00", linewidth = 1.5) +
        scale_x_continuous(transform = x_transform) + scale_y_log10() +
        scale_colour_okabeito(order = c(6)) + theme_classic() +
        labs(y = trait_label, x = environment_label) +
        theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    } else {
      c <- data_fit2 |>
        ggplot(aes(get(standard_name), exp(log_value))) +
        geom_point(alpha = 0.5, col = "grey", data = data_model) +
        geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
                   aes(x = get(standard_name), y = exp(log_value)), alpha = 0.75, col = "#E69F00") +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 1.5, colour = "black") +
        geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = exp(Predicted)), linewidth = 0.9, colour = "grey60") +
        geom_line(data = data_fit2_global, linewidth = 1.5, colour = "black") +
        geom_line(data = data_fit2_global, linewidth = 0.9, colour = "#E69F00") +
        scale_x_continuous(transform = x_transform) + scale_y_log10() + theme_classic() +
        labs(y = trait_label, x = environment_label) +
        theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    }
  } else {
    c <- data_fit2 |>
      ggplot(aes(get(standard_name), log_value)) +
      geom_point(alpha = 0.5, col = "grey", data = data_model) +
      geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE),
                 aes(x = get(standard_name), y = log_value), alpha = 0.75, col = "#E69F00") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 1.5, colour = "black") +
      geom_line(data = data_fit2_re_overall, aes(x = as.numeric(get(ind_variable)), y = Predicted), linewidth = 0.9, colour = "grey60") +
      geom_line(data = data_fit2_global, linewidth = 1.5, colour = "black") +
      geom_line(data = data_fit2_global, linewidth = 0.9, colour = "#E69F00") +
      scale_x_continuous(transform = x_transform) + theme_classic() +
      labs(y = trait_label, x = environment_label) +
      theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
  }

  # ---- Panel d: density of per-species intraspecific slopes ----
  # Vertical lines: grey = OLS overall slope; yellow/orange = fixed β_c (intraspecific effect)
  if (woody_colouring) {
    d <- re_strength |>
      ggplot(aes(x = re_strength)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_density(linewidth = 1.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.1) +
      scale_x_continuous(limits = c(-max(abs(re_strength$re_strength)), max(abs(re_strength$re_strength)))) +
      theme_classic() +
      labs(x = b1_label, y = "Density") +
      theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
  } else {
    d <- ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_density(data = re_strength, aes(x = re_strength), linewidth = 1) +
      geom_vline(xintercept = coef(fit_overall)[2], col = "grey60", linewidth = 1) +
      geom_vline(
        xintercept = if (use_bayes_fitmod)
          fixef(fitmod)[paste0("std_", ind_variable, "_mn"), "Estimate"]
        else
          fixef(fitmod)$cond[paste0("std_", ind_variable, "_mn")],
        col = "#E69F00", linewidth = 1.25
      ) +
      scale_x_continuous(limits = c(-max(abs(re_strength$re_strength)) * 1.1, max(abs(re_strength$re_strength)) * 1.1)) +
      theme_classic() +
      labs(x = b1_label, y = "Density") +
      theme(legend.position = "none", text = element_text(size = 12), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
  }

  return(list(a, c + b + d + plot_layout(nrow = 1)))
}
