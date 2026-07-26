fit_model2d_bayesian_3panel_ratio_only <- function(data_model, ind_variable = "log_annual_precipitation") {

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

  # Helper: extract inter- and intraspecific slope posteriors from a brms model.
  # 'mn_col' is the between-species (interspecific) fixed-effect column name;
  # 'c_col'  is the within-species  (intraspecific) fixed-effect column name.
  # Uses as_draws_df() so columns are selected by name rather than position.
  extract_slopes <- function(fit, mn_col, c_col, env_label) {
    posterior::as_draws_df(fit) %>%
      as_tibble() %>%
      transmute(inter = .data[[mn_col]],
                intra = .data[[c_col]],
                env   = env_label)
  }


  # Returns TRUE when the 95% credible interval spans zero.
  is_not_different_from_zero <- function(fit, hypothesis_string) {
    h <- brms::hypothesis(fit, hypothesis_string, alpha = 0.1)$hypothesis
    h$CI.Lower <= 0 && h$CI.Upper >= 0
  }

  # Group is marked non-significant only when BOTH inter and intra slopes include zero.
  both_slopes_not_different_from_zero <- function(fit, inter_hypothesis, intra_hypothesis) {
    is_not_different_from_zero(fit, inter_hypothesis) &&
      is_not_different_from_zero(fit, intra_hypothesis)
  }

  # Traits with a woody / non-woody split use the growth-form interaction model.
  # The remaining traits (huber_value, wood_density, leaf_P_per_area) are
  # woody-only so a single all-species model is used without a growth-form split.
  woody_split_traits <- c("leaf_mass_per_area", "leaf_N_per_area", "leaf_delta13C", "leaf_area","seed_dry_mass")

  if (trait_name %in% woody_split_traits) {

    # --- Woody-split traits ---
    # Model: log_value ~ (mn + c) * woody + (1 | taxon) + (0 + c | taxon)
    # Saved as list(fit_MAP_stan, fit_MAT_stan, data_model) [length 3].
    fits <- readRDS(paste0("outputs/model_fits/fit2d_bayesian_woody/", trait_name, ".RDS"))
    stopifnot(length(fits) == 3)
    fitmod_prec <- fits[[1]]  # MAP (precipitation) model
    fitmod_temp <- fits[[2]]  # MAT (temperature) model

    # For the growth-form interaction model, the brms fixed-effect columns are:
    #   b_<var>              — non-woody baseline slope
    #   b_<var>:woodywoody   — interaction shift for woody plants
    # The woody slope is therefore the sum of both terms.
    # NOTE: 'woodywoody' assumes the woody level is coded "woody" in the data;
    #       adjust if the factor level name differs.
    posterior_matrix_in_prec <- extract_slopes(
      fitmod_prec,
      mn_col = "b_std_log_annual_precipitation_mn",
      c_col  = "b_std_log_annual_precipitation_c",
      env_label = "MAP"
    ) %>%
      mutate(
        inter_woody = inter + as_tibble(posterior::as_draws_df(fitmod_prec))[["b_std_log_annual_precipitation_mn:woodywoody"]],
        intra_woody = intra + as_tibble(posterior::as_draws_df(fitmod_prec))[["b_std_log_annual_precipitation_c:woodywoody"]]
      )

    posterior_matrix_in_temp <- extract_slopes(
      fitmod_temp,
      mn_col = "b_std_average_temperature_mn",
      c_col  = "b_std_average_temperature_c",
      env_label = "MAT"
    ) %>%
      mutate(
        inter_woody = inter + as_tibble(posterior::as_draws_df(fitmod_temp))[["b_std_average_temperature_mn:woodywoody"]],
        intra_woody = intra + as_tibble(posterior::as_draws_df(fitmod_temp))[["b_std_average_temperature_c:woodywoody"]]
      )

    # Also load the all-species model (no growth-form split) as a combined reference.
    fits_combined <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))

    # Flag panels/groups where both inter and intra slopes are not different from zero.
    no_fill_flags <- bind_rows(
      tibble(
        env = "MAP",
        name = "diff",
        no_fill = both_slopes_not_different_from_zero(
          fitmod_prec,
          "std_log_annual_precipitation_mn = 0",
          "std_log_annual_precipitation_c = 0"
        )
      ),
      tibble(
        env = "MAT",
        name = "diff",
        no_fill = both_slopes_not_different_from_zero(
          fitmod_temp,
          "std_average_temperature_mn = 0",
          "std_average_temperature_c = 0"
        )
      ),
      tibble(
        env = "MAP",
        name = "diff_woody",
        no_fill = both_slopes_not_different_from_zero(
          fitmod_prec,
          "std_log_annual_precipitation_mn + std_log_annual_precipitation_mn:woodywoody = 0",
          "std_log_annual_precipitation_c + std_log_annual_precipitation_c:woodywoody = 0"
        )
      ),
      tibble(
        env = "MAT",
        name = "diff_woody",
        no_fill = both_slopes_not_different_from_zero(
          fitmod_temp,
          "std_average_temperature_mn + std_average_temperature_mn:woodywoody = 0",
          "std_average_temperature_c + std_average_temperature_c:woodywoody = 0"
        )
      ),
      tibble(
        env = "MAP",
        name = "diff_combined",
        no_fill = both_slopes_not_different_from_zero(
          fits_combined[[1]],
          "std_log_annual_precipitation_mn = 0",
          "std_log_annual_precipitation_c = 0"
        )
      ),
      tibble(
        env = "MAT",
        name = "diff_combined",
        no_fill = both_slopes_not_different_from_zero(
          fits_combined[[2]],
          "std_average_temperature_mn = 0",
          "std_average_temperature_c = 0"
        )
      )
    )

    # Compute diff_combined: pointwise intraspecific ratio from the all-species model.
    # Appended after the pivot so it can share the same env/name/value structure.
    posterior_matrix_combined <- bind_rows(
      extract_slopes(fits_combined[[1]], "b_std_log_annual_precipitation_mn",
                     "b_std_log_annual_precipitation_c", "MAP"),
      extract_slopes(fits_combined[[2]], "b_std_average_temperature_mn",
                     "b_std_average_temperature_c", "MAT")
    ) %>%
      transmute(
        env,
        diff_combined = abs(intra) / (abs(intra) + abs(inter))
      ) %>%
      pivot_longer(cols = -env)

    # Compute the intraspecific ratio for each draw and each growth-form group:
    #   diff       = |intra|       / (|intra|       + |inter|)       — non-woody
    #   diff_woody = |intra_woody| / (|intra_woody| + |inter_woody|) — woody
    # Computing pointwise across draws correctly propagates posterior uncertainty
    # through this nonlinear transformation.
    data_for_plot <- bind_rows(posterior_matrix_in_prec, posterior_matrix_in_temp) %>%
      mutate(
        diff       = abs(intra)       / (abs(intra)       + abs(inter)),
        diff_woody = abs(intra_woody) / (abs(intra_woody) + abs(inter_woody))
      ) %>%
      pivot_longer(cols = -env) %>%
      bind_rows(posterior_matrix_combined)

    cat("extracting fitted vals -> ")

    # Summarise each quantity's posterior: median + 95% credible interval.
    # The ratio is bounded [0, 1] by construction, so no capping is needed.
    data_for_plot_group <- data_for_plot %>%
      group_by(env, name) %>%
      summarise(
        mean_value   = median(value),
        upp_quantile = quantile(value, 0.975),
        low_quantile = quantile(value, 0.025),
        .groups = "drop"
      ) %>%
      mutate(woody = case_when(
        str_detect(name, "combined") ~ "c",
        str_detect(name, "woody")    ~ "a",
        TRUE                         ~ "b"
      ))

    # x-positions for each quantity on the coord_flip axis
    spacing <- c(
      inter         = 1.0, intra         = 1.3, diff         = 1.6,
      inter_woody   = 1.1, intra_woody   = 1.4, diff_woody   = 1.7,
      diff_combined = 1.8
    )

    data_for_plot       <- data_for_plot       %>% mutate(xpos = spacing[name])
    data_for_plot_group <- data_for_plot_group %>% mutate(xpos = spacing[name])

    selected_names <- c("diff", "diff_woody", "diff_combined")
    y_axis_label <- bquote(B[r])
    ref_line <- 0.5

    # Keep only the selected posterior quantity for plotting.
    data_for_plot       <- data_for_plot       %>% filter(name %in% selected_names)
    data_for_plot_group <- data_for_plot_group %>% filter(name %in% selected_names)
    data_for_plot <- data_for_plot %>%
      left_join(no_fill_flags, by = c("env", "name")) %>%
      mutate(no_fill = if_else(is.na(no_fill), FALSE, no_fill),
             fill_status = if_else(no_fill, "no_fill", "fill"))

    data_for_plot %>%
      mutate(woody = case_when(
        str_detect(name, "combined") ~ "c",
        str_detect(name, "woody")    ~ "a",
        TRUE                         ~ "b"
      )) %>%
      ggplot(aes(y = value, x = xpos,
                 group = interaction(name, woody), colour = woody, fill = woody)) +
      geom_hline(yintercept = ref_line, linetype = "dashed") +
      geom_violinhalf(aes(alpha = fill_status), flip = FALSE) +
      geom_segment(data = data_for_plot_group,
                   aes(y = low_quantile, yend = upp_quantile, x = xpos),
                   colour = "black", size = 1.2) +
      geom_point(data = data_for_plot_group, aes(y = mean_value), colour = "black") +
      # Median label placed above the violin half (positive x nudge in coord_flip = upward)
      geom_text(data = data_for_plot_group,
                aes(y = mean_value, label = sprintf("%.2f", mean_value)),
                colour = "black", size = 3, position = position_nudge(x = 0.08)) +
      facet_grid(~env) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12),
            axis.text = element_text(size = 10.5),
            axis.text.y = element_blank(),
            strip.background = element_rect(fill = "white", colour = NA),
            strip.text = element_text(size = 16)) +
          labs(y = y_axis_label) +
      xlab(trait_label) +
          scale_y_continuous(limits = c(0, 1)) +
      scale_colour_okabeito(order = c(6, 3, 5), breaks = c("a", "b", "c")) +
      scale_fill_okabeito(order = c(6, 3, 5),   breaks = c("a", "b", "c")) +
        scale_alpha_manual(values = c("fill" = 0.5, "no_fill" = 0)) -> a

    rm(list = c(
      "fits", "fits_combined", "fitmod_prec", "fitmod_temp",
      "posterior_matrix_in_prec", "posterior_matrix_in_temp", "posterior_matrix_combined"
    ), inherits = FALSE)
    gc(verbose = FALSE)
    return(a)

  } else {

    # --- Woody-only traits (huber_value, wood_density, leaf_P_per_area) ---
    # These traits occur only in woody species, so no growth-form split is needed.
    # Uses the all-species Bayesian model from fit2d_bayesian/.
    # Saved as list(fit_MAP_stan, fit_MAT_stan, data_model) [length 3].
    fits <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))

    fitmod_prec <- fits[[1]]  # MAP (precipitation) model
    fitmod_temp <- fits[[2]]  # MAT (temperature) model

    # Flag panels where both inter and intra slopes are not different from zero.
    no_fill_flags <- bind_rows(
      tibble(
        env = "MAP",
        name = "diff",
        no_fill = both_slopes_not_different_from_zero(
          fitmod_prec,
          "std_log_annual_precipitation_mn = 0",
          "std_log_annual_precipitation_c = 0"
        )
      ),
      tibble(
        env = "MAT",
        name = "diff",
        no_fill = both_slopes_not_different_from_zero(
          fitmod_temp,
          "std_average_temperature_mn = 0",
          "std_average_temperature_c = 0"
        )
      )
    )

    # Extract inter- and intraspecific slope posteriors from each climate model.
    posterior_matrix_in_prec <- extract_slopes(
      fitmod_prec,
      mn_col    = "b_std_log_annual_precipitation_mn",
      c_col     = "b_std_log_annual_precipitation_c",
      env_label = "MAP"
    )

    posterior_matrix_in_temp <- extract_slopes(
      fitmod_temp,
      mn_col    = "b_std_average_temperature_mn",
      c_col     = "b_std_average_temperature_c",
      env_label = "MAT"
    )

    # Compute the intraspecific ratio for each posterior draw:
    #   diff = |intra| / (|intra| + |inter|)  ∈ [0, 1]
    # Computing pointwise on draws correctly propagates posterior uncertainty
    # through this nonlinear transformation.
    data_for_plot <- bind_rows(posterior_matrix_in_prec, posterior_matrix_in_temp) %>%
      mutate(
        diff = abs(intra) / (abs(intra) + abs(inter))
      ) %>%
      pivot_longer(cols = -env)

    cat("extracting fitted vals -> ")

    # Summarise the posterior: median + 95% credible interval.
    # The ratio is bounded [0, 1] by construction, so no capping is needed.
    data_for_plot_group <- data_for_plot %>%
      group_by(env, name) %>%
      summarise(
        mean_value   = median(value),
        upp_quantile = quantile(value, 0.975),
        low_quantile = quantile(value, 0.025),
        .groups = "drop"
      )

    # x-positions for each quantity on the coord_flip axis
    spacing <- c(inter = 1.0, intra = 1.3, diff = 1.6)

    data_for_plot       <- data_for_plot       %>% mutate(xpos = spacing[name])
    data_for_plot_group <- data_for_plot_group %>% mutate(xpos = spacing[name])

    selected_name <- "diff"
    y_axis_label <- bquote(B[r])
    ref_line <- 0.5

    # Keep only the selected posterior quantity for plotting; drop raw inter/intra draws.
    data_for_plot       <- data_for_plot       %>% filter(name == selected_name)
    data_for_plot_group <- data_for_plot_group %>% filter(name == selected_name)
    data_for_plot <- data_for_plot %>%
      left_join(no_fill_flags, by = c("env", "name")) %>%
      mutate(no_fill = if_else(is.na(no_fill), FALSE, no_fill),
             fill_status = if_else(no_fill, "no_fill", "fill"))

    data_for_plot %>%
      ggplot(aes(y = value, x = xpos,
                 group = interaction(name), colour = "a", fill = "a")) +
      geom_hline(yintercept = ref_line, linetype = "dashed") +
      geom_violinhalf(aes(alpha = fill_status), flip = FALSE, width = 0.09,
                      position = position_nudge(x = 0.08)) +
      geom_segment(data = data_for_plot_group,
                   aes(y = low_quantile, yend = upp_quantile, x = xpos),
                   colour = "black", size = 1.2,
                   position = position_nudge(x = 0.08)) +
      geom_point(data = data_for_plot_group, aes(y = mean_value),
                 colour = "black", position = position_nudge(x = 0.08)) +
      # Median label placed above the violin half; nudge further than the violin (0.08) to avoid overlap
      geom_text(data = data_for_plot_group,
                aes(y = mean_value, label = sprintf("%.2f", mean_value)),
                colour = "black", size = 3, position = position_nudge(x = 0.17)) +
      facet_grid(~env) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12),
            axis.text = element_text(size = 10.5),
            axis.text.y = element_blank(),
            strip.background = element_rect(fill = "white", colour = NA),
            strip.text = element_text(size = 16)) +
      scale_x_continuous(limits = c(1.55, 1.85)) +
          labs(y = y_axis_label) +
      xlab(trait_label) +
          scale_y_continuous(limits = c(0, 1)) +
        scale_alpha_manual(values = c("fill" = 0.5, "no_fill" = 0)) +
      scale_colour_okabeito(order = c(6), breaks = c("a")) +
      scale_fill_okabeito(order = c(6),   breaks = c("a")) -> a

    rm(list = c(
      "fits", "fitmod_prec", "fitmod_temp",
      "posterior_matrix_in_prec", "posterior_matrix_in_temp"
    ), inherits = FALSE)
    gc(verbose = FALSE)
    return(a)
  }
}

