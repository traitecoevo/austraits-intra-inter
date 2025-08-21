
fit_model_1D <- function(data_model, seed=123) {
  cat("fitting model -> ")
  # model across species means
  fit2 <- glmmTMB(log_value ~ log_annual_precipitation_mn + log_annual_precipitation_c + (log_annual_precipitation_c | taxon_name), data = data_model)

  cat("extracting fitted vals -> ")

  # overall fitted - with REs
  data_fit2_re <- estimate_expectation(fit2) |>
    as_tibble()



  # Global prediction - change in species mean trait with climate
  data_fit2_global <- estimate_means(fit2, by = "log_annual_precipitation_mn") |>
    as_tibble() |>
    # change variable names so these will plot ok
    select(log_annual_precipitation = log_annual_precipitation_mn, log_value = Mean) |>
    # backtransform so we can plot on raw scale
    mutate(
      annual_precipitation = exp(log_annual_precipitation),
      value = exp(log_value),
      # need a taxon name for plotting
      taxon_name = "global average"
    )

  # caluclate species means by removing devistations from center
  data_fit2_mn <- estimate_expectation(fit2,
    data = data_fit2_re |> mutate(log_annual_precipitation_c = 0) |> select(-Predicted) |> distinct()
  ) |>
    as_tibble() |>
    distinct()

  # combine
  data_fit2 <-
    data_fit2_re |>
    rename(log_value = Predicted) |>
    select(-SE, -Residuals) |>
    left_join(
      by = c("log_annual_precipitation_mn", "taxon_name"),
      data_fit2_mn |> select(log_annual_precipitation_mn, taxon_name, log_value_mn = Predicted)
    ) |>
    mutate(
      log_annual_precipitation = log_annual_precipitation_mn + log_annual_precipitation_c,
      annual_precipitation = exp(log_annual_precipitation),
      annual_precipitation_mn = exp(log_annual_precipitation_mn),
      value = exp(log_value),
      value_mn = exp(log_value_mn),
    )

  cat("coooking up a sweet plot -> ")

  # Choose a subset of taxa to plot. We'll choose randomly, but restrict to those with multiple sites and records
  set.seed(seed)
  taxa <- data_model |>
    group_by(taxon_name) |>
    summarise(
      n = n(),
      n_sites = n_distinct(annual_precipitation)
    ) |>
    ungroup() |>
    filter(n > 4, n_sites > 4) |>
    pull(taxon_name) |>
    unique()

  taxa <- taxa |> sample(min(30, length(taxa)))

  data_model_taxa <- data_model |>
    filter(taxon_name %in% taxa) |>
    mutate(taxon2 = taxon_name)

  data_fit2_taxa <- data_fit2 |>
    filter(taxon_name %in% taxa) |>
    mutate(taxon2 = taxon_name)

  p_fit <-
    data_fit2 |>
    ggplot(aes(annual_precipitation, value, group = taxon_name)) +
    # raw data
    geom_point(alpha = 0.5, col = "grey", data = data_model, stroke = 0) +
    # fitted lines
    geom_line(alpha = 0.5, col = "black", linewidth = 0.25) +
    # highlight certain species
    geom_point(aes(group = taxon2), col = "orange", show.legend = FALSE, data = data_model_taxa, size = 0.9, stroke = 0) +
    geom_line(aes(group = taxon2), col = "#d95f02", alpha = 1, show.legend = FALSE, data = data_fit2_taxa) +
    # global line
    geom_line(data = data_fit2_global, col = "red") +
    # species centers
    geom_point(aes(x = annual_precipitation_mn, y = value_mn, group = taxon2), col = "red", alpha = 1, size = 0.9, stroke = 0, data = data_fit2_taxa) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(y = "trait", x = "annual rainfall") +
    theme(legend.position = "none") +
    facet_wrap(~taxon2)
  cat("done\n")

  list(
    data = data_model,
    fit = fit2,
    data_fit = data_fit2,
    plot_fitted = p_fit
  )

}
