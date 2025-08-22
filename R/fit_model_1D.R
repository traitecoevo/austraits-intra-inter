
fit_model_1D <- function(data_model, seed=123, ind_variable = "log_annual_precipitation") {
  cat("fitting model -> ")
  # model across species means
  
  if(str_detect(ind_variable,"log")){f = exp}else{f = I}
  if(str_detect(ind_variable,"log")){standard_name = str_remove(ind_variable, "log_")}else{standard_name = ind_variable}
  if(str_detect(ind_variable,"log")){x_transform = "log10"}else{x_transform = "identity"}
  
  
  form <- as.formula(paste0("log_value ~ ",ind_variable,"_mn + ",ind_variable,"_c + (",ind_variable,"_c | taxon_name)"))
  
  fit2 <- glmmTMB(form, data = data_model)

  cat("extracting fitted vals -> ")

  # overall fitted - with REs
  data_fit2_re <- estimate_expectation(fit2) |>
    as_tibble()



  # Global prediction - change in species mean trait with climate
  data_fit2_global <- estimate_means(fit2, by = paste0(ind_variable,"_mn")) |>
    as_tibble() |>
    # change variable names so these will plot ok
    select(!!ind_variable := paste0(ind_variable,"_mn"), log_value = Mean) |>
    # backtransform so we can plot on raw scale
    mutate(
      !!standard_name := as.numeric(f(get(ind_variable))),
      value = exp(log_value),
      # need a taxon name for plotting
      taxon_name = "global average"
    )

  # caluclate species means by removing devistations from center
  data_fit2_mn <- estimate_expectation(fit2,
    data = data_fit2_re |> mutate(!!(paste0(ind_variable,"_c")) := 0) |> select(-Predicted) |> distinct()
  ) |>
    as_tibble() |>
    distinct()

  # combine
  data_fit2 <-
    data_fit2_re |>
    rename(log_value = Predicted) |>
    select(-SE, -Residuals) |>
    left_join(
      by = c(paste0(ind_variable,"_mn"), "taxon_name"),
      data_fit2_mn |> select(paste0(ind_variable,"_mn"), taxon_name, log_value_mn = Predicted)
    ) |>
    mutate(
      !!ind_variable := get(paste0(ind_variable,"_mn")) + get(paste0(ind_variable,"_c")),
      !!standard_name := as.numeric(f(get(ind_variable))),
      !!paste0(standard_name,"_mn") := as.numeric(f(get(paste0(ind_variable,"_mn")))),
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
      n_sites = n_distinct(get(standard_name))
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
  
  # data_fit2 |>
  #   mutate(!!standard_name := as.numeric(get(standard_name))) -> data_fit2
  # 
  # data_fit2_global |>
  #   mutate(!!standard_name := as.numeric(get(standard_name))) -> data_fit2_global
  # 
  # data_fit2_taxa |>
  #   mutate(!!standard_name := as.numeric(get(standard_name))) -> data_fit2_taxa
  # 
  # data_fit2_taxa |>
  #   mutate(!!standard_name := as.numeric(get(standard_name))) -> data_fit2_taxa
  # 
  p_fit <-
    data_fit2 |>
    ggplot(aes(get(standard_name), value, group = taxon_name)) +
    geom_line() +
    # raw data
    geom_point(alpha = 0.5, col = "grey", data = data_model, stroke = 0) +
    # fitted lines
    geom_line(alpha = 0.5, col = "black", linewidth = 0.25) +
    # highlight certain species
    geom_point(aes(group = taxon2), col = "orange", show.legend = FALSE, data = data_model_taxa, size = 0.9, stroke = 0) +
    geom_line(aes(group = taxon2), data = data_fit2_taxa, col = "#d95f02", alpha = 1, show.legend = FALSE) +
    # global line
    geom_line(data = data_fit2_global, col = "red") +
    # species centers
    geom_point(aes(x = get(paste0(standard_name,"_mn")), y = value_mn, group = taxon2), col = "red", alpha = 1, size = 0.9, stroke = 0, data = data_fit2_taxa) +
    scale_x_continuous(transform = x_transform) +
    scale_y_log10() +
    theme_classic() +
    labs(y = "trait", x = standard_name) +
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
