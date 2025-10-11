fit_model2d_PFT <- function(data_model,ind_variable = "log_annual_precipitation") {
  
  
  if(str_detect(ind_variable,"log")){f = exp}else{f = I}
  if(str_detect(ind_variable,"log")){standard_name = str_remove(ind_variable, "log_")}else{standard_name = ind_variable}
  if(str_detect(ind_variable,"log")){x_transform = "log10"}else{x_transform = "identity"}
  if(ind_variable == "log_annual_precipitation"){
    off_variable = "average_temperature"
  } else{
    off_variable = "log_annual_precipitation"
  }
  
  fit2 <- glmmTMB(log_value ~ (log_annual_precipitation_mn + log_annual_precipitation_c + average_temperature_mn + average_temperature_c)*growth_form_simpler + (log_annual_precipitation_c | taxon_name) + (average_temperature_c | taxon_name), data = data_model)
  
  cat("extracting fitted vals -> ")
  
  data_model %>%
    select(log_annual_precipitation_mn,
           log_annual_precipitation_c,
           average_temperature_mn,
           average_temperature_c,
           growth_form_simpler,
           taxon_name) %>%
    group_by(taxon_name) %>%
    mutate(!!paste0(off_variable,"_c") := mean(get(paste0(off_variable,"_c")), na.rm =T),
           !!paste0(off_variable,"_mn") := mean(get(paste0(off_variable,"_mn")), na.rm =T)) %>%
    drop_na()-> data_for_fit
  
  # overall fitted - with REs
  data_fit2_re <- estimate_expectation(fit2, data = data_for_fit) |>
    as_tibble()
  
  
  
  # Global prediction - change in species mean trait with climate
  data_fit2_global <- estimate_means(fit2, by = c(paste0(ind_variable,"_mn"), "growth_form_simpler")) |>
    as_tibble() |>
    # change variable names so these will plot ok
    select(!!ind_variable := paste0(ind_variable,"_mn"),growth_form_simpler, log_value = Mean) |>
    # backtransform so we can plot on raw scale
    mutate(
      !!standard_name := as.numeric(f(get(ind_variable))),
      value = exp(log_value),
      # need a taxon name for plotting
      taxon_name = "global average"
    )
  
  # caluclate species means by removing devistations from center
  data_fit2_mn <- estimate_expectation(fit2,
                                       data = data_fit2_re |> mutate(log_annual_precipitation_c = 0) |> mutate(average_temperature_c = 0) |> select(-Predicted) |> distinct()) |>
    as_tibble() |>
    distinct()
  
  # combine
  data_fit2 <-
    data_fit2_re |>
    rename(log_value = Predicted) |>
    # select(-SE, -Residuals) |> 
    left_join(
      # by = c(paste0(ind_variable,"_mn"), "taxon_name", "growth_form_simpler"),
      data_fit2_mn |> select(paste0(ind_variable,"_mn"), taxon_name, log_value_mn = Predicted, growth_form_simpler)
    ) |>
    mutate(
      !!ind_variable := get(paste0(ind_variable,"_mn")) + get(paste0(ind_variable,"_c")),
      !!standard_name := as.numeric(f(get(ind_variable))),
      !!paste0(standard_name,"_mn") := as.numeric(f(get(paste0(ind_variable,"_mn")))),
      value = exp(log_value),
      value_mn = exp(log_value_mn),
    )
  
  cat("coooking up a sweet plot -> ")
  
  p_fit <-
    data_fit2 |>
    ggplot(aes(get(standard_name), value, group = taxon_name)) +
    geom_line() +
    # raw data
    geom_point(alpha = 0.5, col = "grey", data = data_model, stroke = 0) +
    # fitted lines
    geom_line(alpha = 0.5, col = "black", linewidth = 0.25) +
    geom_line(data = data_fit2_global, col = "red", linewidth = 0.75) +
    facet_wrap(~growth_form_simpler) +
    # geom_point(aes(x = get(paste0(standard_name,"_mn")), y = value_mn, group = taxon2), col = "red", alpha = 1, size = 0.9, stroke = 0, data = data_fit2_taxa) +
    scale_x_continuous(transform = x_transform) +
    scale_y_log10() +
    theme_classic() +
    labs(y = "leaf_mass_per_area", x = standard_name) +
    theme(legend.position = "none") 
  cat("done\n")
  
  list(
    data = data_model,
    fit = fit2,
    data_fit = data_fit2,
    plot_fitted = p_fit
  )
  
}
