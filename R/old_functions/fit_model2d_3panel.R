fit_model2d_3panel <- function(data_model,ind_variable = "log_annual_precipitation") {
  
  browser()
  
  if(str_detect(ind_variable,"log")){f = exp}else{f = I}
  if(str_detect(ind_variable,"log")){standard_name = str_remove(ind_variable, "log_")}else{standard_name = ind_variable}
  if(str_detect(ind_variable,"log")){x_transform = "log10"}else{x_transform = "identity"}
  if(ind_variable == "log_annual_precipitation"){
    off_variable = "average_temperature"
  } else{
    off_variable = "log_annual_precipitation"
  }
  
  if(ind_variable == "log_annual_precipitation"){
    environment_label = "MAP (mm)"
  } else{
    environment_label = expression(paste(MAT, " (",degree*C,")"))
  }
  
  if(ind_variable == "log_annual_precipitation"){
    b1_label = expression(paste(B["1[Prec,i]"]))
  } else{
    b1_label = expression(paste(B["1[Temp,i]"]))
  }
  
  trait_name <- data_model$trait_name[1]
  if(trait_name == "leaf_mass_per_area"){trait_label = expression(paste(LMA," (g ",m^{-2},")"))}
  if(trait_name == "huber_value"){trait_label = expression(paste(SA:LA))}
  if(trait_name == "leaf_area"){trait_label = expression(paste(Leaf~area," (",mm^{-2},")"))}
  if(trait_name == "leaf_delta13C"){trait_label = expression(paste(delta^{13},C))}
  if(trait_name == "leaf_N_per_area"){trait_label = expression(paste(N[area]," (g ",m^{-2},")"))}
  if(trait_name == "wood_density"){trait_label = expression(paste(Wood~density," (g ",cm^{-3},")"))}
  if(trait_name == "leaf_P_per_area"){trait_label = expression(paste(P[area]," (g ",m^{-2},")"))}
  
  fits_bayesian <- readRDS(paste0("outputs/model_fits/",trait_name,".RDS"))
  fit2 <-fits_bayesian[[1]]
  
  # fit_lme <- glmmTMB(log_value ~ log_annual_precipitation_mn + log_annual_precipitation_c + average_temperature_mn + average_temperature_c + (log_annual_precipitation_c | taxon_name) + (average_temperature_c | taxon_name), data = data_model)
  
  # fit2 <- fit_lme
  cat("extracting fitted vals -> ")
  data_model %>%
    select(log_annual_precipitation_mn,
           log_annual_precipitation_c,
           average_temperature_mn,
           average_temperature_c,
           taxon_name) %>%
    group_by(taxon_name) %>%
    mutate(!!paste0(off_variable,"_c") := mean(get(paste0(off_variable,"_c")), na.rm =T),
           !!paste0(off_variable,"_mn") := mean(get(paste0(off_variable,"_mn")), na.rm =T)) %>%
    drop_na()-> data_for_fit
  
  # overall fitted - with REs
  data_fit2_re <- estimate_expectation(fit2, data = data_for_fit) |>
    as_tibble()
  
  # Global prediction - change in species mean trait with climate
  data_fit2_global <- estimate_means(fit2, by = c(paste0(ind_variable,"_mn"))) |>
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
  
  data_fit2_global_c <- estimate_means(fit2, by = c(paste0(ind_variable,"_c"))) |>
    as_tibble() |>
    # change variable names so these will plot ok
    select(!!ind_variable := paste0(ind_variable,"_c"), log_value = Mean) |>
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
  
  # caluclate species means by removing devistations from center
  data_fit2_c <- estimate_expectation(fit2,
                                       data = data_fit2_re |> mutate(log_annual_precipitation_c = 0) |> mutate(average_temperature_c = 0) |> select(-Predicted) |> distinct()) |>
    as_tibble() |>
    distinct()
  
  data_c_line <- data_fit2_re %>% 
    summarise(!!paste0("max_05_",standard_name,"_mn") := max(get(paste0(ind_variable,"_mn"))), 
              !!paste0("min_05_",standard_name,"_mn") := min(get(paste0(ind_variable,"_mn"))),
              !!paste0("med_",standard_name,"_mn") := median(get(paste0(ind_variable,"_mn"))),
              !!paste0(off_variable,"_c") := 0, 
              !!paste0(off_variable,"_mn") := mean(get(paste0(off_variable,"_mn")))) %>%
    pivot_longer(cols = c(paste0("min_05_",standard_name,"_mn"),
                          paste0("max_05_",standard_name,"_mn")),
                 values_to = paste0(standard_name,"_c_plot")) %>%
    rename(!!paste0(standard_name,"_mn") := paste0("med_",standard_name,"_mn")) %>%
    mutate(!!paste0(standard_name,"_c") := get(paste0(standard_name,"_c_plot")) - get(paste0(standard_name,"_mn"))) %>%
    mutate(taxon_name = "global")
  
  data_c_line %>%
    mutate(value = as.numeric(exp(fixef(fit2)$cond[1] + fixef(fit2)$cond[2]*data_c_line[,"log_annual_precipitation_mn"] +  fixef(fit2)$cond[3]*data_c_line[,"log_annual_precipitation_c"] + fixef(fit2)$cond[4]*data_c_line[,"average_temperature_mn"] + fixef(fit2)$cond[5]*data_c_line[,"average_temperature_c"])[1:2,1])) -> data_c_line
  
  # combine
  data_fit2 <-
    data_fit2_re |>
    rename(log_value = Predicted) |>
    # select(-SE, -Residuals) |> 
    left_join(
      # by = c(paste0(ind_variable,"_mn"), "taxon_name", "growth_form_simpler"),
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
  
  data_model %>%
    select(log_annual_precipitation_mn,
           log_annual_precipitation_c,
           average_temperature_mn,
           average_temperature_c,
           taxon_name) %>%
    group_by(taxon_name) %>%
    mutate(!!paste0(off_variable,"_c") := mean(get(paste0(off_variable,"_c")), na.rm =T),
           !!paste0(off_variable,"_mn") := mean(get(paste0(off_variable,"_mn")), na.rm =T)) %>%
    drop_na()-> data_for_fit_global
  
  fit_overall <- lm(log(value)~log_annual_precipitation + average_temperature, data =data_model, na.action = na.omit)
  data_fit2_re_overall <- estimate_expectation(fit_overall, by = ind_variable) |>
    as_tibble() |>
    mutate(!!ind_variable := f(get(ind_variable)))
  a <-
    data_model |>
    mutate(!!ind_variable := as.numeric(f(get(ind_variable)))) |>
    # mutate(predicted_value = data_fit2_re$Predicted) |> 
    ggplot(aes(x = get(standard_name), y = value)) +
    geom_point(alpha = 0.75, col = "grey") +
    # raw data
    geom_line(data = data_fit2_re_overall, aes(x = (as.numeric(get(ind_variable))), y = exp(Predicted)), size = 1.5) +
    scale_x_continuous(transform = x_transform) +
    scale_y_log10() +
    theme_classic()  +
    labs(y = trait_label, x = environment_label) +
    theme(legend.position = "none",
          text = element_text(size = 12),
          axis.text = element_text(size = 12)) 
  
  # add red line to this plot
  # hist(coef(fit2)$cond$taxon_name$log_annual_precipitation_c)
  re_strength <- tibble(taxon_name = rownames(coef(fit2)$cond$taxon_name),
                        re_strength = coef(fit2)$cond$taxon_name[[paste0(ind_variable,"_c")]])
  
  data_c_line$re_strength <- as.numeric(fixef(fit2)$cond[paste0(standard_name,"_c")])
  data_c_line |>
    mutate(!!paste0(standard_name, "_c_plot") := f(get(paste0(standard_name, "_c_plot")))) -> data_c_line
  b <-
  data_fit2 |>
    left_join(re_strength) %>%
    ggplot(aes(get(standard_name), value, group = taxon_name)) +
    # geom_line() +
    # raw data
    geom_point(alpha = 0.75, col = "grey", data = data_model) +
    # fitted lines
    geom_line(aes(colour = re_strength, group = taxon_name, alpha = abs(re_strength)), linewidth = 0.75) + 
    geom_line(data = data_c_line, aes(x = as.numeric((get(paste0(standard_name,"_c_plot")))), y = (value)), colour = "black", linewidth = 1.5) +
    # geom_point(aes(group = taxon2), col = "orange", show.legend = FALSE, data = data_model_taxa, size = 1, stroke = 0) +
    # geom_line(aes(group = taxon2), data = data_fit2_taxa, col = "#d95f02", alpha = 1, show.legend = FALSE, size = 1) +
    # geom_point(aes(x = get(paste0(standard_name,"_mn")), y = value_mn, group = taxon2), col = "red", alpha = 1, size = 1.25, stroke = 0, data = data_fit2_taxa) +
    scale_x_continuous(transform = x_transform) +
    scale_y_log10() +
    theme_classic() +
    scale_colour_gradient2(
      low = "blue",
      high = "red",
      mid = "purple",
      midpoint = 0
    ) +
    # theme(legend.position = "none") +
    labs(y = "", x = environment_label) +
    theme(legend.position = "none",
          text = element_text(size = 12),
          axis.text = element_text(size = 12)) 
  

  # red dots for species averages
  c <-
    data_fit2 |>
    ggplot(aes(get(standard_name), value, group = taxon_name)) +
    # geom_line() +
    # raw data
    geom_point(alpha = 0.5, col = "grey", data = data_model) +
    geom_point(data = data_fit2 %>% distinct(taxon_name, .keep_all = TRUE), aes(x = get(standard_name), y = value), alpha = 0.75, col = "orange") +
    # fitted lines
    # geom_line(alpha = 0.5, col = "black", linewidth = 0.25) + 
    # geom_point(aes(group = taxon2), col = "orange", show.legend = FALSE, data = data_model_taxa, size = 1, stroke = 0) +
    # geom_line(aes(group = taxon2), data = data_fit2_taxa, col = "#d95f02", alpha = 1, show.legend = FALSE, size = 1) +
    # geom_point(aes(x = get(paste0(standard_name,"_mn")), y = value_mn, group = taxon2), col = "red", alpha = 1, size = 1.25, stroke = 0, data = data_fit2_taxa) +
    geom_line(data = data_fit2_global, col = "black", linewidth = 1.5) +
    # geom_line(data = data_c_line, aes(x = f(annual_precipitation_c_plot), y = value), col = "black") +
    scale_x_continuous(transform = x_transform) +
    scale_y_log10() +
    theme_classic() +
    labs(y = "leaf_mass_per_area", x = standard_name) +
    theme(legend.position = "none") +
    labs(y = "", x = environment_label) +
    theme(legend.position = "none",
          text = element_text(size = 12),
          axis.text = element_text(size = 12)) 
  d <-
    re_strength |>
    ggplot(aes(x = re_strength)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # geom_line() +
    # raw data
    geom_density(linewidth = 1.5) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = b1_label, y= "Density") +
    theme(legend.position = "none",
          text = element_text(size = 12),
          axis.text = element_text(size = 12))
  
  browser()
  return((a + c + b + d) + plot_layout(nrow = 1))
}
