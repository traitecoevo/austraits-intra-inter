fit_model2d_bayesian_3panel_diff_only <- function(data_model,ind_variable = "log_annual_precipitation") {
  
  
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
  if(trait_name == "leaf_mass_per_area"){trait_label = expression(paste(LMA))}
  if(trait_name == "huber_value"){trait_label = expression(paste(SA:LA))}
  if(trait_name == "leaf_area"){trait_label = expression(paste(Leaf~area))}
  if(trait_name == "leaf_delta13C"){trait_label = expression(paste(delta^{13},C))}
  if(trait_name == "leaf_N_per_area"){trait_label = expression(paste(N[area]))}
  if(trait_name == "wood_density"){trait_label = expression(paste(Wood~density,))}
  if(trait_name == "leaf_P_per_area"){trait_label = expression(paste(P[area]))}
  fits <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/",trait_name,".RDS"))
  
  if(trait_name %in% c("leaf_mass_per_area","leaf_N_per_area","leaf_delta13C","leaf_area")){
    if(length(fits) == 2){
      fitmod <-fits[[1]]
      
      posterior_matrix_in <- posterior::as_draws_matrix(fitmod)
      
      
      posterior_matrix_in_prec <- posterior_matrix_in %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[2]`, intra =`beta[3]`,inter_woody =`beta[7]`, intra_woody =`beta[8]`) %>%
        mutate(inter_woody = inter + inter_woody,
               intra_woody = intra + intra_woody) %>%
        as_tibble() %>%
        mutate(env = "MAP")
      
      posterior_matrix_in_temp <- posterior_matrix_in %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[4]`, intra =`beta[5]`,inter_woody =`beta[9]`, intra_woody =`beta[10]`) %>%
        as_tibble() %>%
        mutate(env = "MAT")
      posterior_matrix_in_prec %>%
        bind_rows(posterior_matrix_in_temp) %>%
        mutate(diff = inter - intra,
               diff_woody = inter_woody - intra_woody) %>%
        pivot_longer(cols = -env) -> data_for_plot
      
      data_model <- fits[[2]]
    }
    
    if(length(fits) == 3){
      fitmod_prec <-fits[[1]]
      fitmod_temp <-fits[[2]]
      
      
      posterior_matrix_in_prec <- posterior::as_draws_matrix(fitmod_prec)
      posterior_matrix_in_temp <- posterior::as_draws_matrix(fitmod_temp)
      
      
      posterior_matrix_in_prec <- posterior_matrix_in_prec %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[2]`, intra =`beta[3]`,inter_woody =`beta[5]`, intra_woody =`beta[6]`) %>%
        mutate(inter_woody = inter + inter_woody,
               intra_woody = intra + intra_woody) %>%
        as_tibble() %>%
        mutate(env = "MAP")
      
      posterior_matrix_in_temp <- posterior_matrix_in_temp %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[2]`, intra =`beta[3]`,inter_woody =`beta[5]`, intra_woody =`beta[6]`) %>%
        mutate(inter_woody = inter + inter_woody,
               intra_woody = intra + intra_woody) %>%
        as_tibble() %>%
        mutate(env = "MAT")
      posterior_matrix_in_prec %>%
        bind_rows(posterior_matrix_in_temp) %>%
        # mutate(intra = ifelse(abs(intra) < 0.001, intra/abs(intra)*0.001, intra),
        #        intra_woody = ifelse(abs(intra_woody) < 0.001, intra/abs(intra)*0.001, intra_woody)) %>%
        mutate(intra = as.numeric(intra),
               intra_woody = as.numeric(intra_woody)) %>%
        mutate(inter = as.numeric(inter),
               inter_woody = as.numeric(inter)) %>%
        mutate(diff = inter - intra,
               diff_woody = inter_woody - intra_woody) %>%
        
        pivot_longer(cols = -env) -> data_for_plot
      
      data_model <- fits[[3]]
    }
    cat("extracting fitted vals -> ")
    
    data_for_plot %>%
      group_by(env, name) %>%
      summarise(mean_value = median(value),
                sd_value = sd(value),
                upp_quantile = quantile(value, 0.95),
                low_quantile = quantile(value, 0.05)) %>%
      mutate(upper_sd = mean_value + 1.96*sd_value,
             lower_sd = mean_value - 1.96*sd_value) %>%
      mutate(woody = ifelse(str_detect(name, "woody"), "a",
                            "b")) -> data_for_plot_group
    
    spacing <- c(
      inter = 1,
      intra = 1.3,
      diff = 1.6,
      inter_woody = 1.1,
      intra_woody = 1.4,
      diff_woody = 1.7
    )
    
    data_for_plot$xpos <- spacing[data_for_plot$name]
    data_for_plot_group$xpos <- spacing[data_for_plot_group$name]
    
    spacing_labels <- c(
      Inter = 1.05,
      Intra = 1.35,
      Diff = 1.65
    )
    
    data_for_plot %>%
      filter(grepl("diff",name)) -> data_for_plot
    
    data_for_plot_group %>%
      filter(grepl("diff",name)) -> data_for_plot_group
    
    data_for_plot_group %>%
      mutate(upp_quantile = ifelse(upp_quantile > 2, 2, upp_quantile),
             low_quantile = ifelse(low_quantile < -2, -2, low_quantile)) -> data_for_plot_group
    
    data_for_plot %>%
      mutate(woody = ifelse(str_detect(name, "woody"), "a",
                            "b")) %>%
      ggplot(aes(y = value, x =  xpos, group = interaction(name, woody), colour = woody, fill = woody)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_violinhalf(flip = FALSE, alpha = 0.5) +
      geom_segment(data = data_for_plot_group, aes(y = low_quantile, yend = upp_quantile, x = xpos), colour = "black", size = 1.2) +
      geom_point(data = data_for_plot_group, aes(y = mean_value), colour = "black") +
      facet_grid(~env) +
      coord_flip() +
      # scale_x_continuous(
      #     breaks = spacing_labels,
      #     labels = c(expression(B[1], B[2], B[1]-B[2]))
      #   ) +
      theme_bw() + 
      theme(legend.position = "none",
            text = element_text(size = 12),
            axis.text = element_text(size = 10.5),
            axis.text.y = element_blank(),
            strip.background = element_rect(fill = "white", colour = NA),
            strip.text = element_text(size = 16)) +
      labs(y = bquote(B[1] - B[2]))+
      xlab(trait_label) +
      scale_y_continuous(limits = c(-2,2)) +
      scale_colour_okabeito(order = c(3,6), breaks = c("b","a")) +
      scale_fill_okabeito(order = c(3,6), breaks = c("b","a")) -> a
    return(a)
  } else{
    
    
    if(length(fits) == 2){
      fitmod <-fits[[1]]
      
      posterior_matrix_in <- posterior::as_draws_matrix(fitmod)
      
      
      posterior_matrix_in_prec <- posterior_matrix_in %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[2]`, intra =`beta[3]`) %>%
        as_tibble() %>%
        mutate(env = "MAP")
      
      posterior_matrix_in_temp <- posterior_matrix_in %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[4]`, intra =`beta[5]`) %>%
        as_tibble() %>%
        mutate(env = "MAT")
      posterior_matrix_in_prec %>% 
        bind_rows(posterior_matrix_in_temp) %>%
        mutate(intra = as.numeric(intra)) %>%
        mutate(inter = as.numeric(inter)) %>%
        mutate(diff = inter - intra) %>%
        pivot_longer(cols = -env) -> data_for_plot
      
      
      data_model <- fits[[2]]
    }
    
    if(length(fits) == 3){
      fitmod_prec <-fits[[1]]
      fitmod_temp <-fits[[2]]
      
      
      posterior_matrix_in_prec <- posterior::as_draws_matrix(fitmod_prec)
      posterior_matrix_in_temp <- posterior::as_draws_matrix(fitmod_temp)
      
      
      posterior_matrix_in_prec <- posterior_matrix_in_prec %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[2]`, intra =`beta[3]`) %>%
        as_tibble() %>%
        mutate(env = "MAP")
      
      posterior_matrix_in_temp <- posterior_matrix_in_temp %>%
        as_tibble() %>%
        dplyr::select(inter =`beta[2]`, intra =`beta[3]`) %>%
        as_tibble() %>%
        mutate(env = "MAT")
      posterior_matrix_in_prec %>%
        bind_rows(posterior_matrix_in_temp) %>%
        mutate(intra = as.numeric(intra)) %>%
        mutate(inter = as.numeric(inter)) %>%
        mutate(diff = inter - intra) %>%
        pivot_longer(cols = -env) -> data_for_plot
      
      data_model <- fits[[3]]
    }
    cat("extracting fitted vals -> ")
    
    data_for_plot %>%
      group_by(env, name) %>%
      summarise(mean_value = median(value),
                sd_value = sd(value),
                upp_quantile = quantile(value, 0.95),
                low_quantile = quantile(value, 0.05)) %>%
      mutate(upper_sd = mean_value + 1.96*sd_value,
             lower_sd = mean_value - 1.96*sd_value) %>%
      mutate(woody = ifelse(str_detect(name, "woody"), "a",
                            "b")) -> data_for_plot_group
    
    spacing <- c(
      inter = 1,
      intra = 1.3,
      diff = 1.6
    )
    
    data_for_plot$xpos <- spacing[data_for_plot$name]
    data_for_plot_group$xpos <- spacing[data_for_plot_group$name]
    
    spacing_labels <- c(
      Inter = 1,
      Intra = 1.3,
      Diff = 1.6
    )
    
    data_for_plot %>%
      filter(grepl("diff",name)) -> data_for_plot
    
    data_for_plot_group %>%
      filter(grepl("diff",name)) -> data_for_plot_group
    
    data_for_plot_group %>%
      mutate(upp_quantile = ifelse(upp_quantile > 3.5, 3.5, upp_quantile),
             low_quantile = ifelse(low_quantile < -3.5, -3.5, low_quantile)) -> data_for_plot_group
    
    
    browser()
    
    data_for_plot %>%
      ggplot(aes(y = value, x =  xpos, group = interaction(name), colour = "grey80", fill = "grey80")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_violinhalf(flip = FALSE, alpha = 0.5) +
      geom_segment(data = data_for_plot_group, aes(y = low_quantile, yend = upp_quantile, x = xpos), colour = "black", size = 1.2) +
      geom_point(data = data_for_plot_group, aes(y = mean_value), colour = "black") +
      facet_grid(~env) +
      coord_flip() +
      # scale_x_continuous(
      #     breaks = spacing_labels,
      #     labels = c(expression(B[1], B[2], B[1]-B[2]))
      #   ) +
      theme_bw() + 
      theme(legend.position = "none",
            text = element_text(size = 12),
            axis.text = element_text(size = 10.5),
            axis.text.y = element_blank(),
            strip.background = element_rect(fill = "white", colour = NA),
            strip.text = element_text(size = 16)) +
      labs(y = bquote(B[1] - B[2]))+
      xlab(trait_label) +
      scale_y_continuous(limits = c(-2,2))-> a
    return(a)
    
  }
}

