rel_intra_2d_PFT <- function(){
  
  trait2d_PFT <- list.files("outputs/model_fits/fit2dPFT", full.names = TRUE)
  
  file_in <- list()
  
  for(i in 1:length(trait2d_PFT)){
    file_in[[i]] <- readRDS(trait2d_PFT[i])
    }
  
  trait_names <- str_sub(str_sub(trait2d_PFT, nchar("outputs/model_fits/fit2dPFT/")+1), end = -(nchar(".RDS")+1))
  
  data_for_plot <- list()
  
  for(i in 1:length(trait_names)){
  tibble(names = names((fixef(file_in[[i]]))$cond), values = as.vector(fixef(file_in[[i]])$cond)) |>
    filter(names != "(Intercept)") %>%
    filter(str_detect(names, "average|annual")) |>
    separate(names, ":", into = c("climate","growth_form")) %>%
    mutate(climate = ifelse(str_sub(climate, start = -2) == "_c", paste0(climate,"e"),climate)) %>%
    separate(climate, -2, into = c("climate","type")) %>%
    mutate(growth_form = str_remove(growth_form, "growth_form_simpler")) %>%
    mutate(growth_form = ifelse(is.na(growth_form), "climber", growth_form)) |>
      pivot_wider(names_from = "growth_form", values_from = "values") |>
      pivot_longer(cols = -c(climate, type, climber)) |>
      mutate(value = climber+value) |>
      pivot_wider(names_from = "name", values_from = "value") |>
      pivot_longer(cols = -c(climate, type)) |>
    mutate(trait_name = trait_names[i]) -> data_for_plot[[i]]
  }
  data_for_plot |>
    bind_rows() |>
    pivot_wider(names_from = "type", values_from = "value") |>
    mutate(ratio = ce/mn) |>
    ggplot(aes(x = ratio, y = trait_name, fill = climate, group = climate)) +
    geom_col(position = "dodge") +
    facet_wrap(~name) +
    theme_bw() +
    ylab("")
    
    
    pivot_wider(names_from = "climate", values_from = "values")
  
  pivot_wider(names_from = "names", values_from ="values")
  
  if(!file.exists(paste0("outputs/model_fits/fit2dPFT/",trait_name,".RDS"))|rerun_model){
    fit2dPFT <- glmmTMB(log_value ~ (log_annual_precipitation_mn + log_annual_precipitation_c + average_temperature_mn + average_temperature_c)*growth_form_simpler + (log_annual_precipitation_c | taxon_name) + (average_temperature_c | taxon_name), data = data_model)
    saveRDS(fit2dPFT, paste0("outputs/model_fits/fit2dPFT/",trait_name,".RDS"))
  } else{
    print("model already run")
  }
}
