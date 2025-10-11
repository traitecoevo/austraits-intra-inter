fit_model2d_PFT <- function(data_model, rerun_model = FALSE){
  
  trait_name <- data_model$trait_name[1]
  
  if(!file.exists(paste0("outputs/model_fits/fit2dPFT/",trait_name,".RDS"))|rerun_model){
  fit2dPFT <- glmmTMB(log_value ~ (log_annual_precipitation_mn + log_annual_precipitation_c + average_temperature_mn + average_temperature_c)*growth_form_simpler + (log_annual_precipitation_c | taxon_name) + (average_temperature_c | taxon_name), data = data_model)
  saveRDS(fit2dPFT, paste0("outputs/model_fits/fit2dPFT/",trait_name,".RDS"))
  } else{
    print("model already run")
  }
}
