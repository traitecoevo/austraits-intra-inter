fit_model2d_frequentist_woody <- function(data_model) {
  
  if(data_model$trait_name[1] %in% c()){
    fit_lme <- glmmTMB(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c + std_average_temperature_mn + std_average_temperature_c)*woody  + (1 | taxon_name) +  (0 + std_log_annual_precipitation_c | taxon_name)+ (0 +std_average_temperature_c | taxon_name), data = data_model)
    print(summary(fit_lme))
    
    list_to_save <- list(fit_lme, data_model)
  } else{
    fit_MAP <- glmmTMB(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c)*woody  + (1 | taxon_name) +  (0 + std_log_annual_precipitation_c | taxon_name), data = data_model)
    fit_MAT <- glmmTMB(log_value ~ (std_average_temperature_mn + std_average_temperature_c)*woody + (1 | taxon_name) +  (0 + std_average_temperature_c | taxon_name), data = data_model)
    print(summary(fit_MAP))
    print(summary(fit_MAT))

    
    list_to_save <- list(fit_MAP,fit_MAT, data_model)
  }
  
  saveRDS(list_to_save,paste0("outputs/model_fits/fit2d_frequentist_woody/",data_model$trait_name[1],".RDS"))
}

