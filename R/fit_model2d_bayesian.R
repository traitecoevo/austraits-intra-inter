fit_model2d_bayesian <- function(data_model) {
  if(data_model$trait_name[1] %in% c()){
    fit_lme <- glmmTMB(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c + std_average_temperature_mn + std_average_temperature_c)  + (1 | taxon_name) +  (0 + std_log_annual_precipitation_c | taxon_name)+ (0 +std_average_temperature_c | taxon_name), data = data_model)
    fit_stan <- tmbstan(fit_lme$obj, chains=4, iter=5000, control=list(adapt_delta=0.95), cores = 4)
    # huber value n area
    print(summary(fit_lme))
    list_to_save <- list(fit_stan, data_model)
  } else{
    fit_MAP_stan <- brm(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c) + (1 | taxon_name) + (0 + std_log_annual_precipitation_c | taxon_name), data = data_model, chains = 4, iter = 7500, control = list(adapt_delta = 0.98), cores = 4)
    fit_MAT_stan <- brm(log_value ~ (std_average_temperature_mn + std_average_temperature_c) + (1 | taxon_name) + (0 + std_average_temperature_c | taxon_name), data = data_model, chains = 4, iter = 7500, control = list(adapt_delta = 0.98), cores = 4)
    print(summary(fit_MAP_stan))
    print(summary(fit_MAT_stan))

    list_to_save <- list(fit_MAP_stan, fit_MAT_stan, data_model)
  }
  saveRDS(list_to_save,paste0("outputs/model_fits/fit2d_bayesian/",data_model$trait_name[1],".RDS"))
}
