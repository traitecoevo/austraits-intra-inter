fit_model2d_bayesian_PFT <- function(data_model, variable = "log_annual_precipitation") {
  if(variable == "log_annual_precipitation"){
  fit_lme <- glmmTMB(log_value ~ (log_annual_precipitation_mn + log_annual_precipitation_c)*growth_form_simpler + (log_annual_precipitation_c | taxon_name), data = data_model)
  }
  if(variable == "average_temperature"){
  fit_lme <- glmmTMB(log_value ~ (average_temperature_mn + average_temperature_c)*growth_form_simpler + (average_temperature_c | taxon_name), data = data_model)
  }
  fit_stan <- tmbstan(fit_lme$obj, chains=4, iter=5000, control=list(adapt_delta=0.95))
  latest_warnings <- warnings()
  fit_both <- list(fit_lme, fit_stan, latest_warnings)
  saveRDS(fit_both,paste0("outputs/model_fits/",data_model$trait_name[1],"_",variable,"PFT.RDS"))
  return(fit_both)
}
