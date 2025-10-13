fit_model2d_bayesian <- function(data_model) {
  browser()
  fit_lme <- glmmTMB(log_value ~ log_annual_precipitation_mn + log_annual_precipitation_c + average_temperature_mn + average_temperature_c + (log_annual_precipitation_c | taxon_name) + (average_temperature_c | taxon_name), data = data_model)
  fit_stan <- tmbstan(fit_lme$obj, chains=4, iter=500, control=list(adapt_delta=0.9))
  # fit_stan <- fit_lme
  fit_both <- list(fit_lme, fit_stan)
  saveRDS(fit_both,paste0("outputs/model_fits/",data_model$trait_name[1],".RDS"))
  return(fit_both)
}
