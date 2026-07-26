fit_model2d_frequentist <- function(data_model) {
  if(data_model$trait_name[1] %in% c()){
  fit_lme <- glmmTMB(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c + std_average_temperature_mn + std_average_temperature_c) + (1 | taxon_name) +  (0 + std_log_annual_precipitation_c | taxon_name)+ (0 +std_average_temperature_c | taxon_name), data = data_model)
  print(summary(fit_lme))

  list_to_save <- list(fit_lme, data_model)
  } else{
  hessian_warn_MAP <- FALSE
  fit_MAP <- withCallingHandlers(
    glmmTMB(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c) + (1 | taxon_name) + (0 + std_log_annual_precipitation_c | taxon_name), data = data_model, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
    warning = function(w) {
      if (grepl("non-positive-definite Hessian", conditionMessage(w))) {
        hessian_warn_MAP <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  if (hessian_warn_MAP) {
    message("MAP: non-positive-definite Hessian with BFGS; retrying with default optimizer.")
    fit_MAP <- glmmTMB(log_value ~ (std_log_annual_precipitation_mn + std_log_annual_precipitation_c) + (1 | taxon_name) + (0 + std_log_annual_precipitation_c | taxon_name), data = data_model)
  }

  hessian_warn_MAT <- FALSE
  fit_MAT <- withCallingHandlers(
    glmmTMB(log_value ~ (std_average_temperature_mn + std_average_temperature_c) + (1 | taxon_name) + (0 + std_average_temperature_c | taxon_name), data = data_model, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
    warning = function(w) {
      if (grepl("non-positive-definite Hessian", conditionMessage(w))) {
        hessian_warn_MAT <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  if (hessian_warn_MAT) {
    message("MAT: non-positive-definite Hessian with BFGS; retrying with default optimizer.")
    fit_MAT <- glmmTMB(log_value ~ (std_average_temperature_mn + std_average_temperature_c) + (1 | taxon_name) + (0 + std_average_temperature_c | taxon_name), data = data_model)
  }
  print(summary(fit_MAP))
  print(summary(fit_MAT))
  list_to_save <- list(fit_MAP,fit_MAT, data_model)
  }

  saveRDS(list_to_save,paste0("outputs/model_fits/fit2d_frequentist/",data_model$trait_name[1],".RDS"))
}
