find_ITV <- function(data_model) {
    fitted <- glmmTMB(log_value ~ 1 + (1 | taxon_name), data = data_model, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

    taxon_var <- as.numeric(VarCorr(fitted)$cond$taxon_name[1])
    residual_sd <- sigma(fitted)
    taxon_sd <- sqrt(taxon_var)

    variation_tbl <- tibble::tibble(
        variation_type = c("residual", "taxon_name"),
        standard_deviation = c(residual_sd, taxon_sd)
    )

    variation_tbl %>% 
    mutate(variation = standard_deviation^2) %>%
    select(-standard_deviation) %>%
    pivot_wider(names_from = "variation_type", values_from = "variation") %>%
    mutate(ITV_frac = residual/(residual+taxon_name)) %>%
    mutate(trait_name = data_model$trait_name[1]) -> out

    print(out)
}
