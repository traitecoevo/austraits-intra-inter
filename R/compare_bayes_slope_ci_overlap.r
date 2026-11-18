
# Compare the posterior credible intervals for inter- and intraspecific slopes.
# Woody-split traits are evaluated for non-woody and woody groups;
# huber_value and wood_density are evaluated only for the combined grouping.
compare_bayes_slope_ci_overlap <- function() {

  split_traits <- c("leaf_mass_per_area", "leaf_N_per_area", "leaf_delta13C", "leaf_area")
  combined_only_traits <- c("huber_value", "wood_density")

  summarise_overlap <- function(draws) {
    inter_ci <- quantile(draws$inter_draws, probs = c(0.025, 0.975), na.rm = TRUE)
    intra_ci <- quantile(draws$intra_draws, probs = c(0.025, 0.975), na.rm = TRUE)

    tibble(
      inter_low  = unname(inter_ci[[1]]),
      inter_high = unname(inter_ci[[2]]),
      intra_low  = unname(intra_ci[[1]]),
      intra_high = unname(intra_ci[[2]]),
      ci_overlap  = max(inter_ci[[1]], intra_ci[[1]]) <= min(inter_ci[[2]], intra_ci[[2]])
    )
  }

  extract_draws <- function(fit, inter_term, intra_term, inter_shift = NULL, intra_shift = NULL) {
    draws <- posterior::as_draws_df(fit) %>% as_tibble()

    inter_draws <- draws[[inter_term]]
    intra_draws <- draws[[intra_term]]

    if (!is.null(inter_shift)) {
      inter_draws <- inter_draws + draws[[inter_shift]]
    }
    if (!is.null(intra_shift)) {
      intra_draws <- intra_draws + draws[[intra_shift]]
    }

    tibble(inter_draws = inter_draws, intra_draws = intra_draws)
  }

  map_dfr(c(split_traits, combined_only_traits), function(trait_name) {

    if (trait_name %in% split_traits) {
      fits <- readRDS(paste0("outputs/model_fits/fit2d_bayesian_woody/", trait_name, ".RDS"))
      fit_map <- fits[[1]]
      fit_mat <- fits[[2]]

      bind_rows(
        extract_draws(
          fit_map,
          inter_term = "b_std_log_annual_precipitation_mn",
          intra_term = "b_std_log_annual_precipitation_c"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAP",
            woody_grouping = "non_woody"
          ),
        extract_draws(
          fit_map,
          inter_term = "b_std_log_annual_precipitation_mn",
          intra_term = "b_std_log_annual_precipitation_c"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAP",
            woody_grouping = "combined"
          ),
        extract_draws(
          fit_map,
          inter_term = "b_std_log_annual_precipitation_mn",
          intra_term = "b_std_log_annual_precipitation_c",
          inter_shift = "b_std_log_annual_precipitation_mn:woodywoody",
          intra_shift = "b_std_log_annual_precipitation_c:woodywoody"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAP",
            woody_grouping = "woody"
          ),
        extract_draws(
          fit_mat,
          inter_term = "b_std_average_temperature_mn",
          intra_term = "b_std_average_temperature_c"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAT",
            woody_grouping = "non_woody"
          ),
        extract_draws(
          fit_mat,
          inter_term = "b_std_average_temperature_mn",
          intra_term = "b_std_average_temperature_c"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAT",
            woody_grouping = "combined"
          ),
        extract_draws(
          fit_mat,
          inter_term = "b_std_average_temperature_mn",
          intra_term = "b_std_average_temperature_c",
          inter_shift = "b_std_average_temperature_mn:woodywoody",
          intra_shift = "b_std_average_temperature_c:woodywoody"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAT",
            woody_grouping = "woody"
          )
      )

    } else if (trait_name %in% combined_only_traits) {
      fits <- readRDS(paste0("outputs/model_fits/fit2d_bayesian/", trait_name, ".RDS"))
      fit_map <- fits[[1]]
      fit_mat <- fits[[2]]

      bind_rows(
        extract_draws(
          fit_map,
          inter_term = "b_std_log_annual_precipitation_mn",
          intra_term = "b_std_log_annual_precipitation_c"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAP",
            woody_grouping = "combined"
          ),
        extract_draws(
          fit_mat,
          inter_term = "b_std_average_temperature_mn",
          intra_term = "b_std_average_temperature_c"
        ) %>%
          summarise_overlap() %>%
          mutate(
            trait = trait_name,
            env = "MAT",
            woody_grouping = "combined"
          )
      )
    }
  }) %>%
    select(trait, env, woody_grouping, everything())
}
