# Finds the name of the base_term:woody interaction coefficient in a fitted
# model's term names, trying both interaction orderings. Returns NA if absent.
find_interaction_term <- function(term_names, base_term) {
  cand1 <- paste0(base_term, ":woodywoody")
  cand2 <- paste0("woodywoody:", base_term)
  if (cand1 %in% term_names) return(cand1)
  if (cand2 %in% term_names) return(cand2)
  NA_character_
}

# Delta-method z-test of whether a single growth-form group's slope (for
# base_term, from a model fit with a base_term * woody interaction) differs
# from zero. group_label is "non-woody" (base_term alone) or "woody"
# (base_term + its interaction with woody). Works with glmmTMB fits (fixed
# effects under $cond) and plain lm/glm fits (coef()/vcov() directly).
test_group_slope_zero <- function(fit, base_term, group_label) {
  if (inherits(fit, "glmmTMB")) {
    coefs <- fixef(fit)$cond
    vc    <- as.matrix(vcov(fit)$cond)
  } else {
    coefs <- coef(fit)
    vc    <- as.matrix(vcov(fit))
  }
  term_names <- names(coefs)

  if (!(base_term %in% term_names)) {
    return(tibble::tibble(
      group = group_label,
      estimate = NA_real_,
      std.error = NA_real_,
      z.value = NA_real_,
      p.value = NA_real_,
      significant = NA
    ))
  }

  base_est <- unname(coefs[[base_term]])

  if (group_label == "non-woody") {
    est <- base_est
    se  <- sqrt(vc[base_term, base_term])
  } else {
    int_term <- find_interaction_term(term_names, base_term)
    if (is.na(int_term)) {
      est <- base_est
      se  <- sqrt(vc[base_term, base_term])
    } else {
      int_est <- unname(coefs[[int_term]])
      est <- base_est + int_est
      se  <- sqrt(vc[base_term, base_term] + vc[int_term, int_term] + 2 * vc[base_term, int_term])
    }
  }

  z <- est / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)

  tibble::tibble(
    group = group_label,
    estimate = est,
    std.error = se,
    z.value = z,
    p.value = p,
    significant = p < 0.05
  )
}
