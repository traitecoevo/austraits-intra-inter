# Per-species within-species slopes from a frequentist woody-split glmmTMB fit,
# built explicitly as fixed + random effects rather than via coef(), which matches
# fixed and random terms by name. Under the by-group model the woody random slope
# is named std_<var>_c:woodywoody — the same name as the fixed woody interaction —
# so coef() would wrongly add the two.
#
# Returns one row per taxon:
#   re_strength  — non-woody fixed slope + the taxon's random slope deviation
#   woody_effect — fixed woody interaction on the slope (add it for woody taxa)
# so downstream code can keep using re_strength + woody_effect for woody taxa.
# Works for both the pooled ((0 + std_<var>_c | taxon_name), one random-slope
# column) and by-group (diag(0 + std_<var>_c:woody | taxon_name), one column per
# growth form; a taxon's random slope sits in its own group's column, the other
# being 0) random-effect structures. Only the slope columns (std_<var>_c...) are
# used, so group-specific random intercepts (diag(0 + woody | taxon_name)) are ignored.
species_within_slopes <- function(fit, ind_variable) {
  c_term  <- paste0("std_", ind_variable, "_c")
  fe      <- fixef(fit)$cond
  re      <- ranef(fit)$cond$taxon_name
  re_cols <- grep(paste0("^", c_term), colnames(re), value = TRUE)

  tibble(
    taxon_name   = rownames(re),
    re_strength  = unname(fe[[c_term]] + rowSums(re[, re_cols, drop = FALSE])),
    woody_effect = unname(fe[[paste0(c_term, ":woodywoody")]])
  )
}
