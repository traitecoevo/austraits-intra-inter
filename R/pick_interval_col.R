# Normalises interval column names returned by different modelbased/brms
# prediction methods (estimate_means(), estimate_expectation(), fitted()).
# Returns NA when the prediction object does not include intervals.
pick_interval_col <- function(df, candidates) {
  col_name <- candidates[candidates %in% names(df)][1]
  if (is.na(col_name)) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[col_name]])
}
