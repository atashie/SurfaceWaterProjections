#' @importFrom stats aggregate approx coef cor.test lm median quantile sd
NULL

# Global variable declarations for data.table NSE (non-standard evaluation)
# This avoids R CMD check NOTEs about "no visible binding for global variable"
utils::globalVariables(c(
  ".", "water_year", "month", "dowy", "Q", "PPT",
  "Q_annual", "P_annual", "n_valid", "expected_days",
  "dQ", "dP", "elasticity", "dQ_w", "dP_w", "e_w",
  "cum_Q", "cum_P", "slope", "qp_bimodality",
  "dS", "S", "avg_storage",
  "CLASS", "human_interference_class"
))
