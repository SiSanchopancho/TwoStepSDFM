#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @import utils
NULL

utils::globalVariables(c(
  "Value", "Lower 95%-CI", "Upper 95%-CI",
  "Factor", "Variable", "Loading", "(Co-)Variable", "(Co-)Variance",
  "Ridge Penalty", "CV Errors", "Avg. Lasso Penalty",
  "Sparsity Ratio", "# of LARS Steps", "Series",
  "Variance Explained", "Component", "radial_position", "pct_label",
  "Eigen Value", "cut_off"
))