test_that("", {
  
  # Load data
  data("factor_model")
  data_zoo <- scale(factor_model$data)
  data_mat <- t(coredata(scale(factor_model$data)))
  no_of_factors <- dim(factor_model$factors)[2]
  no_of_vars <- dim(data_zoo)[2]
  selected <- rep(floor(no_of_vars * 0.5), no_of_factors)

  # Parsing zoo object
  expect_silent(sdfm_fit_zoo <- twoStepSDFM(data = data_zoo, delay = factor_model$delay, selected = selected,
                                            no_of_factors = no_of_factors))
  
  # Parsing a matrix object
  expect_silent(sdfm_fit_mat <- twoStepSDFM(data = data_mat, delay = factor_model$delay, selected = selected,
                                            no_of_factors = no_of_factors))
  
  # Basic checks 
  expect_equal(t(coredata(sdfm_fit_zoo$data)), sdfm_fit_mat$data)
  expect_equal(sdfm_fit_zoo$loading_matrix_estimate, sdfm_fit_mat$loading_matrix_estimate)
  expect_equal(t(coredata(sdfm_fit_zoo$smoothed_factors)), sdfm_fit_mat$smoothed_factors)
  expect_equal(sdfm_fit_zoo$smoothed_state_variance, sdfm_fit_mat$smoothed_state_variance)
  expect_equal(sdfm_fit_zoo$factor_var_lag_order, sdfm_fit_mat$factor_var_lag_order)
  expect_equal(sdfm_fit_zoo$error_var_cov_cholesky_factor, sdfm_fit_mat$error_var_cov_cholesky_factor)
  expect_equal(sdfm_fit_zoo$llt_success_code, sdfm_fit_mat$llt_success_code)
  expect_equal(sdfm_fit_zoo$factor_fcast_horizon, sdfm_fit_mat$factor_fcast_horizon)
  expect_equal(sdfm_fit_zoo$data_delay, sdfm_fit_mat$data_delay)
  expect_true(is.zoo(sdfm_fit_zoo$smoothed_factors))
  expect_true(!is.zoo(sdfm_fit_mat$smoothed_factors))
  
})

