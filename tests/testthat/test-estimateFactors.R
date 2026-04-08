test_that("", {

  # Load data
  data("factor_model")
  data <- scale(factor_model$data)

  # Parsing zoo object
  expect_no_error(factor_fit_zoo <- noOfFactors(data = data, min_no_factors = 1, max_no_factors = 5))
  
  # Parsing matrix object
  expect_no_error(factor_fit_matrix <- noOfFactors(data = t(coredata(data)), min_no_factors = 1, max_no_factors = 5))
  
  # Check whether the results are equal
  expect_equal(factor_fit_zoo, factor_fit_matrix)
  
})
