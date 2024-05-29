
test_that("Test matrix denoising functions", {
  #We expect an error when the aspect ration of the matrix is equal to 1 since the function is designed to work with rectangular matrices
  expect_true(optimal_SVHT_coef_gamma_known(10/10000) > 0)
  #We expect the mu_beta value for any rectangular  matrix aspect ratio (which m < n) to be less than 1.
  expect_true(optimal_SVHT_coef_gamma_unknown(10/10000) > 0)
})


