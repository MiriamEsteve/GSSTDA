
test_that("Test matrix denoising functions", {
  #We expect an error when the aspect ration of the matrix is equal to 1 since the function is designed to work with rectangular matrices
  expect_error(get_omega(1))
  #We expect the mu_beta value for any rectangular  matrix aspect ratio (which m < n) to be less than 1.
  expect_true(get_mu_beta(0.3) < 1)
})


