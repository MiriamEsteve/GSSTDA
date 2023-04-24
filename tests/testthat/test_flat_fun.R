
test_that("Test matrix denoising functions", {
  #We expect an error when the aspect ration of the matrix is equal to 1 since the functin is designed to work with rectangular matrices
  expect_error(get_omega(1))
  #We expect the mu_beta value for any rectangular  matrix aspect ratio (which m < n) to be less than 1.
  expect_true(get_mu_beta(0.3) < 1)
  #We expect a message informing that data must be provided in the appropiate format when data is not a numeric matrix or dataframe.
  expect_output(denoise_rectangular_matrix(data.frame(matrix(c("a","b","c","d","e","f"),ncol = 2))), "Data must be provided in numeric matrix or data.frame formats...")
  #We expect the ouput of denoise_rectangular_matrix to be a function when a rectangular matrix function is incuded.
  expect_true(base::is.matrix(denoise_rectangular_matrix(data.frame(matrix(c(1,2,3,4,5,1,3,2,4,5,1,7,7,7,2),ncol = 3)))))
})


