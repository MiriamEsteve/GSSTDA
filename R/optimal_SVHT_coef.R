#' @title Optimal Singular Value Hard Threshold (SVHT) Coefficient
#'
#' @description Computes the optimal SVHT coefficient for matrix denoising based on the noise level.
#'
#' @param beta A numeric vector representing the aspect ratio \( m/n \) of the matrix to be denoised, where \( 0 < \beta \leq 1 \).
#' @param sigma_known A logical value indicating if the noise level is known (TRUE) or unknown (FALSE).
#'
#' @return A numeric vector containing the optimal SVHT coefficient for each aspect ratio in \code{beta}.
#' @export
#' @examples
#' beta <- 0.5
#' sigma_known <- TRUE
#' optimal_SVHT_coef(beta, sigma_known)
optimal_SVHT_coef <- function(beta, sigma_known) {
  if (sigma_known) {
    coef <- optimal_SVHT_coef_sigma_known(beta)
  } else {
    coef <- optimal_SVHT_coef_sigma_unknown(beta)
  }
  return(coef)
}


#' @title Optimal SVHT Coefficient with Known Noise Level
#'
#' @description Computes the optimal SVHT coefficient when the noise level is known.
#'
#' @param beta A numeric vector representing the aspect ratio \( m/n \), where \( 0 < \beta \leq 1 \).
#'
#' @return A numeric vector of optimal SVHT coefficients corresponding to each aspect ratio in \code{beta}.
#' @export
#' @examples
#' beta <- 0.5
#' optimal_SVHT_coef_sigma_known(beta)
optimal_SVHT_coef_sigma_known <- function(beta) {
  stopifnot(all(beta > 0))
  stopifnot(all(beta <= 1))
  stopifnot(length(beta) == prod(dim(beta)))  # beta must be a vector

  w <- (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1))
  lambda_star <- sqrt(2 * (beta + 1) + w)
  return(lambda_star)
}


#' @title Optimal SVHT Coefficient with Unknown Noise Level
#'
#' @description Computes the optimal SVHT coefficient when the noise level is unknown.
#'
#' @param beta A numeric vector representing the aspect ratio \( m/n \), where \( 0 < \beta \leq 1 \).
#'
#' @return A numeric vector of optimal SVHT coefficients corresponding to each aspect ratio in \code{beta}.
#' @export
#' @examples
#' beta <- 0.5
#' optimal_SVHT_coef_sigma_unknown(beta)
optimal_SVHT_coef_sigma_unknown <- function(beta) {
  options(warn = -1)  # Suppress warnings
  stopifnot(all(beta > 0))
  stopifnot(all(beta <= 1))
  stopifnot(length(beta) == prod(dim(beta)))  # beta must be a vector

  coef <- optimal_SVHT_coef_sigma_known(beta)

  MPmedian <- sapply(beta, function(b) MedianMarcenkoPastur(b))
  omega <- coef / sqrt(MPmedian)
  return(omega)
}


#' @title Marcenko-Pastur Integral
#'
#' @description Calculates the integral of the Marcenko-Pastur distribution from the lower bound to a specified value.
#'
#' @param x A numeric value representing the upper limit of the integral.
#' @param beta A numeric value representing the aspect ratio \( m/n \), where \( 0 < \beta \leq 1 \).
#'
#' @return A numeric value representing the integral of the Marcenko-Pastur distribution from the lower bound to \code{x}.
#' @export
#' @examples
#' x <- 0.5
#' beta <- 0.5
#' MarcenkoPasturIntegral(x, beta)
MarcenkoPasturIntegral <- function(x, beta) {
  if (beta <= 0 || beta > 1) {
    stop("beta beyond")
  }
  lobnd <- (1 - sqrt(beta))^2
  hibnd <- (1 + sqrt(beta))^2
  if (x < lobnd || x > hibnd) {
    stop("x beyond")
  }
  dens <- function(t) sqrt((hibnd - t) * (t - lobnd)) / (2 * pi * beta * t)
  I <- integrate(dens, lobnd, x)$value
  cat(sprintf("x=%.3f, beta=%.3f, I=%.3f\n", x, beta, I))
  return(I)
}


#' @title Median of the Marcenko-Pastur Distribution
#'
#' @description Calculates the median of the Marcenko-Pastur distribution for a given aspect ratio.
#'
#' @param beta A numeric value representing the aspect ratio \( m/n \), where \( 0 < \beta \leq 1 \).
#'
#' @return A numeric value representing the median of the Marcenko-Pastur distribution for the specified \code{beta}.
#' @export
#' @examples
#' beta <- 0.5
#' MedianMarcenkoPastur(beta)
MedianMarcenkoPastur <- function(beta) {
  MarPas <- function(x) 1 - incMarPas(x, beta, 0)
  lobnd <- (1 - sqrt(beta))^2
  hibnd <- (1 + sqrt(beta))^2
  change <- TRUE
  while (change && (hibnd - lobnd > 0.001)) {
    change <- FALSE
    x <- seq(lobnd, hibnd, length.out = 5)
    y <- sapply(x, MarPas)
    if (any(y < 0.5)) {
      lobnd <- max(x[y < 0.5])
      change <- TRUE
    }
    if (any(y > 0.5)) {
      hibnd <- min(x[y > 0.5])
      change <- TRUE
    }
  }
  med <- (hibnd + lobnd) / 2
  return(med)
}


#' @title Incomplete Marcenko-Pastur Integral
#'
#' @description Calculates the incomplete Marcenko-Pastur integral from a lower limit to the upper bound of the distribution's support.
#'
#' @param x0 A numeric value representing the lower limit of the integral.
#' @param beta A numeric value representing the aspect ratio \( m/n \), where \( 0 < \beta \leq 1 \).
#' @param gamma A numeric value representing an exponent parameter.
#'
#' @return A numeric value representing the incomplete Marcenko-Pastur integral.
#' @export
#' @examples
#' x0 <- 0.5
#' beta <- 0.5
#' gamma <- 1
#' incMarPas(x0, beta, gamma)
incMarPas <- function(x0, beta, gamma) {
  if (beta > 1) {
    stop("beta beyond")
  }
  topSpec <- (1 + sqrt(beta))^2
  botSpec <- (1 - sqrt(beta))^2
  MarPas <- function(x) {
    ifelse((topSpec - x) * (x - botSpec) > 0,
           sqrt((topSpec - x) * (x - botSpec)) / (beta * x) / (2 * pi),
           0)
  }
  if (gamma != 0) {
    fun <- function(x) x^gamma * MarPas(x)
  } else {
    fun <- MarPas
  }
  I <- integrate(fun, x0, topSpec)$value
  return(I)
}


#' @title Example of Using Optimal SVHT Coefficient for Matrix Denoising
#'
#' @description Demonstrates how to use the optimal SVHT coefficient to denoise a matrix with known noise level.
#'
#' @examples
#' # Define dimensions
#' m <- 50
#' n <- 100
#' # Create a noisy matrix
#' Y <- matrix(rnorm(m * n), m, n)
#' # Aspect ratio
#' beta <- m / n
#' # Known noise level
#' sigma <- 1
#' # Compute the threshold
#' coef <- optimal_SVHT_coef(beta, TRUE) * sqrt(n) * sigma
#' # Perform SVD
#' svd_result <- svd(Y)
#' U <- svd_result$u
#' D <- svd_result$d
#' V <- svd_result$v
#' # Apply hard thresholding
#' D[D < coef] <- 0
#' # Reconstruct the denoised matrix
#' Xhat <- U %*% diag(D) %*% t(V)

