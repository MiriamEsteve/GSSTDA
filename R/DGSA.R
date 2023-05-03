#' @title Flatten normal tissues
#' @description Given a matrix containing the expression values of
#' \code{n} healthy tissue samples, it produces the flattened vector matrix
#' as reported in "Disease-specific genomic analysis: identifying
#' the signature of pathology biology".
#' @param normal_tiss A normal tissue data gene expression matrix.
#' The columns should be the samples and the rows should be the genes.
#' @return A gene expression matrix containing the flattened
#' version of the vectors.
#' @export
#' @examples
#' \dontrun{
#' normal_tissue_matrix <- matrix(stats::rnorm(36), nrow=6)
#' flatten_normal_tiss(normal_tissue_matrix)
#' }
flatten_normal_tiss <- function(normal_tiss){
  matrix_flatten_normal_tiss <- normal_tiss
  for(i in 1:ncol(normal_tiss)){
    matrix_flatten_normal_tiss[,i] <- stats::fitted(stats::lm(normal_tiss[,i] ~ 0 + ., data = data.frame(normal_tiss)[,-i]))
  }
  return(matrix_flatten_normal_tiss)
}


#' @title Computes lambda
#' @description Computes the value of lambda as defined in: "The Optimal
#' Hard Threshold for Singular Values is \eqn{\sqrt(4/ 3)}".
#' @param bet Beta value. Aspect ratio of the input matrix.  \deqn{\frac{m}{n}},
#' were m is the number of rows of the input matrix and n the number of columns.
#' @return Numeric. Lambda value.
#' @export
#' @examples
#' get_lambda(0.3)
get_lambda <- function(bet){
  w = (8*bet) / ((bet + 1) + sqrt(bet^2 + 14*bet + 1))
  lambda_star <- sqrt(2*(bet + 1) + w)
  return(lambda_star)
}


#' @title Marcenko-Pastur distribution to integrate.
#' @description This function is an auxiliary function that includes
#' the Marcenco-Pastur function that has to be integrated to find the
#' upper bound of integration that produces and area of 1/2.
#' @param t Parameter t.
#' @param bet Beta value. Aspect ratio of the input matrix,
#' \eqn{\frac{m}{n}}, were m is the number of rows of the input
#' matrix and n the number of columns.
#' @export
#' @return It returns the function value for a specific t and a particular
#' aspect ration m/n.
#' @examples fun_to_int(1,0.3)
fun_to_int <- function(t,bet){
  b_p <- (1 + sqrt(bet))^2
  b_m <- (1 - sqrt(bet))^2
  numerator <- sqrt((b_p - t)*(t - b_m))
  denominator <- 2*pi*t*bet
  res <- numerator/denominator
  return(res)
}


#' @title Get mu sub beta
#' @description This function identifies the upper bound of integration of the
#' Marcenko-Pastur distribution, as described in "The Optimal Hard
#' Threshold for Singular Values is \eqn{\sqrt(4/ 3)}". It explores 100
#' values in a given interval. Then it selects the values closest to 1/2
#' on the left and on the right. As the upper bound of integration,
#' if the distance between the left and right approximations is lower than
#' a given threshold (1e-10), it converges and the upper bound that produces
#' an area of 1/2 is defined as the mean of the left and right approximations.
#' @param bet Beta value. Aspect ratio of the input matrix, \deqn{\frac{m}{n}},
#' were m is the number of rows of the input matrix and n the number of columns.
#' @return It returns the mu beta value. This is the upper limit of integration
#' where the Marcenko-Pastur distribution is equal to 1/2.
#' @export
#' @examples
#' get_mu_beta(0.3)
get_mu_beta <- function(bet){
  thresh_diff <- 1e-10
  lbond <- (1 - sqrt(bet))^2
  hbond <- (1 + sqrt(bet))^2
  seqvals <- seq(lbond,hbond,length.out = 100)
  while(T){
    values_int <- c()
    for (i in 1:length(seqvals)){
      res <- stats::integrate(fun_to_int,lbond,seqvals[i],bet)$value
      values_int <- c(values_int,res)
    }
    if(abs(max(values_int[values_int < 0.5]) - 0.5) < thresh_diff & abs(min(values_int[values_int > 0.5]) - 0.5) < thresh_diff){
      final_seqval <- (max(seqvals[values_int < 0.5]) + min(seqvals[values_int > 0.5]))/2
      break
    }else{
      seqvals <- seq(seqvals[max(which(values_int < 0.5))],seqvals[min(which(values_int > 0.5))],length.out = 100)
    }
  }
  return(final_seqval)
}


#' @title Compute the omega value
#' @description It computes the omega value as described in "The Optimal Hard
#' Threshold for Singular Values is \eqn{\sqrt(4/ 3)}".
#' @param bet Beta value. Aspect ratio of the input matrix, \deqn{\frac{m}{n}},
#' were m is the number of rows of the input matrix and n the number of columns.
#' @return Numeric. Omega value.
#' @export
#' @examples
#' get_omega(0.3)
get_omega <- function(bet){
  lamb <- get_lambda(bet)
  mu_beta <- get_mu_beta(bet)
  omega <- lamb/sqrt(mu_beta)
  return(omega)
}


#' @title Rectangular Matrix Denoiser.
#' @description It takes a rectangular matrix composed by the addition of
#' a signal matrix and a Gaussian noise matrix and returns a matrix of the same
#' dimension that is denoised through a Singular Value Decomposition
#' truncation process. The selection of the number of singular values is
#' chosen following the proposal by "The optimal hard threshold
#' for singular values is \eqn{\sqrt(4/ 3)}". It should be used after
#' the function \code{flatten_normal_tiss}.
#' @param matrix_flatten_normal_tiss A rectangular noisy matrix to denoise. It is return by
#' \code{flatten_normal_tiss} function.
#' @return A the normal space which has the same dimension denoised version of the matrix
#' returned by \code{flatten_normal_tiss}.
#'
#' @export
#' @examples
#' \dontrun{
#' denoise_rectangular_matrix(matrix(c(1,2,3,4,5,2,3,1,2,3),ncol = 2))
#' }
denoise_rectangular_matrix <- function(matrix_flatten_normal_tiss){
  ### Model reduction ##
  # Find omega value
  omega_found <- get_omega(ncol(matrix_flatten_normal_tiss)/nrow(matrix_flatten_normal_tiss))

  # Calculate Singular Value Decomposition of matrix_flatten_normal_tiss
  svd <- base::svd(matrix_flatten_normal_tiss)
  # Obtain values D, U, V.
  ## d = a vector containing the singular values of x
  ## u = a matrix whose columns contain the left singular vectors of x, present if nu > 0
  ## v = a matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension c(p, nv).
  d <- svd$d
  u <- svd$u
  v <- svd$v

  # Threshold of the singular matrix
  threshold_singular <- stats::median(d)*omega_found

  # Count values up to the threshold
  up_to_threshold <- length(d[d > threshold_singular])
  d[(up_to_threshold + 1):length(d)] <- 0

  # Build rectangular renoised matrix of the normal tiss space
  normal_space <- u %*% diag(d) %*% t(v)
  rownames(normal_space) <- rownames(matrix_flatten_normal_tiss)
  colnames(normal_space) <- colnames(matrix_flatten_normal_tiss)

  return(normal_space)
}


#' @title Generate disease component matrix.
#' @description This function produces a disease component matrix
#' from an expression matrix and the denoised flattened matrix constructed
#' from "healthy tissue data".
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the gens. Both tumour and healthy samples should be included.
#' @param normal_space Denoised flattened matrix constructed from
#' "healthy tissue data". Output of the function \code{denoise_rectangular_matrix}.
#' @return Disease component matrix that contains the disease component
#' of the provided normal space
#'
#' @export
#' @examples
#' full_data <- matrix(stats::rnorm(120),ncol=20)
#' normal_tissue <- full_data[,11:20]
#' normal_tissue_f <- flatten_normal_tiss(normal_tissue)
#' normal_tissue_f_d <- denoise_rectangular_matrix(normal_tissue_f)
#' disease_component <- generate_disease_component(full_data,normal_tissue_f_d)
generate_disease_component <- function(full_data, normal_space){
  # calculate the distance of all points to the normal space by calculating the regression residuals
  disease_component <- full_data
  for(i in 1:ncol(full_data)){
    disease_component[,i] <- stats::resid(stats::lm(full_data[,i] ~ 0 + ., data = data.frame(normal_space)))
  }
  return(disease_component)
}
