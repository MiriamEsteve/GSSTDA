#' @title Computes lambda
#'
#' @description Computes the value of lambda as defined in: The Optimal Hard Threshold for Singular Values is \eqn{\sqrt(4/ 3)}
#'
#' @param bet numeric
#'
#' @return A numeric lambda value.
#' @export
#'
#' @examples
#' get_lambda(0.3)
get_lambda <- function(bet){
  w = (8*bet) / ((bet + 1) + sqrt(bet^2 + 14*bet + 1))
  lambda_star <- sqrt(2*(bet + 1) + w)
  return(lambda_star)
}


#' @title Marcenko-Pastur distribution to integrate.
#'
#' @description This function is an auxiliary function includes  the marcenco-pastur function that was to be integrated to find the upper integration bound that produces and area of 1/2.
#'
#' @param t Parameter t
#' @param bet Beta value. Aspect ratio of the input matrix.  \eqn{\frac{m}{n}}, were m is the number  of rows of the input matrix and n the number of columns.
#'
#' @return returns the function value for a specific t and a particular aspect ration m/n.
#' @export
#'
#' @examples fun_to_int(1,0.3)
fun_to_int <- function(t,bet){
  b_p <- (1 + sqrt(bet))^2
  b_m <- (1 - sqrt(bet))^2
  numerator <- sqrt((b_p - t)*(t - b_m))
  denominator <- 2*pi*t*bet
  res <- numerator/denominator
  return(res)
}



#' @title Get mu_beta
#'
#' @description This function identifies the upper integration bound of the Marcenko-Pastur distribution. It explores 100 values in a given interval. The selects the values clusest to 1/2 on the left and right. As the upper integration boued
#' if the distance between the left and right approximations is lower than a given threshold 1e-10 it converges and the upper bound producing an area of 1/2 is defined as the average of the left and right
#' approximations.
#'
#' @param bet Beta value. Aspect ratio of the input matrix.  \deqn{\frac{m}{n}}, were m is the number of rows of the input matrix and n the number of columns.
#'
#' @return Returns the mu beta value. This is the upper limit of integration where the Marcenko-Pastur distribution is equal to 1/2.
#' @export
#'
#' @examples
#' get_mu_beta(0.3)
get_mu_beta <- function(bet){
  thresh_diff <- 1e-10
  lbond <- (1 - sqrt(bet))^2
  hbond <- (1 + sqrt(bet))^2
  seqvals <- seq(lbond,hbond,length.out = 100)
  #end_loop <- FALSE
  counter <- 1
  while(T){ #end_loop == FALSE
    #print(counter)
    values_int <- c()
    for (i in 1:length(seqvals)){
      res <- stats::integrate(fun_to_int,lbond,seqvals[i],bet)$value
      values_int <- c(values_int,res)
    }
    if(abs(max(values_int[values_int < 0.5]) - 0.5) < thresh_diff & abs(min(values_int[values_int > 0.5]) - 0.5) < thresh_diff){
      #print("Done")
      final_seqval <- (max(seqvals[values_int < 0.5]) + min(seqvals[values_int > 0.5]))/2
      #end_loop <- TRUE
      break
    }else{
      seqvals <- seq(seqvals[max(which(values_int < 0.5))],seqvals[min(which(values_int > 0.5))],length.out = 100)
    }
    counter <- counter + 1
  }
  return(final_seqval)
}


#' @title Compute omega
#'
#' @description Computes the omega value as described in Add reference.
#'
#' @param bet Beta value. Aspect ratio of the input matrix.  \deqn{\frac{m}{n}}, were m is the number of rows of the input matrix and n the number of columns.
#'
#' @return omega. Returns the omega value.
#' @export
#'
#' @examples
#' get_omega(0.3)
get_omega <- function(bet){
  lamb <- get_lambda(bet)
  mu_beta <- get_mu_beta(bet)
  omega <- lamb/sqrt(mu_beta)
  return(omega)
}


#' @title Rectangular matrix denoiser
#'
#' @description Takes a rectangular matrix composed by the addition of a signal and a gaussian noise matrix and returns a matrix of the same dimension that is denoised through a singular value truncation proccess.
#'
#' @param input_mat A noisy matrix to denoise.
#'
#' @return Returns a same dimension denoised version of the matrix.
#' @export
#'
#' @examples
#' denoise_rectangular_matrix(matrix(c(1,2,3,4,5,2,3,1,2,3),ncol = 2))
denoise_rectangular_matrix <- function(input_mat){
  if(base::is.matrix(input_mat) | (base::is.data.frame(input_mat) & all(base::sapply(input_mat, is.numeric)))){
    if(base::is.data.frame(input_mat)){
      input_mat <- base::as.matrix(input_mat)
    }
    omega_found <- get_omega(ncol(input_mat)/nrow(input_mat))
    svd_input_mat <- base::svd(input_mat)
    D <- diag(svd_input_mat$d)
    U <- svd_input_mat$u
    V <- svd_input_mat$v
    threshold_singular <- stats::median(svd_input_mat$d)*omega_found
    up_to_sv <- length(svd_input_mat$d[svd_input_mat$d > threshold_singular])
    diag(D)[(up_to_sv + 1):length(diag(D))] <- 0
    mat_denoised <- U %*% D %*% t(V)
    rownames(mat_denoised) <- rownames(input_mat)
    colnames(mat_denoised) <- colnames(input_mat)
    return(mat_denoised)
  }else{
    print("Data must be provided in numeric matrix or data.frame formats...")
  }
}
