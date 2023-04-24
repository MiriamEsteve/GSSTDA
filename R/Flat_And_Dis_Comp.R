#' @title Flatten normal tissues
#'
#' @description Given a matrix which contains the expression values of \code{n} healthy tissue samples produces the flattened vector matrix
#' as reported in Disease-specific genomic analysis: identifying the signature of pathological biology.
#'
#' @param normal_tiss A normal tissue data gene expression matrix.
#'
#' @return A gene expression matrix containing the flattened version of the vectors.
#' @export
#'
#' @examples
#' normal_tissue_matrix <- matrix(stats::rnorm(36),nrow=6)
#' flatten_normal_tiss(normal_tissue_matrix)
flatten_normal_tiss <- function(normal_tiss){
  df_out <- normal_tiss
  for(i in 1:ncol(normal_tiss)){
    df_out[,i] <- stats::fitted(stats::lm(normal_tiss[,i] ~ 0 + ., data = data.frame(normal_tiss)[,-i]))
  }
  return(df_out)
}


#' @title Generate disease component matrix.
#'
#' @description This functions produces a disease component matrix from the complete dataset and the denoised flattened normal data.
#'
#' @param full_data Matrix containing the full dataset.
#' @param normal_space Denoised flattened healthy tissue data.
#'
#' @return A matrix containing the disease component of the complete dataset.
#' @export
#'
#' @examples
#' full_data <- matrix(stats::rnorm(120),ncol=20)
#' normal_tissue <- full_data[,11:20]
#' normal_tissue_f <- flatten_normal_tiss(normal_tissue)
#' normal_tissue_f_d <- denoise_rectangular_matrix(normal_tissue_f)
#' disease_component <- generate_disease_component(full_data,normal_tissue_f_d)
generate_disease_component <- function(full_data, normal_space){
  disease_component <- full_data
  pb <- utils::txtProgressBar(min = 0, max = ncol(full_data), style = 3)
  for(i in 1:ncol(full_data)){
    utils::setTxtProgressBar(pb, i)
    disease_component[,i] <- stats::resid(stats::lm(full_data[,i] ~ 0 + ., data = data.frame(normal_space)))
  }
  return(disease_component)
}

