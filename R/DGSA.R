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
#' \donttest{
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
#' @param gamma A parameter that indicates the magnitude of the noise.
#' By default gamma is unknown.
#' @return A the normal space which has the same dimension denoised version of the matrix
#' returned by \code{flatten_normal_tiss}.
#' @export
denoise_rectangular_matrix <- function(matrix_flatten_normal_tiss, gamma){
  #Transpose matrix_flatten_normal_tiss
  matrix_flatten_normal_tiss <- t(matrix_flatten_normal_tiss)
  R <- nrow(matrix_flatten_normal_tiss)
  l <- ncol(matrix_flatten_normal_tiss)
  beta <- R/l

  # Calculate Singular Value Decomposition of matrix_flatten_normal_tiss
  svd <- base::svd(matrix_flatten_normal_tiss)
  # Obtain values D, U, V.
  ## d = a vector containing the singular values of x
  ## u = a matrix whose columns contain the left singular vectors of x, present if nu > 0
  ## v = a matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension c(p, nv).
  d <- svd$d
  u <- svd$u
  v <- svd$v

  ## Gamma is observed
  if(!(is.na(gamma))){
      tau <- optimal_SVHT_coef_gamma_known(beta) * sqrt(l) * gamma
  }else{
      tau <- optimal_SVHT_coef_gamma_unknown(beta) * stats::median(d)
  }

  # Apply hard thresholding
  d[d < tau] <- 0

  # Build rectangular renoised matrix of the normal tiss space
  normal_space <- u %*% diag(d) %*% t(v)
  rownames(normal_space) <- rownames(matrix_flatten_normal_tiss)
  colnames(normal_space) <- colnames(matrix_flatten_normal_tiss)

  normal_space <- t(normal_space)
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
#' \donttest{
#' full_data <- matrix(stats::rnorm(120),ncol=20)
#' normal_tissue <- full_data[,11:20]
#' normal_tissue_f <- flatten_normal_tiss(normal_tissue)
#' normal_tissue_f_d <- denoise_rectangular_matrix(normal_tissue_f, gamma=NA)
#' disease_component <- generate_disease_component(full_data,normal_tissue_f_d)}
generate_disease_component <- function(full_data, normal_space){
  # calculate the distance of all points to the normal space by calculating the regression residuals
  disease_component <- full_data
  for(i in 1:ncol(full_data)){
    disease_component[,i] <- stats::resid(stats::lm(full_data[,i] ~ 0 + ., data = data.frame(normal_space)))
  }
  return(disease_component)
}


#' @title plot dgsa
#' @description
#' It draws the heatmap of the DGSA result by selecting the 100 genes with
#' the highest variability between samples.
#' @param selected_matrix_disease_component Disease component matrix of
#' the selected genes that contains the disease component of all patients.
#' Output of the function \code{generate_disease_component}.
#' @param case_tag Character vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' It must be indicated for each patient whether he/she is healthy or not.
#' One value should be used to indicate whether the patient is healthy and
#' another value should be used to indicate whether the patient's sample is
#' tumourous. The user will then be asked which one indicates whether
#' the patient is healthy. Only two values are valid in the vector in total.
#' @import ComplexHeatmap
#' @import circlize
#' @export
#' @return The heatmap of the DGSA result.
plot_dgsa <- function(selected_matrix_disease_component, case_tag){
  col_fun = circlize::colorRamp2(c(-4, 0,4),
                                 c("red", "black", "green"))
  row_text_size = 10
  ha = ComplexHeatmap::HeatmapAnnotation(Group = case_tag,
                                         annotation_legend_param = list(
                                           Group = list(title = "Group",
                                                        at = c(unique(case_tag)))
                                         ))
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(selected_matrix_disease_component,
                                               cluster_columns = T,col = col_fun,
                                               cluster_rows = F,
                                               heatmap_legend_param = list(
                                                 title = "Expression",
                                                 at = c(-4, 0, 4)),
                                               row_names_gp = grid::gpar(fontsize = 5.25),
                                               column_names_gp = grid::gpar(fontsize = 0),
                                               top_annotation = ha))
}

#' @title results dgsa
#' @description
#' It calculates the 100 genes with the highest variability in the matrix
#' disease component between samples and use them to draw the heat map.
#' @param matrix_disease_component Disease component matrix that contains
#' the disease component of all patients. Output of the function
#' \code{generate_disease_component}.
#' @param case_tag Character vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' It must be indicated for each patient whether he/she is healthy or not.
#' One value should be used to indicate whether the patient is healthy and
#' another value should be used to indicate whether the patient's sample is
#' tumourous. The user will then be asked which one indicates whether
#' the patient is healthy. Only two values are valid in the vector in total.
#' @export
#' @return A heatmap of the 100 genes with the highest variability in the matrix
#' disease component.
results_dgsa <- function(matrix_disease_component, case_tag){
  genes_sd <- apply(matrix_disease_component,1,stats::sd)
  selected_genes_sd <- names(genes_sd[order(genes_sd,decreasing = T)])[1:100]
  selected_matrix_disease_component_sd <- matrix_disease_component[selected_genes_sd,]

  # DT::datatable(selected_matrix_disease_component_sd)

  plot_dgsa(selected_matrix_disease_component_sd, case_tag)

  return(selected_genes_sd)
}
