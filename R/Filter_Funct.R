#' @title PAD-S Filtering function
#' @description A filtering function for mapper that projects $$R$^n$ into $R$.
#' It calculates for each colum of the matrix (each patient), its value
#' of the filtering function. Specifically, it computes
#' the vector magnitude in the \[L_{p}\] norm (as well
#' as k powers of this magnitude) of the vector resulting of
#' weighting each element of the column vector by the Z score obtained
#' by fitting a cox proportional hazard model to the level of each gene.
#' For further information see "Progression Analysis of Disease with Survival
#' (PAD-S) by SurvMap identifies different prognostic subgroups of breast
#' cancer in a large combined set of transcriptomics and methylation studies"
#' @param exp_matrix Disease component matrix (output of the function of
#' \code{generate_disease_component}), after having selected the rows
#' corresponding to the selected genes.
#' @param p integer. It indicates the p norm to be calculated.
#' If k = 1 and p = 2, the function computes the standard
#' (Euclidean) vector magnitude of each column. For larger values of p the
#' weight of genes with larger levels is greater.
#' @param k integer. Powers of the vector magnitude. If k = 1 and p = 2,
#' the function computes the standard (Euclidean) vector magnitude
#' of each column.
#' @param cox_all A matrix with the output of the
#' \code{cox_all_genes} function that stores the information of all cox
#' proportional hazard model tests for each gene in the dataset.
#' @return A numeric vector including the values produced by the function
#' for each sample in the dataset.
lp_norm_k_powers_surv <- function(exp_matrix, p, k, cox_all){

  #Redundante rownames(exp-MATRIX). z == evento
  cox_vector <- cox_all[rownames(exp_matrix),"z"]

  exp_matrix <- exp_matrix * cox_vector

  #Prepare exp_matrix and cox_vector
  exp_matrix[exp_matrix < 0] <- ifelse(!is.na(exp_matrix[exp_matrix < 0]), exp_matrix[exp_matrix < 0] - 1, exp_matrix[exp_matrix < 0])
  exp_matrix[exp_matrix >= 0] <- ifelse(!is.na(exp_matrix[exp_matrix >= 0]), exp_matrix[exp_matrix >= 0] + 1, exp_matrix[exp_matrix >= 0])

  lp_norm <- apply(exp_matrix, 2, function(x) (sum(abs(x)^p)^(k/p)))
  return(lp_norm)

}
