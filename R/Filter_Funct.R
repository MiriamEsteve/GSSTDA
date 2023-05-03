#' @title PAD-S Filtering function
#' @description A filtering function for mapper that projects $$R$^n$ into $R$.
#' It calculates for each column of the matrix (each patient), its value
#' of the filtering function. Specifically, it computes
#' the vector magnitude in the \[L_{p}\] norm (as well
#' as k powers of this magnitude) of the vector resulting of
#' weighting each element of the column vector by the Z score obtained
#' by fitting a cox proportional hazard model to the level of each gene.
#' For further information see "Progression Analysis of Disease with Survival
#' (PAD-S) by SurvMap identifies different prognostic subgroups of breast
#' cancer in a large combined set of transcriptomics and methylation studies"
#' @param genes_disease_component Disease component matrix (output of the function of
#' \code{generate_disease_component}), after having selected the rows
#' corresponding to the selected genes.
#' @param p integer. It indicates the p norm to be calculated.
#' If k = 1 and p = 2, the function computes the standard
#' (Euclidean) vector magnitude of each column. For larger values of p the
#' weight of genes with larger levels is greater.
#' @param k integer. Powers of the vector magnitude. If k = 1 and p = 2,
#' the function computes the standard (Euclidean) vector magnitude
#' of each column.
#' @param cox_all_matrix A matrix with the output of the
#' \code{cox_all_genes} function that stores the information of all cox
#' proportional hazard model tests for each gene in the dataset.
#' @return A numeric vector including the values produced by the function
#' for each sample in the dataset.
#' @export
lp_norm_k_powers_surv <- function(genes_disease_component, p, k, cox_all_matrix){
  #Select "z" value (p-value) of cox_all_matrix
  cox_vector <- cox_all_matrix[rownames(genes_disease_component),"z"]

  #Prepare genes_disease_component and cox_vector
  genes_disease_component[genes_disease_component < 0] <- ifelse(!is.na(genes_disease_component[genes_disease_component < 0]),
                                                                 genes_disease_component[genes_disease_component < 0] - 1,
                                                                 genes_disease_component[genes_disease_component < 0])
  genes_disease_component[genes_disease_component > 0] <- ifelse(!is.na(genes_disease_component[genes_disease_component > 0]),
                                                                 genes_disease_component[genes_disease_component > 0] + 1,
                                                                 genes_disease_component[genes_disease_component > 0])
  cox_vector[cox_vector < 0] <- ifelse(!is.na(cox_vector[cox_vector < 0]), cox_vector[cox_vector < 0] - 1,
                                       cox_vector[cox_vector < 0])
  cox_vector[cox_vector > 0] <- ifelse(!is.na(cox_vector[cox_vector > 0]), cox_vector[cox_vector > 0] + 1,
                                       cox_vector[cox_vector > 0])

  genes_disease_component <- genes_disease_component * cox_vector

  #Calculate kp-norm
  lp_norm <- apply(genes_disease_component, 2, function(x) (sum(abs(x)^p)^(k/p)))
  return(lp_norm)

}
