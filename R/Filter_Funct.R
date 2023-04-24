#' @title SurvPAD Filtering function
#'
#' @description A filtering function for mapper that projects $$R$^n$ into $R$ that
#'
#' @param exp_matrix Matrix including the fit residual of the dataset and the healthy state model.
#' @param p Integer
#' @param k Integer
#' @param cox_all A matrix with the output of the \code{cox_all_genes} function that stores the information of all cox proportional hazard model tests for each gene in the dataset.
#'
#' @return A numeric vector including the values produced by the function for each sample in the dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' lp_norm_k_powers(disease_state,2,1, cox_data)}
#'
lp_norm_k_powers_surv <- function(exp_matrix, p, k, cox_all){
  if(is.matrix(exp_matrix)){
    cox_vector <- cox_all[rownames(exp_matrix),"z"]

    #Prepare exp_matrix and cox_vector
    exp_matrix[exp_matrix < 0] <- ifelse(!is.na(exp_matrix[exp_matrix < 0]), exp_matrix[exp_matrix < 0] - 1, exp_matrix[exp_matrix < 0])
    exp_matrix[exp_matrix > 0] <- ifelse(!is.na(exp_matrix[exp_matrix > 0]), exp_matrix[exp_matrix > 0] + 1, exp_matrix[exp_matrix > 0])
    cox_vector[cox_vector < 0] <- ifelse(!is.na(cox_vector[cox_vector < 0]), cox_vector[cox_vector < 0] - 1, cox_vector[cox_vector < 0])
    cox_vector[cox_vector > 0] <- ifelse(!is.na(cox_vector[cox_vector > 0]), cox_vector[cox_vector > 0] + 1, cox_vector[cox_vector > 0])

    #en caso de = 0????

    exp_matrix <- exp_matrix * cox_vector
    lp_norm <- apply(exp_matrix, 2, function(x) (sum(abs(x)^p)^(k/p)))
    return(lp_norm)
  }else{
    stop("A matrix object should be provided.")
  }
}
