#' @title Survival analysis based on gene expression levels.
#' @description It carries out univariate cox proportional hazard models for
#' the expression levels of each gene included in the provided dataset (matrix_disease_component)
#' and their link with relapse-free or overall survival.
#' @param case_disease_component Expression data for disease samples.
#' @param survival_time Numeric vector that includes time to the event information
#' @param survival_event Numeric vector that indicates if relapse or death
#' have been produced (0 and 1s).
#' @return A matrix with the results of the application of proportional
#' hazard models using the expression levels of each gene as covariate.
#' The \code{coef} column corresponds to the regression coefficient; the
#' \code{exp_coef} column corresponds to the value of e^coef  (which is
#' interpreted as the odds ratio); the \code{se_coef} column corresponds
#' to the standard error of each coefficient; the \code{z} column corresponds
#' to the value of coef/se_coef (the higher the Z value, the higher the
#' significance of the variable) and the \code{Pr_z} column corresponds to
#' the p-value for each Z value.
#' @export
#' @import survival
#' @examples
#' \dontrun{
#' cox_all_genes(case_disease_component,survival_time,survival_event)
#' }
cox_all_genes <- function(case_disease_component, survival_time, survival_event){
  print("Calculating the matrix of Zcox")
  pb <- utils::txtProgressBar(min = 0, max = nrow(case_disease_component), style = 3)

  list_out <- list()
  for(i in 1:nrow(case_disease_component)){
    utils::setTxtProgressBar(pb, i)

    temp <- summary(survival::coxph(survival::Surv(survival_time,survival_event)~case_disease_component[i,]))$coefficients[1,]
    list_out[[i]] <- temp
  }
  cox_all_matrix <- data.frame(do.call("rbind",list_out))
  colnames(cox_all_matrix) <-  c("coef","exp_coef","se_coef","z","Pr_z")
  rownames(cox_all_matrix) <- rownames(case_disease_component)
  cox_all_matrix <- as.matrix(cox_all_matrix)
  return(cox_all_matrix)
}


#' @title Gene selection based on variability and the relationship
#' to survival.
#' @description
#' It selects genes for mapper based on the product of standard deviation
#' of the rows (genes) in the disease component matrix
#' plus one times the Z score obtained by fitting a cox proportional
#' hazard model to the level of each gene. For further information see
#' "Topology based data analysis identifies a subgroup of breast cancers
#' with a unique mutational profile and excellent survival"
#' @param case_disease_component Disease component matrix (output of the function
#' \code{generate_disease_component}) having selected only the columns
#' belonging to disease samples. The names of the rows must be the names of the genes.
#' @param cox_all_matrix Output from the \code{cox_all_genes} function. Data.frame with
#' information on the relationship between genes and survival.
#' @param gen_select_type Option. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis).
#' @param num_gen_select Number of genes to be selected (those with the highest
#' product value).
#' @return Character vector with the names of the selected genes.
#' @export
#' @examples
#' \dontrun{
#' gene_selection_surv(case_disease_component, cox_all_matrix, gen_select_type, num_gen_select)
#' }
gene_selection_surv <- function(case_disease_component, cox_all_matrix, gen_select_type, num_gen_select){
  # Same operation to both methods
  probes_test <- apply(case_disease_component, 1,stats::sd)+1

  if(gen_select_type == "abs"){
    probes_test <- probes_test * abs(cox_all_matrix[,4])
    genes_selected <- names(probes_test[order(probes_test,decreasing = T)])[1:num_gen_select]
  }else{
    probes_test <- probes_test * cox_all_matrix[,4]
    if(num_gen_select %% 2 == 0){ num_gen_select <- num_gen_select/2}
    else{ num_gen_select <- (num_gen_select + 1)/2}

    genes_selected <- names(c(probes_test[order(probes_test,decreasing = T)][1:num_gen_select],
                              probes_test[order(probes_test,decreasing = F)][1:num_gen_select]))
  }
  return(genes_selected)
}
