#' @title Survival analysis based on gene expression levels.
#' @description It carries out univariate cox proportional hazard models for
#' the expression levels of each gene included in the provided dataset (eData)
#' and their link with relapse-free or overall survival.
#' @param expression_vector_disease Expression data for disease samples
#' @param time_vector Numeric vector that includes time to the event information
#' @param event_vector Numeric vector that indicates if relapse or death
#' have been produced (0 and 1s).
#' @return A matrix with the results of the application of proportional
#' hazard models using the expression levels of each gene as covariate.
#' @import survival
#' @examples
#' \dontrun{
#' cox_all_genes(expression_vector_disease,time_vector,event_vector)
#' }
cox_all_genes <- function(expression_vector_disease, time_vector, event_vector){
  pb <- utils::txtProgressBar(min = 0, max = nrow(expression_vector_disease), style = 3)
  list_out <- list()
  for(i in 1:nrow(expression_vector_disease)){
    utils::setTxtProgressBar(pb, i)
    temp <- summary(survival::coxph(survival::Surv(time_vector,as.numeric(event_vector))~expression_vector_disease[i,]))$coefficients[1,]
    list_out[[i]] <- temp
  }
  df_out <- data.frame(do.call("rbind",list_out))
  colnames(df_out) <-  c("coef","exp_coef","se_coef","z","Pr_z")
  rownames(df_out) <- rownames(expression_vector_disease)
  df_out <- as.matrix(df_out)
  return(df_out)
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
#' @param disease_component Disease component matrix (output of the function
#' \code{generate_disease_component}). The names of the rows must be the names
#' of the genes.
#' @param p_Data Data.frame with the phenotype data.
#' @param status_Col_Name Column of p_Data for sample filtering.
#' @param status_Value Value of Status_Col_Name for the sample filtering
#' co-variate.
#' @param cox_all Output from the \code{cox_all_genes} function. Data.frame with
#' information on the relationship between genes and survival.
#' @param n_top Number of genes to be selected (those with the highest
#' product value).
#' @param type_sel Option. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis).
#' @return Character vector with the names of the selected genes.
#' @examples
#' \dontrun{
#' gene_selection_surv(disease_component, p_Data, status_Col_Name, status_Value, cox_all, n_top)
#' }
gene_selection_surv <- function(disease_component, p_Data, status_Col_Name, status_Value, cox_all, n_top, type_sel = c("Top_Bot","Abs")){
  cox_all <- cox_all[rownames(disease_component),]
  if(type_sel == "Top_Bot"){
    print("Is Top_Bot")
    probes_test <- (apply(disease_component[,p_Data[,status_Col_Name] == status_Value],
                          1,stats::sd)+1) * cox_all[,4]
    if(n_top %% 2 == 0){
      n_top <- n_top/2
      print(n_top)
    }else{
      n_top <- (n_top + 1)/2
      print(n_top)
    }
    selected_probes <- names(c(probes_test[order(probes_test,decreasing = T)][1:n_top],probes_test[order(probes_test,decreasing = F)][1:n_top]))
  }else if(type_sel == "Abs"){
    print("Is Abs")
    probes_test <- (apply(disease_component[, p_Data[,status_Col_Name] == status_Value],1,stats::sd)+1) * abs(cox_all[,4])
    selected_probes <- names(probes_test[order(probes_test,decreasing = T)])[1:n_top]
  }
  return(selected_probes)
}


#' @title Gene selector based on association to survival.
#'
#' @param cox_all A matrix output from the \code{cox_all_genes} function
#' @param quantile_threshold A two element vector indicating the bottom and top quantile for gene selection.
#'
#' @return Returns a list of genes that present z-values above or below the selected quantile thresholds.
#' @export
#'
#' @examples
#' \dontrun{
#' get_survival_related_genes(cox_all)
#' }
get_survival_related_genes <- function(cox_all, quantile_threshold = c(0.05,0.95)){
  genes_asso_surv <- rownames(cox_all[cox_all[,"z"] < stats::quantile(cox_all[,"z"],probs = quantile_threshold[1]) | cox_all[,"z"] > stats::quantile(cox_all[,"z"],probs = quantile_threshold[2]),])
  return(genes_asso_surv)

}
