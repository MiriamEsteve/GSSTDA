#' @title Gene selection based on variability
#'
#' @param disease_component A matrix with the disease component data of the entire dataset obtained by \code{generate_disease_component} function.
#' @param percentile Percentile for gene selection.
#'
#' @return A vector containing the names of the selected genes.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_selection(disease_component, 0.99)
#' }
gene_selection <- function(disease_component, percentile = 0.85){
  max_abs_5_95 <- apply(disease_component, 1, function(x) base::max(base::abs(stats::quantile(x, c(0.05,0.95)))))
  selected_genes <- names(max_abs_5_95[max_abs_5_95 > stats::quantile(max_abs_5_95, percentile)])
  return(selected_genes)
}



#' @title gene_selection_surv
#'
#' @description Selects markers for mapper based on the product of standard deviation plus one times the Z score obtained through fitting a cox proportional hazard model to the expression of each gene.
#'
#' @param disease_component A matrix with the disease component data of the entire dataset obtained by \code{generate_disease_component} function.
#' @param p_Data Data.frame with the phenotype data.
#' @param status_Col_Name Column for sample filtering.
#' @param status_Value Value for the sample filtering co-variate.
#' @param cox_all Output from the \code{cox_all_genes} function.
#' @param n_top Number of markers presenting the maximum values for the product.
#' @param type_sel Select top from bottom and top or from absolute value.
#'
#' @return vector de los genes seleccionados
#' @export
#'
#' @examples
#' \dontrun{
#' gene_selection_surv(disease_component, p_Data, status_Col_Name, status_Value, cox_all, n_top)
#' }
gene_selection_surv <- function(disease_component, p_Data, status_Col_Name, status_Value, cox_all, n_top, type_sel = c("Top_Bot","Abs")){
  cox_all <- cox_all[rownames(disease_component),]
  if(type_sel == "Top_Bot"){
    print("Is Top_Bot")
    probes_test <- (apply(disease_component[,p_Data[,status_Col_Name] == status_Value],1,stats::sd)+1) * cox_all[,4]
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


#' @title Survival analysis based on gene expression levels.
#'
#' @description Carries out univariate cox protportional hazard models for the expression levels of each gene included in the dataset and its link with relapse-free or overall survival.
#'
#' @param expression_vector_disease Expression data for disease samples
#' @param time_vector Vector including time to relapse or time to death information
#' @param event_vector Numeric vector indicating if relapse or death have been produced.
#'
#' @return A matrix with the results of the application of proportional hazard models using the expression levels of eahc gene as covariate.
#' @export
#'
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
