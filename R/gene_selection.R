#' @title Survival analysis based on gene expression levels.
#' @description It carries out univariate cox proportional hazard models for
#' the expression levels of each gene included in the provided dataset (matrix_disease_component)
#' and their link with relapse-free or overall survival.
#' @param matrix_disease_component Expression data for disease samples
#' @param survival_time Numeric vector that includes time to the event information
#' @param survival_event Numeric vector that indicates if relapse or death
#' have been produced (0 and 1s).
#' @return A matrix with the results of the application of proportional
#' hazard models using the expression levels of each gene as covariate.
#' @import survival
#' @examples
#' \dontrun{
#' cox_all_genes(matrix_disease_component,survival_time,survival_event)
#' }
cox_all_genes <- function(matrix_disease_component, survival_time, survival_event){
  list_out <- list()
  for(i in 1:nrow(matrix_disease_component)){
    temp <- summary(survival::coxph(survival::Surv(survival_time,survival_event)~matrix_disease_component[i,]))$coefficients[1,]
    list_out[[i]] <- temp
  }
  df_out <- data.frame(do.call("rbind",list_out))
  colnames(df_out) <-  c("coef","exp_coef","se_coef","z","Pr_z")
  rownames(df_out) <- rownames(matrix_disease_component)
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
#' @param matrix_disease_component Disease component matrix (output of the function
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
gene_selection_surv <- function(matrix_disease_component, p_Data, status_Col_Name, status_Value, cox_all, n_top, type_sel = c("Top_Bot","Abs")){
  cox_all <- cox_all[rownames(disease_component),]
  if(type_sel == "Top_Bot"){
    print("Is Top_Bot")
    probes_test <- (apply(disease_component[,p_Data[,status_Col_Name] == status_Value],
                          1,stats::sd)+1) * cox_all[,4]
    if(n_top %% 2 == 0){
      n_top <- n_top/2
    }else{
      n_top <- (n_top + 1)/2
    }
    selected_probes <- names(c(probes_test[order(probes_test,decreasing = T)][1:n_top],probes_test[order(probes_test,decreasing = F)][1:n_top]))
  }else if(type_sel == "Abs"){
    print("Is Abs")
    probes_test <- (apply(disease_component[, p_Data[,status_Col_Name] == status_Value],1,stats::sd)+1) * abs(cox_all[,4])
    selected_probes <- names(probes_test[order(probes_test,decreasing = T)])[1:n_top]
  }
  return(selected_probes)
}
