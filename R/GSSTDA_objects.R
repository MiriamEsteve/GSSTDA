#' @title Disease-Specific Genomic Analysis
#' @description Disease-Specific Genomic Analysis (dsga).
#' This analysis, developed by Nicolau *et al.*, allows the calculation of
#' the "disease component" of a expression matrix which consists of, through
#' linear models, eliminating the part of the data  that is considered normal
#' or healthy and keeping only the component that is due to the disease. It
#' is intended to precede other techniques like classification or clustering.
#' For more information see *Disease-specific genomic analysis: identifying
#' the signature of pathologic biology* (doi: 10.1093/bioinformatics/btm033).
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' @param survival_time Numerical vector of the same length as the number of
#' columns of \code{full_data}. In addition, the patients must be in the same
#' order as in \code{full_data}. For the patients whose sample is pathological
#' should be indicated the time between the disease diagnosis and event
#' (death, relapse or other). If the event has not occurred, it should be
#' indicated the time until the end of follow-up. Patients whose sample is
#' from healthy tissue must have an NA value
#' @param survival_event Numerical vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. For the the patients with pathological sample should
#' be indicated whether the event has occurred (1) or not (0). Only these
#' values are valid and healthy patients must have an NA value.
#' @param case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
#' @param control_tag Tag of the healthy sample.E.g. "T"
#' @param gamma A parameter that indicates the magnitude of the noise assumed in
#' the flat data matrix for the generation of the Healthy State Model. If it
#' takes the value `NA` the magnitude of the noise is assumed to be unknown.
#' By default gamma is unknown.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{dsga} object. It contains: the \code{full_data} without
#' NAN's values, the label designated for healthy samples (\code{control_tag}),
#' the \code{case_tag} vector without NAN's values, the \code{survival_event},
#' the the \code{survival_time} the matrix with the normal space (linear space
#' generated from normal tissue samples) and the matrix of the disease
#' components (the transformed full_data matrix from which the normal component
#' has been removed).
#' @export
#' @examples
#' \donttest{
#' dsga_obj <- dsga(full_data,  survival_time, survival_event, case_tag)}
dsga <- function(full_data, survival_time, survival_event, case_tag, control_tag = NA, gamma = NA, na.rm = TRUE){
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)

  #Select the control_tag
  return_check <- check_vectors(full_data, survival_time, survival_event, case_tag, control_tag, na.rm)
  control_tag <- return_check[[1]]
  full_data <- return_check[[2]]
  survival_event <- return_check[[3]]
  survival_time <- return_check[[4]]
  case_tag <- return_check[[5]]

  ################### BLOCK I: Pre-process. dsga (using "NT" control_tag) ##############################
  message("\nBLOCK I: The pre-process dsga is started")
  #   Select the normal tissue data gene expression matrix.
  normal_tiss <- full_data[,which(case_tag == control_tag)]

  #   Obtain the gene expression matrix containing the flattened version of the vectors.
  matrix_flatten_normal_tiss <- flatten_normal_tiss(normal_tiss)
  #   Obtain the normal space
  normal_space <- denoise_rectangular_matrix(matrix_flatten_normal_tiss, gamma)
  #   Obtain the disease component of the normal_space
  matrix_disease_component <- generate_disease_component(full_data, normal_space)

  message("\nBLOCK I: The pre-process dsga is finished\n")

  ############################################  Create the object #########################################
  dsga_object <- list("full_data" = full_data,
                      "control_tag" = control_tag,
                      "case_tag" = case_tag,
                      "survival_event" = survival_event,
                      "survival_time" = survival_time,
                      "normal_space" = normal_space,
                      "matrix_disease_component" = matrix_disease_component)

  class(dsga_object) <- "dsga_object"

  return(dsga_object)
}


#' @title Gene selection and filter function
#' @description Gene selection and calculation of filter function values.
#' After fitting a Cox proportional hazard model to each gene, this function
#' makes a selection of genes according to both their variability within
#' the database and their relationship with survival. Subsequently, with the
#' genes selected, the values of the filtering functions are calculated for
#' each patient. The filter function allows to summarise each vector of each
#' individual in a single data. This function takes into account the survival
#' associated with each gene. In particular, the implemented filter function
#' performs the vector magnitude in the Lp norm (as well as k powers
#' of this magnitude) of the vector resulting of weighting each element of
#' the column vector by the Z score obtained in the cox proportional
#' hazard model.
#' @param data_object Object with:
#' - full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' - survival_time Numerical vector of the same length as the number of columns
#' of \code{full_data}. In addition, the patients must be in the same
#' order as in \code{full_data}. For the patients whose sample is pathological
#' should be indicated the time between the disease diagnosis and event
#' (death, relapse or other). If the event has not occurred, it should be
#' indicated the time until the end of follow-up. Patients whose sample is
#' from healthy tissue must have an NA value
#' - survival_event Numerical vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. For the the patients with pathological sample should
#' be indicated whether the event has occurred (1) or not (0). Only these
#' values are valid and healthy patients must have an NA value.
#' - case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
#' @param gen_select_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param percent_gen_select Percentage (from zero to one hundred) of genes
#' to be selected to be used in mapper. 10 default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gene_selection_object}. It contains:
#' - the \code{full_data} without NAN's values (\code{data})
#' - the \code{cox_all_matrix} (a matrix with the results of the application of
#' proportional hazard models: with the regression coefficients, the odds ratios,
#' the standard errors of each coefficient, the Z values (coef/se_coef) and
#' the p-values for each Z value)
#' - a vector with the name of the selected genes
#' - the matrix of disease components with only the rows of the selected genes
#' (\code{genes_disease_component})
#' - and the vector of the values of the filter function.
#' @export
#' @examples
#' \donttest{
#' data_object <- list("full_data" = full_data, "survival_time" = survival_time,
#' "survival_event" = survival_event, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#' gene_selection_obj <- gene_selection(data_object,
#' gen_select_type ="top_bot", percent_gen_select=10)}
gene_selection <- function(data_object, gen_select_type,
                          percent_gen_select, na.rm = TRUE){
  UseMethod("gene_selection")
}

#' @title Private gene_selection_
#' @description Private function to gene selection
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' @param survival_time Numerical vector of the same length as the number of
#' columns of \code{full_data}. In addition, the patients must be in the same
#' order as in \code{full_data}. For the patients whose sample is pathological
#' should be indicated the time between the disease diagnosis and event
#' (death, relapse or other). If the event has not occurred, it should be
#' indicated the time until the end of follow-up. Patients whose sample is
#' from healthy tissue must have an NA value
#' @param survival_event Numerical vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. For the the patients with pathological sample should
#' be indicated whether the event has occurred (1) or not (0). Only these
#' values are valid and healthy patients must have an NA value.
#' @param control_tag_cases Numeric vector with the indices of the columns
#' corresponding to the healthy sample patients.
#' @param gen_select_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param num_gen_select Number of genes to be selected to be used in mapper.
#' @param matrix_disease_component Optional, only necessary in case of gene
#' selection after dsga has been performed. Matrix of the disease components
#' (the transformed \code{full_data} matrix from which the normal component has
#' been removed) from the \code{dsga_function}.
#' @return A \code{gene_selection_object}. It contains:
#' - the matrix with which the gene selection has been performed without NAN's
#' values (\code{data}). It is the \code{matrix_disease_component} in case it has been
#' performed from a \code{dsga_object} or \code{full_data} in the opposite case.
#' - the \code{cox_all_matrix} (a matrix with the results of the application of
#' proportional hazard models: with the regression coefficients, the odds ratios,
#' the standard errors of each coefficient, the Z values (coef/se_coef) and
#' the p-values for each Z value)
#' - a vector with the name of the selected genes
#' - the matrix of disease components with only the rows of the selected genes
#' (\code{genes_disease_component})
#' - and the vector of the values of the filter function.
#' @export
#' @examples
#' \donttest{
#' gen_select_type <- "Top_Bot"
#' percent_gen_select <- 10
#' control_tag_cases <- which(case_tag == "NT")
#' gene_selection_obj <- gene_selection_(full_data, survival_time, survival_event,
#' control_tag_cases, gen_select_type ="top_bot", num_gen_select = 10)}
gene_selection_ <- function(full_data, survival_time, survival_event,
                           control_tag_cases, gen_select_type, num_gen_select,
                           matrix_disease_component = NULL){

  message("\nBLOCK II: The gene selection is started\n")

  if(is.null(matrix_disease_component)) {
    matrix_disease_component <- full_data
  }
  #Remove NAN's values (case_tag == control_tag) of survival_time and survival_event
  survival_time <- survival_time[-control_tag_cases]
  survival_event <- survival_event[-control_tag_cases]
  #Select the disease component of the "T" control_tag
  case_full_data <- full_data[,-control_tag_cases]
  case_disease_component <- matrix_disease_component[,-control_tag_cases]

  # Univariate cox proportional hazard models for the expression levels of each gene included in the
  #provided dataset
  cox_all_matrix <- cox_all_genes(case_full_data, survival_time, survival_event)

  #Selects genes for mapper
  genes_selected <- gene_selection_surv(case_disease_component, cox_all_matrix, gen_select_type,
                                        num_gen_select)

  # Select genes in matrix_disease_component or full_data (if don't apply Block I)
  genes_disease_component <- matrix_disease_component[genes_selected,]

  # Filter the genes_disease_component
  filter_values <- lp_norm_k_powers_surv(genes_disease_component, 2, 1, cox_all_matrix)

  message("\nBLOCK II: The gene selection is finished\n")

  gene_selection_object <- list( "data" = matrix_disease_component,
                                "cox_all_matrix" = cox_all_matrix,
                                "genes_selected" = genes_selected,
                                "genes_disease_component" = genes_disease_component,
                                "filter_values" = filter_values
  )
  class(gene_selection_object) <- "gene_selection_object"

  return(gene_selection_object)
}

#' @title gene_selection_classes.dsga_object
#'
#' @description Private function to select Gene with dsga object
#' @param data_object dsga object information
#' @param gen_select_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param percent_gen_select Percentage (from zero to one hundred) of genes
#' to be selected to be used in mapper. 10 default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gene_selection} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#' @export
#' @examples
#' \donttest{
#' dsga_obj <- dsga(full_data, survival_time, survival_event, case_tag, na.rm = "checked")
#'
#' gene_selection_object <- gene_selection(dsga_obj, gen_select_type ="top_bot",
#'                                       percent_gen_select = 10)}
gene_selection.dsga_object <- function(data_object, gen_select_type, percent_gen_select, na.rm = TRUE){
  print(class(data_object))

  matrix_disease_component <- data_object[["matrix_disease_component"]]
  #Check and obtain gene selection (we use in the gene_select_surv)
  num_gen_select <- check_gene_selection(nrow(matrix_disease_component),
                                         gen_select_type, percent_gen_select)

  full_data <- data_object[["full_data"]]
  control_tag <- data_object[["control_tag"]]
  survival_event <- data_object[["survival_event"]]
  survival_time <- data_object[["survival_time"]]
  case_tag <- data_object[["case_tag"]]

  control_tag_cases <- which(case_tag == control_tag)
  gene_selection_object <- gene_selection_(full_data, survival_time, survival_event,
                                         control_tag_cases, gen_select_type, num_gen_select,
                                         matrix_disease_component)

  return(gene_selection_object)
}

#' @title gene_selection_classes.default
#'
#' @description Private function to select Gene without dsga process
#' @param data_object Object with:
#' - full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' - survival_time Numerical vector of the same length as the number of
#' columns of \code{full_data}. In addition, the patients must be in the same
#' order as in \code{full_data}. For the patients whose sample is pathological
#' should be indicated the time between the disease diagnosis and event
#' (death, relapse or other). If the event has not occurred, it should be
#' indicated the time until the end of follow-up. Patients whose sample is
#' from healthy tissue must have an NA value
#' - survival_event Numerical vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. For the the patients with pathological sample should
#' be indicated whether the event has occurred (1) or not (0). Only these
#' values are valid and healthy patients must have an NA value.
#' - case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
#' @param gen_select_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param percent_gen_select Percentage (from zero to one hundred) of genes
#' to be selected to be used in mapper. 10 default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gene_selection} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#' @export
#' @examples
#' \donttest{
#' data_object <- list("full_data" = full_data, "survival_time" = survival_time,
#' "survival_event" = survival_event, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#' gene_selection_object <- gene_selection(data_object, gen_select_type ="top_bot",
#'                                       percent_gen_select = 10)}
gene_selection.default <- function(data_object, gen_select_type, percent_gen_select, na.rm = TRUE){
  full_data <- data_object[["full_data"]]
  survival_event <- data_object[["survival_event"]]
  survival_time <- data_object[["survival_time"]]
  case_tag <- data_object[["case_tag"]]

  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)
  #Select the control_tag. This do it inside of the dsga function
  #Check and obtain gene selection (we use in the gene_select_surv)
  num_gen_select <- check_gene_selection(nrow(full_data), gen_select_type, percent_gen_select)

  #Select the control_tag
  return_check <- check_vectors(full_data, survival_time, survival_event, case_tag, na.rm)
  control_tag <- return_check[[1]]
  full_data <- return_check[[2]]
  survival_event <- return_check[[3]]
  survival_time <- return_check[[4]]
  case_tag <- return_check[[5]]

  control_tag_cases <- which(case_tag == control_tag)

  gene_selection_object <- gene_selection_(full_data, survival_time, survival_event,
                                         control_tag_cases, gen_select_type, num_gen_select)

  return(gene_selection_object)
}


#' @title Mapper object
#'
#' @description TDA are persistent homology and mapper. Persistent homology
#' borrows ideas from abstract algebra to identify particular aspects
#' related to the shape of the data such as the number of connected
#' components and the presence of higher-dimensional holes, whereas
#' mapper condenses the information of high-dimensional datasets into
#' a combinatory graph or simplicial complex that is referred to as
#' the skeleton of the dataset. This implementation is the mapper of one
#' dimension, i.e. using only one filter function value.
#' @param data Input matrix whose columns correspond to the individuals
#' and rows to the features.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param num_intervals Number of intervals used to create the first sample
#' partition based on filtering values. 5 default option.
#' @param percent_overlap Percentage of overlap between intervals. Expressed
#' as a percentage. 40 default option.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("correlation") and euclidean ("euclidean"). "correlation"
#' default option.
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' "hierarchical" default option.
#' @param num_bins_when_clustering Number of bins to generate the
#' histogram employed by the standard optimal number of cluster finder
#' method. Parameter not necessary if the "optimal_clustering_mode" option
#' is "silhouette" or the "clustering_type" is "PAM". 10 default option.
#' @param linkage_type Linkage criteria used in hierarchical clustering.
#' Choose between "single" for single-linkage clustering, "complete" for
#' complete-linkage clustering or "average" for average linkage clustering
#' (or UPGMA). Only necessary for hierarchical clustering.
#' "single" default option.
#' @param optimal_clustering_mode Method for selection optimal number of
#' clusters. It is only necessary if the chosen type of algorithm is
#' hierarchical. In this case, choose between "standard" (the method used
#' in the original mapper article) or "silhouette". In the case of the PAM
#' algorithm, the method will always be "silhouette".
#' @param silhouette_threshold Minimum value of \eqn{\overline{s}}{s-bar} that a set of
#' clusters must have to be chosen as optimal. Within each interval of the
#' filter function, the average silhouette values \eqn{\overline{s}}{s-bar} are computed
#' for all possible partitions from $2$ to $n-1$, where $n$ is the number of
#' samples within a specific interval. The $n$ that produces the highest value
#' of \eqn{\overline{s}}{s-bar} and that exceeds a specific threshold is selected as the
#' optimum number of clusters. If no partition produces an \eqn{\overline{s}}{s-bar}
#' exceeding the chosen threshold, all samples are then assigned to a unique
#' cluster. The default value is $0.25$. The threshold of $0.25$ for
#' \eqn{\overline{s}}{s-bar} has been chosen based on standard practice, recognizing it
#' as a moderate value that reflects adequate separation and cohesion within
#' clusters.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{mapper_obj} object. It contains the values of the intervals
#' (interval_data), the samples included in each interval (sample_in_level),
#' information about the cluster to which the individuals in each interval
#' belong (clustering_all_levels), a list including the individuals contained
#' in each detected node (node_samples), their size (node_sizes), the
#' average of the filter function values of the individuals of each node
#' (node_average_filt) and the adjacency matrix linking the nodes (adj_matrix).
#' @export
#' @examples
#' \donttest{
#' control_tag_cases <- which(case_tag == "NT")
#' gene_selection_object <- gene_selection_(full_data, survival_time, survival_event,
#' control_tag_cases, gen_select_type ="top_bot", num_gen_select = 10)
#'
#' mapper_object <- mapper(data = gene_selection_object[["genes_disease_component"]],
#' filter_values = gene_selection_object[["filter_values"]],
#' num_intervals = 5,
#' percent_overlap = 40, distance_type = "correlation",
#' clustering_type = "hierarchical",
#' linkage_type = "single")}
mapper <- function(data, filter_values, num_intervals = 5, percent_overlap = 40,
                   distance_type = "correlation", clustering_type = "hierarchical",
                   num_bins_when_clustering = 10, linkage_type = "single",
                   optimal_clustering_mode = NA, silhouette_threshold = 0.25, na.rm=TRUE){
  # Don't call by GSSTDA function
  if (na.rm != "checked"){
    # Check the data introduces
    data <- check_full_data(data)
    # Check mapper arguments
    check_return <- check_arg_mapper(data, filter_values, distance_type, clustering_type,
                                     linkage_type, optimal_clustering_mode, silhouette_threshold)

    data <- check_return[[1]]
    filter_values <- check_return[[2]]
    optimal_clustering_mode <- check_return[[3]]
  }

  mapper_object_ini <- list("data" = data,
                            "filter_values" = filter_values,
                            "num_intervals" = num_intervals,
                            "percent_overlap" = percent_overlap/100,
                            "distance_type" = distance_type,
                            "optimal_clustering_mode" = optimal_clustering_mode,
                            "num_bins_when_clustering" = num_bins_when_clustering,
                            "clustering_type" = clustering_type,
                            "linkage_type" = linkage_type,
                            "optimal_clustering_mode" = optimal_clustering_mode,
                            "silhouette_threshold" = silhouette_threshold)

  class(mapper_object_ini) <- "mapper_initialization"

  mapper_object <- one_D_Mapper(mapper_object_ini)

  return(mapper_object)
}


