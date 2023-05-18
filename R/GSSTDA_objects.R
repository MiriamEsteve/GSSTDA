#' @title DGSA
#'
#' @description Disease-Specific Genomic Analysis (DGSA)
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' @param survival_time Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated the time between
#' disease diagnosis and death (if not dead until the end of follow-up)
#' and healthy patients must have an NA value.
#' @param survival_event Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated whether
#' the patient has died (1) or not (0). Only these values are valid
#' and healthy patients must have an NA value.
#' @param case_tag Character vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' It must be indicated for each patient whether he/she is healthy or not.
#' One value should be used to indicate whether the patient is healthy and
#' another value should be used to indicate whether the patient's sample is
#' tumourous. The user will then be asked which one indicates whether
#' the patient is healthy. Only two values are valid in the vector in total.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{DGSA} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#' @export
#' @examples
#' \dontrun{
#' DGSA_obj <- DGSA(full_data,  survival_time, survival_event, case_tag)}
DGSA <- function(full_data,  survival_time, survival_event, case_tag, na.rm = TRUE){
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)

  #Select the control_tag
  return_check <- check_vectors(full_data, survival_time, survival_event, case_tag, na.rm)
  control_tag <- return_check[[1]]
  full_data <- return_check[[2]]
  survival_event <- return_check[[3]]
  survival_time <- return_check[[4]]
  case_tag <- return_check[[5]]

  ################### BLOCK I: Pre-process. DGSA (using "NT" control_tag) ##############################
  print("\nBLOCK I: The pre-process DGSA is started")
  #   Select the normal tissue data gene expression matrix.
  normal_tiss <- full_data[,which(case_tag == control_tag)]

  #   Obtain the gene expression matrix containing the flattened version of the vectors.
  matrix_flatten_normal_tiss <- flatten_normal_tiss(normal_tiss)
  #   Obtain the normal space
  normal_space <- denoise_rectangular_matrix(matrix_flatten_normal_tiss)
  #   Obtain the disease component of the normal_space
  matrix_disease_component <- generate_disease_component(full_data, normal_space)

  print("\nBLOCK I: The pre-process DGSA is finished")

  ############################################  Create the object #########################################
  DGSA_object <- list("full_data" = full_data,
                      "control_tag" = control_tag,
                      "case_tag" = case_tag,
                      "survival_event" = survival_event,
                      "survival_time" = survival_time,
                      "normal_space" = normal_space,
                      "matrix_disease_component" = matrix_disease_component)

  class(DGSA_object) <- "DGSA_object"

  return(DGSA_object)
}


#' @title geneSelection
#'
#' @description Gene selection
#' @param data_object Object with:
#' - full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' - survival_time Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated the time between
#' disease diagnosis and death (if not dead until the end of follow-up)
#' and healthy patients must have an NA value.
#' - survival_event Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated whether
#' the patient has died (1) or not (0). Only these values are valid
#' and healthy patients must have an NA value.
#' - case_tag Character vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' It must be indicated for each patient whether he/she is healthy or not.
#' One value should be used to indicate whether the patient is healthy and
#' another value should be used to indicate whether the patient's sample is
#' tumourous. The user will then be asked which one indicates whether
#' the patient is healthy. Only two values are valid in the vector in total.
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
#' @return A \code{DGSA} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#' @export
#' @examples
#' \dontrun{
#' geneSelection_obj <- geneSelection(data_object,
#' gen_select_type, percent_gen_select)}
geneSelection <- function(data_object, gen_select_type,
                          percent_gen_select, na.rm = TRUE){
  UseMethod("geneSelection")
}

#' @title gene_selection
#'
#' @description Private function to select Gene
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' @param survival_time Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated the time between
#' disease diagnosis and death (if not dead until the end of follow-up)
#' and healthy patients must have an NA value.
#' @param survival_event Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated whether
#' the patient has died (1) or not (0). Only these values are valid
#' and healthy patients must have an NA value.
#' @param control_tag_cases Character vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' It must be indicated for each patient whether he/she is healthy or not.
#' One value should be used to indicate whether the patient is healthy and
#' another value should be used to indicate whether the patient's sample is
#' tumourous. The user will then be asked which one indicates whether
#' the patient is healthy. Only two values are valid in the vector in total.
#' @param gen_select_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param num_gen_select Number of genes to be selected to be used in mapper.
#' @return A \code{geneSelection} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#'
#' @export
#' @examples
#' \dontrun{
#' geneSelection_obj <- gene_selection(full_data, survival_time, survival_event, control_tag_cases,
#' gen_select_type, num_gen_select)}
gene_selection <- function(full_data, survival_time, survival_event, control_tag_cases,
                           gen_select_type, num_gen_select){

  print("\nBLOCK II: The gene selection is started")
  #Remove NAN's values (case_tag == control_tag) of survival_time and survival_event
  survival_time <- survival_time[-control_tag_cases]
  survival_event <- survival_event[-control_tag_cases]
  #Select the disease component of the "T" control_tag
  case_disease_component <- full_data[,-control_tag_cases]

  # Univariate cox proportional hazard models for the expression levels of each gene included in the
  #provided dataset
  cox_all_matrix <- cox_all_genes(case_disease_component, survival_time, survival_event)

  #Selects genes for mapper
  genes_selected <- gene_selection_surv(case_disease_component, cox_all_matrix, gen_select_type,
                                        num_gen_select)

  # Select genes in matrix_disease_component or full_data (if don't apply Block I)
  genes_disease_component <- full_data[genes_selected,]

  # Filter the genes_disease_component
  filter_values <- lp_norm_k_powers_surv(genes_disease_component, 2, 1, cox_all_matrix)

  print("\nBLOCK II: The gene selection is finished")

  geneSelection_object <- list( "data" = full_data,
                                "cox_all_matrix" = cox_all_matrix,
                                "genes_selected" = genes_selected,
                                "genes_disease_component" = genes_disease_component,
                                "filter_values" = filter_values
  )
  class(geneSelection_object) <- "geneSelection_object"

  return(geneSelection_object)
}

#' @title gene_selection_classes.DGSA_object
#'
#' @description Private function to select Gene with DGSA object
#' @param data_object DGSA object information
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
#' @return A \code{geneSelection} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#'
#' @export
#' @examples
#' \dontrun{
#' geneSelection_obj <- geneSelection.DGSA_object(data_object, gen_select_type,
#'                                                        percent_gen_select)}
geneSelection.DGSA_object <- function(data_object, gen_select_type, percent_gen_select, na.rm = TRUE){
  print(class(data_object))

  matrix_disease_component <- data_object[["matrix_disease_component"]]
  #Check and obtain gene selection (we use in the gene_select_surv)
  num_gen_select <- check_gene_selection(nrow(matrix_disease_component), gen_select_type, percent_gen_select)

  control_tag <- data_object[["control_tag"]]
  survival_event <- data_object[["survival_event"]]
  survival_time <- data_object[["survival_time"]]
  case_tag <- data_object[["case_tag"]]

  control_tag_cases <- which(case_tag == control_tag)
  geneSelection_object <- gene_selection(matrix_disease_component, survival_time, survival_event,
                                         control_tag_cases, gen_select_type, num_gen_select)

  return(geneSelection_object)
}

#' @title gene_selection_classes.matrix
#'
#' @description Private function to select Gene without DGSA process
#' @param data_object Object with:
#' - full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' - survival_time Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated the time between
#' disease diagnosis and death (if not dead until the end of follow-up)
#' and healthy patients must have an NA value.
#' - survival_event Numerical vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' For the patients with tumour sample should be indicated whether
#' the patient has died (1) or not (0). Only these values are valid
#' and healthy patients must have an NA value.
#' - case_tag Character vector of the same length as the number of
#' columns of full_data. Patients must be in the same order as in full_data.
#' It must be indicated for each patient whether he/she is healthy or not.
#' One value should be used to indicate whether the patient is healthy and
#' another value should be used to indicate whether the patient's sample is
#' tumourous. The user will then be asked which one indicates whether
#' the patient is healthy. Only two values are valid in the vector in total.
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
#' @return A \code{geneSelection} object. It contains: the full_data without NAN's values,
#' the control tag of the healthy patient, the matrix with the normal space and
#' the matrix of the disease components.
#'
#' @export
#' @examples
#' \dontrun{
#' geneSelection_obj <- geneSelection.default(data_object, gen_select_type, percent_gen_select)}
geneSelection.default <- function(data_object, gen_select_type, percent_gen_select, na.rm = TRUE){
  full_data <- data_object[["full_data"]]
  survival_event <- data_object[["survival_event"]]
  survival_time <- data_object[["survival_time"]]
  case_tag <- data_object[["case_tag"]]

  print("gene_selection_classes")
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)
  #Select the control_tag. This do it inside of the DGSA function
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

  geneSelection_object <- gene_selection(full_data, survival_time, survival_event,
                                         control_tag_cases, gen_select_type, num_gen_select)

  return(geneSelection_object)
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
#' @param full_data Input matrix whose columns correspond to the individuals
#' and rows to the features.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param num_intervals Number of intervals used to create the first sample
#' partition based on filtering values. 5 default option.
#' @param percent_overlap Percentage of overlap between intervals. Expressed
#' as a percentage. 40 default option.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("cor") and euclidean ("euclidean"). "cor"
#' default option.
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' "hierarchical" default option.
#' @param num_bins_when_clustering Number of bins to generate the
#' histogram employed by the standard optimal number of cluster finder
#' method. Parameter not necessary if the "optimal_clust_mode" option
#' is "silhouette" or the "clust_type" is "PAM". 10 default option.
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
#' \dontrun{
#' num_rows <- 100
#' full_data <- data.frame( x=2*cos(1:num_rows), y=sin(1:num_rows) )
#' filter_values <- list(2*cos(1:num_rows))
#' mapper_obj <- mapper(full_data, filter_values, num_intervals = 4,
#'                      percent_overlap = 0.5, distance_type = "euclidean",
#'                      num_bins_when_clustering = 8,
#'                      clustering_type = "hierarchical",
#'                      linkage_type = "single")}
mapper <- function(full_data, filter_values, num_intervals = 5, percent_overlap = 40,
                   distance_type = "cor", clustering_type = "hierarchical",
                   num_bins_when_clustering = 10, linkage_type = "single",
                   optimal_clustering_mode="", na.rm=TRUE){
  # Don't call by GSSTDA function
  if (na.rm != "checked"){
    # Check the full_data introduces
    full_data <- check_full_data(full_data)
    # Check mapper arguments
    check_return <- check_arg_mapper(full_data, filter_values, distance_type, clustering_type,
                                     linkage_type)

    full_data <- check_return[[1]]
    filter_values <- check_return[[2]]
    optimal_clustering_mode <- check_return[[3]]
  }

  mapper_object_ini <- list("full_data" = full_data,
                            "filter_values" = filter_values,
                            "num_intervals" = num_intervals,
                            "percent_overlap" = percent_overlap/100,
                            "distance_type" = distance_type,
                            "optimal_clustering_mode" = optimal_clustering_mode,
                            "num_bins_when_clustering" = num_bins_when_clustering,
                            "clustering_type" = clustering_type,
                            "linkage_type" = linkage_type,
                            "optimal_clustering_mode" = optimal_clustering_mode)

  class(mapper_object_ini) <- "mapper_initialization"

  mapper_object <- one_D_Mapper(mapper_object_ini)

  return(mapper_object)
}


