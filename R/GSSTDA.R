#' @title GSSTDA
#'
#' @description Gene Structure Survival using Topological Data Analysis
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
#' @param survival_time Time between disease diagnosis and death (if not dead until the end of follow-up).
#' @param survival_event \code{logical}. Whether the patient has died or not.
#' @param case_tag The tag of the healthy patient (healthy or not).
#' @param num_intervals Number of intervals used to create the first sample
#' partition based on filtering values.
#' @param percent_overlap Percentage of overlap between intervals. Expressed
#' as a fraction from zero to one.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("cor") and euclidean ("euclidean"). "cor"
#' default option.
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' "hierarchical" default option.
#' @param num_bins_when_clustering Number of bins to generate the
#' histogram employed by the standard optimal number of cluster finder
#' method. Parameter not necessary if the "optimal_clust_mode" option
#' is "silhouette" or the "clust_type" is "PAM".
#' @param linkage_type Linkage criteria used in hierarchical clustering.
#' Choose between "single" for single-linkage clustering, "complete" for
#' complete-linkage clustering or "average" for average linkage clustering
#' (or UPGMA). Only necessary for hierarchical clustering.
#' "single" default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#' @return A \code{GSSTDA} object.
#' @export
#' @examples
#' \dontrun{
#' num_rows <- 100
#' full_data <- data.frame( x=2*cos(1:num_rows), y=sin(1:num_rows) )
#' filter_values <- list(2*cos(1:num_rows))
#' GSSTDA <- GSSTDA(full_data,  survival_time, survival_event, case_tag, num_intervals = 4,
#'                      percent_overlap = 0.5, distance_type = "euclidean",
#'                      num_bins_when_clustering = 8,
#'                      clustering_type = "hierarchical",
#'                      linkage_type = "single")}
GSSTDA <- function(full_data, survival_time, survival_event, case_tag, num_intervals, percent_overlap,
                   distance_type, clustering_type, num_bins_when_clustering, linkage_type, na.rm=TRUE){
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data)
  # Select the control_tag
  control_tag <- check_vectors(ncol(full_data), survival_time, survival_event, case_tag)
  #Don't check filter_values because it is not created.
  filter_values <- ""
  check_return <- check_arg_mapper(full_data, filter_values, distance_type, clustering_type,
                                              linkage_type)
  full_data <- check_return[[1]]
  filter_values <- check_return[[2]]
  optimal_clustering_mode <- check_return[[3]]


  ################### BLOCK I: Pre-process. DGSA (using "NT" control_tag) ##############################
  #   Select the normal tissue data gene expression matrix.
  control_tag_cases <- which(case_tag == control_tag)
  normal_tiss <- full_data[,control_tag_cases]

  #   Obtain the gene expression matrix containing the flattened version of the vectors.
  matrix_flatten_normal_tiss <- flatten_normal_tiss(normal_tiss)
  #   Obtain the normal space
  normal_space <- denoise_rectangular_matrix(matrix_flatten_normal_tiss)
  #   Obtain the disease component of the normal_space
  matrix_disease_component <- generate_disease_component(full_data, normal_space)

  print("The pre-process DGSA is finished")

  ################### BLOCK II: Gene selection (using "T" control_tag) ##################################
  #Remove NAN's values (case_tag == control_tag) of survival_time and survival_event
  survival_time <- survival_time[-control_tag_cases]
  survival_event <- survival_event[-control_tag_cases]
  # Univariate cox proportional hazard models for the expression levels of each gene included in the provided dataset
  cox_all_matrix <- cox_all_genes(matrix_disease_component, survival_time, survival_event)

  print("The gene selection is finished")

  ################### BLOCK III: Create mapper object where the arguments are checked ###################
  # Transpose full_data: rows = patient, columns = genes
  full_data <- t(full_data)

  #   Check filter_values
  #full_data_and_filter_values <- check_filter_values(filter_values)
  #full_data <- full_data_and_filter_values[[1]]
  #filter_values <- full_data_and_filter_values[[2]]
  filter_values <- c(2*cos(1:ncol(full_data)))
  names(filter_values) <- colnames(full_data)

  mapper_obj <- mapper(full_data, filter_values, num_intervals, percent_overlap, distance_type,
                       clustering_type, num_bins_when_clustering, linkage_type, na.rm = "checked")

  print("The mapper process is finished")

  # Create the object
  GSSTDA_object <- list("normal_space" = normal_space,
                        "matrix_disease_component" = matrix_disease_component
                        )

  class(GSSTDA_object) <- "GSSTDA_obj"
  return(GSSTDA_object)
}


