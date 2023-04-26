#' @title GSSTDA_obj
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
#' GSSTDA_obj <- GSSTDA_obj(full_data,  survival_time, survival_event, case_tag, num_intervals = 4,
#'                      percent_overlap = 0.5, distance_type = "euclidean",
#'                      num_bins_when_clustering = 8,
#'                      clustering_type = "hierarchical",
#'                      linkage_type = "single")}
GSSTDA_obj <- function(full_data, survival_time, survival_event, case_tag, num_intervals, percent_overlap, distance_type, clustering_type, num_bins_when_clustering, linkage_type, na.rm=TRUE){
  #Check the arguments introduces in the function
  check_full_data(full_data)
  check_vectors()
  optimal_clustering_mode <- check_arg_mapper(full_data, filter_values, distance_type, clustering_type, linkage_type)

  # Pre-process. DGSA


  # Create mapper object where the arguments are checked
  mapper_obj <- mapper(full_data, filter_values, num_intervals, percent_overlap, distance_type, clustering_type, num_bins_when_clustering, linkage_type, na.rm = "checked")


  # Create the object
  GSSTDA_object_ini <- list( unlist(mapper_obj),
                             "survival_time" = survival_time,
                             "survival_event" = survival_event,
                             "case_tag" = case_tag
                             )

  class(GSSTDA_object_ini) <- "GSSTDA_initialization"
  return(GSSTDA_object_ini)
}


#' @title GSSTDA
#'
#' @description Gene Structure Survival using Topological Data Analysis
#' @param GSSTDA_obj G-SS-TDA object return by \code{GSSTDA_obj} function.
#' @return A \code{GSSTDA} output.
#' @export
#' @examples
#' \dontrun{
#' GSSTDA_obj <- GSSTDA_obj(full_data, num_intervals = 4,
#'                      percent_overlap = 0.5, distance_type = "euclidean",
#'                      num_bins_when_clustering = 8,
#'                      clustering_type = "hierarchical",
#'                      linkage_type = "single")}
#' GSSTDA <- GSSTDA(GSSTDA_obj)
GSSTDA <- function(GSSTDA_obj){



}
