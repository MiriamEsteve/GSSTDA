#' @title GSSTDA
#'
#' @description Gene Structure Survival using Topological Data Analysis
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
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
#' GSSTDA_obj <- GSSTDA(full_data, num_intervals = 4,
#'                      percent_overlap = 0.5, distance_type = "euclidean",
#'                      num_bins_when_clustering = 8,
#'                      clustering_type = "hierarchical",
#'                      linkage_type = "single")}
GSSTDA <- function(full_data, num_intervals, percent_overlap, distance_type, clustering_type, num_bins_when_clustering, linkage_type, na.rm=TRUE){
  #Read the data set
  yes_no <- readline(prompt="Are the columns of the data set the subjects and the rows the features?: yes/no ")
  if(yes_no == "no" | yes_no == "n" | yes_no == ""){
    #Transpose the data set
    full_data <- t(full_data)
  }
  #Convert full_data to matrix type
  full_data <- as.matrix(full_data)

  #Omit NAN's values
  if (na.rm == TRUE){
    full_data <- stats::na.omit(full_data)
    print("Missing values and NaN's are omitted")
  }

  # Control tag
  control_tag <- readline(prompt="What is the tag of control patient?")


  #mapper_object_ini <- list("full_data" = full_data,
  #                             "filter_values" = filter_values,
  #                             "num_intervals" = num_intervals,
  #                             "percent_overlap" = percent_overlap,
  #                             "distance_type" = distance_type,
  #                             "optimal_clustering_mode" = optimal_clustering_mode,
  #                             "num_bins_when_clustering" = num_bins_when_clustering,
  #                             "clustering_type" = clustering_type,
  #                             "linkage_type" = linkage_type)

  #class(mapper_object_ini) <- "mapper_initialization"
  #return(mapper_object)
}
