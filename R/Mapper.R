#' @title Mapper
#'
#' @description TDA are persistent homology and mapper. Persistent homology
#' borrows ideas from abstract algebra to identify particular aspects
#' related to the shape of the data such as the number of connected
#' components and the presence of higher-dimensional holes, whereas
#' mapper condenses the information of high-dimensional datasets into
#' a combinatory graph or simplicial complex that is referred to as
#' the skeleton of the dataset. This implementation is the mapper of one
#' dimension, i.e. using only one filter function value.
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
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
mapper <- function(full_data, filter_values, num_intervals, percent_overlap, distance_type, clustering_type, num_bins_when_clustering, linkage_type, na.rm=TRUE){
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

  #Check if filter_values is a vector
  if(!is.vector(filter_values)){
    stop("filter_values must be a valid values vector")
  }

  #Check if there are more of two rows in the full_data

  #Check if the names of the filter_values are the same as the columns of full_data.
  if(!setequal(names(filter_values), colnames(full_data))){
    stop("The name of the filter_values must be the same as the subject name of the full_data.")
  }

  #Check distance_type
  distances <- c("cor","euclidean")
  if(!distance_type %in% distances){
    stop(paste("Invalid distance selected. Choose one of the folowing:", paste(distances, collapse = ", ")))
  }

  #Check clustering_type
  clust_types <- c("hierarchical","PAM")
  if(!clustering_type %in% clust_types){
    stop(paste("Invalid clustering method selected. Choose one of the folowing:", paste(clust_types,collapse = ", ")))
  }

  optimal_clustering_mode <- "silhouette"

  if(clustering_type == "hierarchical"){
    opt_clust_modes <- c("standard","silhouette")
    option <- readline(prompt="Choose one of the folowing optimal cluster number method: standard/silhouette ")
    if(option == "standard" | option == ""){
      optimal_clustering_mode <- "standard"
    }
  }

  #Check linkage_type
  link_types <- c("single","average","complete")
  if(!linkage_type %in% link_types){
    stop(paste("Invalid linkage method selected. Choose one of the folowing:", paste(link_types,collapse = ", ")))
  }


  mapper_object_ini <- list("full_data" = full_data,
                               "filter_values" = filter_values,
                               "num_intervals" = num_intervals,
                               "percent_overlap" = percent_overlap,
                               "distance_type" = distance_type,
                               "optimal_clustering_mode" = optimal_clustering_mode,
                               "num_bins_when_clustering" = num_bins_when_clustering,
                               "clustering_type" = clustering_type,
                               "linkage_type" = linkage_type)

  class(mapper_object_ini) <- "mapper_initialization"

  #mapper_object <- one_D_Mapper(mapper_object_ini)

  return(mapper_object)
}


#' @title one_D_Mapper
#'
#' @description Wrapping function to carry out the complete process.
#'
#' @param mapper_object_ini Mapper TDA initializated object generated
#' by \code{mapper} function.
#' @return A \code{mapper_obj} object. It contains the values of the intervals
#' (interval_data), the samples included in each interval (sample_in_level),
#' information about the cluster to which the individuals in each interval
#' belong (clustering_all_levels), a list including the individuals contained
#' in each detected node (node_samples), their size (node_sizes), the
#' average of the filter function values of the individuals of each node
#' (node_average_filt) and the adjacency matrix linking the nodes (adj_matrix).
#' @export
#'
#' @examples
#' \dontrun{
#' one_D_Mapper(mapper_object_ini)}
one_D_Mapper <- function(mapper_object_ini){

  full_data <- mapper_object_ini[["full_data"]]
  filter_values <- mapper_object_ini[["filter_values"]]

  #Getting intervals.
  interval_data <- get_intervals_One_D(filter_values, mapper_object_ini[["num_intervals"]], mapper_object_ini[["percent_overlap"]])

  #Getting samples on each interval.
  samp_in_lev <- samples_in_levels(interval_data, filter_values)

  #Clustering all levels.
  test_clust_all_levels <- clust_all_levels(full_data, samp_in_lev, mapper_object_ini[["distance_type"]], mapper_object_ini[["clustering_type"]],
                                            mapper_object_ini[["linkage_type"]], mapper_object_ini[["optimal_clustering_mode"]],  mapper_object_ini[["num_bins_when_clustering"]])
  #Transforming levels into nodes.
  node_samples <- levels_to_nodes(test_clust_all_levels)

  #Computing adjacency matrix.
  adj_matrix_out <- compute_node_adjacency(node_samples)

  #Generating the object of the output data
  mapper_object <- list("interval_data" = interval_data,
                           "sample_in_level" = samp_in_lev,
                           "clustering_all_levels" = test_clust_all_levels,
                           "node_samples" = node_samples,
                           "node_sizes" = unlist(lapply(node_samples,length)),
                           "node_average_filt" = lapply(node_samples,function(x,y) mean(y[x]),filter_values),
                           "adj_matrix" = adj_matrix_out)

  class(mapper_object) <- "mapper_object"
  return(mapper_object)
}


#' @title Get information from results
#' @description This functions retrieves additional information about the
#' mapper results.
#' @param mapper_object Result object from a \code{mapper} function.
#' @return A mapper_information object (a list) with informative parameters
#' of the mapper result, in this order: the number of nodes, the average node
#' size, the standard deviation of the node size, the number of connections
#' between nodes, the proportion of connections to all possible connections
#' and the number of ramifications.
#' @export
#' @examples
#' \dontrun{
#' get_information_from_results(mapper_object)
#' }
get_information_from_results <- function(mapper_object){
  n_nodes <- length(mapper_object[["node_sizes"]])
  av_node_size <- mean(mapper_object[["node_sizes"]])
  sd_node_size <- stats::sd(mapper_object[["node_sizes"]])

  adj_mat <- mapper_object[["adj_matrix"]]
  upper_tri <- upper.tri(adj_mat)
  lower_tri <- lower.tri(adj_mat)

  n_connections <- sum(adj_mat[upper_tri] == 1)
  prop_connections <- n_connections/length(adj_mat[upper_tri])

  adj_mat[lower_tri] <- t(adj_mat)[lower_tri]
  diag(adj_mat) <- 0
  n_ramifications <- colSums(adj_mat)-2
  n_ramifications[n_ramifications %in% c(-1,-2)] <- 0
  n_ramifications <- sum(n_ramifications)

  #Generating the object of the output data
  mapper_information <- list("node_sizes" = n_nodes,
                                "average_nodes"= av_node_size,
                                "standard_desviation_nodes " = sd_node_size,
                                "number_connections" = n_connections,
                                "proportion_connections" = prop_connections,
                                "number_ramifications" = n_ramifications)

  class(mapper_information) <- "mapper_information"

  return(mapper_information)
}


#' @title Plot mapper
#' @description This function produces an interactive network plot using
#' the \code{visNetork} function from the mapper results.
#' @param mapper_object A list produced as an output of the \code{one_D_Mapper}
#' function.
#' @param trans_node_size Logical, it indicates whether you want to resize
#' the size of the nodes. \code{TRUE} default option.
#' @param exp_to_res Only necessary if trans_node_size is \code{TRUE}. An
#' exponent of the form 1/n to which the node sizes must be raised in order
#' to resize them.
#' @return Plots an interactive network using the \code{visNetwork} function.
#' @export
#' @import visNetwork
#' @examples
#' \dontrun{
#' plot_mapper(mapper_object)}
plot_mapper <- function(mapper_object,trans_node_size = TRUE,exp_to_res = 1/2){
  arr_ind <- base::which(arr.ind = TRUE,mapper_object[["adj_matrix"]] == 1)
  df_out <- base::data.frame(base::rownames(mapper_object[["adj_matrix"]])[arr_ind[,1]],base::colnames(mapper_object[["adj_matrix"]])[arr_ind[,2]])
  df_out <- base::cbind(arr_ind,df_out)
  base::rownames(df_out) <- 1:base::nrow(df_out)
  base::colnames(df_out) <- c("from","to","from_name","to_name")
  nodes_to_net <- base::unique(base::data.frame(c(df_out[,1]-1,df_out[,2]-1),c(df_out[,3],df_out[,4])))
  nodes_to_net$node_size <- mapper_object[["node_sizes"]]
  if(trans_node_size){
    nodes_to_net$node_size <- (nodes_to_net$node_size)^exp_to_res
  }
  base::colnames(nodes_to_net) <- c("id","label","size")
  nodes_to_net$color <- map_to_color(base::log2(base::unlist(mapper_object[["node_average_filt"]]) + 2))
  edges_to_net <- df_out[,c(1,2)]-1
  base::colnames(edges_to_net) <- c("from","to")
  visNetwork::visNetwork(nodes_to_net,edges_to_net[!edges_to_net$from == edges_to_net$to,],)
}

