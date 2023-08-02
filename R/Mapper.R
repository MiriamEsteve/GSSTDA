#' @title one_D_Mapper
#'
#' @description Wrapping function to carry out the complete process.
#'
#' @param mapper_object_ini Mapper TDA initializated object generated
#' by \code{mapper} function.
#'
#' @export
#' @return A \code{mapper_obj} object. It contains the values of the intervals
#' (interval_data), the samples included in each interval (sample_in_level),
#' information about the cluster to which the individuals in each interval
#' belong (clustering_all_levels), a list including the individuals contained
#' in each detected node (node_samples), their size (node_sizes), the
#' average of the filter function values of the individuals of each node
#' (node_average_filt) and the adjacency matrix linking the nodes (adj_matrix).
#' Moreover, information is provided on the number of nodes, the average node
#' size, the standard deviation of the node size, the number of connections
#' between nodes, the proportion of connections to all possible connections
#' and the number of ramifications.
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

  node_sizes = unlist(lapply(node_samples,length))
  # average of the filter function of each node
  node_average_filt = lapply(node_samples,function(x,y) mean(y[x]),filter_values)

  # additional parameters
  n_nodes <- length(node_sizes)
  av_node_size <- mean(node_sizes)
  sd_node_size <- stats::sd(node_sizes)

  adj_mat <- adj_matrix_out
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
  mapper_object <- list("interval_data" = interval_data,
                        "sample_in_level" = samp_in_lev,
                        "clustering_all_levels" = test_clust_all_levels,
                        "node_samples" = node_samples,
                        "node_sizes" = node_sizes,
                        "node_average_filt" = node_average_filt,
                        "adj_matrix" = adj_matrix_out,
                        "n_sizes" = n_nodes,
                        "average_nodes"= av_node_size,
                        "standard_desviation_nodes" = sd_node_size,
                        "number_connections" = n_connections,
                        "proportion_connections" = prop_connections,
                        "number_ramifications" = n_ramifications)

  class(mapper_object) <- "mapper_object"
  return(mapper_object)
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
#' \donttest{
#' # Create data object
#' data_object <- list("full_data" = full_data, "survival_time" = survival_time,
#'                    "survival_event" = survival_event, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#'
#' #Select gene from data object
#' geneSelection_object <- geneSelection(data_object, gen_select_type="top_bot",
#'  percent_gen_select=10)
#'
#' mapper_object <- mapper(full_data = geneSelection_object[["genes_disease_component"]],
#' filter_values = geneSelection_object[["filter_values"]],
#' num_intervals = 5,
#' percent_overlap = 40, distance_type = "cor",
#' clustering_type = "hierarchical",
#' linkage_type = "single")
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

