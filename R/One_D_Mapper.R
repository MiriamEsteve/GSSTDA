#' @title Extract intervals from filter function output values.
#' @description
#' It calculates the intervals of the given output values of a filter function
#' according to the given number and percentage of overlap.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param num_intervals Number of intervals to divide the filtering function
#' values in.
#' @param percent_overlap Percentage of overlap between intervals.
#' @return Returns a list with the set of intervals for the filtering function
#' values.
#' @examples
#' \dontrun{
#' get_intervals_One_D(filter_values,num_intervals,percent_overlap)}
get_intervals_One_D <- function(filter_values,num_intervals,percent_overlap){
  range_filt <- max(filter_values) - min(filter_values)
  n_ov <- num_intervals -1
  l_int <- range_filt/(num_intervals - (n_ov*percent_overlap))
  p_int <- percent_overlap*l_int
  list_int <- list()
  list_int[[1]] <- c(min(filter_values)-0.1,min(filter_values + l_int))
  names(list_int)[1] <- "Level_1"
  for(i in 2:num_intervals){
    if(i < num_intervals){
      list_int[[i]] <- c(list_int[[i-1]][2]-p_int,list_int[[i-1]][2]-p_int + l_int)
      names(list_int)[i] <- paste("Level_",i,sep ="")
    }
    else if(i == num_intervals){
      list_int[[i]] <- c(list_int[[i-1]][2]-p_int,list_int[[i-1]][2]-p_int + l_int + 0.1)
      names(list_int)[i] <- paste("Level_",i,sep ="")
    }
  }
  return(list_int)
}


#' @title Samples in levels
#' @description This function returns a list of vectors containing the
#' individuals included at each level, i.e. the vectors of individuals with
#' a value of the filter function within each of the intervals.
#' @param interval_data Filter function intervals. List with the set of
#' intervals for the filtering function values produced by the
#' \code{get_intervals_One_D} function.
#' @param filter_values Vector obtained after applying the filtering
#' function to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @return A list of character vectors with the samples included
#' in each of the levels (i.e. each of the intervals of the values of the
#' filter functions).
#' @examples
#' \dontrun{
#' samples_in_levels(interval_data,filter_values)}
samples_in_levels <- function(interval_data,filter_values){
  return(lapply(interval_data,function(x,y) names(which(y >= x[[1]] & y < x[[2]])),filter_values))
}


#' @title Get clusters for a particular data level
#' @description
#' It performs clustering of the samples belonging to a particular level (to a
#' particular interval of the filter function) with the proposed clustering
#' algorithm and the proposed method to determine the optimal number
#' of clusters.
#' @param full_data_i Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("cor") and euclidean ("euclidean").
#' @param clustering_type Type of clustering method.
#' Choose between "hierarchical" and "PAM" (“partition around medoids”)
#' options.
#' @param linkage_type Linkage criteria used in hierarchical clustering.
#' Choose between "single" for single-linkage clustering, "complete" for
#' complete-linkage clustering or "average" for average linkage clustering
#' (or UPGMA). Only necessary for hierarchical clustering. The value provided
#' if the type of clustering chosen is hierarchical will be ignored
#' @param optimal_clustering_mode Method for selection optimal number of
#' clusters. It is only necessary if the chosen type of algorithm is
#' hierarchical. In this case, choose between "standard" (the method used
#' in the original mapper article) or "silhouette". In the case of the PAM
#' algorithm, the method will always be "silhouette".
#' @param num_bins_when_clustering Number of bins to generate the histogram
#' employed by the standard optimal number of cluster finder method.
#' Parameter not necessary if the "optimal_clust_mode" option is "silhouette"
#' or the "clust_type" is "PAM".
#' @param level_name Name of the studied level. # ERROR No usado
#' @return Returns a interger vector with the samples included in each cluster
#' for the specific level analyzed. The names of the vector values are the
#' names of the samples and the vector values are the node number
#' to which the individual belongs.
#' @import cluster
#' @examples
#' \dontrun{
#' clust_lev(full_data_i, distance_type, clustering_type, linkage_type,
#'           optimal_clustering_mode, num_bins_when_clustering,level_name)}
clust_lev <- function(full_data_i, distance_type, clustering_type, linkage_type,
                      optimal_clustering_mode, num_bins_when_clustering, level_name){
  #Distance type
  if(distance_type == "cor"){
    level_dist <- stats::as.dist(1-stats::cor(full_data_i))
  }else{
    level_dist <- stats::dist(base::t(full_data_i),method = distance_type)
  }

  #Clustering type
  max_dist_lev <- base::max(level_dist)

  if(clustering_type == "PAM"){
    av_sil <- c()
    n_clust <- c()
    for(i in 1:(ncol(full_data_i)-1)){
      temp_clust <- cluster::pam(x =level_dist,diss = TRUE,k = i)
      if(i == 1){
        av_sil <- c(av_sil,0)
        n_clust <- c(n_clust,1)
      }else{
        av_sil <-c(av_sil,mean(cluster::silhouette(temp_clust$clustering,level_dist)[,3]))
        n_clust <- c(n_clust,i)
      }
    }
    if(max(av_sil) >= 0.25){
      op_clust <- n_clust[which.max(av_sil)]
      cluster_indices_level <- cluster::pam(x =level_dist,diss = TRUE,k = op_clust)$clustering
      return(cluster_indices_level)
    }else{
      cluster_indices_level <- rep(1,ncol(full_data_i))
      names(cluster_indices_level) <- colnames(full_data_i)
      return(cluster_indices_level)
    }
  } else if(clustering_type == "hierarchical"){
    level_hclust_out <- stats::hclust(level_dist,method = linkage_type)

    #Optimal clustering mode
    if(optimal_clustering_mode == "standard"){
      heights <- level_hclust_out$height
      breaks_for_bins <- base::seq(from=min(heights), to=max_dist_lev, by=(max_dist_lev - base::min(heights))/num_bins_when_clustering)
      histogram <- graphics::hist(c(heights,max_dist_lev), breaks=breaks_for_bins, plot=FALSE)
      hist_gap <- (histogram$counts == 0)
      if(all(!hist_gap)){
        print("There is no gap... therefore only one cluster...")
        cluster_indices_level <- base::rep(1,ncol(full_data_i))
        names(cluster_indices_level) <- base::colnames(full_data_i)
        return(cluster_indices_level)
      }else{
        threshold_value <- histogram$mids[min(which(hist_gap == TRUE))]
        cluster_indices_level <- base::as.vector(stats::cutree(level_hclust_out, h=threshold_value))
        base::names(cluster_indices_level) <- base::colnames(full_data_i)
        return(cluster_indices_level)
      }
    }
    else if(optimal_clustering_mode == "silhouette"){
      max_dist_lev <- base::max(level_dist)
      level_hclust_out <- stats::hclust(level_dist,method=linkage_type)
      n_clust <- c()
      av_sil <- c()
      for(i in 2:(length(level_hclust_out$order)-1)){
        n_clust <- c(n_clust,i)
        test <- cluster::silhouette(stats::cutree(level_hclust_out,i),level_dist)
        av_sil <- c(av_sil,mean(test[,3]))
      }
      if(max(av_sil) >= 0.25){
        op_clust <- n_clust[which.max(av_sil)]
        cluster_indices_level <- stats::cutree(level_hclust_out,op_clust)
        return(cluster_indices_level)
      }else{
        cluster_indices_level <- stats::cutree(level_hclust_out,1)
        return(cluster_indices_level)
      }
    }
  }
}


#' @title Get clusters for all data level
#' @description It performs the clustering of the samples in each
#' of the levels. That is to say, in each interval of values of the
#' filtering function, the samples with a value within that interval
#' are clustered using the proposed clustering algorithm and the
#' proposed method to determine the optimal number of clusters.
#' @param full_data Input data matrix whose columns are the individuals
#' and rows are the features.BR cambiar nombre.
#' @param samp_in_lev A list of character vectors with the individuals
#' included in each of the levels (i.e. each of the intervals of the values
#' of the filter functions). It is the output of the \code{samples_in_levels}
#' function.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("cor") and euclidean ("euclidean").
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' @param linkage_type Linkage criteria used in hierarchical clustering.
#' Choose between "single" for single-linkage clustering, "complete" for
#' complete-linkage clustering or "average" for average linkage clustering
#' (or UPGMA). Only necessary for hierarchical clustering.
#' @param optimal_clustering_mode Method for selection optimal number of
#' clusters. It is only necessary if the chosen type of algorithm is
#' hierarchical. In this case, choose between "standard" (the method used
#' in the original mapper article) or "silhouette". In the case of the
#' PAM algorithm, the method will always be "silhouette". "silhouette"
#' @param num_bins_when_clustering Number of bins to generate the
#' histogram employed by the standard optimal number of cluster finder
#' method. Parameter not necessary if the "optimal_clust_mode" option
#' is "silhouette" or the "clust_type" is "PAM".
#' @return List of interger vectors. Each of the vectors contains information
#' about the nodes at each level and the individuals contained in them. The
#' names of the vector values are the names of the samples and the vector
#' values are the node number of that level to which the individual belongs.
#' @examples
#' \dontrun{
#' clust_all_levels(full_data, samp_in_lev, distance_type, clustering_type,
#'                  linkage_type, optimal_clustering_mode, num_bins_when_clustering)}
clust_all_levels <- function(full_data, samp_in_lev, distance_type, clustering_type,
                             linkage_type, optimal_clustering_mode, num_bins_when_clustering){

  list_out <- base::list()
  for(i in 1:base::length(samp_in_lev)){
    if(length(samp_in_lev[[i]]) > 2){
      clust_level_temp <- clust_lev(full_data[,samp_in_lev[[i]]], distance_type, clustering_type,
                                    linkage_type, optimal_clustering_mode, num_bins_when_clustering,
                                    base::paste("Level",i,sep="_"))
    }else{
      if(base::length(samp_in_lev[[i]]) < 3 & base::length(samp_in_lev[[i]]) > 0){
        clust_level_temp <- base::rep(1,length(samp_in_lev[[i]]))
        base::names(clust_level_temp) <- samp_in_lev[[i]]
      }else if(length(samp_in_lev[[i]]) == 0){
        clust_level_temp <- NA
      }
    }
    list_out[[i]] <- clust_level_temp
  }
  base::names(list_out) <- base::names(samp_in_lev)
  return(list_out)
}


#' @title Extract Information about Nodes
#' @description Extract the nodes information based on information about
#' clustering. The individuals who are part of each node are identified
#' @param clust_all_levels_list A list with information on the levels
#' obtained from the \code{clust_all_levels} function.
#' @return A list including the individuals content of each detected node.
#' List of character vectors. Each of the vectors contains the names
#' of the individuals at each node.
#' @examples
#' \dontrun{
#' levels_to_nodes(clust_all_levels_list)}
levels_to_nodes <- function(clust_all_levels_list){
  nodes_list <- list()
  node_counter <- 1
  for(i in 1:base::length(clust_all_levels_list)){
    if(!base::all(base::is.na(clust_all_levels_list[[i]]))){
      clusters <- base::unique(clust_all_levels_list[[i]])
      for(j in 1:base::length(clusters)){
        nodes_list[[base::paste("Node",node_counter,sep="_")]] <- base::names(clust_all_levels_list[[i]][clust_all_levels_list[[i]] == clusters[j]])
        node_counter <- node_counter + 1
      }
    }
  }
  return(nodes_list)
}


#' @title Computes the adjacency matrix.
#' @description It computes the adjacency matrix between nodes. Two nodes
#' are considered connected if they share at least one individual.
#' @param nodes_list Output of the \code{levels_to_nodes} function. List of
#' character vectors. Each of the vectors contains the names of the
#' individuals at each node.
#' @return It returns a matrix of magnitude n nodes x n nodes that stores a
#' 1 if there are shared samples in two given nodes and a 0 otherwise.
#' @examples
#' \dontrun{
#' compute_node_adjacency(nodes_list)}
compute_node_adjacency <- function(nodes_list){
  adj_matrix <- base::matrix(0,nrow = base::length(nodes_list),ncol = base::length(nodes_list))
  for(i in 1:(base::length(nodes_list))){
    for(j in i:(base::length(nodes_list))){
      if(length(base::intersect(nodes_list[[i]],nodes_list[[j]])) > 0){
        adj_matrix[i,j] <- 1
      }
    }
  }
  base::colnames(adj_matrix) <- base::names(nodes_list)
  base::rownames(adj_matrix) <- base::names(nodes_list)
  return(adj_matrix)
}


#' @title Map to color
#' @description
#' Auxiliary function that maps a numeric vector, the average node
#' filtering function values, to a color vector.
#' @param x A vector of numeric values storing the average filtering
#' function values found in the samples placed into a specific node.
#' @param limits A two element numeric vector including the range of values.
#' This is optional.
#' @return A vector of the same length of x with colors ranging from blue to
#' red.
#' @import grDevices
#' @examples
#' \dontrun{
#' map_to_color(base::log2(2))}
map_to_color <- function(x,limits=NULL){
  pallette_ob <-  grDevices::colorRampPalette(colors = c("blue","red"))(100)
  if(is.null(limits)){
    limits=range(x)}
  map_to_col <- pallette_ob[base::findInterval(x,base::seq(limits[1],limits[2],length.out=length(pallette_ob)+1), all.inside=TRUE)]
  return(map_to_col)
}

