#' @title check_full_data
#'
#' @description Checking the full_data introduces in the package
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#' @return Return \code{full_data} without NAN's and as a matrix
#' @examples
#' \dontrun{
#' check_full_data(full_data, na.rm = TRUE)}
check_full_data <- function(full_data, na.rm = TRUE){
  #Read the data set
  yes_no <- readline(prompt="Are the columns of the data set the patient and the rows the genes?: yes/no ")
  if(yes_no == "no" | yes_no == "n" | yes_no == ""){
    #Transpose the data set. Columns = patient and rows = genes
    full_data <- t(full_data)
  }
  #Convert full_data to matrix type
  full_data <- as.matrix(full_data)

  #Omit NAN's values
  #if (na.rm == TRUE){
    # Remove rows (genes) with NA's values
  #  full_data <- full_data[rowSums(is.na(full_data))==0,]
  #  print("Missing values and NaN's are omitted")
  #}
  return(full_data)
}


#' @title check_vectors
#' @description Checking the \code{survival_time}, \code{survival_event} and \code{case_tag} introduces in the \code{GSSTDA} object.
#'
#' @param col_full_data Column names of the genes of the full_data (maybe remove by na.rm = TRUE)
#' @param survival_time Time between disease diagnosis and death (if not dead until the end of follow-up).
#' @param survival_event \code{logical}. Whether the patient has died or not.
#' @param case_tag The tag of the healthy patient (healthy or not).
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#' @return control_tag Return the tag of the healthy patient
#' @examples
#' \dontrun{control_tag <- check_vectors(col_full_data, survival_time, survival_event, case_tag)}
check_vectors <- function(col_full_data, survival_time, survival_event, case_tag, na.rm = TRUE){
  ncol_full_data <- length(col_full_data)

  # Check if the arguments are vectors; a valid type of data; and the vectors are the same dimension as a full_data
  if(!is.vector(survival_time) | !is.numeric(survival_time) | length(survival_time) != ncol_full_data){
    stop("survival_time must be a valid values vector and its length must be the same as the number of patients (columns) of the full_data.")
  }
  if(!is.vector(survival_event) | !length(unique(survival_event)) | length(survival_event) != ncol_full_data){
    stop("survival_event must be a valid values vector. Only two type of event. Also, its length must be the same as the number of patients (columns) of the full_data.")
  }
  if(!is.vector(case_tag) | !length(unique(case_tag))){
    stop("case_tag must be a valid values vector. Only two type of tags.")
  }

  #if(na.rm == TRUE){
    # Select the survival_time, survival_event and case_tag with the genes selected in the full_data
  #  survival_time <- survival_time[col_full_data]
  #  survival_event <- survival_event[col_full_data]
  #  case_tag <- case_tag[col_full_data]
  #}

  control_tag_opt <- unique(case_tag)
  control_tag <- readline(prompt=paste("What is the tag of the healthy patient? (", control_tag_opt[1], " or ", control_tag_opt[2], "): " , sep="") )
  if(!(control_tag %in% control_tag_opt)){
    stop("Provide one of the specified case tag")
  }else if(control_tag == ""){
    control_tag = "NT"
  }
  return(control_tag)
}

#' @title check_filter_values
#'
#' @description Checking the filter_values introduces in the \code{mapper} object.
#'
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#'
#' @return \code{filter_value} and \code{full_data} without NAN's
#' @examples
#' \dontrun{check_arg_mapper(filter_values, distance_type, clustering_type, linkage_type)}
check_filter_values <- function(filter_values, na.rm = TRUE){
  # Check if filter_values is a vector
  if(!is.vector(filter_values)){
    stop("filter_values must be a valid values vector")
  }

  #Omit NAN's values
  #if (na.rm == TRUE){
    # Remove colums (subjects) and their filter values with NA's values
  #  filter_values <- filter_values[colnames(full_data)]
    # Remove filter values and respective columns with NA's values
  #  filter_values <- stats::na.omit(filter_values)
  #  full_data <- full_data[,names(filter_values)]
  #  print("Missing values and NaN's are omitted")
  #}
  return(c(full_data, filter_values))
}

#' @title check_arg_mapper
#'
#' @description Checking the arguments introduces in the \code{mapper} object.
#'
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("cor") and euclidean ("euclidean"). "cor"
#' default option.
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' "hierarchical" default option.
#' @param linkage_type Linkage criteria used in hierarchical clustering.
#' Choose between "single" for single-linkage clustering, "complete" for
#' complete-linkage clustering or "average" for average linkage clustering
#' (or UPGMA). Only necessary for hierarchical clustering.
#' "single" default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#'
#' @return \code{optimal_clustering_mode}
#' @examples
#' \dontrun{check_arg_mapper(filter_values, distance_type, clustering_type, linkage_type)}
check_arg_mapper <- function(full_data, filter_values, distance_type, clustering_type, linkage_type, na.rm = TRUE){
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

  # Check if filter_values == "" the filter_values is not calculated yet. So, we checked only the others args
  if(filter_values != ""){
    full_data_and_filter_values <- check_filter_values(filter_values)
    full_data <- full_data_and_filter_values[[1]]
    filter_values <- full_data_and_filter_values[[2]]
  }

  return(c(full_data, filter_values, optimal_clustering_mode))
}
