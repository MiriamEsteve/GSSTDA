#' @title check_full_data
#'
#' @description Checking the full_data introduces in the package
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#' @return Return \code{full_data} without NAN's and as a matrix
check_full_data <- function(full_data, na.rm = TRUE){
  # If this function has been executed don't do nothing
  if(na.rm == "checked"){
    return(full_data)
  }

  #Read the data set
  yes_no <- readline(prompt="Are the columns of the data set the patient and the rows the genes?: yes/no ")

  if(yes_no == "no" | yes_no == "n"){
    #Transpose the data set. Columns = patient and rows = genes
    full_data <- t(full_data)
  }
  #Convert full_data to matrix type
  full_data <- as.matrix(full_data)

  #Omit NAN's values
  if (na.rm == TRUE){
    nrow_ini = nrow(full_data)
    # Remove rows (genes) with NA's values
    full_data <- full_data[rowSums(is.na(full_data))==0,]

    message(paste(nrow(full_data) - nrow_ini, " missing values and NaN's are omitted in the genes (rows)"))
  }
  return(full_data)
}


#' @title check_vectors
#' @description Checking the \code{survival_time}, \code{survival_event} and \code{case_tag} introduces in the \code{GSSTDA} object.
#'
#' @param full_data The genes of the full_data (maybe remove by na.rm = TRUE)
#' @param survival_time Time between disease diagnosis and event (if there was
#' no event until the end of follow-up).
#' @param survival_event \code{logical}. Whether or not the event has occurred.
#' @param case_tag The tag of the healthy sample (healthy or not).
#' @param control_tag Tag of the healthy sample.E.g. "T"
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#' @return control_tag Return the tag of the healthy sample.
check_vectors <- function(full_data, survival_time, survival_event, case_tag, control_tag, na.rm = TRUE){
  ncol_full_data <- ncol(full_data)
  # Check if the arguments are vectors; a valid type of data; and the vectors are the same dimension as a full_data
  if(!is.vector(survival_time) | !is.numeric(survival_time) | length(survival_time) != ncol_full_data){
    stop("survival_time must be a valid values vector and its length must be the same as the number of patients (columns) of the full_data.")
  }

  # Omit NAN's values in checking
  if(!is.vector(survival_event) | !(length(unique(stats::na.omit(survival_event))) == 2 & is.numeric(stats::na.omit(survival_event))) | length(survival_event) != ncol_full_data){
    stop("survival_event must be a valid values vector. Only two type of event (0 or 1). Also, its length must be the same as the number of patients (columns) of the full_data.")
  }

  #If exits NAN's values remove it and check if it contain only two cases and it has the same dimension (columns) as full_data
  if(!is.vector(case_tag)){
    stop("case_tag must be a valid values vector.")
  }
  if( any(is.na(case_tag))){
    without_nan_patient <- which(!is.na(case_tag))
    case_tag <- case_tag[without_nan_patient]
    full_data <- full_data[,without_nan_patient]
    survival_event <- survival_event[without_nan_patient]
    survival_time <- survival_time[without_nan_patient]
    ncol_full_data <- ncol(full_data)
    message("NAN's values in patient was removed in case_tag, full_data, survival_time and survival_event")
  }
  if(length(unique(case_tag)) != 2){
    stop("case_tag must has only two type of tags.")
  }
  if(length(case_tag) != ncol_full_data){
    stop("The length of case_tag must be the same as the number of patients (columns) of the full_data.")
  }

  control_tag_opt <- unique(case_tag)
  if(is.na(control_tag)){
    control_tag <- readline(prompt=paste("What is the tag of the healthy patient (value in the case_tag)? (", control_tag_opt[1], " or ", control_tag_opt[2], "): " , sep="") )
  }

  if(!(control_tag %in% control_tag_opt)){
    print(paste("The case tag is '", control_tag_opt[1], "' by default"))
    control_tag <- control_tag_opt[1]
  }
  return(list(control_tag, full_data, survival_event, survival_time, case_tag))
}

#' @title check_filter_values
#'
#' @description Checking the filter_values introduces in the \code{mapper} object.
#'
#' @param full_data Matrix with the columns of the input matrix
#' corresponding to the individuals belonging to the level. This matrix could be the
#' genes_disease_component.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#'
#' @return \code{filter_value} and \code{full_data} without NAN's
check_filter_values <- function(full_data, filter_values, na.rm = TRUE){
  # Check if filter_values is a vector
  if(!is.vector(filter_values)){
    stop("filter_values must be a valid values vector")
  }

  #Check if the names of the filter_values are the same as the cols of full_data.
  if(!setequal(names(filter_values), colnames(full_data))){
    stop("The name of the filter_values must be the same as the patient name of the full_data (or genes_disease_component).")
  }

  #Omit NAN's values
  if (na.rm == TRUE){
    # Remove rows (subjects) and their filter values with NA's values
    filter_values <- filter_values[colnames(full_data)]
    # Remove filter values and respective rows with NA's values
    filter_values <- stats::na.omit(filter_values)
    full_data <- full_data[,names(filter_values)]
  }
  return(list(full_data, filter_values))
}

#' @title check_gene_selection
#'
#' @description Checking the arguments introduces in the gene selection process.
#'
#' @param num_genes Number of genes in the full_data
#' @param gen_select_type Type of gene selection to be used. Choose between "top_bot" (top-botton)
#' and "abs" (absolute)
#' @param percent_gen_select Percentage of genes to be selected
#'
#' @return num_gen_select Number of genes to be selected according to the percent_gen_select value
check_gene_selection <- function(num_genes, gen_select_type, percent_gen_select){
  #Convert text to lowercase
  gen_select_type <- tolower(gen_select_type)
  #Check gen_select_type
  gen <- c("top_bot","abs")
  if(!gen_select_type %in% gen){
    stop(paste("Invalid gene selection type selected. Choose one of the folowing: ", paste(gen, collapse = ", ")))
  }

  #Number of genes to be selected in gene_selection_surv function
  num_gen_select <- trunc((percent_gen_select/100) * num_genes)

  return(num_gen_select)
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
#' Choose between correlation ("correlation") and euclidean ("euclidean"). "correlation"
#' default option.
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' "hierarchical" default option.
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
#' If \code{FALSE}, an error occurs in case of \code{NA} rows.
#'
#' @return \code{optimal_clustering_mode}
check_arg_mapper <- function(full_data, filter_values, distance_type, clustering_type, linkage_type, optimal_clustering_mode, na.rm = TRUE){
  #Check distance_type
  distances <- c("correlation","euclidean")
  if(!distance_type %in% distances){
    stop(paste("Invalid distance selected. Choose one of the folowing: ", paste(distances, collapse = ", ")))
  }

  #Check clustering_type
  clust_types <- c("hierarchical","PAM")
  if(!clustering_type %in% clust_types){
    stop(paste("Invalid clustering method selected. Choose one of the folowing: ", paste(clust_types,collapse = ", ")))
  }

  if(is.na(optimal_clustering_mode)){
    optimal_clustering_mode <- "silhouette"

    if(clustering_type == "hierarchical"){
      option <- readline(prompt="Choose one of the following optimal cluster number method: standard/silhouette: ")

      if(option != "standard"){
        optimal_clustering_mode <- "silhouette"
      }
      else{
        optimal_clustering_mode <- "standard"
      }
    }
  }else{
    #Check optimal_clustering_mode
    optimal_clustering <- c("silhouette","standard")
    if(!optimal_clustering_mode %in% optimal_clustering){
      stop(paste("Invalid optimal_clustering selected. Choose one of the folowing: ", paste(optimal_clustering, collapse = ", ")))
    }
  }
  message("The optimal clustering mode is '", optimal_clustering_mode, " '")

  #Check linkage_type
  link_types <- c("single","average","complete")
  if(!linkage_type %in% link_types){
    stop(paste("Invalid linkage method selected. Choose one of the folowing: ", paste(link_types,collapse = ", ")))
  }

  # Check if filter_values == [] the filter_values is not calculated yet. So, we checked only the others args
  if(length(filter_values) != 0 & na.rm != "checked"){
    full_data_and_filter_values <- check_filter_values(full_data, filter_values)
    full_data <- full_data_and_filter_values[[1]]
    filter_values <- full_data_and_filter_values[[2]]
  }

  return(list(full_data, filter_values, optimal_clustering_mode))
}
