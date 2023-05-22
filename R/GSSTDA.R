#' @title Gene Structure Survival using Topological Data Analysis (GSSTDA).
#'
#' @description Gene Structure Survival using Topological Data Analysis.
#' This function implements an analysis for expression array data
#' based on the *Progression Analysis of Disease* developed by Nicolau
#' *et al.* (doi: 10.1073/pnas.1102826108) that allows the information
#' contained in an expression matrix to be condensed into a combinatory graph.
#' The novelty is that information on survival is integrated into the analysis.
#'
#' The analysis consists of 3 parts: a preprocessing of the data, the gene
#' selection and the filter function, and the mapper algorithm. The
#' preprocessing is specifically the Disease Specific Genomic Analysis (proposed
#' by Nicolau *et al.*) that consists of, through linear models, eliminating the
#' part of the data that is considered "healthy" and keeping only the component
#' that is due to the disease. The genes are then selected according to their
#' variability and whether they are related to survival and the values of the
#' filtering function for each patient are calculated taking into account the
#' survival associated with each gene. Finally, the mapper algorithm is applied
#' from the disease component matrix and the values of the filter function
#' obtaining a combinatory graph.
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
#' @param gen_select_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param percent_gen_select Percentage (from zero to one hundred) of genes
#' to be selected to be used in mapper. 10 default option.
#' @param num_intervals Parameter for the mapper algorithm. Number of
#' intervals used to create the first sample partition based on
#' filtering values. 5 default option.
#' @param percent_overlap Parameter for the mapper algorithm. Percentage
#' of overlap between intervals. Expressed as a percentage. 40 default option.
#' @param distance_type Parameter for the mapper algorithm.
#' Type of distance to be used for clustering. Choose between correlation
#' ("cor") and euclidean ("euclidean"). "cor" default option.
#' @param clustering_type Parameter for the mapper algorithm. Type of
#' clustering method. Choose between "hierarchical" and "PAM"
#' (“partition around medoids”) options. "hierarchical" default option.
#' @param num_bins_when_clustering Parameter for the mapper algorithm.
#' Number of bins to generate the histogram employed by the standard
#' optimal number of cluster finder method. Parameter not necessary if the
#' "optimal_clust_mode" option is "silhouette" or the "clust_type" is "PAM".
#' 10 default option.
#' @param linkage_type Parameter for the mapper algorithm. Linkage criteria
#' used in hierarchical clustering. Choose between "single" for single-linkage
#' clustering, "complete" for complete-linkage clustering or "average" for
#' average linkage clustering (or UPGMA). Only necessary for hierarchical
#' clustering. "single" default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{GSSTDA} object. It contains:
#' - the matrix with the normal space \code{normal_space},
#' - the matrix of the disease components normal_space \code{matrix_disease_component},
#' - a matrix with the results of the application of proportional hazard models
#' for each gene (\code{cox_all_matrix)},
#' - the genes selected for mapper \code{genes_disease_componen},
#' - the matrix of the disease components with information from these genes only
#' \code{genes_disease_component}
#' - and a \code{mapper_obj} object. This \code{mapper_obj} object contains the
#' values of the intervals (interval_data), the samples included in each
#' interval (sample_in_level), information about the cluster to which the
#' individuals in each interval belong (clustering_all_levels), a list including
#' the individuals contained in each detected node (node_samples), their size
#' (node_sizes), the average of the filter function values of the individuals
#' of each node (node_average_filt) and the adjacency matrix linking the nodes
#' (adj_matrix). Moreover, information is provided on the number of nodes,
#' the average node size, the standard deviation of the node size, the number
#' of connections between nodes, the proportion of connections to all possible
#' connections and the number of ramifications.
#' @export
#' @examples
#' \dontrun{
#' GSSTDA <- GSSTDA(full_data,  survival_time, survival_event, case_tag,
#'                  gen_select_type="Top_Bot", percent_gen_select=10,
#'                  num_intervals = 4, percent_overlap = 50,
#'                  distance_type = "euclidean", num_bins_when_clustering = 8,
#'                  clustering_type = "hierarchical", linkage_type = "single")}
GSSTDA <- function(full_data, survival_time, survival_event, case_tag, gen_select_type="Top_Bot",
                   percent_gen_select=10, num_intervals=5, percent_overlap=40, distance_type="cor",
                   clustering_type="hierarchical", num_bins_when_clustering=10, linkage_type="single",
                   na.rm=TRUE){
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)
  #Select the control_tag. This do it inside of the DGSA function
  #Check and obtain gene selection (we use in the gene_select_surv). It execute in Block II
  #num_gen_select <- check_gene_selection(nrow(full_data), gen_select_type, percent_gen_select)

  #Don't check filter_values because it is not created.
  filter_values <- c()
  check_return <- check_arg_mapper(full_data, filter_values, distance_type, clustering_type,
                                              linkage_type, na.rm)

  full_data <- check_return[[1]]
  filter_values <- check_return[[2]]
  optimal_clustering_mode <- check_return[[3]]


  ################### BLOCK I: Pre-process. DGSA (using "NT" control_tag) ##############################
  DGSA_obj <- DGSA(full_data, survival_time, survival_event, case_tag, na.rm = "checked")
  matrix_disease_component <- DGSA_obj[["matrix_disease_component"]]
  control_tag <- DGSA_obj[["control_tag"]]
  full_data <- DGSA_obj[["full_data"]]
  survival_event <- DGSA_obj[["survival_event"]]
  survival_time <- DGSA_obj[["survival_time"]]
  case_tag <- DGSA_obj[["case_tag"]]

  ################### BLOCK II: Gene selection (using "T" control_tag) ##################################
  geneSelection_object <- geneSelection(DGSA_obj, gen_select_type, percent_gen_select)
  cox_all_matrix <- geneSelection_object[["cox_all_matrix"]]
  genes_selected <- geneSelection_object[["genes_selected"]]
  genes_disease_component <- geneSelection_object[["genes_disease_component"]]
  filter_values <- geneSelection_object[["filter_values"]]

  ################### BLOCK III: Create mapper object where the arguments are checked ###################
  cat("\nBLOCK III: The mapper process is started")

  # Transpose genes_disease_component: rows = patient, columns = genes
  #genes_disease_component <- t(genes_disease_component)

  #   Check filter_values
  check_filter <- check_filter_values(genes_disease_component, filter_values)
  genes_disease_component <- check_filter[[1]]
  filter_values <- check_filter[[2]]

  mapper_obj <- mapper(genes_disease_component, filter_values, num_intervals, percent_overlap, distance_type,
                       clustering_type, num_bins_when_clustering, linkage_type, optimal_clustering_mode,
                       na.rm = "checked")

  cat("\nBLOCK III: The mapper process is finished")


  ############################################  Create the object #########################################
  GSSTDA_object <- list("normal_space" = DGSA_obj[["normal_space"]],
                        "matrix_disease_component" = matrix_disease_component,
                        "cox_all_matrix" = cox_all_matrix,
                        "genes_selected" = genes_selected,
                        "genes_disease_component" = genes_disease_component,
                        "mapper_obj" = mapper_obj
                        )

  class(GSSTDA_object) <- "GSSTDA_obj"
  return(GSSTDA_object)
}
