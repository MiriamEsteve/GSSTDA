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
  print("BLOCK I: The pre-process DGSA is started")
  #   Select the normal tissue data gene expression matrix.
  normal_tiss <- full_data[,which(case_tag == control_tag)]

  #   Obtain the gene expression matrix containing the flattened version of the vectors.
  matrix_flatten_normal_tiss <- flatten_normal_tiss(normal_tiss)
  #   Obtain the normal space
  normal_space <- denoise_rectangular_matrix(matrix_flatten_normal_tiss)
  #   Obtain the disease component of the normal_space
  matrix_disease_component <- generate_disease_component(full_data, normal_space)

  print("BLOCK I: The pre-process DGSA is finished")

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


gene_selection <- function(data, control_tag_cases, survival_time, survival_event){
  print("BLOCK II: The gene selection is started")
  #Remove NAN's values (case_tag == control_tag) of survival_time and survival_event
  survival_time <- survival_time[-control_tag_cases]
  survival_event <- survival_event[-control_tag_cases]
  #Select the disease component of the "T" control_tag
  case_disease_component <- data[,-control_tag_cases]

  # Univariate cox proportional hazard models for the expression levels of each gene included in the
  #provided dataset
  cox_all_matrix <- cox_all_genes(case_disease_component, survival_time, survival_event)

  #Selects genes for mapper
  genes_selected <- gene_selection_surv(case_disease_component, cox_all_matrix, gen_select_type,
                                        num_gen_select)

  # Select genes in matrix_disease_component or full_data (if don't apply Block I)
  genes_disease_component <- data[genes_selected,]

  # Filter the genes_disease_component
  filter_values <- lp_norm_k_powers_surv(genes_disease_component, 2, 1, cox_all_matrix)

  print("BLOCK II: The gene selection is finished")

  geneSelection_object <- list( "data" = data,
                                "cox_all_matrix" = cox_all_matrix,
                                "genes_selected" = genes_selected,
                                "genes_disease_component" = genes_disease_component,
                                "filter_values" = filter_values
  )
  class(geneSelection_object) <- "geneSelection_obj"
}

#Generic function
geneSelection <- function(full_data,  survival_time, survival_event, case_tag, gen_select_type,
                          percent_gen_select, na.rm = TRUE){
  UseMethod("gene_selection_classes")
}

gene_selection_classes.DGSA_object <- function(x, gen_select_type, percent_gen_select){
  print(class(x))

  matrix_disease_component <- x[["matrix_disease_component"]]
  #Check and obtain gene selection (we use in the gene_select_surv)
  num_gen_select <- check_gene_selection(nrow(matrix_disease_component), gen_select_type, percent_gen_select)

  control_tag <- x[["control_tag"]]
  survival_event <- x[["survival_event"]]
  survival_time <- x[["survival_time"]]
  case_tag <- x[["case_tag"]]

  control_tag_cases <- which(case_tag == control_tag)
  geneSelection_object <- gene_selection(matrix_disease_component, control_tag_cases, survival_time, survival_event)

  return(geneSelection_object)
}

gene_selection_classes.matrix <- function(data, survival_time, survival_event, case_tag, gen_select_type,
                                          percent_gen_select, na.rm = TRUE){
  print("gene_selection_classes")
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(data, na.rm)
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

  geneSelection_object <- gene_selection(full_data, control_tag_cases, survival_time, survival_event)

  return(geneSelection_object)
}

