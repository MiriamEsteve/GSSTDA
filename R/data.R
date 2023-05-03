#' Gene expression matrix
#'
#' Matrix containing gene expression profiling of 104 breast cancer and
#' 17 normal breast biopsies. Expression profiling data by array
#'   (platform HG-U133_Plus_2).
#'
#' @name full_data
#' @format Gene expression matrix with  20825 rows and 121 columns.
#' \describe{The columns correspond to the patients and the rows to
#'   the genes. The column names correspond to the patient identifier
#'   in GEO. The row names correspond to the gene names.
#' }
#' @source The data are from the study GSE42568 available in
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
#' The data were processed as explained in the \code{details} section.
#' @details
#' Normalized gene expression data GSE42568. Background correction,
#' summarization, and quantile normalization were carried out using the
#' fRMA method implemented in the \code{fRMA} package. Filtered probes that did
#' not target genes with valid gene id were filtered and  those probes
#' targeting the same gene were collapse by thanking those presenting the
#' highest row variance using the \code{WGCNA::collapseRows} function.
#' @usage
#' data(full_data, package = "GSSTDA")
"full_data"

#' Case-control vector
#'
#' Character vector of length 121 containing the group to which each
#' sample belongs.
#'
#' @name case_tag
#' @format  Character vector length 121.
#' \describe{"NT": control, sample from healthy tissue;
#'   "T": case, sample from neoplastic tissue.
#' }
#' @usage
#' data(case_tag, package = "GSSTDA")
#' @source The data are from the study GSE42568. Information extracted from
#' the file GSE42568_family.soft.gz available at
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
"case_tag"

#' Survival event vector
#'
#' Character vector of length 121 containing whether or not the patient
#' is deceased.
#'
#' @name survival_event
#' @format  Character vector of length 121.
#' \describe{A value of "0" indicates that the patient did not pass away
#'  during follow-up, a value of "1" indicates that the patient did. Samples
#'  from healthy tissue contain a value of \code{NA}.
#' }
#' @usage
#' data(survival_event, package = "GSSTDA")
#' @source The data are from the study GSE42568. Information extracted from
#' the file GSE42568_family.soft.gz available at
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
"survival_event"

#' Survival time vector
#'
#' Numeric vector of length 121 containing the time in months until
#' the death of the patient or until the end of the follow-up in case the
#' patient has not passed away.
#'
#' @name survival_time
#' @format  Numeric vector of length 121.
#' \describe{Time in months. Samples from healthy tissue contain a
#'   value of \code{NA}.
#' }
#' @usage
#' data(survival_time, package = "GSSTDA")
#' @source The data are from the study GSE42568. Information extracted from
#' the file GSE42568_family.soft.gz available at
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}
"survival_time"


