#' Gene expression matrix
#'
#' A matrix containing gene expression profiling of 104 breast cancer and
#' 17 normal breast biopsies.
#' @name full_data
#' @format Gene expression matrix.
#' \describe{
#'   \item{full_data}{Gene expression matrix of 104 breast cancer samples
#'   and 17 healthy tissue samples. Expression profiling data by array
#'   (HG-U133_Plus_2).The columns correspond to the patients and the rows to
#'   the genes. The column names correspond to the patient identifier
#'   in GEO. The row names correspond to the gene names.}
#' }
#' @source The data are from the study GSE42568 available in
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
#' The data were processed as explained in the \code{details} section.
#' @details
#' Normalized gene expression data GSE42568. Background correction,
#' summarization, and quantile normalization were carried out using the
#' fRMA method implemented in the fRMA package. Filtered probes that did
#' not target genes with valid gene id were filtered and  those probes
#' targeting the same gene were collapse by thanking those presenting the
#' highest row variance using the \code{WGCNA::collapseRows} function.
"full_data"

#' Case-control vector
#' Character vector of length 121 containing the group to which each
#' sample belongs.
#' @name case_tag
#' @format  Character vector
#' \describe{
#'   \item{code_tag}{Character vector of length 121 containing the group to
#'   which each sample belongs ("NT": control, sample from healthy tissue;
#'   "T": case, sample from neoplastic tissue).}
#' }
#' @source The data are from the study GSE42568. Information extracted from
#' the file GSE42568_family.soft.gz available at
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
"case_tag"

#' Survival event vector
#' Character vector of length 121 containing whether or not the patient
#' is deceased.
#' @name survival_event
#' @format  Character vector
#' \describe{
#'   \item{survival_event}{Character vector of length 121 containing whether
#'   or not the patient is deceased ("0": no, "1": yes). Samples from healthy
#'   tissue contain a value of \code{NA}.}
#' }
#' @source The data are from the study GSE42568. Information extracted from
#' the file GSE42568_family.soft.gz available at
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568}.
"survival_event"

#' Survival time vector
#' Numeric vector of length 121 containing the time in months until
#' the death of the patient or until the end of the follow-up in case the
#' patient has not died.
#' @name survival_time
#' @format  Numeric vector
#' \describe{
#'   \item{survival_time}{Numeric vector of length 121 containing the time
#'   in months until the death of the patient or until the end of the follow-up
#'   in case the patient has not died. Samples from healthy tissue contain a
#'   value of \code{NA}.}
#' }
"survival_time"


