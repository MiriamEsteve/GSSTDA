% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DSGA.R
\name{plot_dsga}
\alias{plot_dsga}
\title{plot dsga}
\usage{
plot_dsga(selected_matrix_disease_component, case_tag)
}
\arguments{
\item{selected_matrix_disease_component}{Disease component matrix of
the selected genes that contains the disease component of all patients.
Output of the function \code{generate_disease_component}.}

\item{case_tag}{Character vector of the same length as the number of
columns of full_data. Patients must be in the same order as in full_data.
It must be indicated for each patient whether he/she is healthy or not.
One value should be used to indicate whether the patient is healthy and
another value should be used to indicate whether the patient's sample is
tumourous. The user will then be asked which one indicates whether
the patient is healthy. Only two values are valid in the vector in total.}
}
\value{
The heatmap of the dsga result.
}
\description{
It draws the heatmap of the dsga result by selecting the 100 genes with
the highest variability between samples.
}
