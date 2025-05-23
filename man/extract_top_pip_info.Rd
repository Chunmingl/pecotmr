% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_rss_diagnostics.R
\name{extract_top_pip_info}
\alias{extract_top_pip_info}
\title{Extract Information for Top Variant from Finemapping Results}
\usage{
extract_top_pip_info(con_data)
}
\arguments{
\item{con_data}{List. The method layer data from a finemapping RDS file.}
}
\value{
A data frame with one row containing the following columns:
  \item{cs_name}{NA (as no CS is identified)}
  \item{variants_per_cs}{NA (as no CS is identified)}
  \item{top_variant}{ID of the variant with the highest PIP}
  \item{top_variant_index}{Index of the top variant in the original data}
  \item{top_pip}{Highest Posterior Inclusion Probability (PIP)}
  \item{top_z}{Z-score of the top variant}
  \item{p_value}{P-value calculated from the top Z-score}
  \item{cs_corr}{NA (as no CS correlation is available)}
}
\description{
This function extracts information about the variant with the highest Posterior 
Inclusion Probability (PIP) from finemapping results, typically used when no 
Credible Sets (CS) are identified in the analysis.
}
\details{
This function is designed to be used when no Credible Sets are identified in 
the finemapping results, but information about the most significant variant 
is still desired. It identifies the variant with the highest PIP and extracts 
relevant statistical information.
}
\note{
This function is particularly useful for capturing information about potentially 
important variants that might be included in Credible Sets under different 
analysis parameters or lower coverage. It maintains a structure similar to 
the output of `extract_cs_info()` for consistency in downstream analyses.
}
\seealso{
\code{\link{extract_cs_info}} for processing when Credible Sets are present.
}
