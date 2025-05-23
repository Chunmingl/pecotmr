% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocboost_pipeline.R
\name{load_multitask_regional_data}
\alias{load_multitask_regional_data}
\title{This function loads a mixture data sets for a specific region, including individual-level data (genotype, phenotype, covariate data)
or summary statistics (sumstats, LD). Run \code{load_regional_univariate_data} and \code{load_rss_data} multiple times for different datasets}
\usage{
load_multitask_regional_data(
  region,
  genotype_list = NULL,
  phenotype_list = NULL,
  covariate_list = NULL,
  conditions_list_individual = NULL,
  match_geno_pheno = NULL,
  maf_cutoff = 0,
  mac_cutoff = 0,
  xvar_cutoff = 0,
  imiss_cutoff = 0,
  association_window = NULL,
  extract_region_name = NULL,
  region_name_col = NULL,
  keep_indel = TRUE,
  keep_samples = NULL,
  keep_variants = NULL,
  phenotype_header = 4,
  scale_residuals = FALSE,
  tabix_header = TRUE,
  sumstat_path_list = NULL,
  column_file_path_list = NULL,
  LD_meta_file_path_list = NULL,
  match_LD_sumstat = NULL,
  conditions_list_sumstat = NULL,
  n_samples = 0,
  n_cases = 0,
  n_controls = 0,
  extract_sumstats_region_name = NULL,
  sumstats_region_name_col = NULL,
  comment_string = "#",
  extract_coordinates = NULL
)
}
\arguments{
\item{region}{The region where tabix use to subset the input dataset.}

\item{genotype_list}{a vector of PLINK bed file containing genotype data.}

\item{phenotype_list}{A vector of phenotype file names.}

\item{covariate_list}{A vector of covariate file names corresponding to the phenotype file vector.}

\item{conditions_list_individual}{A vector of strings representing different conditions or groups.}

\item{match_geno_pheno}{A vector of index of phentoypes matched to genotype if mulitple genotype PLINK files}

\item{maf_cutoff}{Minimum minor allele frequency (MAF) cutoff. Default is 0.}

\item{mac_cutoff}{Minimum minor allele count (MAC) cutoff. Default is 0.}

\item{xvar_cutoff}{Minimum variance cutoff. Default is 0.}

\item{imiss_cutoff}{Maximum individual missingness cutoff. Default is 0.}

\item{association_window}{A string of chr:start-end for the association analysis window (cis or trans). If not provided, all genotype data will be loaded.}

\item{extract_region_name}{A list of vectors of strings (e.g., gene ID ENSG00000269699) to subset the information when there are multiple regions available. Default is NULL.}

\item{region_name_col}{Column name containing the region name. Default is NULL.}

\item{keep_indel}{Logical indicating whether to keep insertions/deletions (INDELs). Default is TRUE.}

\item{keep_samples}{A vector of sample names to keep. Default is NULL.}

\item{phenotype_header}{Number of rows to skip at the beginning of the transposed phenotype file (default is 4 for chr, start, end, and ID).}

\item{scale_residuals}{Logical indicating whether to scale residuals. Default is FALSE.}

\item{tabix_header}{Logical indicating whether the tabix file has a header. Default is TRUE.}

\item{sumstat_path_list}{A vector of file path to the summary statistics.}

\item{column_file_path_list}{A vector of file path to the column file for mapping.}

\item{LD_meta_file_path_list}{A vector of path of LD_metadata, LD_metadata is a data frame specifying LD blocks with columns "chrom", "start", "end", and "path". "start" and "end" denote the positions of LD blocks. "path" is the path of each LD block, optionally including bim file paths.}

\item{match_LD_sumstat}{A vector of index of sumstat matched to LD if mulitple sumstat files}

\item{conditions_list_sumstat}{A vector of strings representing different sumstats.}

\item{n_samples}{User-specified sample size. If unknown, set as 0 to retrieve from the sumstat file.}

\item{n_cases}{User-specified number of cases.}

\item{n_controls}{User-specified number of controls.}

\item{extract_sumstats_region_name}{User-specified gene/phenotype name used to further subset the phenotype data.}

\item{sumstats_region_name_col}{Filter this specific column for the extract_sumstats_region_name.}

\item{comment_string}{comment sign in the column_mapping file, default is #}

\item{extract_coordinates}{Optional data frame with columns "chrom" and "pos" for specific coordinates extraction.}
}
\value{
A list containing the individual_data and sumstat_data:
individual_data contains the following components if exist
\itemize{
  \item residual_Y: A list of residualized phenotype values (either a vector or a matrix).
  \item residual_X: A list of residualized genotype matrices for each condition.
  \item residual_Y_scalar: Scaling factor for residualized phenotype values.
  \item residual_X_scalar: Scaling factor for residualized genotype values.
  \item dropped_sample: A list of dropped samples for X, Y, and covariates.
  \item covar: Covariate data.
  \item Y: Original phenotype data.
  \item X_data: Original genotype data.
  \item X: Filtered genotype matrix.
  \item maf: Minor allele frequency (MAF) for each variant.
  \item chrom: Chromosome of the region.
  \item grange: Genomic range of the region (start and end positions).
  \item Y_coordinates: Phenotype coordinates if a region is specified.
}
sumstat_data contains the following components if exist
\itemize{
  \item sumstats: A list of summary statistics for the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
  \item LD_info: A list of LD information, each sublist contains combined_LD_variants, combined_LD_matrix, ref_panel  \code{load_LD_matrix}.
}
}
\description{
This function loads a mixture data sets for a specific region, including individual-level data (genotype, phenotype, covariate data)
or summary statistics (sumstats, LD). Run \code{load_regional_univariate_data} and \code{load_rss_data} multiple times for different datasets
}
\section{Loading individual level data from multiple corhorts}{
NA
}

\section{Loading summary statistics from multiple corhorts or data set}{
NA
}

