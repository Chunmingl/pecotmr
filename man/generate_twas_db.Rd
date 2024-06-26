% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctwas_wrapper.R
\name{generate_twas_db}
\alias{generate_twas_db}
\title{cTWAS (causal-TWAS) wrapper
Determine Imputability and Variant Selection

This function load TWAS weights and assesses the imputability of genes across different contexts, 
and select top variants for imputable gene-context pair, and output the extracted variant weights 
from best performing method model. Imputability of a gene-context pair are determined by having 
least one of the twas methods' model being imputable based on cross livadtion metrics. The imputable 
model is if the model surpasses the r-square and p-value threashold from cross validation metrics. 
If non of the contexts of a gene has at least one method being imputable, then this gene will be 
considered as unimputable, and do not return any weight results. For imputable gene-context pair, 
we select a subset of variants with high Posterior Inclusion Probabilities (PIPs) from SuSiE 
fine-mapping results for the imputable gene-context pair. After selecting variants, we extracts 
variant weights from the best performing model among imputable models with highest r-square or 
smallest p-value from cross validation metrics. This function is essential for preparing data 
for complex trait weighted analysis with sparse weight.}
\usage{
generate_twas_db(
  weight_db_file,
  contexts = NULL,
  variable_name_obj = c("preset_variants_result", "variant_names"),
  twas_weights_table = "twas_weights",
  max_var_selection,
  min_rsq_threshold = 0.01,
  p_val_cutoff = 0.05,
  data_type = NULL
)
}
\arguments{
\item{variable_name_obj}{Name of the element in SuSiE-TWAS output that lists variant names aligned with TWAS weights.}

\item{max_var_selection}{Maximum number of SNPs to be selected from each gene-context pair.}

\item{min_rsq_threshold}{Minimum \(R^2\) threshold for determining imputability of a gene-context pair.}

\item{p_val_cutoff}{Maximum allowable p-value for a gene to be considered imputable in a given context.}

\item{weight_db_files}{File paths to `.rds` files containing SuSiE-TWAS weights.}
}
\value{
A list of elements containing fine-mapping and TWAS related data:
\itemize{
  \item{refined_twas_weights}{Context imputability, model selection, and selected variants, and extracted 
weights for the selected variants across all contexts for a imputable gene. }
  \item{susie_results}{SuSiE fine-mapping results, including variant selection from the SuSiE-TWAS pipeline.}
  \item{weights}{TWAS weights from the SuSiE-TWAS pipeline, to be harmonized for further analysis.}
  \item{gene}{Gene name of the twas db weight.}
  \item{cv_performance}{Cross validation performance metrics for all contexts.}
}
}
\description{
Causal inference for TWAS with Summary Statistics with Multigroup Data
}
\examples{
results <- generate_twas_db(weight_db_file = "path/to/weights.rds", conditions=c("Mic", "Oli", "Exc"),
                                  max_var_selection = 10, min_rsq_threshold = 0.01, p_val_cutoff = 0.05)
}
