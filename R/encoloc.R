#' xQTL GWAS Enrichment Analysis
#'
#' This function processes GWAS and xQTL finemapped data files and then computes QTL enrichment.
#' For details on the parameters `pi_gwas`, `pi_qtl`, `lambda`, `ImpN`, and `num_threads`,
#' refer to the documentation of the `compute_qtl_enrichment` function.
#'
#' @param xqtl_files Vector of xQTL RDS file paths.
#' @param gwas_files Vector of GWAS RDS file paths.
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_varname_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_varname_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param pi_gwas Optional parameter for GWAS enrichment estimation (see `compute_qtl_enrichment`).
#' @param pi_qtl Optional parameter for xQTL enrichment estimation (see `compute_qtl_enrichment`).
#' @param lambda Shrinkage parameter for enrichment computation (see `compute_qtl_enrichment`).
#' @param ImpN Importance parameter for enrichment computation (see `compute_qtl_enrichment`).
#' @param num_threads Number of threads for parallel processing (see `compute_qtl_enrichment`).
#' @return The output from the compute_qtl_enrichment function.
#' @examples
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' xqtl_files <- c("xqtl_file1.rds", "xqtl_file2.rds")
#' result <- xqtl_enrichment_wrapper(gwas_files, xqtl_files)
#' @export
xqtl_enrichment_wrapper <- function(xqtl_files, gwas_files,
                                    xqtl_finemapping_obj = NULL, gwas_finemapping_obj = NULL,
                                    xqtl_varname_obj = NULL, gwas_varname_obj = NULL,
                                    num_gwas = NULL, pi_qtl = NULL,
                                    lambda = 1.0, ImpN = 25,
                                    num_threads = 1) {
  process_finemapped_data <- function(xqtl_files, gwas_files,
                                      xqtl_finemapping_obj = NULL, gwas_finemapping_obj = NULL,
                                      xqtl_varname_obj = NULL, gwas_varname_obj = NULL) {
    # Load and process GWAS data
    gwas_pip <- list()
    for (file in gwas_files) {
      raw_data <- readRDS(file)[[1]] # changed
      gwas_data <- if (!is.null(gwas_finemapping_obj)) get_nested_element(raw_data, gwas_finemapping_obj) else raw_data
      pip <- gwas_data$pip
      if (!is.null(gwas_varname_obj)) names(pip) <- get_nested_element(raw_data, gwas_varname_obj)
      gwas_pip <- c(gwas_pip, list(pip))
    }

    # Check for unique variant names in GWAS pip vectors
    all_variant_names <- unique(unlist(lapply(gwas_pip, names)))
    if (length(unique(all_variant_names)) != length(all_variant_names)) {
      stop("Non-unique variant names found in GWAS data with different pip values.")
    }
    gwas_pip <- unlist(gwas_pip)

    # Process xQTL data
    xqtl_data <- lapply(xqtl_files, function(file) {
      raw_data <- readRDS(file)[[1]]
      xqtl_data <- tryCatch(
        {
          if (!is.null(xqtl_finemapping_obj)) get_nested_element(raw_data, xqtl_finemapping_obj) else raw_data
        },
        error = function(e) {
          return(NULL)
        }
      )
      if (!is.null(xqtl_data)) {
        list(
          alpha = xqtl_data$alpha,
          pip = setNames(xqtl_data$pip, get_nested_element(raw_data, xqtl_varname_obj)),
          prior_variance = xqtl_data$V
        )
      } else {
        NULL
      }
    })

    # Return results as a list
    return(list(gwas_pip = gwas_pip, xqtl_data = xqtl_data))
  }

  # Load data
  dat <- process_finemapped_data(xqtl_files, gwas_files, xqtl_finemapping_obj, gwas_finemapping_obj, xqtl_varname_obj, gwas_varname_obj)
  # Compute QTL enrichment
  return(compute_qtl_enrichment(
    gwas_pip = dat$gwas_pip, susie_qtl_regions = dat$xqtl_data,
    num_gwas = num_gwas, pi_qtl = pi_qtl,
    lambda = lambda, ImpN = ImpN,
    num_threads = num_threads
  ))
}

#' Function to filter and order colocalization results
#' @noRd
filter_and_order_coloc_results <- function(coloc_results_fil) {
  # Ensure the input has more than one column
  if (ncol(coloc_results_fil) <= 1) {
    stop("Insufficient number of columns in colocalization results")
  }

  cs_num <- ncol(coloc_results_fil) - 1
  ordered_results <- list()

  for (n in 1:cs_num) {
    # Selecting relevant columns and ordering
    tmp_coloc_results_fil <- coloc_results_fil[, c(1, n + 1)] %>% .[order(.[, 2], decreasing = TRUE), ]
    ordered_results[[n]] <- tmp_coloc_results_fil
  }

  return(ordered_results)
}

#' Function to calculate cumulative sum
#' @noRd
calculate_cumsum <- function(coloc_results) {
  cumsum(coloc_results[, 2])
}

#' Function to load and extract LD matrix
#' @noRd
load_and_extract_ld_matrix <- function(ld_meta_file_path, analysis_region, variants, ld_ref = NULL, in_sample = NULL) {
  # Extract variant positions
  var_pos <- str_split(variants, ":", simplify = TRUE)[,2] %>% as.numeric()
  min_var_pos <- min(var_pos)
  max_var_pos <- max(var_pos)
  chr <- str_split(analysis_region, ":", simplify = TRUE)[,1]
  analysis_region_narrow <- paste0(chr, ":", min_var_pos, "-", max_var_pos)

  # --- Determine mode if not explicitly provided ---
  if (is.null(ld_ref) && is.null(in_sample)) {
    ld_ref <- TRUE  # Default assumption
    if (grepl("plink|genotype|\\.bed$|\\.bim$|\\.fam$", ld_meta_file_path)) {
      ld_ref <- FALSE
      in_sample <- TRUE
    } else {
      in_sample <- FALSE
    }
  }

  # --- Enforce exclusivity ---
  if (isTRUE(ld_ref)) in_sample <- FALSE
  if (isTRUE(in_sample)) ld_ref <- FALSE

  # --- LD reference mode ---
  if (ld_ref) {
    message("Using LD reference mode")
    ld_ref_data <- load_LD_matrix(
      LD_meta_file_path = ld_meta_file_path,
      region = analysis_region_narrow
    )
    ext_ld <- ld_ref_data$combined_LD_matrix[variants, variants]
    return(ext_ld)
  }

  # --- In-sample genotype mode ---
  if (in_sample) {
    message("Using in-sample genotype mode")

    if (grepl("\\.txt$", ld_meta_file_path)) {
      geno_meta <- read_tsv(ld_meta_file_path, comment = "#", col_names = c("id", "path"), show_col_types = FALSE)
      chr_num <- gsub("^chr", "", chr)
      geno_path <- geno_meta %>% filter(id == chr_num) %>% pull(path) %>% basename %>% paste0(dirname(ld_meta_file_path), "/", .)

      if (length(geno_path) != 1) stop("No matching entry found in metadata for chromosome ", chr)

      geno_prefix <- str_remove(geno_path, "\\.bed$")
    } else if (grepl("\\.bed$", ld_meta_file_path)) {
      geno_prefix <- str_remove(ld_meta_file_path, "\\.bed$")
    } else {
      stop("In in-sample mode, expected plink file or .txt genotype metadata file.")
    }

    # Load original genotype data 
    ld_mat <- load_genotype_region(geno_prefix, region = analysis_region_narrow, keep_indel = TRUE, keep_variants_path = NULL )
      
    # Change colname format of genotype data
    colnames(ld_mat) <- align_variant_names(colnames(ld_mat), variants)$aligned_variants
    ld_mat <- ld_mat[, variants]  # subset target variants then get imputation and correlation
    if(length(variants) > 1) {
        # Mean imputation
        ld_mat_imputed <- apply(ld_mat, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
        ext_ld <- get_cormat(ld_mat_imputed)
    } else {
        ext_ld = as.matrix(1)
    }
    return(ext_ld)
  }

  stop("Neither LD mode was activated — check inputs.")
}

#' Function to calculate purity
#' @noRd
calculate_purity <- function(variants, ext_ld, squared = FALSE) {
  # This is a placeholder for calculating purity, adjust as per your actual function
  purity <- matrix(susieR:::get_purity(variants, Xcorr = ext_ld, squared), 1, 3)
  purity
}

#' Main processing function
#' This function is designed to summarize coloc results based on the following criteria:
#' 1. Among the colocalized variant pairs, PPH4 has the highest value compared to PPH0-PPH3.
#' 2. PPH4 exceeds threshold, default as 0 since we advocate not using PPH4 concept but rather use CoS
#' 3. We aggregate variants and cumulatively sum their PPH4 values to form a credible set until the threshold, default as 0.95.
#' 4. The cs's purity is computed with the `get_purity` function from the `gaow/susieR` package, and the same purity criteria are employed to filter the credibility set.
#' @noRd
process_coloc_results <- function(coloc_result, LD_meta_file_path, analysis_region, PPH4_thres = 0, coverage = 0.95, min_abs_corr = 0.8, null_index = 0, coloc_index = "PP.H4.abf") {
  # Extract PIP values from coloc_result summary
  coloc_summary <- as.data.frame(coloc_result$summary)
  coloc_pip <- coloc_summary[, grepl("PP", colnames(coloc_summary))]

  # Filter and extract relevant columns from coloc_result results
  # PP.H4 is highest and > 0.8
  coloc_results_df <- as.data.frame(coloc_result$results)
  coloc_filter <- apply(coloc_pip, 1, function(row) {
    max_index <- which.max(row)
    max_value <- row[max_index]
    return(max_value > PPH4_thres && colnames(coloc_pip)[max_index] == coloc_index)
  })

  coloc_res <- list()

  if (sum(coloc_filter) > 0) {
    coloc_results_fil <- coloc_results_df[, c(1, which(coloc_filter) + 1), drop = FALSE]
    coloc_summary_fil <- coloc_summary[which(coloc_filter), , drop = FALSE]

    # prepare to calculate purity
    ordered_results <- filter_and_order_coloc_results(coloc_results_fil)
    cs <- list()
    purity <- NULL

    for (n in 1:length(ordered_results)) {
      tmp_coloc_results_fil <- ordered_results[[n]]
      tmp_coloc_results_fil_csm <- calculate_cumsum(tmp_coloc_results_fil)
      cs[[n]] <- tmp_coloc_results_fil[, 1][1:(which(tmp_coloc_results_fil_csm > coverage) %>% min())]
      variants <- cs[[n]] %>% gsub("chr", "", .)

      # Load and extract LD matrix
      ext_ld <- load_and_extract_ld_matrix(LD_meta_file_path, analysis_region, variants)

      # Calculate purity
      if (null_index > 0 && null_index %in% variants) {
        purity <- rbind(purity, c(-9, -9, -9))
      } else {
        current_purity <- calculate_purity(variants, ext_ld)
        purity <- rbind(purity, current_purity)
      }
    }

    # Process purity data
    purity <- as.data.frame(purity)
    colnames(purity) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
    is_pure <- which(purity[, 1] >= min_abs_corr)

    # Finalize the result
    if (length(is_pure) > 0) {
      cs <- cs[is_pure]
      purity <- purity[is_pure, ]
      true_summary <- coloc_summary_fil[is_pure, ]
      coloc_res$sets <- list(cs = cs, purity = purity, true_summary = true_summary)
    }
  } else {
    message("Coloc results did not find any variants that satisfy the condition of PP.H4 being the highest value and > ", PPH4_thres)
    coloc_res$sets <- list(cs = NULL)
  }

  return(coloc_res)
}

#' Colocalization Analysis Wrapper
#'
#' This function processes xQTL and multiple GWAS finemapped data files for colocalization analysis.
#'
#' @param xqtl_file Path to the xQTL RDS file.
#' @param gwas_files Vector of paths to GWAS RDS files.
#' @param xqtl_finemapping_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_finemapping_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_varname_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_varname_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param xqtl_region_obj Optional table name in xQTL RDS files (default 'susie_fit').
#' @param gwas_region_obj Optional table name in GWAS RDS files (default 'susie_fit').
#' @param region_obj Optional table name of region info in susie_twas output filess (default 'region_info').
#' @param p1, p2, and p12 are results from xqtl_enrichment_wrapper (default 'p1=1e-4, p2=1e-4, p12=5e-6', same as coloc.bf_bf).
#' @param prior_tol When the prior variance is estimated, compare the estimated value to \code{prior_tol} at the end of the computation,
#'   and exclude a single effect from PIP computation if the estimated prior variance is smaller than this tolerance value.
#' @return A list containing the coloc results and the summarized sets.
#' @examples
#' xqtl_file <- "xqtl_file.rds"
#' gwas_files <- c("gwas_file1.rds", "gwas_file2.rds")
#' result <- coloc_wrapper(xqtl_file, gwas_files, LD_meta_file_path)
#' @importFrom dplyr bind_rows
#' @importFrom tidyr replace_na
#' @importFrom coloc coloc.bf_bf
#' @export
coloc_wrapper <- function(xqtl_file, gwas_files,
                          xqtl_finemapping_obj = NULL, xqtl_varname_obj = NULL, xqtl_region_obj = NULL,
                          gwas_finemapping_obj = NULL, gwas_varname_obj = NULL, gwas_region_obj = NULL,
                          filter_lbf_cs = FALSE, filter_lbf_cs_secondary = NULL,
                          prior_tol = 1e-9, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6, ...) {
  region <- NULL # define first to avoid element can not be found
  # Load and process GWAS data
  gwas_lbf_matrices <- lapply(gwas_files, function(file) {
    raw_data <- readRDS(file)[[1]]
    gwas_data <- if (!is.null(gwas_finemapping_obj)) get_nested_element(raw_data, gwas_finemapping_obj) else raw_data
    gwas_lbf_matrix <- as.data.frame(gwas_data$lbf_variable)
      # fsusie has a different structure
    if(is.null(gwas_lbf_matrix) | nrow(gwas_lbf_matrix)==0){
        gwas_lbf_matrix <- do.call(rbind, raw_data[[1]]$fsusie_result$lBF) %>% as.data.frame
        if(nrow(gwas_lbf_matrix) > 0) message("This is a fSuSiE case")
    }
    if (filter_lbf_cs & is.null(filter_lbf_cs_secondary)) {
        gwas_lbf_matrix <- gwas_lbf_matrix[gwas_data$sets$cs_index,, drop = F]
    } else if (!is.null(filter_lbf_cs_secondary)) {
        gwas_lbf_filter_index <- get_filter_lbf_index(gwas_data, coverage = filter_lbf_cs_secondary) 
        gwas_lbf_matrix <- gwas_lbf_matrix[gwas_lbf_filter_index,, drop = F]
    } else {
      gwas_lbf_matrix <- gwas_lbf_matrix[gwas_data$V > prior_tol,, drop = F]
    }
    if (!is.null(gwas_varname_obj)) colnames(gwas_lbf_matrix) <- get_nested_element(raw_data, gwas_varname_obj)
    # fsusie could have NA in variant name
    gwas_lbf_matrix <- gwas_lbf_matrix[,!is.na(colnames(gwas_lbf_matrix))] 
    return(gwas_lbf_matrix)
  })

  # changed: need to remove this check, as the last variant of former block could overlapped with first variant of later block
  # # Validate uniqueness of column names across GWAS matrices
  # all_gwas_colnames <- unique(unlist(lapply(gwas_lbf_matrices, colnames)))
  # if (length(all_gwas_colnames) != sum(sapply(gwas_lbf_matrices, ncol))) {
  #   stop("Duplicate variant names found across GWAS regions analyzed. This is not expected.")
  # }
  # Combine GWAS matrices and replace NAs with zeros
  combined_gwas_lbf_matrix <- bind_rows(gwas_lbf_matrices) %>%
    mutate(across(everything(), ~ replace_na(., 0)))

  # Process xQTL data
  xqtl_raw_data <- readRDS(xqtl_file)[[1]]
  xqtl_data <- if (!is.null(xqtl_finemapping_obj)) {
    tryCatch(
      {
        get_nested_element(xqtl_raw_data, xqtl_finemapping_obj)
      },
      error = function(e) {
        message(paste("no", xqtl_finemapping_obj[2], "in", xqtl_finemapping_obj[1]))
        NULL
      }
    )
  } else {
    xqtl_raw_data
  }
  if (!is.null(xqtl_data)) {
    xqtl_lbf_matrix <- as.data.frame(xqtl_data$lbf_variable)
    # fsusie has a different structure
    if(is.null(xqtl_lbf_matrix) | nrow(xqtl_lbf_matrix)==0){
        xqtl_lbf_matrix <- do.call(rbind, xqtl_raw_data[[1]]$fsusie_result$lBF) %>% as.data.frame
        if(nrow(xqtl_lbf_matrix) > 0) message("This is a fSuSiE case")
    }
    # fsusie data does not have V element in results
    if (filter_lbf_cs & is.null(filter_lbf_cs_secondary)) {
        xqtl_lbf_matrix <- xqtl_lbf_matrix[xqtl_data$sets$cs_index, , drop = F]
    } else if (!is.null(filter_lbf_cs_secondary)) {
        xqtl_lbf_filter_index <- get_filter_lbf_index(xqtl_data, coverage = filter_lbf_cs_secondary) 
        xqtl_lbf_matrix <- xqtl_lbf_matrix[xqtl_lbf_filter_index,, drop = F]
    } else {
      if ("V" %in% names(xqtl_data)) xqtl_lbf_matrix <- xqtl_lbf_matrix[xqtl_data$V > prior_tol, , drop = F] else (message("No V found in orginal data."))
    }
    if (nrow(combined_gwas_lbf_matrix) > 0 && nrow(xqtl_lbf_matrix) > 0) {
      if (!is.null(xqtl_varname_obj)) colnames(xqtl_lbf_matrix) <- get_nested_element(xqtl_raw_data, xqtl_varname_obj)
      # fsusie could have NA in variant name
      xqtl_lbf_matrix <- xqtl_lbf_matrix[,!is.na(colnames(xqtl_lbf_matrix))] 
        
      colnames(xqtl_lbf_matrix) <- align_variant_names(colnames(xqtl_lbf_matrix), colnames(combined_gwas_lbf_matrix))$aligned_variants
      common_colnames <- intersect(colnames(xqtl_lbf_matrix), colnames(combined_gwas_lbf_matrix))
      xqtl_lbf_matrix <- xqtl_lbf_matrix[, common_colnames, drop = FALSE] %>% as.matrix()
      combined_gwas_lbf_matrix <- combined_gwas_lbf_matrix[, common_colnames, drop = FALSE] %>% as.matrix()

      # Report the number of dropped columns from xQTL matrix
      num_dropped_cols <- length(setdiff(colnames(xqtl_lbf_matrix), common_colnames))
      message("Number of columns dropped from xQTL matrix: ", num_dropped_cols)

      # Function to convert region df to str
      convert_to_string <- function(df) paste0("chr", df$chrom, ":", df$start, "-", df$end)
      region <- if (!is.null(xqtl_region_obj)) get_nested_element(xqtl_raw_data, xqtl_region_obj) %>% convert_to_string() else NULL

      # COLOC function
      coloc_res <- coloc.bf_bf(xqtl_lbf_matrix, combined_gwas_lbf_matrix, p1 = p1, p2 = p2, p12 = p12)
    } else {
      coloc_res <- list("No coloc results due to the absence of a GWAS log Bayes factor matrix filtered by prior tolerance.")
    }
  } else {
    coloc_res <- list(paste("no", xqtl_finemapping_obj[2], "in", xqtl_finemapping_obj[1]))
  }
  return(c(coloc_res, analysis_region = region))
}

#' coloc_post_processor function
#' @param coloc_res coloc results from coloc.susie.
#' @param LD_meta_file_path Path to the metadata of LD reference.
#' @param analysis_region Path to the analysis region of coloc result.
#' @return A list containing the coloc results and post processed coloc sets.
#' @export
coloc_post_processor <- function(coloc_res, LD_meta_file_path = NULL, analysis_region = NULL, ...) {
  if (!is.null(LD_meta_file_path)) {
    if (is.null(analysis_region)) {
      stop("LD_meta_file_path is provided but analysis_region is not provided. Please provide analysis_region for purity filter.")
    }
    # Perform purity filter using LD_meta_file_path and analysis_region
    coloc_res <- c(coloc_res, process_coloc_results(coloc_res, LD_meta_file_path, analysis_region = analysis_region))
  } else {
    if (!is.null(analysis_region)) {
      warning("Analysis_region is provided but will not be used as LD_meta_file_path is not provided.")
    }
    warning("LD_meta_file_path not provided. Purity filter cannot be applied.")
  }
  return(coloc_res)
}

# In practice, analysis will contain two lines:
# res <- coloc_wrapper(...)
# post_processed_res <- coloc_post_processor
