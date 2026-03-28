#' Function to Check if Regions are in increasing order and remove duplicated rows
#' @importFrom dplyr distinct arrange group_by mutate ungroup
#' @importFrom magrittr %>%
#' @noRd
order_dedup_regions <- function(df) {
  # Ensure that 'chrom' values are integers, df can be genomic_data or regions_of_interest
  df$chrom <- as.integer(strip_chr_prefix(df$chrom))

  # Remove duplicated rows based on 'chrom' and 'start' columns
  df <- distinct(df, chrom, start, .keep_all = TRUE) %>%
    arrange(chrom, start)

  # for (chr in unique(df$chrom)) {
  #  chr_rows <- which(df$chrom == chr)
  #  if (length(chr_rows) > 1) {
  #    starts <- df$start[chr_rows]
  #    for (i in 2:length(starts)) {
  #      if (starts[i] < starts[i - 1]) {
  #        stop("The input list of regions is not in increasing order within each chromosome.")
  #      }
  #    }
  #  }
  # }

  return(df)
}

#' Function to Find Start and End Rows of Genomic Data for Region of Interest
#' @importFrom dplyr filter arrange slice
#' @noRd
find_intersection_rows <- function(genomic_data, region_chrom, region_start, region_end) {
  # Filter for the specific chromosome
  chrom_data <- genomic_data %>% filter(chrom == region_chrom)

  min_start <- if (nrow(chrom_data) > 0) min(chrom_data$start) else NA
  max_end <- if (nrow(chrom_data) > 0) max(chrom_data$end) else NA

  # Adjust region bounds if they're outside available data range
  if (!is.na(min_start) && region_start < min_start) {
    region_start <- min_start
  }

  if (!is.na(max_end) && region_end > max_end) {
    region_end <- max_end
  }

  # Try to find rows that cover the region start and end
  start_row <- genomic_data %>%
    filter(chrom == region_chrom, start <= region_start, end > region_start) %>%
    slice(1)

  end_row <- genomic_data %>%
    filter(chrom == region_chrom, start < region_end, end >= region_end) %>%
    arrange(desc(end)) %>%
    slice(1)

  if (nrow(start_row) == 0 || nrow(end_row) == 0) {
    stop("Region of interest is not covered by any rows in the data frame.")
  }

  list(start_row = start_row, end_row = end_row)
}

#' Function to Validate Selected Region
#' @noRd
validate_selected_region <- function(start_row, end_row, region_start, region_end) {
  if (!(start_row$start <= region_start && end_row$end >= region_end)) {
    stop("The selected region is not fully covered by the merged region.")
  }
}

#' Extract File Paths Based on Intersection Criteria
#' @noRd
extract_file_paths <- function(genomic_data, intersection_rows, column_to_extract) {
  # Ensure the file_path_column exists in genomic_data
  if (!column_to_extract %in% names(genomic_data)) {
    stop(paste("Column", column_to_extract, "not found in genomic data"))
  }
  # Extract rows based on intersection criteria
  extracted_paths <- genomic_data[[column_to_extract]][which(
    genomic_data$chrom == intersection_rows$start_row$chrom &
      genomic_data$start >= intersection_rows$start_row$start &
      genomic_data$start <= intersection_rows$end_row$start
  )]

  # Return the extracted paths
  return(extracted_paths)
}

#' Intersect LD reference with Regions of Interest
#'
#' @param ld_reference_meta_file A file of data frame with columns "chrom", "start", "end", and "path" representing genomic regions.
#' "chrom" is the chromosome, "start" and "end" are the positions of the LD block, and "path" is the file path for the LD block.
#' @param region A data frame with columns "chrom", "start", and "end" specifying regions of interest.
#' "start" and "end" are the positions of these regions. Or it can take the form of `chr:start-end`
#'
#' @return A list containing processed genomic data, region of interest, and a detailed result list.
#' The result list contains, for each region:
#' - The start and end row indices in the genomic data
#' - File paths from the genomic data corresponding to intersected regions
#' - Optionally, bim file paths if available
#' @importFrom stringr str_split
#' @importFrom dplyr select
#' @importFrom vroom vroom
#' @noRd
get_regional_ld_meta <- function(ld_reference_meta_file, region, complete_coverage_required = FALSE) {
  genomic_data <- vroom(ld_reference_meta_file)
  region <- parse_region(region)
  # Set column names
  names(genomic_data) <- c("chrom", "start", "end", "path")
  names(region) <- c("chrom", "start", "end")

  # Order and deduplicate regions
  genomic_data <- order_dedup_regions(genomic_data)
  region <- order_dedup_regions(region)

  # Process file paths
  file_path <- genomic_data$path %>%
    str_split(",", simplify = TRUE) %>%
    data.frame() %>%
    `colnames<-`(if (ncol(.) == 2) c("LD_file_path", "bim_file_path") else c("LD_file_path"))

  genomic_data <- cbind(genomic_data, file_path) %>% select(-path)

  # Find intersection rows
  intersection_rows <- find_intersection_rows(genomic_data, region$chrom, region$start, region$end)

  # Validate region
  if (complete_coverage_required) {
    validate_selected_region(intersection_rows$start_row, intersection_rows$end_row, region$start, region$end)
  }

  # Extract file paths
  LD_paths <- find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "LD_file_path"))
  bim_paths <- if ("bim_file_path" %in% names(genomic_data)) {
    find_valid_file_paths(ld_reference_meta_file, extract_file_paths(genomic_data, intersection_rows, "bim_file_path"))
  } else {
    NULL
  }

  return(list(
    intersections = list(
      start_index = intersection_rows$start_row,
      end_index = intersection_rows$end_row,
      LD_file_paths = LD_paths,
      bim_file_paths = bim_paths
    ),
    ld_meta_data = genomic_data,
    region = region
  ))
}

#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom stats setNames
# Process an LD matrix from a file path
process_LD_matrix <- function(LD_file_path, bim_file_path) {
  # Order the variants by position
  order_variants_by_position <- function(strings) {
    # Apply the function to each variant to get a vector of positions
    positions <- sapply(strings, function(variant) as.integer(strsplit(variant, ":")[[1]][2]))
    # Check whether the merged variants is orderd
    # when diff() returns 0, it is multiallelic at the position (same position but >1 variations)
    if (!all(diff(positions[order(positions)]) >= 0)) {
      stop("The positions are not in non-decreasing order")
    }
    # Order the variants by position
    strings_ordered <- strings[order(positions)]
    return(strings_ordered)
  }
  # Read the LD matrix
  LD_file_con <- xzfile(LD_file_path)
  LD_matrix <- scan(LD_file_con, quiet = TRUE)
  close(LD_file_con)
  LD_matrix <- matrix(LD_matrix, ncol = sqrt(length(LD_matrix)), byrow = TRUE)

  # Process the bim file to extract variant information
  bim_file_name <- if (!is.null(bim_file_path)) {
    bim_file_path
  } else {
    paste0(LD_file_path, ".bim", sep = "")
  }

  # Process variant names from file paths
  LD_variants <- read.table(bim_file_name)
  if (ncol(LD_variants) == 9) {
    LD_variants <- LD_variants %>%
      setNames(c("chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss")) %>%
      mutate(chrom = as.character(as.integer(strip_chr_prefix(chrom)))) %>%
      mutate(variants = normalize_variant_id(variants))
  } else if (ncol(LD_variants) == 6) {
    LD_variants <- LD_variants %>%
      setNames(c("chrom", "variants", "GD", "pos", "A1", "A2")) %>%
      mutate(chrom = as.character(as.integer(strip_chr_prefix(chrom)))) %>%
      mutate(variants = normalize_variant_id(variants))
  } else {
    stop("Unexpected number of columns in the input file.")
  }


  # Set column and row names of the LD matrix
  colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants
  # Check if the matrix is upper diagonal
  # We assume a matrix is upper diagonal if all elements below the main diagonal are zero
  is_upper_diagonal <- all(LD_matrix[lower.tri(LD_matrix)] == 0)
  if (is_upper_diagonal) {
    # If the matrix is upper diagonal, transpose the upper triangle to the lower triangle
    LD_matrix[lower.tri(LD_matrix)] <- t(LD_matrix)[lower.tri(LD_matrix)]
  } else {
    # If the matrix is lower diagonal, transpose the lower triangle to the upper triangle
    LD_matrix[upper.tri(LD_matrix)] <- t(LD_matrix)[upper.tri(LD_matrix)]
  }
  LD_variants_ordered <- LD_variants[match(order_variants_by_position(LD_variants$variants), LD_variants$variants), ]
  LD_matrix <- LD_matrix[match(LD_variants_ordered$variants, rownames(LD_matrix)), match(LD_variants_ordered$variants, rownames(LD_matrix))]
  list(LD_matrix = LD_matrix, LD_variants = LD_variants_ordered)
}

#' Extract LD matrix and variants for a specific region
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @importFrom utils tail
extract_LD_for_region <- function(LD_matrix, variants, region, extract_coordinates) {
  # Filter variants based on region
  extracted_LD_variants <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)
  if (!is.null(extract_coordinates)) {
    # Preprocess 'extract_coordinate' to ensure 'chrom' is numeric and without 'chr'
    extract_coordinates <- extract_coordinates %>%
      mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
      select(chrom, pos)
    # Now merge with 'LD_variants_region_selected'
    extracted_LD_variants <- extracted_LD_variants %>%
      mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
      # Merge with 'extract_coordinate' after 'chrom' adjustment
      merge(extract_coordinates, by = c("chrom", "pos"))
    # Select the desired columns, assuming 'variants' column is equivalent to the 'variants' in 'LD_variants_region_selected'
    # Select columns dynamically based on the presence of 'variance'
    cols_to_select <- c("chrom", "variants", "pos", "GD", "A1", "A2") # select(chrom, variants, pos, GD, A1, A2)
    if ("variance" %in% names(extracted_LD_variants)) {
      cols_to_select <- c(cols_to_select, "variance")
    }
    extracted_LD_variants <- select(extracted_LD_variants, all_of(cols_to_select))
  }
  # Extract LD matrix
  extracted_LD_matrix <- LD_matrix[extracted_LD_variants$variants, extracted_LD_variants$variants, drop = FALSE]
  list(extracted_LD_matrix = extracted_LD_matrix, extracted_LD_variants = extracted_LD_variants)
}

# Create a combined LD matrix from multiple matrices
create_LD_matrix <- function(LD_matrices, variants) {
  # Updated mergeVariants helper that checks the structure of each element.
  mergeVariants <- function(LD_variants_list) {
    mergedVariants <- character(0)
    # Loop over the list of variant information
    for (LD_variants in LD_variants_list) {
      currentVariants <- if (is.list(LD_variants) && !is.null(LD_variants$variants)) {
        LD_variants$variants
      } else {
        LD_variants
      }

      if (length(currentVariants) == 0) next

      # If the last variant in mergedVariants is the same as the first of currentVariants, skip the duplicate
      if (length(mergedVariants) > 0 && tail(mergedVariants, 1) == currentVariants[1]) {
        mergedVariants <- c(mergedVariants, currentVariants[-1])
      } else {
        mergedVariants <- c(mergedVariants, currentVariants)
      }
    }
    return(mergedVariants)
  }

  unique_variants <- mergeVariants(variants)
  # Initialize an empty combined LD matrix with the unique variants
  LD_matrix <- matrix(0, nrow = length(unique_variants), ncol = length(unique_variants))
  rownames(LD_matrix) <- unique_variants
  colnames(LD_matrix) <- unique_variants

  # Function to align the values from each LD matrix to the combined matrix
  align_matrix <- function(ld_matrix, combined_matrix, variant_names) {
    indices <- match(variant_names, rownames(combined_matrix))
    combined_matrix[indices, indices] <- ld_matrix
    return(combined_matrix)
  }

  # Use Map to pair each LD matrix with its rownames, then Reduce to fill the combined matrix.
  LD_matrix <- Reduce(
    function(x, y) align_matrix(y[[1]], x, y[[2]]),
    Map(list, LD_matrices, lapply(LD_matrices, rownames)),
    LD_matrix
  )

  return(LD_matrix)
}

#' Load and Process Linkage Disequilibrium (LD) Matrix
#'
#' Unified entry point for loading LD data. Auto-detects the source type:
#' \itemize{
#'   \item Pre-computed LD blocks (metadata file with chrom/start/end/path columns)
#'   \item PLINK2 genotype files (.pgen/.pvar[.zst]/.psam + .afreq)
#'   \item PLINK1 genotype files (.bed/.bim/.fam)
#' }
#' For PLINK genotype sources, LD is computed on the fly via \code{compute_LD()}.
#'
#' @param LD_meta_file_path Path to LD metadata file, or prefix for PLINK genotype files.
#' @param region Region of interest: "chr:start-end" string or data.frame with chrom/start/end.
#' @param extract_coordinates Optional data.frame with columns "chrom" and "pos" for
#'   specific coordinates extraction (only for pre-computed LD blocks).
#' @param return_genotype If TRUE, return the genotype matrix X instead of the LD
#'   correlation matrix R. Only valid for PLINK genotype sources.
#' @param n_sample Optional sample size for computing variance (= 2*p*(1-p)*n/(n-1)).
#'   If NULL, ref_panel will not include variance or n_nomiss columns.
#'   Only used for PLINK genotype sources.
#'
#' @return A list with:
#' \describe{
#'   \item{LD_variants}{Character vector of variant IDs (canonical format).}
#'   \item{LD_matrix}{Symmetric LD correlation matrix (or genotype matrix if return_genotype=TRUE).}
#'   \item{ref_panel}{Data.frame with variant metadata (chrom, pos, A2, A1, variant_id,
#'     and optionally allele_freq, variance, n_nomiss).}
#'   \item{block_metadata}{Data.frame with block info (block_id, chrom, block_start, block_end, size, start_idx, end_idx).}
#' }
#' @export
load_LD_matrix <- function(LD_meta_file_path, region, extract_coordinates = NULL,
                           return_genotype = FALSE, n_sample = NULL) {
  source <- resolve_ld_source(LD_meta_file_path)

  if (source$type %in% c("plink2", "plink1")) {
    return(load_LD_from_genotype(source$data_path, region, source$type,
                                 return_genotype = return_genotype,
                                 n_sample = n_sample))
  }

  # Pre-computed LD blocks (.cor.xz)
  if (return_genotype) {
    stop("return_genotype=TRUE requires PLINK genotype files, not pre-computed LD matrices.")
  }
  load_LD_from_blocks(source$meta_path, region, extract_coordinates, n_sample = n_sample)
}

# ---------- Internal: resolve LD source type ----------

#' @noRd
has_plink2_files <- function(prefix) {
  file.exists(paste0(prefix, ".pgen")) &&
    (file.exists(paste0(prefix, ".pvar")) || file.exists(paste0(prefix, ".pvar.zst"))) &&
    file.exists(paste0(prefix, ".psam"))
}

#' @noRd
has_plink1_files <- function(prefix) {
  file.exists(paste0(prefix, ".bed")) &&
    file.exists(paste0(prefix, ".bim")) &&
    file.exists(paste0(prefix, ".fam"))
}

#' Resolve an LD source path to its actual data type.
#'
#' Handles both direct PLINK prefixes and metadata TSV files (which may
#' themselves point to PLINK files or pre-computed .cor.xz matrices).
#'
#' @param path Direct PLINK prefix or metadata TSV file path.
#' @return A list with:
#'   \item{type}{"plink2", "plink1", or "precomputed"}
#'   \item{data_path}{Resolved PLINK prefix (for plink types)}
#'   \item{meta_path}{Metadata TSV path (for precomputed, or when input was a TSV)}
#' @importFrom vroom vroom
#' @noRd
resolve_ld_source <- function(path) {
  # Direct PLINK prefix?
  if (has_plink2_files(path)) return(list(type = "plink2", data_path = path))
  if (has_plink1_files(path)) return(list(type = "plink1", data_path = path))

  # Must be a metadata file
  if (!file.exists(path)) {
    stop("Cannot determine LD source from path: ", path,
         "\n  Expected: PLINK2 prefix (.pgen/.pvar[.zst]/.psam), ",
         "PLINK1 prefix (.bed/.bim/.fam), or LD metadata file.")
  }

  # Peek at first row's path column to determine underlying data type
  meta <- as.data.frame(vroom(path, show_col_types = FALSE, n_max = 1))
  colnames(meta)[1] <- "chrom"
  raw_path <- gsub(",.*$", "", meta$path[1])  # strip comma-separated bim path
  resolved <- file.path(dirname(path), raw_path)

  if (has_plink2_files(resolved)) return(list(type = "plink2", data_path = resolved, meta_path = path))
  if (has_plink1_files(resolved)) return(list(type = "plink1", data_path = resolved, meta_path = path))

  # Pre-computed .cor.xz blocks
  list(type = "precomputed", meta_path = path)
}

# ---------- Internal: load LD from genotype files ----------

#' Load genotype data from PLINK files and compute LD or return genotype matrix.
#' @param source_type Character, "plink1" or "plink2" (from resolve_ld_source).
#' @noRd
load_LD_from_genotype <- function(prefix, region, source_type,
                                  return_genotype = FALSE, n_sample = NULL) {
  # Load genotype matrix and variant info
  result <- if (source_type == "plink2") {
    load_plink2_data(prefix, region = region)
  } else {
    load_plink1_data(prefix, region = region)
  }
  X <- result$X
  variant_info <- result$variant_info

  # Normalize variant IDs to canonical format (chr:pos:A2:A1)
  variant_ids <- normalize_variant_id(
    format_variant_id(variant_info$chrom, variant_info$pos, variant_info$A2, variant_info$A1)
  )
  colnames(X) <- variant_ids

  # Build ref_panel
  ref_panel <- parse_variant_id(variant_ids)
  ref_panel$variant_id <- variant_ids

  # Load allele frequency from .afreq file (required for PLINK sources)
  afreq <- read_afreq(prefix)
  if (is.null(afreq)) {
    stop("Allele frequency file (.afreq or .afreq.zst) not found at prefix: ", prefix,
         "\n  The .afreq file is required for PLINK genotype LD sources.")
  }
  freq_match <- match(variant_info$id, afreq$id)
  n_unmatched <- sum(is.na(freq_match))
  if (n_unmatched > 0) {
    warning(n_unmatched, " out of ", length(freq_match),
            " variants have no allele frequency in .afreq file.")
  }
  ref_panel$allele_freq <- afreq$alt_freq[freq_match]

  # Compute variance if sample size provided
  if (!is.null(n_sample)) {
    p <- ref_panel$allele_freq
    ref_panel$variance <- 2 * p * (1 - p) * n_sample / (n_sample - 1)
    ref_panel$n_nomiss <- n_sample
  }

  # Block metadata (single block spanning the loaded region)
  positions <- variant_info$pos
  block_metadata <- data.frame(
    block_id = 1L,
    chrom = as.character(variant_info$chrom[1]),
    block_start = min(positions),
    block_end = max(positions),
    size = length(variant_ids),
    start_idx = 1L,
    end_idx = length(variant_ids),
    stringsAsFactors = FALSE
  )

  if (return_genotype) {
    # Note: LD_matrix holds the genotype matrix X (not LD) when return_genotype=TRUE
    return(list(
      LD_variants = variant_ids,
      LD_matrix = X,
      ref_panel = ref_panel,
      block_metadata = block_metadata
    ))
  }

  # Compute LD correlation matrix
  R <- compute_LD(X, method = "sample")

  list(
    LD_variants = variant_ids,
    LD_matrix = R,
    ref_panel = ref_panel,
    block_metadata = block_metadata
  )
}

# ---------- Internal: load LD from pre-computed blocks ----------

#' Load pre-computed LD from block-based metadata files.
#' @noRd
load_LD_from_blocks <- function(LD_meta_file_path, region, extract_coordinates = NULL, n_sample = NULL) {
  # Intersect LD metadata with specified regions
  intersected_LD_files <- get_regional_ld_meta(LD_meta_file_path, region)

  LD_file_paths <- intersected_LD_files$intersections$LD_file_paths
  bim_file_paths <- intersected_LD_files$intersections$bim_file_paths

  extracted_LD_matrices_list <- list()
  extracted_LD_variants_list <- list()
  block_chroms <- character(length(LD_file_paths))

  for (j in seq_along(LD_file_paths)) {
    LD_matrix_processed <- process_LD_matrix(LD_file_paths[j], bim_file_paths[j])
    extracted_LD_list <- extract_LD_for_region(
      LD_matrix = LD_matrix_processed$LD_matrix,
      variants = LD_matrix_processed$LD_variants,
      region = intersected_LD_files$region,
      extract_coordinates = extract_coordinates
    )
    extracted_LD_matrices_list[[j]] <- extracted_LD_list$extracted_LD_matrix
    extracted_LD_variants_list[[j]] <- extracted_LD_list$extracted_LD_variants
    if (nrow(extracted_LD_variants_list[[j]]) > 0) {
      block_chroms[j] <- as.character(extracted_LD_variants_list[[j]]$chrom[1])
    } else {
      block_chroms[j] <- as.character(intersected_LD_files$region$chrom)
    }
  }

  LD_matrix <- create_LD_matrix(
    LD_matrices = extracted_LD_matrices_list,
    variants = extracted_LD_variants_list
  )
  LD_variants <- rownames(LD_matrix)

  block_variants <- lapply(extracted_LD_variants_list, function(v) v$variants)
  block_positions <- lapply(extracted_LD_variants_list, function(v) v$pos)
  block_metadata <- data.frame(
    block_id = seq_along(LD_file_paths),
    chrom = block_chroms,
    block_start = sapply(block_positions, min),
    block_end = sapply(block_positions, max),
    size = sapply(block_variants, length),
    start_idx = sapply(block_variants, function(v) min(match(v, LD_variants))),
    end_idx = sapply(block_variants, function(v) max(match(v, LD_variants))),
    stringsAsFactors = FALSE
  )

  rm(extracted_LD_matrices_list)

  ref_panel <- parse_variant_id(rownames(LD_matrix))
  merged_variant_list <- do.call(rbind, extracted_LD_variants_list)
  ref_panel$variant_id <- rownames(LD_matrix)

  if ("allele_freq" %in% colnames(merged_variant_list)) {
    ref_panel$allele_freq <- merged_variant_list$allele_freq[match(rownames(LD_matrix), merged_variant_list$variants)]
  }
  if ("variance" %in% colnames(merged_variant_list)) {
    ref_panel$variance <- merged_variant_list$variance[match(rownames(LD_matrix), merged_variant_list$variants)]
  }
  if ("n_nomiss" %in% colnames(merged_variant_list)) {
    ref_panel$n_nomiss <- merged_variant_list$n_nomiss[match(rownames(LD_matrix), merged_variant_list$variants)]
  }

  # Compute variance from n_sample + allele_freq if not already present
  if (!is.null(n_sample) && (!"variance" %in% colnames(ref_panel) || all(is.na(ref_panel$variance)))) {
    if ("allele_freq" %in% colnames(ref_panel)) {
      p <- ref_panel$allele_freq
      ref_panel$variance <- 2 * p * (1 - p) * n_sample / (n_sample - 1)
      ref_panel$n_nomiss <- n_sample
    }
  }

  list(
    LD_variants = LD_variants,
    LD_matrix = LD_matrix,
    ref_panel = ref_panel,
    block_metadata = block_metadata
  )
}

#' Filter variants by LD Reference
#'
#' Filters a vector of variant IDs to those present in the LD reference panel.
#' Auto-detects the reference type (PLINK2, PLINK1, or pre-computed LD metadata).
#'
#' @param variant_ids variant names in the format chr:pos:ref:alt.
#' @param ld_reference_meta_file Path to LD metadata file or PLINK prefix.
#' @param keep_indel Whether to keep indel variants. Default TRUE.
#' @return A list with:
#'   \item{data}{Character vector of filtered variant IDs.}
#'   \item{idx}{Integer vector of indices into the original variant_ids.}
#' @importFrom dplyr group_by summarise
#' @importFrom vroom vroom
#' @importFrom magrittr %>%
#' @export
filter_variants_by_ld_reference <- function(variant_ids, ld_reference_meta_file, keep_indel = TRUE) {
  variants_df <- parse_variant_id(variant_ids)

  # Derive region to scope the reference lookup
  region_df <- variants_df %>%
    group_by(chrom) %>%
    summarise(start = min(pos), end = max(pos))

  # Use shared helper — no genotype loading
  ref_info <- get_ref_variant_info(ld_reference_meta_file, region_df)
  ref_chrom <- as.integer(strip_chr_prefix(ref_info$chrom))
  ref_key <- paste0(ref_chrom, ":", ref_info$pos)

  variant_key <- paste0(variants_df$chrom, ":", variants_df$pos)
  keep_indices <- which(variant_key %in% ref_key)

  if (!keep_indel) {
    snp_idx <- which(is_snp_alleles(variants_df$A1, variants_df$A2))
    keep_indices <- intersect(keep_indices, snp_idx)
  }

  message(length(variant_ids) - length(keep_indices), " out of ", length(variant_ids),
          " total variants dropped due to absence on the reference LD panel.")

  list(data = variant_ids[keep_indices], idx = keep_indices)
}

#' Partition LD Matrix into Block-Specific Matrices
#'
#' This function takes the output from load_LD_matrix and partitions the combined LD matrix
#' into a list of smaller matrices based on the block_indices, making it easier to work with
#' large LD matrices that span multiple blocks.
#'
#' @param ld_data A list as returned by load_LD_matrix, containing LD_matrix,
#'                LD_variants, ref_panel, and block_metadata.
#' @param merge_small_blocks Logical, whether to merge blocks smaller than min_merged_block_size (default: TRUE).
#' @param min_merged_block_size Integer, minimum number of variants for a block after merging (default: 500).
#' @param max_merged_block_size Integer, maximum number of variants in a block after merging (default: 10000).
#'
#' @return returns a list containing:
#' \describe{
#' \item{ld_matrices}{A list of matrices, each representing LD for a specific block.}
#' \item{variant_indices}{A data frame that maps variant IDs to their corresponding block.}
#' \item{block_metadata}{Information about each block including size, chromosome, start and end positions.}
#' }
#' @noRd
partition_LD_matrix <- function(ld_data, merge_small_blocks = TRUE,
                                min_merged_block_size = 500, max_merged_block_size = 10000) {
  # Extract components from ld_data
  combined_matrix <- ld_data$LD_matrix
  block_metadata <- ld_data$block_metadata
  variant_ids <- ld_data$LD_variants

  # Error if matrix is empty
  if (is.null(combined_matrix) || nrow(combined_matrix) == 0 || ncol(combined_matrix) == 0) {
    stop("Empty or NULL LD matrix provided.")
  }

  # Ensure the row and column names of the matrix match the variant_ids
  if (is.null(rownames(combined_matrix)) || is.null(colnames(combined_matrix)) ||
    !identical(rownames(combined_matrix), variant_ids) || !identical(colnames(combined_matrix), variant_ids)) {
    rownames(combined_matrix) <- variant_ids
    colnames(combined_matrix) <- variant_ids
  }

  # Validate the block structure of the matrix (skip if only one block)
  if (nrow(block_metadata) > 1) {
    validate_block_structure(combined_matrix, block_metadata, variant_ids)
  }

  # Optionally merge small blocks
  if (merge_small_blocks && any(block_metadata$size < min_merged_block_size) && nrow(block_metadata) > 1) {
    block_metadata <- merge_blocks(block_metadata, min_merged_block_size, max_merged_block_size)
  }

  # Partition the matrix based on block metadata
  result <- extract_block_matrices(combined_matrix, block_metadata, variant_ids)
  return(result)
}

# Helper function to validate block structure
validate_block_structure <- function(matrix, block_metadata, variant_ids) {
  validation_failed <- FALSE
  validation_message <- character(0)

  # Iterate over all pairs of distinct blocks
  for (i in 1:(nrow(block_metadata) - 1)) {
    for (j in (i + 1):nrow(block_metadata)) {
      # Get indices for each block
      block_i_start <- block_metadata$start_idx[i]
      block_i_end <- block_metadata$end_idx[i]
      block_j_start <- block_metadata$start_idx[j]
      block_j_end <- block_metadata$end_idx[j]

      # Check if indices are valid
      if (block_i_start > length(variant_ids) || block_i_end > length(variant_ids) ||
        block_j_start > length(variant_ids) || block_j_end > length(variant_ids)) {
        validation_failed <- TRUE
        validation_message <- c(
          validation_message,
          paste("Block indices are out of range for blocks", i, "and", j)
        )
        next
      }

      # Get variants for each block, EXCLUDING THE BOUNDARY VARIANTS
      # For block i, exclude the last variant (potential overlap with next block)
      # For block j, exclude the first variant (potential overlap with previous block)
      block_i_vars <- variant_ids[block_i_start:(block_i_end - 1)]
      block_j_vars <- variant_ids[(block_j_start + 1):block_j_end]

      # Only check if there are variants in both blocks
      if (length(block_i_vars) > 0 && length(block_j_vars) > 0) {
        # Extract the cross-block submatrix for these variants
        cross_block_unique <- matrix[block_i_vars, block_j_vars, drop = FALSE]

        # Check if any values are non-zero
        max_value <- max(abs(cross_block_unique))
        if (max_value > 1e-10) {
          validation_failed <- TRUE
          validation_message <- c(
            validation_message,
            paste(
              "Non-zero correlation detected between unique variants of blocks", i, "and", j,
              "- Max value:", max_value
            )
          )
        }
      }
    }
  }

  # If validation failed, stop execution with error message
  if (validation_failed) {
    stop(
      "Matrix does not have the expected block structure:\n",
      paste(validation_message, collapse = "\n")
    )
  }
}
# Helper function to merge small blocks
# Check if two blocks can be merged based on chromosome and size constraints
can_merge <- function(block1, block2, max_size) {
  # Check if blocks are on the same chromosome
  same_chrom <- block1$chrom == block2$chrom

  # Check if combined size is within limits
  combined_size <- block1$size + block2$size
  size_ok <- combined_size <= max_size

  return(same_chrom && size_ok)
}

# Merge two blocks in the metadata dataframe
merge_two_blocks <- function(block_metadata, idx1, idx2) {
  # Ensure idx1 < idx2 for consistent behavior
  if (idx1 > idx2) {
    temp <- idx1
    idx1 <- idx2
    idx2 <- temp
  }

  # Create a copy to avoid modifying the input directly
  result <- block_metadata

  # Update the end index and size of the first block
  result$end_idx[idx1] <- block_metadata$end_idx[idx2]
  result$size[idx1] <- block_metadata$size[idx1] + block_metadata$size[idx2]

  # Remove the second block
  result <- result[-idx2, ]

  # Renumber block IDs
  result$block_id <- seq_len(nrow(result))

  return(result)
}

# Find small blocks and their best merge candidates
find_merge_candidates <- function(block_metadata, min_size, max_size) {
  candidates <- data.frame(
    block_idx = integer(),
    merge_with = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(block_metadata))) {
    # Skip if block size is adequate
    if (block_metadata$size[i] >= min_size) next

    # Check previous block
    prev_idx <- i - 1
    can_merge_prev <- prev_idx >= 1 &&
      can_merge(block_metadata[i, ], block_metadata[prev_idx, ], max_size)

    # Check next block
    next_idx <- i + 1
    can_merge_next <- next_idx <= nrow(block_metadata) &&
      can_merge(block_metadata[i, ], block_metadata[next_idx, ], max_size)

    # Determine best merge option
    if (can_merge_prev && can_merge_next) {
      # Choose the smaller of the two to minimize impact
      if (block_metadata$size[prev_idx] <= block_metadata$size[next_idx]) {
        merge_with <- prev_idx
      } else {
        merge_with <- next_idx
      }
      candidates <- rbind(candidates, data.frame(block_idx = i, merge_with = merge_with))
    } else if (can_merge_prev) {
      candidates <- rbind(candidates, data.frame(block_idx = i, merge_with = prev_idx))
    } else if (can_merge_next) {
      candidates <- rbind(candidates, data.frame(block_idx = i, merge_with = next_idx))
    }
    # If no valid merge candidates, leave the block as is
  }

  return(candidates)
}

# Helper function to merge small blocks
merge_blocks <- function(block_metadata, min_size, max_size) {
  # If there's only one block or empty input, just return it
  if (nrow(block_metadata) <= 1) {
    return(block_metadata)
  }

  new_block_metadata <- block_metadata
  made_changes <- TRUE

  # Iterate until no more merges are possible
  while (made_changes) {
    # Find all current merge candidates
    candidates <- find_merge_candidates(new_block_metadata, min_size, max_size)

    # Stop if no candidates found
    if (nrow(candidates) == 0) {
      made_changes <- FALSE
      break
    }

    # Process the first candidate (we only do one merge per iteration to avoid index issues)
    new_block_metadata <- merge_two_blocks(
      new_block_metadata,
      candidates$block_idx[1],
      candidates$merge_with[1]
    )
  }

  return(new_block_metadata)
}

# Helper function to extract block matrices
extract_block_matrices <- function(matrix, block_metadata, variant_ids) {
  ld_matrices <- list()
  variant_mapping <- data.frame(
    variant_id = character(),
    block_id = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(block_metadata))) {
    start_idx <- block_metadata$start_idx[i]
    end_idx <- block_metadata$end_idx[i]

    # Skip empty blocks
    if (end_idx < start_idx) next

    # Ensure indices are within bounds
    if (start_idx > length(variant_ids) || end_idx > length(variant_ids)) {
      warning(paste("Block", i, "has indices outside the range of variant_ids. Skipping."))
      next
    }

    # Extract variant IDs for this block
    block_variants <- variant_ids[start_idx:end_idx]

    # Extract submatrix for this block - use named indexing
    block_matrix <- matrix[block_variants, block_variants, drop = FALSE]

    # Store in list
    ld_matrices[[i]] <- block_matrix

    # Update variant mapping
    block_mapping <- data.frame(
      variant_id = block_variants,
      block_id = i,
      stringsAsFactors = FALSE
    )
    variant_mapping <- rbind(variant_mapping, block_mapping)

  }

  return(list(
    ld_matrices = ld_matrices,
    variant_indices = variant_mapping,
    block_metadata = block_metadata
  ))
}
