#' Deduplicate and sort genomic regions by chromosome and start position.
#' @importFrom dplyr distinct arrange
#' @importFrom magrittr %>%
#' @noRd
order_dedup_regions <- function(df) {
  df$chrom <- as.integer(strip_chr_prefix(df$chrom))
  df <- distinct(df, chrom, start, .keep_all = TRUE) %>%
    arrange(chrom, start)
  df
}

#' Find the first and last rows of genomic_data that overlap a query region.
#' Clamps the query to the available data range before searching.
#' @importFrom dplyr filter arrange slice
#' @noRd
find_intersection_rows <- function(genomic_data, region_chrom, region_start, region_end) {
  chrom_data <- genomic_data %>% filter(chrom == region_chrom)
  if (nrow(chrom_data) == 0) stop("No data for chromosome ", region_chrom)

  # Clamp query to available range
  region_start <- max(region_start, min(chrom_data$start))
  region_end   <- min(region_end,   max(chrom_data$end))

  start_row <- genomic_data %>%
    filter(chrom == region_chrom, start <= region_start, end > region_start) %>%
    slice(1)
  end_row <- genomic_data %>%
    filter(chrom == region_chrom, start < region_end, end >= region_end) %>%
    arrange(desc(end)) %>%
    slice(1)

  if (nrow(start_row) == 0 || nrow(end_row) == 0) {
    stop("Region ", region_chrom, ":", region_start, "-", region_end,
         " is not covered by any rows in the LD metadata.")
  }
  list(start_row = start_row, end_row = end_row)
}

#' Validate that start_row..end_row fully covers [region_start, region_end].
#' @noRd
validate_selected_region <- function(start_row, end_row, region_start, region_end) {
  if (start_row$start > region_start || end_row$end < region_end) {
    stop("Region ", region_start, "-", region_end, " is not fully covered by the LD metadata ",
         "(available: ", start_row$start, "-", end_row$end, ").")
  }
}

#' Extract values of a column for rows spanning the intersection range.
#' @noRd
extract_file_paths <- function(genomic_data, intersection_rows, column_to_extract) {
  if (!column_to_extract %in% names(genomic_data)) {
    stop("Column '", column_to_extract, "' not found in genomic data.")
  }
  idx <- which(genomic_data$chrom == intersection_rows$start_row$chrom &
               genomic_data$start >= intersection_rows$start_row$start &
               genomic_data$start <= intersection_rows$end_row$start)
  genomic_data[[column_to_extract]][idx]
}

#' Find LD blocks overlapping a query region from a metadata TSV file.
#'
#' @param ld_reference_meta_file TSV with columns chrom, start, end, path.
#'   The path column may be comma-separated: "ld_file,bim_file".
#' @param region "chr:start-end" string or data.frame with chrom/start/end.
#' @param complete_coverage_required If TRUE, error when the region extends
#'   beyond available LD blocks.
#' @return A list with: intersections (LD_file_paths, bim_file_paths),
#'   ld_meta_data, and parsed region.
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

#' Read a pre-computed LD matrix (.cor.xz) and its bim file, returning a
#' symmetric matrix with variants ordered by position.
#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom stats setNames
#' @noRd
process_LD_matrix <- function(LD_file_path, bim_file_path) {
  # Read .cor.xz matrix
  LD_file_con <- xzfile(LD_file_path)
  LD_matrix <- scan(LD_file_con, quiet = TRUE)
  close(LD_file_con)
  LD_matrix <- matrix(LD_matrix, ncol = sqrt(length(LD_matrix)), byrow = TRUE)

  # Read bim (variant metadata)
  bim_file_name <- if (!is.null(bim_file_path)) bim_file_path else paste0(LD_file_path, ".bim")
  LD_variants <- read.table(bim_file_name)
  col_names <- if (ncol(LD_variants) == 9) {
    c("chrom", "variants", "GD", "pos", "A1", "A2", "variance", "allele_freq", "n_nomiss")
  } else if (ncol(LD_variants) == 6) {
    c("chrom", "variants", "GD", "pos", "A1", "A2")
  } else {
    stop("Unexpected number of columns (", ncol(LD_variants), ") in bim file: ", bim_file_name)
  }
  LD_variants <- LD_variants %>%
    setNames(col_names) %>%
    mutate(chrom = as.character(as.integer(strip_chr_prefix(chrom))),
           variants = normalize_variant_id(variants))

  # Label and symmetrize the matrix
  colnames(LD_matrix) <- rownames(LD_matrix) <- LD_variants$variants
  if (all(LD_matrix[lower.tri(LD_matrix)] == 0)) {
    LD_matrix[lower.tri(LD_matrix)] <- t(LD_matrix)[lower.tri(LD_matrix)]
  } else {
    LD_matrix[upper.tri(LD_matrix)] <- t(LD_matrix)[upper.tri(LD_matrix)]
  }

  # Order variants by genomic position
  pos_order <- order(sapply(LD_variants$variants, function(v) as.integer(strsplit(v, ":")[[1]][2])))
  LD_variants <- LD_variants[pos_order, ]
  LD_matrix <- LD_matrix[LD_variants$variants, LD_variants$variants]

  list(LD_matrix = LD_matrix, LD_variants = LD_variants)
}

#' Subset an LD matrix and variant info to a genomic region, optionally
#' further restricted to specific coordinates.
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @noRd
extract_LD_for_region <- function(LD_matrix, variants, region, extract_coordinates) {
  extracted <- subset(variants, chrom == region$chrom & pos >= region$start & pos <= region$end)

  if (!is.null(extract_coordinates)) {
    extract_coordinates <- extract_coordinates %>%
      mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
      select(chrom, pos)
    extracted <- extracted %>%
      mutate(chrom = as.integer(strip_chr_prefix(chrom))) %>%
      merge(extract_coordinates, by = c("chrom", "pos"))
    keep_cols <- intersect(c("chrom", "variants", "pos", "GD", "A1", "A2",
                             "variance", "allele_freq", "n_nomiss"), names(extracted))
    extracted <- select(extracted, all_of(keep_cols))
  }

  mat <- LD_matrix[extracted$variants, extracted$variants, drop = FALSE]
  list(extracted_LD_matrix = mat, extracted_LD_variants = extracted)
}

#' Combine multiple block-level LD matrices into one, handling boundary overlaps.
#' @importFrom utils tail
#' @noRd
create_LD_matrix <- function(LD_matrices, variants) {
  # Merge variant lists, deduplicating boundary overlaps
  merge_variants <- function(variant_list) {
    merged <- character(0)
    for (v in variant_list) {
      ids <- if (is.list(v) && !is.null(v$variants)) v$variants else v
      if (length(ids) == 0) next
      if (length(merged) > 0 && tail(merged, 1) == ids[1]) ids <- ids[-1]
      merged <- c(merged, ids)
    }
    merged
  }

  all_variants <- merge_variants(variants)
  combined <- matrix(0, nrow = length(all_variants), ncol = length(all_variants),
                     dimnames = list(all_variants, all_variants))

  # Place each block into the combined matrix
  for (i in seq_along(LD_matrices)) {
    v <- rownames(LD_matrices[[i]])
    idx <- match(v, all_variants)
    combined[idx, idx] <- LD_matrices[[i]]
  }
  combined
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

  # Use shared helper -- no genotype loading
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

#' Validate that cross-block entries are zero (excluding boundary variants).
#' @noRd
validate_block_structure <- function(matrix, block_metadata, variant_ids) {
  msgs <- character(0)
  n <- length(variant_ids)

  for (i in 1:(nrow(block_metadata) - 1)) {
    for (j in (i + 1):nrow(block_metadata)) {
      si <- block_metadata$start_idx[i]; ei <- block_metadata$end_idx[i]
      sj <- block_metadata$start_idx[j]; ej <- block_metadata$end_idx[j]
      if (si > n || ei > n || sj > n || ej > n) {
        msgs <- c(msgs, paste("Block indices out of range for blocks", i, "and", j))
        next
      }
      # Exclude boundary variants (potential overlaps)
      vi <- variant_ids[si:(ei - 1)]
      vj <- variant_ids[(sj + 1):ej]
      if (length(vi) > 0 && length(vj) > 0) {
        max_val <- max(abs(matrix[vi, vj, drop = FALSE]))
        if (max_val > 1e-10) {
          msgs <- c(msgs, paste("Non-zero correlation between blocks", i, "and", j,
                                "- max:", max_val))
        }
      }
    }
  }
  if (length(msgs) > 0) stop("Matrix lacks expected block structure:\n", paste(msgs, collapse = "\n"))
}

#' @noRd
can_merge <- function(block1, block2, max_size) {
  block1$chrom == block2$chrom && (block1$size + block2$size) <= max_size
}

#' @noRd
merge_two_blocks <- function(block_metadata, idx1, idx2) {
  if (idx1 > idx2) { tmp <- idx1; idx1 <- idx2; idx2 <- tmp }
  result <- block_metadata
  result$end_idx[idx1] <- block_metadata$end_idx[idx2]
  result$size[idx1] <- block_metadata$size[idx1] + block_metadata$size[idx2]
  result <- result[-idx2, ]
  result$block_id <- seq_len(nrow(result))
  result
}

#' Find blocks below min_size and identify the best neighbor to merge with.
#' @noRd
find_merge_candidates <- function(block_metadata, min_size, max_size) {
  candidates <- data.frame(block_idx = integer(), merge_with = integer(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(block_metadata))) {
    if (block_metadata$size[i] >= min_size) next
    prev_ok <- i > 1 && can_merge(block_metadata[i, ], block_metadata[i - 1, ], max_size)
    next_ok <- i < nrow(block_metadata) && can_merge(block_metadata[i, ], block_metadata[i + 1, ], max_size)
    merge_with <- if (prev_ok && next_ok) {
      if (block_metadata$size[i - 1] <= block_metadata$size[i + 1]) i - 1 else i + 1
    } else if (prev_ok) i - 1
      else if (next_ok) i + 1
      else next
    candidates <- rbind(candidates, data.frame(block_idx = i, merge_with = merge_with))
  }
  candidates
}

#' Iteratively merge blocks below min_size with their smallest neighbor.
#' @noRd
merge_blocks <- function(block_metadata, min_size, max_size) {
  if (nrow(block_metadata) <= 1) return(block_metadata)
  repeat {
    candidates <- find_merge_candidates(block_metadata, min_size, max_size)
    if (nrow(candidates) == 0) break
    block_metadata <- merge_two_blocks(block_metadata, candidates$block_idx[1], candidates$merge_with[1])
  }
  block_metadata
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
