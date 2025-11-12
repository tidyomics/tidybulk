#' Scale the counts of transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description scale_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and Scales transcript abundance compansating for sequencing depth (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo quo_name quo_is_null
#' @importFrom stats median
#' @importFrom SummarizedExperiment assays colnames
#' @importFrom lifecycle is_present deprecate_warn
#'
#' @name scale_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param method A character string. The scaling method passed to the back-end function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_sample A character string. The name of the reference sample. If NULL the sample with highest total read count will be selected as reference.
#' @param .subset_for_scaling A gene-wise quosure condition. This will be used to filter rows (features/genes) of the dataset. For example
#' @param suffix A character string to append to the scaled abundance column name. Default is "_scaled".
#' @param chunk_sample_size An integer indicating how many samples to process per chunk. Default is `Inf` (no chunking). For HDF5-backed data or large datasets, set to a finite value (e.g., 50) to enable memory-efficient chunked processing with BiocParallel parallelization.
#'
#' @param reference_selection_function DEPRECATED. please use reference_sample.
#' @param ... Further arguments.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#' @details Scales transcript abundance compensating for sequencing depth
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with minimum_counts and minimum_proportion parameters)
#' are filtered out from the scaling procedure.
#' The scaling inference is then applied back to all unfiltered data.
#'
#' Underlying method
#' edgeR::calcNormFactors(.data, method = c("TMM","TMMwsp","RLE","upperquartile"))
#'
#'
#'
#' @return A tbl object with additional columns with scaled data as `<NAME OF COUNT COLUMN>_scaled`
#'
#'
#' @examples
#' ## Load airway dataset for examples
#'
#'   data('airway', package = 'airway')
#'   # Ensure a 'condition' column exists for examples expecting it
#'
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#'
#'
#'  airway |>
#'    identify_abundant() |>
#'    scale_abundance()
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology, 11(3), R25. doi:10.1186/gb-2010-11-3-r25
#'
#' @docType methods
#' @rdname scale_abundance-methods
#' @export

setGeneric("scale_abundance", function(.data,
                                       
                                       
                                       abundance = SummarizedExperiment::assayNames(.data)[1],
                                       method = "TMM",
                                       reference_sample = NULL,
                                       .subset_for_scaling = NULL,
                                       suffix = "_scaled",
                                       chunk_sample_size = Inf,
                                       # DEPRECATED
                                       reference_selection_function = NULL,
                                       ...,
                                       .abundance = NULL)
  standardGeneric("scale_abundance"))



#' @importFrom magrittr multiply_by
#' @importFrom magrittr divide_by
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment colData
#' @importFrom utils tail
#' @importFrom stats na.omit
#' @importFrom stringr str_c
#' @importFrom dplyr mutate select pull arrange slice n_distinct across
#' @importFrom tidyr nest unnest pivot_longer pivot_wider drop_na
#' @importFrom purrr map when
#' @importFrom stringr str_subset str_remove str_replace str_replace_all
#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom purrr map2_dfr
#' @importFrom SummarizedExperiment "colData<-"
#' @importFrom SummarizedExperiment "rowData<-"
#' @importFrom magrittr "%$%"
#' @importFrom SummarizedExperiment "assays<-"
#' @importFrom rlang enquo quo_is_null quo_get_expr
#' @importFrom lifecycle is_present deprecate_warn
#' @importFrom methods is
#'
.scale_abundance_se = function(.data,
                               
                               
                               abundance = SummarizedExperiment::assayNames(.data)[1],
                               method = "TMM",
                               reference_sample = NULL,
                               .subset_for_scaling = NULL,
                               suffix = "_scaled",
                               chunk_sample_size = Inf,
                               # DEPRECATED
                               reference_selection_function = NULL,
                               ...,
                               .abundance = NULL) {
  
  # Fix NOTEs
  . = NULL
  
  # Soft-deprecate .abundance, prefer abundance (character)
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "scale_abundance(.abundance)", "scale_abundance(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  my_assay <- abundance
  
  # Check if package is installed, otherwise install
  check_and_install_packages("edgeR")
  
  
  # DEPRECATION OF reference function
  if (is_present(reference_selection_function) && !is.null(reference_selection_function)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.1.8", "tidybulk::scale_abundance(reference_selection_function = )", details = "The argument reference_selection_function is now deprecated please use reference_sample. By default the reference selection function is max()")
    
  }
  
  # Check that reference sample exists
  if(!is.null(reference_sample) && !reference_sample %in% (.data |> colnames()))
    stop("tidybulk says: your reference sample is not among the samples in your data frame")
  
  # Force evaluation of .subset_for_scaling to avoid S4 dispatch issues with default arguments
  .subset_for_scaling <- force(.subset_for_scaling)
  .subset_for_scaling = enquo(.subset_for_scaling)
  
  .data_filtered =
    filter_if_abundant_were_identified(.data)
  
  # Only apply filtering if a non-NULL subset condition was provided  
  subset_expr <- rlang::quo_get_expr(.subset_for_scaling)
  has_subset <- !is.null(subset_expr) && !identical(subset_expr, quote(NULL))
  
  if (has_subset)
    .data_filtered = filter_genes_on_condition(.data_filtered, !!.subset_for_scaling)
  
  # Filter based on user condition
  
  # Check I have genes left
  if (nrow(.data_filtered) == 0)
    stop("tidybulk says: there are 0 genes that passes the filters (.abundant and/or .subset_for_scaling). Please check your filtering or your data.")
  
  # Get the rownames of features that passed filtering
  features_to_use <- rownames(.data_filtered)

  # We just carry the gene names, no need to stress the memory with the full data frame
  rm(.data_filtered)


  # If not enough genes, warning
  if(length(features_to_use)<100) warning(warning_for_scaling_with_few_genes)
  
  # Set column name for value scaled
  value_scaled = paste0(my_assay, suffix)
  
  # Get reference
  reference <-
    reference_sample
  
  if (is.null(reference)){
      # Get filtered counts for reference selection
  library_size_filtered = assays(.data[features_to_use, ])[[my_assay]] |> colSums(na.rm  = TRUE)
    reference = library_size_filtered |>
    sort() |>
    tail(1) |>
    names()
  }

  
  
  # Communicate the reference if chosen by default
  if(is.null(reference_sample)) message(sprintf("tidybulk says: the sample with largest library size %s was chosen as reference for scaling", reference))
  

  # Calculate TMM normalization factors once for all samples
  chunk_counts_filtered <- assays(.data[features_to_use, ])[[my_assay]] |> na.omit()
  chunk_library_size <- chunk_counts_filtered |> colSums(na.rm = TRUE)
  
  nf <- edgeR::calcNormFactors(
    chunk_counts_filtered,
    refColumn = reference,
    method = method
  )
  
  # Calculate multiplier for all samples
  multiplier <- 
    (chunk_library_size[reference] * nf[reference]) |>
    divide_by(chunk_library_size * nf)
  
  # Add to sample info
  colData(.data)$TMM <- nf
  colData(.data)$multiplier <- multiplier
  
  # If chunk_sample_size is finite, process in chunks with parallelization
  if (!is.finite(chunk_sample_size) || chunk_sample_size >= ncol(.data)) {
    
    # Standard processing without chunking - just apply the multipliers
    .data <- apply_scaling_only(.data, my_assay, multiplier, value_scaled)
    
  } else {

    sample_names <- colnames(.data)
    sample_indices <- seq_along(sample_names)
    sample_groups <- split(sample_names, ceiling(sample_indices / chunk_sample_size))
    
    # Check if BiocParallel is available, otherwise install
    check_and_install_packages("BiocParallel")
    
    # Get the current BiocParallel backend
    bp_param <- BiocParallel::bpparam()
    
    # Inform user about parallelization settings
    if (is(bp_param, "SerialParam")) {
      message("tidybulk says: Processing chunks serially. For parallel processing, configure BiocParallel with BiocParallel::register() before calling this function. For example: BiocParallel::register(BiocParallel::MulticoreParam(workers = 4, progressbar = TRUE))")
    } else {
      message(sprintf("tidybulk says: Processing %d chunks in parallel using %s with %d workers", 
                      length(sample_groups), class(bp_param)[1], BiocParallel::bpnworkers(bp_param)))
    }
    
    chunk_results <-
      BiocParallel::bplapply(
        sample_groups,
        function(current_samples) {
          # Extract chunk
          chunk_se <- .data[, current_samples, drop = FALSE]
          
          # Get multipliers for this chunk
          chunk_multiplier <- multiplier[current_samples]
          
          # Apply scaling to chunk (TMM already calculated)
          chunk_scaled <- apply_scaling_only(
            chunk_se,
            my_assay, 
            chunk_multiplier, 
            value_scaled
          )
          
          chunk_scaled
        },
        BPPARAM = bp_param
      )
    
    # Combine chunks - use SummarizedExperiment::cbind for proper S4 dispatch
    message("tidybulk says: Combining chunks")
    .data <- do.call(SummarizedExperiment::cbind, args = chunk_results)

  }

    scaled_symbol <- rlang::sym(value_scaled)
    scaled_quosure <- drop_enquo_env(rlang::new_quosure(scaled_symbol))
    
   .data |>
      memorise_methods_used(c("edger", "tmm")) |>
      add_tt_columns(.abundance_scaled = !!scaled_quosure)
  
}

#' scale_abundance
#'
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("scale_abundance",
          "SummarizedExperiment",
          .scale_abundance_se)

#' scale_abundance
#'
#' @docType methods
#' @rdname scale_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("scale_abundance",
          "RangedSummarizedExperiment",
          .scale_abundance_se)


  # Core scaling function - applies pre-calculated multipliers to create scaled assay
  apply_scaling_only <- function(data_obj, assay_name, multipliers, scaled_name) {
    
    # Get the assay data and apply scaling
    chunk_assay <- assay(data_obj, assay_name)
    is_delayed <- is(chunk_assay, "DelayedArray")
    
    if (is_delayed) {
      # For DelayedArray, use sweep to create a delayed operation
      check_and_install_packages("DelayedArray")
      my_counts_scaled <- list(
        DelayedArray::sweep(chunk_assay, 2, multipliers, "*")
      ) |> setNames(scaled_name)
    } else {
      # For regular matrices, use matrix multiplication
      my_counts_scaled <- list(
        chunk_assay %*% diag(multipliers)
      ) |> setNames(scaled_name)
    }
    
    colnames(my_counts_scaled[[1]]) <- colnames(chunk_assay)
    
    # Add the scaled assay
    assays(data_obj, withDimnames = FALSE) <- assays(data_obj) |> c(my_counts_scaled)
    data_obj
  }
  