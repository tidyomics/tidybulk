
#' Identify abundant transcripts/genes
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description 
#' Identifies transcripts/genes that are consistently expressed above a threshold across samples. 
#' This function adds a logical column `.abundant` to indicate which features pass the filtering criteria.
#'
#' @param .data A `tbl` or `SummarizedExperiment` object containing transcript/gene abundance data
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param design A design matrix for more complex experimental designs. If provided, this is passed to filterByExpr instead of factor_of_interest.
#' @param formula_design ...
#' @param minimum_counts ...
#' @param minimum_proportion ...
#' @param minimum_count_per_million ...
#' @param ... Further arguments.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#' @param factor_of_interest The name of the column containing groups/conditions for filtering. 
#'        Used by edgeR's filterByExpr to define sample groups. 
#'        \strong{DEPRECATED:} Use 'design' or 'formula_design' instead. This argument will be removed in a future release.
#'
#' @details 
#' This function uses edgeR's filterByExpr() function to identify consistently expressed features.
#' A feature is considered abundant if it has CPM > minimum_counts in at least minimum_proportion 
#' of samples in at least one experimental group (defined by factor_of_interest or design).
#'
#' @return 
#' Returns the input object with an additional logical column `.abundant` indicating which 
#' features passed the abundance threshold criteria.
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
#' # Basic usage
#' airway |> identify_abundant()
#'
#' # With custom thresholds
#' airway |> identify_abundant(
#'   minimum_counts = 5,
#'   minimum_proportion = 0.5
#' )
#'
#' # Using a factor of interest
#' airway |> identify_abundant(factor_of_interest = condition)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @references
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of 
#' multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 
#' 40(10), 4288-4297. DOI: 10.1093/bioinformatics/btp616
#'
#' @importFrom rlang enquo
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#' @importFrom magrittr not
#' @importFrom stats as.formula
#'
#' @docType methods
#' @rdname identify_abundant-methods
#' @export
setGeneric("identify_abundant", function(.data,
                                         abundance = assayNames(.data)[1],
                                         design = NULL,
                                         formula_design = NULL,
                                         minimum_counts = 10,
                                         minimum_proportion = 0.7,
                                         minimum_count_per_million = NULL,
                                         factor_of_interest = NULL,
                                         ..., 
                                         .abundance = NULL
                                         ) # add factor_of_interest
  standardGeneric("identify_abundant"))





.identify_abundant_se = function(.data,
                                 abundance = assayNames(.data)[1],
                                 design = NULL,
                                 formula_design = NULL,
                                 minimum_counts = 10,
                                 minimum_proportion = 0.7,
                                 minimum_count_per_million = NULL,
                                 factor_of_interest = NULL,
                                 ..., 
                                 .abundance = NULL
                                 ) { # <-- add back for backward compatibility
 
  factor_of_interest <- enquo(factor_of_interest)
  
  if (!quo_is_null(factor_of_interest)) {
  # Tidy deprecation warning for factor_of_interest
    lifecycle::deprecate_warn(
      when = "2.0.0",
      what = "identify_abundant(factor_of_interest)",
      with = "identify_abundant(formula_design)",
      details = "The argument 'factor_of_interest' is deprecated and will be removed in a future release. Please use the 'design' or 'formula_design' argument instead."
    )
  }
  
  # Fix NOTEs
  . = NULL
  
  # If formula_design is provided, use it to create the design matrix
  if (!is.null(formula_design)) {
    if (!is.null(design)) warning("Both formula_design and design provided; formula_design will be used.")
    design <- model.matrix(formula_design, data = colData(.data))
  }
  
  # Map factor_of_interest to design if design is NULL and factor_of_interest is provided (for backward compatibility)
  if (!quo_is_null(factor_of_interest) && is.null(design)) {
    # If factor_of_interest is a quosure or symbol, convert to character
    factor_of_interest_chr <- rlang::quo_name(factor_of_interest)

    
    design_formula <- as.formula(paste("~", paste(factor_of_interest_chr, collapse = "+")))
    
      design <- model.matrix(design_formula, data = as.data.frame(colData(.data)))

  }
  
  # Soft-deprecate .abundance, prefer abundance (character)
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "identify_abundant(.abundance)", "identify_abundant(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::quo_name(rlang::ensym(.abundance))
    }
  }
  my_assay <- abundance
  
  if (minimum_counts < 0)
    stop("The parameter minimum_counts must be > 0")
  
  if (minimum_proportion < 0 | minimum_proportion > 1)
    stop("The parameter minimum_proportion must be between 0 and 1")
  
  # If column is present use this instead of doing more work
  if(".abundant" %in% colnames(colData(.data))){
    message("tidybulk says: the column .abundant already exists in colData. Nothing was done")
    return(.data)
  }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("edgeR")

  # Get gene to exclude
  # If minimum_count_per_million is provided, use it and ignore minimum_counts
  if (!is.null(minimum_count_per_million)) {
    gene_to_exclude =
      .data |>
      filterByExpr_SE(
        design = design,
        min.prop = minimum_proportion,
        CPM.Cutoff = minimum_count_per_million,
        assay_name = my_assay
      ) |>
      not() |>
      which() |>
      names()
  } else {
    gene_to_exclude =
      .data |>
      filterByExpr_SE(
        min.count = minimum_counts,
        design = design,
        min.prop = minimum_proportion,
        assay_name = my_assay
      ) |>
      not() |>
      which() |>
      names()
  }
  
  rowData(.data)$.abundant = (rownames(rowData(.data)) %in% gene_to_exclude) |> not()
  
  # Return
  .data
}


#' identify_abundant
#' @param minimum_count_per_million Minimum CPM cutoff to use for filtering (passed to CPM.Cutoff in filterByExpr). If provided, this will override the minimum_counts parameter. Default is NULL (uses edgeR default).
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("identify_abundant",
          "SummarizedExperiment",
          .identify_abundant_se
)

#' identify_abundant
#' @param minimum_count_per_million Minimum CPM cutoff to use for filtering (passed to CPM.Cutoff in filterByExpr). If provided, this will override the minimum_counts parameter. Default is NULL (uses edgeR default).
#'
#' @docType methods
#' @rdname identify_abundant-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("identify_abundant",
          "RangedSummarizedExperiment",
          .identify_abundant_se
)



#' Filter to keep only abundant transcripts/genes
#'
#' \lifecycle{questioning}
#'
#' @description 
#' Filters the data to keep only transcripts/genes that are consistently expressed above 
#' a threshold across samples. This is a filtering version of identify_abundant() that 
#' removes low-abundance features instead of just marking them.
#'
#' @param .data A `tbl` or `SummarizedExperiment` object containing transcript/gene abundance data
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param design A design matrix for more complex experimental designs. If provided, this is passed to filterByExpr instead of factor_of_interest.
#' @param formula_design ...
#' @param minimum_counts ...
#' @param minimum_proportion ...
#' @param minimum_count_per_million ...
#' @param ... Further arguments.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#' @param factor_of_interest The name of the column containing groups/conditions for filtering. 
#'        Used by edgeR's filterByExpr to define sample groups. 
#'        \strong{DEPRECATED:} Use 'design' or 'formula_design' instead. This argument will be removed in a future release.
#'
#' @details 
#' This function uses edgeR's filterByExpr() function to identify and keep consistently expressed features.
#' A feature is kept if it has CPM > minimum_counts in at least minimum_proportion of samples
#' in at least one experimental group (defined by factor_of_interest or design).
#' 
#' This function is similar to identify_abundant() but instead of adding an .abundant column,
#' it filters out the low-abundance features directly.
#'
#' @return 
#' Returns a filtered version of the input object containing only the features that passed
#' the abundance threshold criteria.
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
#' # Basic usage
#' airway |> keep_abundant()
#'
#' # With custom thresholds
#' airway |> keep_abundant(
#'   minimum_counts = 5,
#'   minimum_proportion = 0.5
#' )
#'
#' # Using a factor of interest
#' airway |> keep_abundant(factor_of_interest = condition)
#'
#' @references
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of 
#' multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 
#' 40(10), 4288-4297. DOI: 10.1093/bioinformatics/btp616
#'
#' @importFrom rlang enquo
#' @importFrom dplyr filter
#'
#' @name keep_abundant
#' @title Filter to keep only abundant transcripts/genes
#' @description This function is similar to identify_abundant() but instead of adding an .abundant column,
#' it filters out the low-abundance features directly.
#'
#' @param .data A `tbl` or `SummarizedExperiment` object containing transcript/gene abundance data
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param design A design matrix for more complex experimental designs. If provided, this is passed to filterByExpr instead of factor_of_interest.
#' @param formula_design A formula for creating the design matrix
#' @param minimum_counts The minimum count threshold for a feature to be considered abundant
#' @param minimum_proportion The minimum proportion of samples in which a feature must be abundant
#' @param minimum_count_per_million The minimum count per million threshold
#' @param factor_of_interest The name of the column containing groups/conditions for filtering. DEPRECATED: Use 'design' or 'formula_design' instead.
#' @param ... Further arguments.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#' @return 
#' Returns a filtered version of the input object containing only the features that passed
#' the abundance threshold criteria.
#'
#' @docType methods
#' @export
setGeneric("keep_abundant", function(.data,
                                     abundance = assayNames(.data)[1],
                                     design = NULL,
                                     formula_design = NULL,
                                     minimum_counts = 10,
                                     minimum_proportion = 0.7,
                                     minimum_count_per_million = NULL,
                                     factor_of_interest = NULL,
                                     ..., 
                                     .abundance = NULL
) # add factor_of_interest
  standardGeneric("keep_abundant"))


.keep_abundant_se = function(.data,
                                         abundance = assayNames(.data)[1],
                                         design = NULL,
                                         formula_design = NULL,
                                         minimum_counts = 10,
                                         minimum_proportion = 0.7,
                                         minimum_count_per_million = NULL,
                                         factor_of_interest = NULL,
                                         ..., 
                                         .abundance = NULL
                                         ) # add factor_of_interest
{
  # Fix NOTEs
  . = NULL
  
  # Tidy deprecation warning for factor_of_interest
  factor_of_interest <- enquo(factor_of_interest)
  
  if (!quo_is_null(factor_of_interest)) {
    lifecycle::deprecate_warn(
      when = "2.0.0",
      what = "keep_abundant(factor_of_interest)",
      with = "keep_abundant(formula_design)",
      details = "The argument 'factor_of_interest' is deprecated and will be removed in a future release. Please use the 'design' or 'formula_design' argument instead."
    )
  }
  
  # If formula_design is provided, use it to create the design matrix
  if (!is.null(formula_design)) {
    if (!is.null(design)) warning("Both formula_design and design provided; formula_design will be used.")
    design <- model.matrix(formula_design, data = colData(.data))
  }
  
  # Soft-deprecate .abundance, prefer abundance (character)
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "keep_abundant(.abundance)", "keep_abundant(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::quo_name(.abundance)
    }
  }
  .data =
    .data |>
    identify_abundant(
      minimum_counts = minimum_counts,
      minimum_proportion = minimum_proportion,
      abundance = abundance,
      design = design,
      minimum_count_per_million = minimum_count_per_million,
      factor_of_interest = !!factor_of_interest # pass through
    )
  idx <- rowData(.data)[[".abundant"]]
  .data[idx,]
}

#' keep_abundant
#'
#' @docType methods
#' @inheritParams keep_abundant
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_abundant",
          "SummarizedExperiment",
          .keep_abundant_se)

#' keep_abundant
#'
#' @docType methods
#' @inheritParams keep_abundant
#' @return A `SummarizedExperiment` object
#'
setMethod("keep_abundant",
          "RangedSummarizedExperiment",
          .keep_abundant_se)



