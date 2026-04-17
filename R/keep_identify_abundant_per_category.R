
#' Identify abundant transcripts/genes in a per-category manner
#'
#' \lifecycle{experimental}
#'
#' @description 
#' Identifies transcripts/genes that are consistently expressed above a threshold across samples. 
#' This function adds a logical column `.abundant` to indicate which features pass the filtering criteria.
#'
#' @param .data A `SummarizedExperiment` object containing transcript/gene abundance data
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param design A design matrix for more complex experimental designs.
#' @param formula_design A formula to generate the design matrix. Overrides `design` when both are provided
#' @param minimum_counts The minimum count threshold for a feature to be considered abundant
#' @param minimum_proportion The minimum proportion of samples in which a feature must be abundant
#' @param minimum_count_per_million ...
#' @param minimum_category The minimum number of categories/experimental groups that have sufficient abundance samples.
#' @param force Should existing .abundant column be replaced; defaults to FALSE.
#' @param coerce_design Should non-categorical design matrix be coerced; defaults to FALSE.
#' @param ... Further arguments taken for compatibility, but are not used.
#'
#' @details 
#' This function provides an alternative filtering solution for edgeR's filterByExpr() function to identify consistently expressed features.
#' A feature is considered abundant if it has CPM > minimum_counts in at least minimum_proportion 
#' of samples in at least minimum_category experimental groups (defined by design/formula_design).
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
#' airway |> identify_abundant_per_category()
#'
#' # With custom thresholds
#' airway |> identify_abundant_per_category(
#'   minimum_counts = 5,
#'   minimum_proportion = 0.5
#' )
#'
#' # Using a factor of interest
#' airway |> identify_abundant_per_category(formula_design = ~condition)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @references
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of 
#' multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 
#' 40(10), 4288-4297. DOI: 10.1093/bioinformatics/btp616
#'
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#' @importFrom magrittr not
#' @importFrom stats as.formula
#'
#' @docType methods
#' @rdname identify_abundant_per_category-methods
#' @export
setGeneric("identify_abundant_per_category",
           function(.data,
                    abundance = assayNames(.data)[1],
                    design = NULL,
                    formula_design = NULL,
                    minimum_counts = 10,
                    minimum_proportion = 0.7,
                    minimum_count_per_million = NULL,
                    minimum_category = 1,
                    force = FALSE,
                    coerce_design = FALSE,
                    ...
                    )
             standardGeneric("identify_abundant_per_category"))

.identify_abundant_per_category_se = function(
    .data,
    abundance = assayNames(.data)[1],
    design = NULL,
    formula_design = NULL,
    minimum_counts = 10,
    minimum_proportion = 0.7,
    minimum_count_per_million = NULL,
    minimum_category=1,
    force=FALSE,
    coerce_design = FALSE,
    ...
) {
  # Fix NOTEs
  . = NULL
  
  # If formula_design is provided, use it to create the design matrix
  if (!is.null(formula_design)) {
    if (!is.null(design)) warning("Both formula_design and design provided; formula_design will be used.")
    design <- model.matrix(formula_design, data = colData(.data))
  }else{
    # if neither is provided, populate design from ~ 1 (single-group, intercept only)
    if (is.null(design)){
      warning("Neither formula_design and design provided; build design from formula ~ 1 .")
      design <- model.matrix( ~ 1, data = colData(.data))
    }
  }
  
  # check if design is categorical
  if(any(!design%in%c(0,1))){
    if(coerce_design){
      message("tidybulk says: the design matrix has columns containing elements other than 0 or 1, but user choose to coerce it by setting coerce_design=TRUE.")
      non_cat_cols = apply(design,2,function(x){any(!x%in%c(0,1))})
      design[,non_cat_cols] = 1
    }else{
      stop("The design matrix has columns containing elements other than 0 or 1. Filtering aborted. Check design or use coerce_design=TRUE to coerce those columns to 1.")
    }
  }
  
  # check if parameters are valid
  if (minimum_counts < 0)
    stop("The parameter minimum_counts must be >= 0")
  
  if(!is.null(minimum_count_per_million)){
    if (minimum_count_per_million < 0)
      stop("The parameter minimum_count_per_million must be >= 0")
    if (minimum_counts > 0) message("tidybulk says: both minimum_counts and minimum_count_per_million are provided. Use only minimum_count_per_million.")
  }

  if (minimum_proportion < 0 | minimum_proportion > 1)
    stop("The parameter minimum_proportion must be between 0 and 1")
  
  if (minimum_category < 0)
    stop("The parameter minimum_category must be >= 0")
  
  # If column is present, use this for user awareness
  if(".abundant" %in% colnames(rowData(.data))){
    if(force){
      message("tidybulk says: the column .abundant already exists in colData, but user choose to replace it by setting force=TRUE.")
    }else{
      message("tidybulk says: the column .abundant already exists in colData. Nothing was done. Use force=TRUE to replace existing column.")
      return(.data)
    }
  }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("edgeR")
  
  # Get logical vector of abundant features (TRUE = keep/abundant)
  # Check if abundant in a sample 
  if(is.null(minimum_count_per_million)){
    count_is_abundant = assays(.data)[[abundance]] >= minimum_counts
  }else{
    CPM = edgeR::cpm(assays(.data)[[abundance]],lib.size=colSums(assays(.data)[[abundance]]))
    count_is_abundant = CPM  >= minimum_count_per_million
  }

  # calculate the proportion of abundant samples across categories for each feature
  prop_abundant_per_category = apply(count_is_abundant,1,function(x){
    unname(sapply(split(x,apply(design,1,paste0,collapse="_")),function(xx){
      sum(xx)/length(xx)
    }))
  },simplify = FALSE)
  
  # determine if to keep feature based on the proportions
  if_keep_feature = sapply(prop_abundant_per_category,function(this_prop){
    sum(this_prop >= minimum_proportion) >= minimum_category
  })
  
  # add/replace .abundant column
  rowData(.data)$.abundant = if_keep_feature
  
  # Return
  .data
}

#' identify_abundant_per_category
#' @param minimum_count_per_million Minimum CPM cutoff to use for filtering. If provided, this will override the minimum_counts parameter. Default is NULL.
#'
#' @docType methods
#' @rdname identify_abundant_per_category-methods
#'
#'
setMethod("identify_abundant_per_category",
          "SummarizedExperiment",
          .identify_abundant_per_category_se
)

#' identify_abundant_per_category
#' @param minimum_count_per_million Minimum CPM cutoff to use for filtering. If provided, this will override the minimum_counts parameter. Default is NULL.
#'
#' @docType methods
#' @rdname identify_abundant_per_category-methods
#'
#'
setMethod("identify_abundant_per_category",
          "RangedSummarizedExperiment",
          .identify_abundant_per_category_se
)



#' Filter to keep only abundant transcripts/genes in a per-category manner
#'
#' \lifecycle{experimental}
#'
#' @description 
#' Filters the data to keep only transcripts/genes that are consistently expressed above 
#' a threshold across samples. This is a filtering version of identify_abundant_per_category() that 
#' removes low-abundance features instead of just marking them.
#'
#' @param .data A `SummarizedExperiment` object containing transcript/gene abundance data
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param design A design matrix for more complex experimental designs.
#' @param formula_design A formula to generate the design matrix. Overrides `design` when both are provided
#' @param minimum_counts The minimum count threshold for a feature to be considered abundant
#' @param minimum_proportion The minimum proportion of samples in which a feature must be abundant
#' @param minimum_count_per_million ...
#' @param minimum_category The minimum number of categories/experimental groups that have sufficient abundance samples.
#' @param force Should existing .abundant column be replaced; defaults to FALSE.
#' @param coerce_design Should non-categorical design matrix be coerced; defaults to FALSE.
#' @param ... Further arguments taken for compatibility, but are not used.
#'
#' @details 
#' This function provides an alternative filtering solution for edgeR's filterByExpr() function to identify consistently expressed features.
#' A feature is considered abundant if it has CPM > minimum_counts in at least minimum_proportion 
#' of samples in at least minimum_category experimental groups (defined by design/formula_design).
#' 
#' This function is similar to identify_abundant_per_category() but in addition to writing to an .abundant column,
#' it also filters out the low-abundance features directly.
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
#' airway |> keep_abundant_per_category()
#'
#' # With custom thresholds
#' airway |> keep_abundant_per_category(
#'   minimum_counts = 5,
#'   minimum_proportion = 0.5
#' )
#'
#' # Using a factor of interest
#' airway |> keep_abundant_per_category(formula_design = ~condition)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @references
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of 
#' multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 
#' 40(10), 4288-4297. DOI: 10.1093/bioinformatics/btp616
#'
#' @importFrom dplyr filter
#'
#'
#' @docType methods
#' @rdname keep_abundant_per_category-methods
#' @export
#' 
setGeneric("keep_abundant_per_category", 
           function(.data,
                    abundance = assayNames(.data)[1],
                    design = NULL,
                    formula_design = NULL,
                    minimum_counts = 10,
                    minimum_proportion = 0.7,
                    minimum_count_per_million = NULL,
                    minimum_category = 1,
                    force = FALSE,
                    coerce_design = FALSE,
                    ...
                    )
             standardGeneric("keep_abundant_per_category"))

.keep_abundant_per_category_se = function(.data,
                             abundance = assayNames(.data)[1],
                             design = NULL,
                             formula_design = NULL,
                             minimum_counts = 10,
                             minimum_proportion = 0.7,
                             minimum_count_per_million = NULL,
                             minimum_category = 1,
                             force = FALSE,
                             coerce_design = FALSE,
                             ...
)
{
  # Fix NOTEs
  . = NULL

  # If formula_design is provided, use it to create the design matrix
  if (!is.null(formula_design)) {
    if (!is.null(design)) warning("Both formula_design and design provided; formula_design will be used.")
    design <- model.matrix(formula_design, data = colData(.data))
  }else{
    # if neither is provided, populate design from ~ 1 (single-group, intercept only)
    if (is.null(design)){
      warning("Neither formula_design and design provided; build design from formula ~ 1 .")
      design <- model.matrix( ~ 1, data = colData(.data))
    }
  }
  
  # check if design is categorical
  if(any(!design%in%c(0,1))){
    if(coerce_design){
      message("tidybulk says: the design matrix has columns containing elements other than 0 or 1, but user choose to coerce it by setting coerce_design=TRUE.")
      non_cat_cols = apply(design,2,function(x){any(!x%in%c(0,1))})
      design[,non_cat_cols] = 1
    }else{
      stop("The design matrix has columns containing elements other than 0 or 1. Filtering aborted. Check design or use coerce_design=TRUE to coerce those columns to 1.")
    }
  }
  
  # identify abundant features
  .data =
    .data |>
    identify_abundant_per_category(
      abundance = abundance,
      design = design,
      minimum_counts = minimum_counts,
      minimum_proportion = minimum_proportion,
      minimum_count_per_million = minimum_count_per_million,
      minimum_category = minimum_category,
      force=force,
      coerce_design=coerce_design,
      ...
    )
  
  # filter identified features and return
  idx <- rowData(.data)[[".abundant"]]
  .data[idx,]
}

#' keep_abundant_per_category
#'
#' @docType methods
#' @rdname keep_abundant_per_category-methods
#'
setMethod("keep_abundant_per_category",
          "SummarizedExperiment",
          .keep_abundant_per_category_se)

#' keep_abundant
#'
#' @docType methods
#' @rdname keep_abundant_per_category-methods
#'
setMethod("keep_abundant_per_category",
          "RangedSummarizedExperiment",
          .keep_abundant_per_category_se)



