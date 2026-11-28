#' Adjust transcript abundance for unwanted variation
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description adjust_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with an additional adjusted abundance column. This method uses scaled counts if present.
#'
#' @importFrom rlang enquo
#' @importFrom rlang enquo quo_is_null quo_name quo
#' @importFrom lifecycle deprecate_warn
#' @importFrom stringr str_c
#' @importFrom tibble as_tibble
#' @importFrom dplyr select pull
#' @importFrom stats as.formula
#' @importFrom SummarizedExperiment assay assayNames assays colData
#' @importFrom stats model.matrix rnorm
#'
#' @name adjust_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param .formula DEPRECATED - A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_interest + batch)
#' @param .factor_unwanted A tidy select, e.g. column names without double quotation. c(batch, country) These are the factor that we want to adjust for, including unwanted batcheffect, and unwanted biological effects.
#' @param .factor_of_interest A tidy select, e.g. column names without double quotation. c(treatment) These are the factor that we want to preserve.
#' @param method A character string. Methods include combat_seq (default), combat and limma_remove_batch_effect.
#'
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param transform DEPRECATED - A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param inverse_transform DEPRECATED - A function that is the inverse of transform (e.g. expm1 is inverse of log1p). This is needed to tranform back the counts after analysis.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#' @details This function adjusts the abundance for (known) unwanted variation.
#' At the moment just an unwanted covariate is allowed at a time using Combat (DOI: 10.1093/bioinformatics/bts034)
#'
#' Underlying method:
#' 	sva::ComBat(data, batch = my_batch,	mod = design,	prior.plots = FALSE, ...)
#'
#' @return A consistent object (to the input) with additional columns for the adjusted counts as `<COUNT COLUMN>_adjusted`
#'
#'
#'
#'
#' @examples
#' ## Load airway dataset for examples
#'
#'   data('airway', package = 'airway')
#'   # Ensure a 'condition' column exists for examples expecting it
#'
#'     SummarizedExperiment::colData(airway)$condition <- as.factor(SummarizedExperiment::colData(airway)$dex)
#'
#'
#'
#'
#'
#' cm = airway
#' # Create a balanced two-level batch within each condition to avoid confounding
#' cond <- SummarizedExperiment::colData(cm)$condition
#' cm$batch <- rep(NA_character_, ncol(cm))
#' for (lev in unique(cond)) {
#'   idx <- which(cond == lev)
#'   cm$batch[idx] <- rep(c('A','B'), length.out = length(idx))
#'
#' }
#' cm$batch <- as.factor(cm$batch)
#'
#' # sva is a Suggests package and is not built for every Bioconductor devel.
#' if (requireNamespace("sva", quietly = TRUE)) {
#'   cm |>
#'     identify_abundant() |>
#'     adjust_abundance(.factor_unwanted = batch, .factor_of_interest = condition, method = "combat_seq")
#' }
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics, 2(3), lqaa078. doi:10.1093/nargab/lqaa078
#'
#' Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8(1), 118–127. doi:10.1093/biostatistics/kxj037
#'
#' Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47. doi:10.1093/nar/gkv007
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#' @export
#'
#'
setGeneric("adjust_abundance", function(.data,
                                        
                                        # DEPRECATED
                                        .formula = NULL,
                                        .factor_unwanted =NULL,
                                        .factor_of_interest = NULL,
                                        
                                        
                                        abundance = assayNames(.data)[1],
                                        .abundance = NULL,
                                        method = "combat_seq",
                                        
                                        ...,
                                        
                                        # DEPRECATED
                                        log_transform = NULL,
                                        transform = NULL,
                                        inverse_transform = NULL
                                        
)
standardGeneric("adjust_abundance"))




.adjust_abundance_se = function(.data,
                                
                                # DEPRECATED
                                .formula = NULL,
                                
                                .factor_unwanted = NULL,
                                .factor_of_interest = NULL,
                                
                                abundance = assayNames(.data)[1],
                                .abundance = NULL,
                                
                                method = "combat_seq",
                                
                                ...,
                                
                                # DEPRECATED
                                log_transform = NULL,
                                
                                transform = NULL,
                                inverse_transform = NULL
) {
  
  # Fix NOTEs
  . = NULL
  
  .abundance <- enquo(.abundance)
  
  # Deprecation logic for .abundance
  if (!quo_is_null(.abundance)) {
    
    lifecycle::deprecate_warn("2.0.0", "test_differential_abundance(.abundance)", "test_differential_abundance(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::quo_name(.abundance)
    }
  }
  
  if (is.null(abundance)) {
    stop("You must provide the `abundance` argument (character).")
  }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("sva")
  
  
  # DEPRECATION OF log_transform
  if (
    (is_present(transform) & !is.null(transform)) |
    is_present(inverse_transform) & !is.null(inverse_transform)
  ) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.11.6", "tidybulk::test_differential_abundance(transform = )", details = "The argument transform and inverse_transform is now deprecated, please use method argument instead specifying \"combat\" or \"combat_seq\".")
    
  }
  
  # Set column name for value scaled
  # Use the actual assay being adjusted for naming
  value_adjusted = paste0(abundance, adjusted_string)
  
  
  # DEPRECATION OF .formula
  if (is_present(.formula) & !is.null(.formula)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.11.6", "tidybulk::test_differential_abundance(.formula = )", details = "The argument .formula is now deprecated, please use factor_unwanted and factor_of_interest. Using the formula, the first factor is of interest and the second is unwanted")
    
    # Check that .formula includes at least two covariates
    if (parse_formula(.formula) |> length() |> st(2))
      stop(
        "The .formula must contain two covariates, the first being the factor of interest, the second being the factor of unwanted variation"
      )
    
    # Check that .formula includes no more than two covariates at the moment
    if (parse_formula(.formula) |> length() |> gt(3))
      warning("tidybulk says: Only the second covariate in the .formula is adjusted for")
    
    
    .factor_of_interest = quo(!!as.symbol(parse_formula(.formula)[1]))
    .factor_unwanted = quo(!!as.symbol(parse_formula(.formula)[2]))
    
  } else {
    
    .factor_of_interest = enquo(.factor_of_interest)
    .factor_unwanted = enquo(.factor_unwanted)
  }
  
  # Create design matrix
  design =
    model.matrix(
      object = as.formula(sprintf("~ %s", colData(.data) |> as_tibble() |> select(!!.factor_of_interest) |> colnames() |> str_c(collapse = '+'))),
      # get first argument of the .formula
      data = colData(.data)
    )
  
  my_batch = colData(.data) |> as_tibble() |> select(!!.factor_unwanted)
  
  
  
  # If no assay is specified take first
  my_assay = abundance
 
  if(tolower(method) == "combat"){
    
    my_assay_adjusted =
      .data |>
      assay(my_assay) |> # Check if log transform is needed
      log1p() |>
      # Add little noise to avoid all 0s for a covariate that would error combat code (not statistics that would be fine)
      (\(y) y + rnorm(length(y), 0, 0.000001))()
    
    
    for(i in colnames(my_batch)){
      my_assay_adjusted =
        my_assay_adjusted |>
        
        # Run combat
        sva::ComBat(
          batch = my_batch[,i] |> pull(1),
          mod = design,
          prior.plots = FALSE,
          ...
        )
    }
    
    # Tranfrom back
    my_assay_adjusted =
      my_assay_adjusted |>
      expm1() |>
      apply(2, pmax, 0)
    
  }
  else if(tolower(method) == "combat_seq"){
    
    my_assay_adjusted =
      .data |>
      
      assay(my_assay)
    
    for(i in ncol(my_batch)){
      
      my_assay_adjusted =
        my_assay_adjusted |>
        sva::ComBat_seq(batch = my_batch[,i] |> pull(1),
                        covar_mod = design,
                        full_mod=TRUE,
                        ...)
    }
    
  }
  else if(tolower(method) == "limma_remove_batch_effect") {
    
    unwanted_covariate_matrix =
      model.matrix(
        object = as.formula(sprintf("~ 0 + %s", colData(.data) |> as_tibble() |> select(!!.factor_unwanted) |> colnames() |> str_c(collapse = '+'))),
        # get first argument of the .formula
        data = colData(.data)
      )
    
    my_assay_adjusted =
      .data |>
      assay(my_assay) |>
      edgeR::cpm(log = TRUE) |>
      limma::removeBatchEffect(
        design = design,
        covariates = unwanted_covariate_matrix,
        ...
      ) |>
      expm1() |>
      apply(2, pmax, 0)
    
  } else {
    stop("tidybulk says: the argument \"method\" must be \"combat_seq\", \"combat\", or \"limma_remove_batch_effect\"")
  }
  
  
  # Add the assay
  my_assay_scaled = list(my_assay_adjusted) |> setNames(value_adjusted)
  
  assays(.data) =  assays(.data) |> c(my_assay_scaled)
  
  # Return
  .data |>
    
    # Add methods
    memorise_methods_used("sva") |>
    
    # Attach column internals
    add_tt_columns(.abundance_adjusted = !!(((function(x, v)	 enquo(v))(x,!!as.symbol(value_adjusted))) |> drop_enquo_env()) )
  
}

#' adjust_abundance
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("adjust_abundance",
          "SummarizedExperiment",
          .adjust_abundance_se)

#' adjust_abundance
#'
#' @docType methods
#' @rdname adjust_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("adjust_abundance",
          "RangedSummarizedExperiment",
          .adjust_abundance_se)

