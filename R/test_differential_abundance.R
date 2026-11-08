#' Perform differential transcription testing using edgeR quasi-likelihood (QLT), edgeR likelihood-ratio (LR), limma-voom, limma-voom-with-quality-weights or DESeq2
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_differential_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo quo_name
#' @importFrom magrittr not
#' @importFrom dplyr select mutate
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom stats lm
#'
#'
#' @name test_differential_abundance
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. If there is more than one factor, they should be in the order factor of interest + additional factors.
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param contrasts This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A character vector. Available methods are "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT), "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb". Only one method can be specified at a time. For glmmSeq, pass `formula_dispersion` (a fixed-effects formula for edgeR) to plug in tagwise dispersion, or omit it to let each gene estimate phi.
#' @param test_above_log2_fold_change A positive real value. This works for edgeR and limma_voom methods. It uses the `treat` function, which tests that the difference in abundance is bigger than this threshold rather than zero \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the back-end functions: edgeR and limma-voom (i.e., edgeR::normLibSizes; "TMM","TMMwsp","RLE","upperquartile"). Setting the parameter to \"none\" will skip the compensation for sequencing-depth for the method edgeR or limma_voom.
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param prefix A character string. The prefix you would like to add to the result columns. It is useful if you want to compare several methods.
#' @param formula_dispersion A fixed-effects formula for edgeR tagwise dispersion (glmmSeq). If NULL, each gene estimates phi.
#' @param significance_threshold DEPRECATED - A real between 0 and 1 (usually 0.05).
#' @param fill_missing_values DEPRECATED - A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param .contrasts DEPRECATED - This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param ... Further arguments passed to some of the internal experimental functions. For glmmSeq: `.scaling_factor` skips TMM and uses that column from [scale_abundance()].
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#'
#' @details This function provides the option to use edgeR \url{https://doi.org/10.1093/bioinformatics/btp616}, limma-voom \url{https://doi.org/10.1186/gb-2014-15-2-r29}, limma_voom_sample_weights \url{https://doi.org/10.1093/nar/gkv412} or  DESeq2 \url{https://doi.org/10.1186/s13059-014-0550-8} to perform the testing.
#' All methods use raw counts, irrespective of if scale_abundance or adjust_abundance have been calculated, therefore it is essential to add covariates such as batch effects (if applicable) in the formula.
#'
#' Underlying method for edgeR framework:
#'
#' 	.data |>
#'
#' 	# Filter
#'	keep_abundant(
#'			factor_of_interest = !!(as.symbol(parse_formula(.formula)[1])),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'		) |>
#'
#'			# Format
#'			select(!!.transcript,!!.sample,!!.abundance) |>
#'			spread(!!.sample,!!.abundance) |>
#'			as_matrix(rownames = !!.transcript) |>
#'
#'			# edgeR
#'			edgeR::DGEList(counts = .) |>
#'			edgeR::normLibSizes(method = scaling_method) |>
#'			edgeR::estimateDisp(design) |>
#'
#'			# Fit
#'			edgeR::glmQLFit(design) |> // or glmFit according to choice
#'			edgeR::glmQLFTest(coef = 2, contrast = my_contrasts) // or glmLRT according to choice
#'
#'
#'	Underlying method for DESeq2 framework:
#'
#'	keep_abundant(
#'			factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'	) |>
#'
#'	# DESeq2
#'	DESeq2::DESeqDataSet(design = .formula) |>
#'	DESeq2::DESeq() |>
#'	DESeq2::results()
#'
#'
#'
#' Underlying method for glmmSeq framework:
#'
#' counts = assay(.data, my_assay)
#'
#' # Tagwise phi from edgeR if formula_dispersion is given.
#'
#'   glmmSeq( .formula,
#'            countdata = counts ,
#'            metadata =   metadata |> as.data.frame(),
#'            dispersion = dispersion,
#'            progress = TRUE,
#'            method = method |> str_remove("(?i)^glmmSeq_" ),
#'   )
#'
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the test (e.g.,  log fold change, p-value and false discovery rate).
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
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#'
#'  # edgeR (default method)
#'
#'  airway |>
#'  identify_abundant() |>
#' 	test_differential_abundance( ~ condition, method = "edgeR_quasi_likelihood" )
#'
#'  # You can also explicitly specify the method
#'  airway |>
#'  identify_abundant() |>
#' 	test_differential_abundance( ~ condition, method = "edgeR_quasi_likelihood" )
#'
#' 	# The function `test_differential_abundance` operates with contrasts too
#'
#'  airway |>
#'  identify_abundant(factor_of_interest = condition) |>
#'  test_differential_abundance(
#' 	    ~ 0 + condition,
#' 	    contrasts = c( "conditiontrt - conditionuntrt"),
#'     method = "edgeR_quasi_likelihood"
#'  )
#'
#'  # DESeq2 - equivalent for limma-voom
#'
#' my_se_mini = airway
#' my_se_mini$condition  = factor(my_se_mini$condition)
#'
#' # demontrating with `fitType` that you can access any arguments to DESeq()
#' my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'        test_differential_abundance( ~ condition, method="deseq2", fitType="local")
#'
#' # testing above a log2 threshold, passes along value to lfcThreshold of results()
#' res <- my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'         test_differential_abundance( ~ condition, method="deseq2",
#'             fitType="local",
#'             test_above_log2_fold_change=4 )
#'
#' # Use random intercept and random effect models
#'
#'  airway[1:50,] |>
#'   identify_abundant(factor_of_interest = condition) |>
#'   test_differential_abundance(
#'     ~ condition + (1 + condition | cell),
#'     method = "glmmseq_lme4",
#'     formula_dispersion = ~ condition + cell,
#'     cores = 1
#'   )
#'
#' # confirm that lfcThreshold was used
#' # Not run to keep the example fast
#' # res |>
#' #     mcols() |>
#' #     DESeq2::DESeqResults() |>
#' #     DESeq2::plotMA()
#'
#'
#' # The function `test_differential_abundance` operates with contrasts too
#'
#'  my_se_mini |>
#'  identify_abundant() |>
#'  test_differential_abundance(
#' 	    ~ 0 + condition,
#' 	    contrasts = list(c("condition", "trt", "untrt")),
#' 	    method="deseq2",
#'      fitType="local"
#'  )
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), 4288-4297. doi:10.1093/nar/gks042
#'
#' Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. doi:10.1186/s13059-014-0550-8
#'
#' Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology, 15(2), R29. doi:10.1186/gb-2014-15-2-r29
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#' @export
#'
setGeneric("test_differential_abundance", function(.data,
                                                   .formula,
                                                   
                                                   
                                                   abundance =  assayNames(.data)[1],
                                                   contrasts = NULL,
                                                   method = c("edgeR_quasi_likelihood", "edgeR_likelihood_ratio", "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"),
                                                   test_above_log2_fold_change = NULL,
                                                   scaling_method = "TMM",
                                                   omit_contrast_in_colnames = FALSE,
                                                   prefix = "",
                                                   formula_dispersion = NULL,
                                                   ...,
                                                   
                                                   # DEPRECATED
                                                   significance_threshold = NULL,
                                                   fill_missing_values = NULL,
                                                   .contrasts = NULL,
                                                   .abundance = NULL
)
standardGeneric("test_differential_abundance"))



#' @importFrom rlang inform
.test_differential_abundance_se = function(.data,
                                           .formula,
                                           
                                           
                                           abundance =  assayNames(.data)[1],
                                           contrasts = NULL,
                                           method = c("edgeR_quasi_likelihood", "edgeR_likelihood_ratio", "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"),
                                           test_above_log2_fold_change = NULL,
                                           scaling_method = "TMM",
                                           omit_contrast_in_colnames = FALSE,
                                           prefix = "",
                                           formula_dispersion = NULL,
                                           ...,
                                           
                                           # DEPRECATED
                                           significance_threshold = NULL,
                                           fill_missing_values = NULL,
                                           .contrasts = NULL,
                                           .abundance = NULL)
{

  .abundance <- enquo(.abundance)

  # Deprecation logic for .abundance
  if (!quo_is_null(.abundance)) {
    
    lifecycle::deprecate_warn("2.0.0", "test_differential_abundance(.abundance)", "test_differential_abundance(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::quo_name(.abundance)
    }
  }
  
  # Fix NOTEs
  . = NULL
  
  # DEPRECATION OF .constrasts
  if (is_present(.contrasts) & !is.null(.contrasts)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.7.4", "tidybulk::test_differential_abundance(.contrasts = )", details = "The argument .contrasts is now deprecated please use contrasts (without the dot).")
    
    contrasts = .contrasts
  }
  
  # Clearly state what counts are used
  rlang::inform("=====================================
tidybulk says: All testing methods use raw counts, irrespective of if scale_abundance
or adjust_abundance have been calculated. Therefore, it is essential to add covariates
such as batch effects (if applicable) in the formula.
=====================================", .frequency_id = "All testing methods use raw counts",  .frequency = "once")
  
  
  # Validate method parameter
  if(length(method) > 1) {
    stop("tidybulk says: only one method can be specified at a time. Please choose one method from the available options.")
  }
  
  # Use the first method if multiple are provided (for backward compatibility)
  method = method[1]
  
  # Test test_above_log2_fold_change
  if(!is.null(test_above_log2_fold_change) && test_above_log2_fold_change < 0)
    stop("tidybulk says: test_above_log2_fold_change should be a positive real or NULL")
  
  # Filter abundant if performed
  .data = filter_if_abundant_were_identified(.data)
  
  if(tolower(method) %in% c("edger_quasi_likelihood", "edger_likelihood_ratio", "edger_robust_likelihood_ratio"))
    my_differential_abundance =
    get_differential_transcript_abundance_bulk_SE(
      .data,
      .formula,
      abundance = abundance,
      .contrasts = contrasts,
      sample_annotation = colData(.data),
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  
  else if (grepl("voom", method))
    my_differential_abundance =
    get_differential_transcript_abundance_bulk_voom_SE(
      .data,
      .formula,
      abundance = abundance,
      .contrasts = contrasts,
      sample_annotation = colData(.data),
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  
  else if(tolower(method)=="deseq2")
    my_differential_abundance =
    get_differential_transcript_abundance_deseq2_SE(
      .data,
      .formula,
      abundance = abundance,
      .contrasts = contrasts,
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      ...
    )
  
  
  else if(	tolower(method) %in% c("glmmseq_lme4", "glmmseq_glmmtmb"))
    my_differential_abundance =
    get_differential_transcript_abundance_glmmSeq_SE(
      .data,
      .formula,
      abundance = abundance,
      .contrasts = contrasts,
      sample_annotation = colData(.data),
      method = method,
      test_above_log2_fold_change = test_above_log2_fold_change,
      scaling_method = scaling_method,
      omit_contrast_in_colnames = omit_contrast_in_colnames,
      prefix = prefix,
      formula_dispersion = formula_dispersion,
      ...
    )
  else
    stop("tidybulk says: the only methods supported at the moment are \"edgeR_quasi_likelihood\" (i.e., QLF), \"edgeR_likelihood_ratio\" (i.e., LRT), \"edger_robust_likelihood_ratio\", \"DESeq2\", \"limma_voom\", \"limma_voom_sample_weights\", \"glmmseq_lme4\", \"glmmseq_glmmtmb\"")
  
  # Add results
  stats_matrix <- my_differential_abundance$result |>
    as_matrix(rownames = "transcript")
  stats_matrix <- stats_matrix[match(rownames(rowData(.data)), rownames(stats_matrix)), , drop = FALSE]
  rowData(.data) = cbind(rowData(.data), stats_matrix)
  
  
  # Add bibliography
  data_obj_intermediate <- .data
  method_lower <- tolower(method)

  method_map <- list(
    edger_likelihood_ratio = c("edger", "edgeR_likelihood_ratio"),
    edger_quasi_likelihood = c("edger", "edgeR_quasi_likelihood"),
    edger_robust_likelihood_ratio = c("edger", "edger_robust_likelihood_ratio"),
    limma_voom = "voom",
    limma_voom_sample_weights = "voom_sample_weights",
    deseq2 = "deseq2",
    glmmseq_lme4 = "glmmseq",
    glmmseq_glmmtmb = "glmmseq"
  )

  if (method_lower %in% names(method_map)) {
    data_obj_intermediate <- memorise_methods_used(data_obj_intermediate, method_map[[method_lower]])
  } else {
    stop("tidybulk says: method not supported")
  }

  if (!is.null(test_above_log2_fold_change)) {
    data_obj_intermediate <- memorise_methods_used(data_obj_intermediate, "treat")
  }

    data_obj_intermediate <- attach_to_metadata(data_obj_intermediate, my_differential_abundance$result_raw, paste0(method, "_fit"))
    data_obj_intermediate <- attach_to_metadata(data_obj_intermediate, my_differential_abundance$de_object, paste0(method, "_object"))
   
    rlang::inform(
      sprintf("tidybulk says: to access the DE object do `metadata(.)$tidybulk$%s_object`", method),
      .frequency_id = sprintf("Access DE results %s", method),
      .frequency = "always"
    )

        rlang::inform(
      sprintf("tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$%s_fit`", method),
      .frequency_id = sprintf("Access DE results %s", method),
      .frequency = "always"
    )
    data_obj_intermediate
  
  
  
}

#' test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_abundance",
  "SummarizedExperiment",
  .test_differential_abundance_se
)

#' test_differential_abundance
#'
#' @docType methods
#' @rdname test_differential_abundance-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_abundance",
  "RangedSummarizedExperiment",
  .test_differential_abundance_se
)

#' Perform differential expression testing using edgeR quasi-likelihood (QLT), edgeR likelihood-ratio (LR), limma-voom, limma-voom-with-quality-weights or DESeq2
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_differential_expression() is an alias for test_differential_abundance() that takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo quo_name
#' @importFrom magrittr not
#' @importFrom dplyr select mutate
#' @importFrom SummarizedExperiment rowData colData
#'
#'
#' @name test_differential_expression
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula representing the desired linear model. If there is more than one factor, they should be in the order factor of interest + additional factors.
#' @param abundance The name of the transcript/gene abundance column (character, preferred)
#' @param contrasts This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A character vector. Available methods are "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT), "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb". Only one method can be specified at a time. For glmmSeq, pass `formula_dispersion` (a fixed-effects formula for edgeR) to plug in tagwise dispersion, or omit it to let each gene estimate phi.
#' @param test_above_log2_fold_change A positive real value. This works for edgeR and limma_voom methods. It uses the `treat` function, which tests that the difference in abundance is bigger than this threshold rather than zero \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the back-end functions: edgeR and limma-voom (i.e., edgeR::normLibSizes; "TMM","TMMwsp","RLE","upperquartile"). Setting the parameter to \"none\" will skip the compensation for sequencing-depth for the method edgeR or limma_voom.
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param prefix A character string. The prefix you would like to add to the result columns. It is useful if you want to compare several methods.
#' @param formula_dispersion A fixed-effects formula for edgeR tagwise dispersion (glmmSeq). If NULL, each gene estimates phi.
#' @param significance_threshold DEPRECATED - A real between 0 and 1 (usually 0.05).
#' @param fill_missing_values DEPRECATED - A boolean. Whether to fill missing sample/transcript values with the median of the transcript. This is rarely needed.
#' @param .contrasts DEPRECATED - This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param ... Further arguments passed to some of the internal experimental functions. For glmmSeq: `.scaling_factor` skips TMM and uses that column from [scale_abundance()].
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#'
#' @details This function provides the option to use edgeR \url{https://doi.org/10.1093/bioinformatics/btp616}, limma-voom \url{https://doi.org/10.1186/gb-2014-15-2-r29}, limma_voom_sample_weights \url{https://doi.org/10.1093/nar/gkv412} or  DESeq2 \url{https://doi.org/10.1186/s13059-014-0550-8} to perform the testing.
#' All methods use raw counts, irrespective of if scale_abundance or adjust_abundance have been calculated, therefore it is essential to add covariates such as batch effects (if applicable) in the formula.
#'
#' Underlying method for edgeR framework:
#'
#' 	.data |>
#'
#' 	# Filter
#'	keep_abundant(
#'			factor_of_interest = !!(as.symbol(parse_formula(.formula)[1])),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'		) |>
#'
#'			# Format
#'			select(!!.transcript,!!.sample,!!.abundance) |>
#'			spread(!!.sample,!!.abundance) |>
#'			as_matrix(rownames = !!.transcript) |>
#'
#'			# edgeR
#'			edgeR::DGEList(counts = .) |>
#'			edgeR::normLibSizes(method = scaling_method) |>
#'			edgeR::estimateDisp(design) |>
#'
#'			# Fit
#'			edgeR::glmQLFit(design) |> // or glmFit according to choice
#'			edgeR::glmQLFTest(coef = 2, contrast = my_contrasts) // or glmLRT according to choice
#'
#'
#'	Underlying method for DESeq2 framework:
#'
#'	keep_abundant(
#'			factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#'			minimum_counts = minimum_counts,
#'			minimum_proportion = minimum_proportion
#'	) |>
#'
#'	# DESeq2
#'	DESeq2::DESeqDataSet(design = .formula) |>
#'	DESeq2::DESeq() |>
#'	DESeq2::results()
#'
#'
#'
#' Underlying method for glmmSeq framework:
#'
#' counts = assay(.data, my_assay)
#'
#' # Tagwise phi from edgeR if formula_dispersion is given.
#'
#'   glmmSeq( .formula,
#'            countdata = counts ,
#'            metadata =   metadata |> as.data.frame(),
#'            dispersion = dispersion,
#'            progress = TRUE,
#'            method = method |> str_remove("(?i)^glmmSeq_" ),
#'   )
#'
#'
#' @return A consistent object (to the input) with additional columns for the statistics from the test (e.g.,  log fold change, p-value and false discovery rate).
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
#'     SummarizedExperiment::colData(airway)$condition <- SummarizedExperiment::colData(airway)$dex
#'
#'
#'
#'  # edgeR (default method)
#'
#'  airway |>
#'  identify_abundant() |>
#' 	test_differential_expression( ~ condition, method = "edgeR_quasi_likelihood" )
#'
#'  # You can also explicitly specify the method
#'  airway |>
#'  identify_abundant() |>
#' 	test_differential_expression( ~ condition, method = "edgeR_quasi_likelihood" )
#'
#' 	# The function `test_differential_expression` operates with contrasts too
#'
#'  airway |>
#'  identify_abundant(factor_of_interest = condition) |>
#'  test_differential_expression(
#' 	    ~ 0 + condition,
#' 	    contrasts = c( "conditiontrt - conditionuntrt"),
#'     method = "edgeR_quasi_likelihood"
#'  )
#'
#'  # DESeq2 - equivalent for limma-voom
#'
#' my_se_mini = airway
#' my_se_mini$condition  = factor(my_se_mini$condition)
#'
#' # demontrating with `fitType` that you can access any arguments to DESeq()
#' my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'        test_differential_expression( ~ condition, method="deseq2", fitType="local")
#'
#' # testing above a log2 threshold, passes along value to lfcThreshold of results()
#' res <- my_se_mini  |>
#'    identify_abundant(factor_of_interest = condition) |>
#'         test_differential_expression( ~ condition, method="deseq2",
#'             fitType="local",
#'             test_above_log2_fold_change=4 )
#'
#' # Use random intercept and random effect models
#'
#'  airway[1:50,] |>
#'   identify_abundant(factor_of_interest = condition) |>
#'   test_differential_expression(
#'     ~ condition + (1 + condition | cell),
#'     method = "glmmseq_lme4",
#'     formula_dispersion = ~ condition + cell,
#'     cores = 1
#'   )
#'
#' # confirm that lfcThreshold was used
#'
#' \dontrun{
#'     res |>
#'         mcols() |>
#'         DESeq2::DESeqResults() |>
#'         DESeq2::plotMA()
#' }
#'
#' # The function `test_differential_expression` operates with contrasts too
#'
#'  my_se_mini |>
#'  identify_abundant() |>
#'  test_differential_expression(
#' 	    ~ 0 + condition,
#' 	    contrasts = list(c("condition", "trt", "untrt")),
#' 	    method="deseq2",
#'          fitType="local"
#'  )
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), 4288-4297. doi:10.1093/nar/gks042
#'
#' Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. doi:10.1186/s13059-014-0550-8
#'
#' Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology, 15(2), R29. doi:10.1186/gb-2014-15-2-r29
#'
#' @docType methods
#' @rdname test_differential_expression-methods
#' @export
#'
setGeneric("test_differential_expression", function(.data,
                                                   .formula,
                                                   
                                                   
                                                   abundance =  assayNames(.data)[1],
                                                   contrasts = NULL,
                                                   method = c("edgeR_quasi_likelihood", "edgeR_likelihood_ratio", "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"),
                                                   test_above_log2_fold_change = NULL,
                                                   scaling_method = "TMM",
                                                   omit_contrast_in_colnames = FALSE,
                                                   prefix = "",
                                                   formula_dispersion = NULL,
                                                   ...,
                                                   
                                                   # DEPRECATED
                                                   significance_threshold = NULL,
                                                   fill_missing_values = NULL,
                                                   .contrasts = NULL,
                                                   .abundance = NULL
)
standardGeneric("test_differential_expression"))

#' test_differential_expression
#'
#' @docType methods
#' @rdname test_differential_expression-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_expression",
  "SummarizedExperiment",
  function(.data, .formula, abundance = assayNames(.data)[1], contrasts = NULL, method = c("edgeR_quasi_likelihood", "edgeR_likelihood_ratio", "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"), test_above_log2_fold_change = NULL, scaling_method = "TMM", omit_contrast_in_colnames = FALSE, prefix = "", formula_dispersion = NULL, ..., significance_threshold = NULL, fill_missing_values = NULL, .contrasts = NULL, .abundance = NULL) {
    test_differential_abundance(.data, .formula, abundance = abundance, contrasts = contrasts, method = method, test_above_log2_fold_change = test_above_log2_fold_change, scaling_method = scaling_method, omit_contrast_in_colnames = omit_contrast_in_colnames, prefix = prefix, formula_dispersion = formula_dispersion, ..., significance_threshold = significance_threshold, fill_missing_values = fill_missing_values, .contrasts = .contrasts, .abundance = .abundance)
  }
)

#' test_differential_expression
#'
#' @docType methods
#' @rdname test_differential_expression-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod(
  "test_differential_expression",
  "RangedSummarizedExperiment",
  function(.data, .formula, abundance = assayNames(.data)[1], contrasts = NULL, method = c("edgeR_quasi_likelihood", "edgeR_likelihood_ratio", "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"), test_above_log2_fold_change = NULL, scaling_method = "TMM", omit_contrast_in_colnames = FALSE, prefix = "", formula_dispersion = NULL, ..., significance_threshold = NULL, fill_missing_values = NULL, .contrasts = NULL, .abundance = NULL) {
    test_differential_abundance(.data, .formula, abundance = abundance, contrasts = contrasts, method = method, test_above_log2_fold_change = test_above_log2_fold_change, scaling_method = scaling_method, omit_contrast_in_colnames = omit_contrast_in_colnames, prefix = prefix, formula_dispersion = formula_dispersion, ..., significance_threshold = significance_threshold, fill_missing_values = fill_missing_values, .contrasts = .contrasts, .abundance = .abundance)
  }
)




#' Get differential transcription information to a tibble using edgeR.
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom purrr when
#' @importFrom magrittr extract2
#' @importFrom stringr str_detect str_replace
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param test_above_log2_fold_change A positive real value. This works for edgeR and limma_voom methods. It uses the `treat` function, which tests that the difference in abundance is bigger than this threshold rather than zero \url{https://pubmed.ncbi.nlm.nih.gov/19176553}.
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::normLibSizes; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#' @return A tibble with edgeR results
#'
get_differential_transcript_abundance_bulk_SE <- function(
    .data,
    .formula,
    abundance = assayNames(.data)[1],
    sample_annotation,
    .contrasts = NULL,
    method = c("edgeR_quasi_likelihood", "edgeR_likelihood_ratio", "edger_robust_likelihood_ratio", "DESeq2", "limma_voom", "limma_voom_sample_weights", "glmmseq_lme4", "glmmseq_glmmtmb"),
    test_above_log2_fold_change = NULL,
    scaling_method = "TMM",
    omit_contrast_in_colnames = FALSE,
    prefix = "",
    ...,
    .abundance = NULL
) {
  # Deprecation logic for .abundance
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_differential_transcript_abundance_bulk_SE(.abundance)", "get_differential_transcript_abundance_bulk_SE(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  my_assay <- abundance
  
  # Check if omit_contrast_in_colnames is correctly setup
  if(omit_contrast_in_colnames & length(.contrasts) > 1){
    warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
    omit_contrast_in_colnames = FALSE
  }
  
  # Create design matrix
  design =
    model.matrix(
      object = .formula,
      data = sample_annotation
    )
  
  # Replace `:` with ___ because it creates error with edgeR
  if(design |> colnames() |> str_detect(":") |> any()) {
    message("tidybulk says: the interaction term `:` has been replaced with `___` in the design matrix, in order to work with edgeR.")
    colnames(design) = design |> colnames() |> str_replace(":", "___")
  }
  
  # Print the design column names in case I want contrasts
  message(
    sprintf(
      "tidybulk says: The design column names are \"%s\"",
      design |> colnames() |> paste(collapse = ", ")
    )
  )
  
  if(length(.contrasts) > 0) {
    my_contrasts = limma::makeContrasts(contrasts = .contrasts, levels = design)
  } else {
    my_contrasts = NULL
  }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("edgeR")
  
  # If no assay is specified take first
  # if(abundance |> quo_is_symbol()) {
  #   abundance = quo_name(abundance)
  # } else {
  #   abundance = .data |>
  #     assayNames() |>
  #     extract2(1)
  # }	
  # my_assay = abundance
  
  edgeR_object =
    .data |>
    assay(my_assay) |>
    edgeR::DGEList() 
  
  if(scaling_method != "none")
    edgeR_object = edgeR_object |> edgeR::normLibSizes(method = scaling_method)
  
  

  method_lower <- tolower(method)
  if(method_lower %in% c("edger_likelihood_ratio", "edger_quasi_likelihood")) {
    edgeR_object <- edgeR_object |> edgeR::estimateDisp(design)
    fit_object <- if(method_lower == "edger_likelihood_ratio") {
      edgeR_object |> edgeR::glmFit(design)
    } else if(method_lower == "edger_quasi_likelihood") {
      edgeR_object |> edgeR::glmQLFit(design)
    }
  } else if(method_lower == "edger_robust_likelihood_ratio") {
    edgeR_object <- edgeR_object |> edgeR::estimateGLMRobustDisp(design)
    fit_object <- edgeR_object |> edgeR::glmFit(design)
  }
  # Return
  if(my_contrasts |> is.null() | omit_contrast_in_colnames)	{
    if(!is.null(test_above_log2_fold_change))
      fit_object = fit_object |> edgeR::glmTreat(coef = 2, contrast = my_contrasts, lfc=test_above_log2_fold_change)
    else if(tolower(method) %in%  c("edger_likelihood_ratio", "edger_robust_likelihood_ratio"))
      fit_object = fit_object |> edgeR::glmLRT(coef = 2, contrast = my_contrasts)
    else if(tolower(method) ==  "edger_quasi_likelihood")
      fit_object = fit_object |> edgeR::glmQLFTest(coef = 2, contrast = my_contrasts)
    else
      stop("tidybulk says: method not supported")	
    
    
    # Convert to tibble
    result = 
      fit_object |> 
      edgeR::topTags(n = Inf) %$%
      table |>
      as_tibble(rownames = "transcript") |>
      
      # # Mark DE genes
      # mutate(significant = FDR < significance_threshold) 	|>
      
      # Arrange
      arrange(FDR)
  }  else {
    
    result = 
      
      1:ncol(my_contrasts) |>
      map_dfr(function(contrast_index) {
        # Test one contrast at a time. Passing the whole matrix makes edgeR
        # run a joint F-test and repeat the same statistics for every contrast.
        this_contrast <- my_contrasts[, contrast_index]
        if(!is.null(test_above_log2_fold_change))
          fit_object = fit_object |> edgeR::glmTreat(coef = 2, contrast = this_contrast, lfc=test_above_log2_fold_change)
        else if(tolower(method) %in%  c("edger_likelihood_ratio", "edger_robust_likelihood_ratio"))
          fit_object = fit_object |> edgeR::glmLRT(coef = 2, contrast = this_contrast)
        else if(tolower(method) ==  "edger_quasi_likelihood")
          fit_object = fit_object |> edgeR::glmQLFTest(coef = 2, contrast = this_contrast)
        else
          stop("tidybulk says: method not supported")	
        
        fit_object |>
          
          
          # Convert to tibble
          edgeR::topTags(n = Inf) %$%
          table |>
          as_tibble(rownames = "transcript") |>
          mutate(constrast = colnames(my_contrasts)[contrast_index])
      }) |>
      pivot_wider(values_from = -c(transcript, constrast),
                  names_from = constrast, names_sep = "___")
  }
  
  if (exists("result") && !is.null(result)) {
    result <- result |>
      setNames(c(
        colnames(result)[1],
        sprintf("%s%s", prefix, colnames(result)[2:ncol(result)])
      ))
    
    list(
      result_raw = fit_object, 
      de_object = edgeR_object,
      result = result
    )
  } else {
    stop("tidybulk says: Internal error -- result object was not created in get_differential_transcript_abundance_bulk_SE.")
  }
}


#' Get differential transcription information to a tibble using voom.
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom purrr when
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See voom makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "limma_voom", "limma_voom_sample_weights"
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::normLibSizes; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param .abundance DEPRECATED. The name of the transcript/gene abundance column (symbolic, for backward compatibility)
#'
#' @return A tibble with voom results
#'
get_differential_transcript_abundance_bulk_voom_SE <- function(
    .data,
    .formula,
    abundance = assayNames(.data)[1],
    sample_annotation,
    .contrasts = NULL,
    method = "limma_voom",
    test_above_log2_fold_change = NULL,
    scaling_method = "TMM",
    omit_contrast_in_colnames = FALSE,
    prefix = "",
    ...,
    .abundance = NULL
) {
  # Deprecation logic for .abundance
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_differential_transcript_abundance_bulk_voom_SE(.abundance)", "get_differential_transcript_abundance_bulk_voom_SE(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  my_assay <- abundance
  
  # Check if omit_contrast_in_colnames is correctly setup
  if(omit_contrast_in_colnames & length(.contrasts) > 1){
    warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
    omit_contrast_in_colnames = FALSE
  }
  
  
  # Create design matrix
  design =
    model.matrix(
      object = .formula,
      data = sample_annotation
    )
  
  # Print the design column names in case I want contrasts
  message(
    sprintf(
      "tidybulk says: The design column names are \"%s\"",
      design |> colnames() |> paste(collapse = ", ")
    )
  )
  
  # Check if package is installed, otherwise install
  check_and_install_packages("limma")
  
  if(length(.contrasts) > 0) {
    my_contrasts =
      limma::makeContrasts(contrasts = .contrasts, levels = design)
  } else {
    my_contrasts = NULL
  }
  
  
  
  # If no assay is specified take first
  # my_assay = ifelse(
  #   abundance |> quo_is_symbol(),
  #   quo_name(abundance),
  #   .data |>
  #     assayNames() |>
  #     extract2(1)
  # )
  
  voom_object =
    .data |>
    
    assay(my_assay) |>
    edgeR::DGEList() 
  
  # Scale data if method is not "none"
  # use if else instead of when
  if(scaling_method != "none")
    voom_object = voom_object |> edgeR::normLibSizes(method = scaling_method)
  
  
  if(tolower(method) == "limma_voom")
    voom_object = voom_object |> limma::voom(design, plot=FALSE)
  else if(tolower(method) == "limma_voom_sample_weights")
    voom_object = voom_object |> limma::voomWithQualityWeights(design, plot=FALSE)
  else
    stop("tidybulk says: method not supported")
  
  # select method
  result = 
    voom_object |>
    
    limma::lmFit(design)
  
  # Return
  
  
  if(my_contrasts |> is.null() | omit_contrast_in_colnames) {
    result = result |>
      
      # Contrasts
      limma::contrasts.fit(contrasts=my_contrasts, coefficients =  when(my_contrasts, is.null(.) ~ 2)) |>
      limma::eBayes() 
    
    if(is.null(test_above_log2_fold_change)) {
      result = result |> limma::topTable(n = Inf)
    } else {
      result = result |> limma::treat(lfc=test_above_log2_fold_change) |> limma::topTreat(n = Inf)
    }
    
    
    result = result |>
      # Convert to tibble
      as_tibble(rownames = "transcript") |>
      
      # # Mark DE genes
      # mutate(significant = adj.P.Val < significance_threshold) 	|>
      
      # Arrange
      arrange(adj.P.Val)		
  } else {
    
    
    result = 
      1:ncol(my_contrasts) |>
      map_dfr(
        function(contrast_index) {
          
          
          result = result |>
            
            # Contrasts
            limma::contrasts.fit(contrasts=my_contrasts[, contrast_index]) |>
            limma::eBayes() 
          
          if(is.null(test_above_log2_fold_change)) {
            result = result |> limma::topTable(n = Inf)
          } else {
            result = result |> limma::treat(lfc=test_above_log2_fold_change) |> limma::topTreat(n = Inf)
          }
          
          result = result |>
            
            # Convert to tibble
            as_tibble(rownames = "transcript") |>
            mutate(constrast = colnames(my_contrasts)[contrast_index])

        }) |>
      pivot_wider(values_from = -c(transcript, constrast),
                  names_from = constrast, names_sep = "___")
    
  }
  
  
  result =
    # Attach prefix
    result |> setNames(c(
      colnames(result)[1],
      sprintf("%s%s", prefix, colnames(result)[2:ncol(result)])
    ))
  
  list(
    result_raw = result,
    de_object = voom_object,
    result = result
  )
  
}


#' Get differential transcription information to a tibble using glmmSeq
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom purrr when
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::normLibSizes; "TMM","TMMwsp","RLE","upperquartile")
#' @param .scaling_factor A tidyeval (column name) for the precalculated TMM scaling
#' @param formula_dispersion A fixed-effects formula for edgeR tagwise dispersion. If NULL, each gene estimates phi.
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param ... Additional arguments for glmmSeq
#'
#' @return A tibble with glmmSeq results
#'
get_differential_transcript_abundance_glmmSeq_SE <- function(
    .data,
    .formula,
    abundance = assayNames(.data)[1],
    .contrasts = NULL,
    sample_annotation,
    method,
    test_above_log2_fold_change = NULL,
    scaling_method = "TMM",
    .scaling_factor = NULL,
    omit_contrast_in_colnames = FALSE,
    prefix = "",
    formula_dispersion = NULL,
    ...,
    .abundance = NULL
) {
  
  .scaling_factor = enquo(.scaling_factor)
  
  # Deprecation logic for .abundance
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_differential_transcript_abundance_glmmSeq_SE(.abundance)", "get_differential_transcript_abundance_glmmSeq_SE(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  my_assay <- abundance
  
  # Check if contrasts are of the same form
  if(
    .contrasts |> is.null() |> not() &
    .contrasts |> class() |> equals("list") |> not()
  )
    stop("tidybulk says: for DESeq2 the list of constrasts should be given in the form list(c(\"condition_column\",\"condition1\",\"condition2\")) i.e. list(c(\"genotype\",\"knockout\",\"wildtype\"))")
  
  # Check if omit_contrast_in_colnames is correctly setup
  if(omit_contrast_in_colnames & length(.contrasts) > 1){
    warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
    omit_contrast_in_colnames = FALSE
  }
  
  # # Check if package is installed, otherwise install
  # if (find.package("edgeR", quiet = TRUE) |> length() |> equals(0)) {
  #   message("tidybulk says: Installing edgeR needed for differential transcript abundance analyses")
  #   if (!requireNamespace("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager", repos = "https://cloud.r-project.org")
  #   BiocManager::install("edgeR", ask = FALSE)
  # }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("glmmSeq")
  
  metadata =
    .data |>
    colData()
  
  counts =
    .data |>
    assay(my_assay)
  
  if (!is.null(formula_dispersion)) {
    .data <- estimate_dispersion(.data, formula_dispersion, my_assay)
    message(
      "tidybulk says: plugging rowData column dispersion_shrinked into glmmSeq."
    )
    dispersion <- stats::setNames(
      as.numeric(rowData(.data)$dispersion_shrinked),
      rownames(.data)
    )[rownames(counts)]
  } else {
    message(
      "tidybulk says: no formula_dispersion given; ",
      "glmmSeq will estimate phi for each gene. ",
      "This is only recommended when the degrees of freedom ",
      "(i.e. sample size - parameters) is very high (e.g. > 50)."
    )
    dispersion <- NA
  }
  
  # Scaling
  if(.scaling_factor |> quo_is_symbolic())
    sizeFactors = .data |> pivot_sample() |> pull(!!.scaling_factor)
  else
    sizeFactors <- counts |> edgeR::normLibSizes(method = scaling_method)
  
  
  glmmSeq_object =
    glmmSeq( .formula,
             countdata = counts ,
             metadata =   metadata |> as.data.frame(),
             dispersion = dispersion,
             sizeFactors = sizeFactors,
             progress = TRUE,
             method = method |> str_remove("(?i)^glmmSeq_" ),
             ...
    )
  
  result_df <- glmmSeq_object |>
    summary_lmmSeq() |>
    as_tibble(rownames = "transcript") |>
    mutate(across(starts_with("P_"), list(adjusted = function(x) p.adjust(x, method="BH")), .names = "{.col}_{.fn}")) |>
    
    # Attach attributes
    reattach_metadata(.data) |>
    
    # select method
    memorise_methods_used("glmmSeq")
  
  # Attach prefix
  result_df <- result_df |>
    setNames(c(
      colnames(result_df)[1],
      sprintf("%s%s", prefix, colnames(result_df)[2:ncol(result_df)])
    ))
  
  list(result = result_df, result_raw = glmmSeq_object, de_object = glmmSeq_object)
  
  
}


#' Get differential transcription information to a tibble using DESeq2
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom magrittr set_colnames
#' @importFrom stats model.matrix
#' @importFrom purrr when
#'
#'
#' @param .data A tibble
#' @param .formula a formula with no response variable, referring only to numeric variables
#' @param .contrasts A character vector. See edgeR makeContrasts specification for the parameter `contrasts`. If contrasts are not present the first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param method A string character. Either "edgeR_quasi_likelihood" (i.e., QLF), "edgeR_likelihood_ratio" (i.e., LRT)
#' @param scaling_method A character string. The scaling method passed to the backend function (i.e., edgeR::normLibSizes; "TMM","TMMwsp","RLE","upperquartile")
#' @param omit_contrast_in_colnames If just one contrast is specified you can choose to omit the contrast label in the colnames.
#' @param ... Additional arguments for DESeq2
#'
#' @return A tibble with DESeq2 results
#'
get_differential_transcript_abundance_deseq2_SE <- function(.data,
                                                            .formula,
                                                            abundance = assayNames(.data)[1],
                                                            
                                                            .abundance = NULL,
                                                            .contrasts = NULL,
                                                            method = "deseq2",
                                                            
                                                            test_above_log2_fold_change = NULL,
                                                            
                                                            scaling_method = "TMM",
                                                            omit_contrast_in_colnames = FALSE,
                                                            prefix = "",
                                                            ...) {
  
  # Deprecation logic for .abundance (symbolic)
  if (!is.null(.abundance)) {
    lifecycle::deprecate_warn("2.0.0", "get_differential_transcript_abundance_deseq2_SE(.abundance)", "get_differential_transcript_abundance_deseq2_SE(abundance)")
    if (missing(abundance) || is.null(abundance)) {
      abundance <- rlang::as_name(rlang::ensym(.abundance))
    }
  }
  
  .abundance <- rlang::enquo(.abundance)
  
  # Fix NOTEs
  . = NULL
  pvalue = NULL
  padj = NULL
  
  # Check if contrasts are of the same form
  if(
    .contrasts |> is.null() |> not() &
    .contrasts |> class() |> equals("list") |> not()
  )
    stop("tidybulk says: for DESeq2 the list of constrasts should be given in the form list(c(\"condition_column\",\"condition1\",\"condition2\")) i.e. list(c(\"genotype\",\"knockout\",\"wildtype\"))")
  
  # Check if omit_contrast_in_colnames is correctly setup
  if(omit_contrast_in_colnames & length(.contrasts) > 1){
    warning("tidybulk says: you can omit contrasts in column names only when maximum one contrast is present")
    omit_contrast_in_colnames = FALSE
  }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("DESeq2")
  
  
  if (is.null(test_above_log2_fold_change)) {
    test_above_log2_fold_change <- 0
  }
  
  my_contrasts = .contrasts
  
  # If no assay is specified take first
  my_assay <- abundance
  
  deseq2_object =
    
    # DESeq2
    DESeq2::DESeqDataSetFromMatrix(
      countData = .data |> assay(my_assay),
      colData = colData(.data),
      design = .formula
    ) |>
    DESeq2::DESeq(...)
  

      # Simplified logic, no anonymous function, broken up for clarity

      has_no_contrasts <- is.null(my_contrasts)
      is_continuous <- class(deseq2_object@colData[, parse_formula(.formula)[1]]) %in% c("numeric", "integer", "double")
      should_omit_contrast_names <- !is.null(my_contrasts) & omit_contrast_in_colnames

      if (has_no_contrasts & is_continuous) {
        # Simple comparison continuous
        result <- deseq2_object |>
          DESeq2::results(lfcThreshold = test_above_log2_fold_change) |>
          as_tibble(rownames = "transcript")
      } else if (has_no_contrasts) {
        # Simple comparison discrete
        factor_levels <- deseq2_object@colData[, parse_formula(.formula)[1]] |> as.factor() |> levels()
        result <- deseq2_object |>
          DESeq2::results(
            contrast = c(
              parse_formula(.formula)[1],
              factor_levels[2],
              factor_levels[1]
            ),
            lfcThreshold = test_above_log2_fold_change
          ) |>
          as_tibble(rownames = "transcript")
      } else if (should_omit_contrast_names) {
        # Single contrast, omit contrast names in columns
        result <- deseq2_object |>
          DESeq2::results(contrast = my_contrasts[[1]], lfcThreshold = test_above_log2_fold_change) |>
          as_tibble(rownames = "transcript")
      } else {
        # Multiple contrasts
        result <- 1:length(my_contrasts) |>
          map_dfr(function(contrast_index) {
            deseq2_object |>
              DESeq2::results(contrast = my_contrasts[[contrast_index]], lfcThreshold = test_above_log2_fold_change) |>
              as_tibble(rownames = "transcript") |>
              mutate(constrast = sprintf("%s %s-%s", my_contrasts[[contrast_index]][1], my_contrasts[[contrast_index]][2], my_contrasts[[contrast_index]][3]))
          }) |>
          pivot_wider(
            values_from = -c(transcript, constrast),
            names_from = constrast,
            names_sep = "___"
          )
      }
      
      # Attach prefix without using pipe
      cn <- colnames(result)
      cn_new <- c(cn[1], sprintf("%s%s", prefix, cn[2:length(cn)]))
      colnames(result) <- cn_new

  # Return
  list(
    result_raw = result,
    de_object = deseq2_object
  )
  
  
}
