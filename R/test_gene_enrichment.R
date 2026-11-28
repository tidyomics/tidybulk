#' analyse gene enrichment with EGSEA
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description test_gene_enrichment() takes as input a `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a `tbl` of gene set information
#'
#' @importFrom rlang enquo quo_name
#' @importFrom magrittr not
#' @importFrom dplyr filter arrange mutate pull bind_rows
#' @importFrom SummarizedExperiment colData rowData assays
#' @importFrom tidyselect all_of

#'
#'
#' @name test_gene_enrichment
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .entrez The ENTREZ ID of the transcripts/genes
#' @param .abundance The name of the transcript/gene abundance column
#' @param contrasts This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#' @param methods A character vector. One or 3 or more methods to use in the testing (currently EGSEA errors if 2 are used). Type EGSEA::egsea.base() to see the supported GSE methods.
#' @param gene_sets A character vector or a list. It can take one or more of the following built-in collections as a character vector: c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"), to be used with EGSEA buildIdx. c1 is human specific. Alternatively, a list of user-supplied gene sets can be provided, to be used with EGSEA buildCustomIdx. In that case, each gene set is a character vector of Entrez IDs and the names of the list are the gene set names.
#' @param species A character. It can be human, mouse or rat.
#' @param cores An integer. The number of cores available
#'
#' @param method DEPRECATED. Please use methods.
#' @param .contrasts DEPRECATED - This parameter takes the format of the contrast parameter of the method of choice. For edgeR and limma-voom is a character vector. For DESeq2 is a list including a character vector of length three. The first covariate is the one the model is tested against (e.g., ~ factor_of_interest)
#'
#' @details This wrapper executes ensemble gene enrichment analyses of the dataset using EGSEA (DOI:0.12688/f1000research.12544.1)
#'
#'
#' dge =
#' 	data |>
#' 	keep_abundant(
#' 		factor_of_interest = !!as.symbol(parse_formula(.formula)[[1]]),
#' 		!!.sample, !!.entrez, !!.abundance
#' 	) |>
#'
#' 	# Make sure transcript names are adjacent
#' 	[...] |>
#' 	as_matrix(rownames = !!.entrez) |>
#' 	edgeR::DGEList(counts = .)
#'
#' idx =  buildIdx(entrezIDs = rownames(dge), species = species, msigdb.gsets = msigdb.gsets,
#'	               kegg.exclude = kegg.exclude)
#'
#' dge |>
#'
#' 	# Calculate weights
#' 	limma::voom(design, plot = FALSE) |>
#'
#' 	# Execute EGSEA
#' 	egsea(
#' 		contrasts = my_contrasts,
#' 		baseGSEAs = methods,
#' 		gs.annots = idx,
#' 		sort.by = "med.rank",
#' 		num.threads = cores,
#' 		report = FALSE
#' 	)
#'
#' @return A consistent object (to the input)
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
#' library(tidySummarizedExperiment)
#'
#' \dontrun{
#' library("EGSEA")
#'
#' df_entrez <- airway |>
#'   mutate(entrez = .feature) |>
#'   aggregate_duplicates(.transcript = entrez )
#'
#' 	test_gene_enrichment(
#'			df_entrez,
#'			~ condition,
#'			.entrez = entrez,
#'			.abundance = count,
#'			methods = c("roast" , "safe", "gage"  ,  "padog" , "globaltest", "ora" ),
#'			gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
#'			species="human",
#'			cores = 2
#'		)
#' }
#'
#'
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Alhamdoosh, M., Ng, M., Wilson, N. J., Sheridan, J. M., Huynh, H., Wilson, M. J., & Ritchie, M. E. (2017). Combining multiple tools outperforms individual methods for gene set enrichment analysis in single-cell RNA-seq data. Genome Biology, 18(1), 174. doi:10.1186/s13059-017-1279-y
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#' @export
#'
#'
setGeneric("test_gene_enrichment", function(.data,
                                            .formula,
                                            
                                            .entrez,
                                            .abundance = NULL,
                                            contrasts = NULL,
                                            methods = c("camera" ,    "roast" ,     "safe",       "gage"  ,     "padog" ,     "globaltest",  "ora" ),
                                            gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
                                            species,
                                            cores = parallel::detectCores(),
                                            
                                            # DEPRECATED
                                            method = NULL,
                                            .contrasts = NULL
)
standardGeneric("test_gene_enrichment"))





#' @importFrom lifecycle deprecate_warn
#' @importFrom stringr str_replace
#' @importFrom dplyr everything
#'
#'
#'
.test_gene_enrichment_SE = 		function(.data,
                                      .formula,
                                      
                                      .entrez,
                                      .abundance = NULL,
                                      contrasts = NULL,
                                      methods = c("camera" ,    "roast" ,     "safe",       "gage"  ,     "padog" ,     "globaltest",  "ora" ),
                                      gene_sets = c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "kegg_disease", "kegg_metabolism", "kegg_signaling"),
                                      species,
                                      cores = parallel::detectCores(),
                                      
                                      # DEPRECATED
                                      method = NULL,
                                      .contrasts = NULL
)	{
  
  # Fix NOTEs
  . = NULL
  
  # DEPRECATION OF reference function
  if (is_present(method) & !is.null(method)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.3.2", "tidybulk::test_gene_enrichment(method = )", details = "The argument method is now deprecated please use methods")
    methods = method
  }
  
  # DEPRECATION OF .constrasts
  if (is_present(.contrasts) & !is.null(.contrasts)) {
    
    # Signal the deprecation to the user
    deprecate_warn("1.7.4", "tidybulk::test_differential_abundance(.contrasts = )", details = "The argument .contrasts is now deprecated please use contrasts (without the dot).")
    
    contrasts = .contrasts
  }
  
  .entrez = enquo(.entrez)
  
  # Check that there are no entrez missing
  entrez_col_name <- quo_name(.entrez)
  if (entrez_col_name %in% colnames(rowData(.data))) {
    entrez_values <- rowData(.data)[, entrez_col_name]
    has_na_entrez <- any(is.na(entrez_values))
    
    if (has_na_entrez) {
      warning("tidybulk says: There are NA entrez IDs. Those genes will be filtered")
      # Filter out NAs using row indexing instead of dplyr filter
      na_rows <- is.na(entrez_values)
      .data <- .data[!na_rows, ]
    }
  } else {
    stop("tidybulk says: the .entrez parameter appears to not be set")
  }
  
  # Check if duplicated entrez
  if(rowData(.data)[,quo_name(.entrez)] |> duplicated() |> any())
    stop("tidybulk says: There are duplicated .entrez IDs. Please use aggregate_duplicates(.transcript = entrez).")
  
  # For use within when
  .my_data = .data
  
  
  # Comply with CRAN NOTES
  . = NULL
  
  # Check if at least two samples for each group
  if (.data |>
      pivot_sample() |>
      count(!!as.symbol(parse_formula(.formula))) |>
      distinct(n) |>
      pull(n) |>
      min() |>
      st(2))
    stop("tidybulk says: You need at least two replicates for each condition for EGSEA to work")
  
  
  # Create design matrix
  design =	model.matrix(	object = .formula,	data = .data |> colData() 	)
  
  # Print the design column names in case I want contrasts
  message(
    sprintf(
      "tidybulk says: The design column names are \"%s\"",
      design |> colnames() |> paste(collapse = ", ")
    )
  )
  
  if (length(contrasts) > 0) {
    my_contrasts <- limma::makeContrasts(contrasts = contrasts, levels = design)
  } else {
    my_contrasts <- NULL
  }
  
  # Check if package is installed, otherwise install
  check_and_install_packages("EGSEA")
  
  if (!"EGSEA" %in% (.packages())) {
    stop("EGSEA package not loaded. Please run library(\"EGSEA\"). With this setup, EGSEA require manual loading, for technical reasons.")
  }
  
  # Extract assay data
  assay_data <- .data |>
    assays() |>
    as.list() |>
    _[[1]] |>
    as.matrix()
  
  # Change rownames to entrez
  if (!quo_is_null(.entrez)) {
    rownames(assay_data) <-
      .my_data |>
      pivot_transcript() |>
      pull(!!.entrez)
  }
  
  # Filter missing entrez
  valid_rows <- !is.na(rownames(assay_data))
  assay_data <- assay_data[valid_rows, , drop = FALSE]
  
  # Check if we have data
  if (nrow(assay_data) == 0) {
    stop("tidybulk says: No data remaining after filtering")
  }
  
  dge = edgeR::DGEList(counts = assay_data)
  
  # Add gene ids for Interpret Results tables in report
  dge$genes = rownames(dge$counts)
  
  if (is.list(gene_sets)) {
    
    idx =  buildCustomIdx(geneIDs = rownames(dge), species = species, gsets=gene_sets)
    nonkegg_genesets = idx
    kegg_genesets = NULL
    
  } else {
    
    # Specify gene sets to include
    msig_all <- c("h", "c1", "c2", "c3", "c4", "c5", "c6", "c7")
    kegg_all <- c("kegg_disease", "kegg_metabolism", "kegg_signaling")
    
    # Record which collections used (kegg, msigdb) for bibliography
    collections_bib = c()
    
    # Identify any msigdb sets to be included
    msigdb.gsets <- gene_sets[gene_sets %in% msig_all]
    if (length(msigdb.gsets) >= 1) {
      collections_bib = c(collections_bib, "msigdb")
    }
    
    # Have to identify kegg sets to exclude for EGSEA
    kegg_to_exclude = kegg_all[!(kegg_all %in% gene_sets)]
    
    # If all 3 kegg sets are excluded then set to "all" as specifying the 3 names gives empty kegg object
    if (length(kegg_to_exclude) == 3) {
      kegg.exclude = "all"
    } else {
      kegg.exclude = kegg_to_exclude |> str_replace("kegg_", "")
      collections_bib = c(collections_bib, "kegg")
    }
    
    
    idx =  buildIdx(entrezIDs = rownames(dge), species = species,  msigdb.gsets = msigdb.gsets,
                    kegg.exclude = kegg.exclude)
    
    # Due to a bug with kegg pathview overlays, this collection is run without report
    # https://support.bioconductor.org/p/122172/#122218
    
    kegg_genesets = idx[which(names(idx)=="kegg")]
    nonkegg_genesets = idx[which(names(idx)!="kegg")]
  }
  
  # Specify column to use to sort results in output table
  # If only one method is specified there is no med.rank column
  if (length(methods) == 1) {
    sort_column = "p.value"
  } else {
    sort_column = "med.rank"
  }
  
  
  if (length(nonkegg_genesets) != 0) {
    res =
      dge |>
      
      # Calculate weights
      limma::voom(design, plot = FALSE) |>
      
      # Execute EGSEA
      egsea(
        contrasts = my_contrasts,
        gs.annots = nonkegg_genesets,
        baseGSEAs = methods,
        sort.by = sort_column,
        num.threads = cores
      )
    
    gsea_web_page = "https://www.gsea-msigdb.org/gsea/msigdb/cards/%s.html"
    
    res_formatted_nonkegg =
      res@results |>
      map2_dfr(
        (.) |> names(),
        ~ .x[[1]][[1]] |>
          as_tibble(rownames = "pathway") |>
          mutate(data_base = .y)
      ) |>
      arrange(sort_column) |>
      
      # Ensure required columns exist
      mutate(
        data_base = if("data_base" %in% names(.)) data_base else NA_character_,
        pathway = if("pathway" %in% names(.)) pathway else NA_character_,
        web_page = if("pathway" %in% names(.) && nrow(.) > 0) sprintf(gsea_web_page, pathway) else NA_character_
      ) |>
      # Select columns that exist
      (function(data) {
        available_cols <- names(data)
        select_cols <- c("data_base", "pathway", "web_page", sort_column)
        existing_cols <- select_cols[select_cols %in% available_cols]
        if (length(existing_cols) > 0) {
          select(data, all_of(existing_cols), everything())
        } else {
          data
        }
      })()
  }
  
  if (length(kegg_genesets) != 0) {
    message("tidybulk says: due to a bug in the call to KEGG database (http://supportupgrade.bioconductor.org/p/122172/#122218), the analysis for this database is run without report production.")
    
    res_kegg =
      dge |>
      
      # Calculate weights
      limma::voom(design, plot = FALSE) |>
      
      # Execute EGSEA
      egsea(
        contrasts = my_contrasts,
        gs.annots = kegg_genesets,
        baseGSEAs = methods,
        sort.by = sort_column,
        num.threads = cores,
        report = FALSE
      )
    
    res_formatted_kegg =
      res_kegg@results |>
      map2_dfr(
        (.) |> names(),
        ~ .x[[1]][[1]] |>
          as_tibble(rownames = "pathway") |>
          mutate(data_base = .y)
      ) |>
      arrange(sort_column) |>
      # Ensure required columns exist
      mutate(
        data_base = if("data_base" %in% names(.)) data_base else NA_character_,
        pathway = if("pathway" %in% names(.)) pathway else NA_character_,
        web_page = if("pathway" %in% names(.) && nrow(.) > 0) sprintf(gsea_web_page, pathway) else NA_character_
      ) |>
      # Select columns that exist
      (function(data) {
        available_cols <- names(data)
        select_cols <- c("data_base", "pathway", "web_page", sort_column)
        existing_cols <- select_cols[select_cols %in% available_cols]
        if (length(existing_cols) > 0) {
          select(data, all_of(existing_cols), everything())
        } else {
          data
        }
      })()
    
  }
  
  # output tibble
  if (exists("res_formatted_nonkegg") & exists("res_formatted_kegg")) {
    out = bind_rows(res_formatted_nonkegg, res_formatted_kegg)
  } else if (exists("res_formatted_nonkegg")) {
    out = res_formatted_nonkegg
  } else {
    out = res_formatted_kegg
  }
  
  # add to bibliography
  if (exists("collections_bib")) {
    out |> memorise_methods_used(c("egsea", collections_bib, methods))
  }
  
}

#' test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
          "SummarizedExperiment",
          .test_gene_enrichment_SE)

#' test_gene_enrichment
#'
#' @docType methods
#' @rdname test_gene_enrichment-methods
#'
#' @return A consistent object (to the input)
setMethod("test_gene_enrichment",
          "RangedSummarizedExperiment",
          .test_gene_enrichment_SE)
