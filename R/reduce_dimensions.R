#' Dimension reduction of the transcript abundance data
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description reduce_dimensions() takes as input a `SummarizedExperiment` and calculates the reduced dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo quo_name
#' @importFrom magrittr not
#' @importFrom dplyr filter distinct select mutate rename
#' @importFrom tidyr pivot_wider
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData rowData assays
#' @importFrom stats prcomp
#'
#'
#' @name reduce_dimensions
#'
#' @param .data A `SummarizedExperiment`
#' @param .abundance The name of the assay to use for dimension reduction. This argument must be explicitly specified.
#'
#' @param method A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param .dims An integer. The number of dimensions your are interested in (e.g., 4 for returning the first four principal components).
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE", or uwot::tumap if you choose method="umap"
#'
#' @param log_transform DEPRECATED - A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS; DOI.org/10.1186/gb-2010-11-3-r25),
#' principal component analysis (PCA), or tSNE (Jesse Krijthe et al. 2018)
#'
#' Underlying method for PCA:
#' prcomp(scale = scale, ...)
#'
#' Underlying method for MDS:
#' limma::plotMDS(ndim = .dims, plot = FALSE, top = top)
#'
#' Underlying method for tSNE:
#' Rtsne::Rtsne(data, ...)
#'
#' Underlying method for UMAP:
#'
#'  df_source =
#' .data |>
#'
#'   # Filter NA symbol
#'   filter(!!.feature |> is.na() |> not()) |>
#'
#'   # Prepare data frame
#'   distinct(!!.feature,!!.element,!!.abundance) |>
#'
#'   # Filter most variable genes
#'   keep_variable_transcripts(top) |>
#'   reduce_dimensions(.abundance = counts, method="PCA",  .dims = calculate_for_pca_dimensions ) |>
#'   as_matrix(rownames = quo_name(.element)) |>
#'   uwot::tumap(...)
#'
#'
#' @return A tbl object with additional columns for the reduced dimensions
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
#'
#' counts.MDS =
#'  airway |>
#'  identify_abundant() |>
#'  reduce_dimensions(.abundance = counts, method="MDS", .dims = 3)
#'
#'
#' counts.PCA =
#'  airway |>
#'  identify_abundant() |>
#'  reduce_dimensions(.abundance = counts, method="PCA", .dims = 3)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' Krijthe, J. H. (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation. R package version 0.15.
#'
#' McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv preprint arXiv:1802.03426.
#'
#'
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#' @export
#'
#'
setGeneric("reduce_dimensions", function(.data,
                                         .abundance,
                                         method,
                                         .dims = 2,
                                         
                                         top = 500,
                                         of_samples = TRUE,
                                         transform = log1p,
                                         scale = TRUE,
                                         ...,
                                         
                                         # DEPRECATED
                                         log_transform = NULL
                                         
)
standardGeneric("reduce_dimensions"))


.reduce_dimensions_se = function(.data,
                                 .abundance,
                                 
                                 method,
                                 .dims = 2,
                                 top = 500,
                                 of_samples = TRUE,
                                 transform = log1p,
                                 scale = TRUE,
                                 ...) {
  
  # Fix NOTEs
  . = NULL
  
  if(missing(.abundance))
    stop("tidybulk says: please specify `.abundance` explicitly (e.g. `.abundance = counts_scaled`). If needed, create an assay scaled proportionally to library size with `scale_abundance()` and pass that assay via `.abundance`.")
  
  .abundance = enquo(.abundance)
  
  if(.abundance |> quo_is_symbolic()) my_assay = quo_name(.abundance)
  else stop("tidybulk says: please specify `.abundance` explicitly (e.g. `.abundance = counts_scaled`). If needed, create an assay scaled proportionally to library size with `scale_abundance()` and pass that assay via `.abundance`.")
  
  # adjust top for the max number of features I have
  if(top > nrow(.data)){
    warning(sprintf(
      "tidybulk says: the \"top\" argument %s is higher than the number of features %s",
      top,
      nrow(.data)
    ))
    
    top = min(top, nrow(.data))
  }
  
  my_assay =
    .data |>
    
    # Filter abundant if performed
    filter_if_abundant_were_identified() |>
    
    assay(my_assay) |>
    
    # Filter most variable genes
    keep_variable_transcripts_SE(top = top, transform = transform) |>
    
    # Check if log transform is needed
    transform()
  
  my_reduction_function  =
    if (tolower(method) == tolower("MDS")) {
      get_reduced_dimensions_MDS_bulk_SE
    } else if (tolower(method) == tolower("PCA")) {
      get_reduced_dimensions_PCA_bulk_SE
    } else if (tolower(method) == tolower("tSNE")) {
      get_reduced_dimensions_TSNE_bulk_SE
    } else if (tolower(method) == tolower("UMAP")) {
      get_reduced_dimensions_UMAP_bulk_SE
    } else {
      stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\", or \"UMAP\" ")
    }
  
  # Both dataframe and raw result object are returned
  reduced_dimensions =
    my_reduction_function(
      my_assay,
      .dims = .dims,
      top = top,
      of_samples = of_samples,
      transform = transform,
      scale=scale,
      ...
    )
  
  .data |>
    
    # Add dimensions to metadata
    (\(.) {
      .x = (.)
      if (of_samples) {
        # Ensure reduced_dimensions$result has same number of rows as colData
        result_df = reduced_dimensions$result
        if (!is.null(result_df) && is.data.frame(result_df) && nrow(result_df) > 0) {
          if (nrow(result_df) != nrow(colData(.x))) {
            message(sprintf("Row count mismatch: result_df (%d), colData (%d)", nrow(result_df), nrow(colData(.x))))
            # Align rows by sample names
            sample_names = rownames(colData(.x))
            result_df = result_df[match(sample_names, rownames(result_df)), , drop = FALSE]
          }
          if (nrow(result_df) == nrow(colData(.x))) {
            # Only remove the sample column if it matches rownames or is named 'sample'
            col_names = colnames(result_df)
            if (tolower(col_names[1]) == "sample" || all(rownames(colData(.x)) == result_df[[1]])) {
              dim_cols = col_names[-1]
              colData(.x) = cbind(colData(.x), result_df[, dim_cols, drop = FALSE])
            } else {
              colData(.x) = cbind(colData(.x), result_df)
            }
          } else {
            warning("Skipping cbind: row counts still do not match after alignment.")
          }
        }
      } else {
        rowData(.x) = rowData(.x) |> cbind(reduced_dimensions$result)
      }
      .x
    })() |>
    
    # Add bibliography
    (\(.) {
      if (tolower(method) == tolower("MDS")) {
        memorise_methods_used(., "limma")
      } else if (tolower(method) == tolower("PCA")) {
        memorise_methods_used(., "stats")
      } else if (tolower(method) == tolower("tSNE")) {
        memorise_methods_used(., "rtsne")
      } else if (tolower(method) == tolower("UMAP")) {
        memorise_methods_used(., "uwot")
      } else {
        stop("tidybulk says: method must be either \"MDS\" or \"PCA\" or \"tSNE\", or \"UMAP\" ")
      }
    })() |>
    
    # Attach edgeR for keep variable filtering
    memorise_methods_used(c("edger")) |>
    # Add raw object
    attach_to_metadata(reduced_dimensions$raw_result, method) |>
    
    # Communicate the attribute added
    (\(.) {
      
      rlang::inform(sprintf("tidybulk says: to access the raw results do `metadata(.)$tidybulk$%s`", method), .frequency_id = sprintf("Access %s results", method),  .frequency = "always")
      
      (.)
    })()
  
  
}

#' reduce_dimensions
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("reduce_dimensions",
          "SummarizedExperiment",
          .reduce_dimensions_se)

#' reduce_dimensions
#'
#' @docType methods
#' @rdname reduce_dimensions-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("reduce_dimensions",
          "RangedSummarizedExperiment",
          .reduce_dimensions_se)




#' Get dimensionality information to a tibble using MDS
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_MDS_bulk_SE <-
  function(.data,
           .dims = 2,
           top = 500,
           of_samples = TRUE,
           transform = log1p,
           scale = NULL # This is only a dummy argument for making it compatibble with PCA
  ) {
    # Comply with CRAN NOTES
    . = NULL
    Component = NULL
    `Component value` = NULL
    
    # Get components from dims
    components = 1:.dims
    
    # For MDS, we need to create pairs of dimensions (x,y)
    # If .dims is odd, we need to add one more to make pairs
    if((length(components) %% 2) != 0 ) {
      components = components |> append(components[length(components)])
    }
    # Ensure we only create exactly .dims dimensions
    if (.dims == 2) {
      components = c(1, 2)
      components_list = list(c(1, 2))
    } else {
      # For MDS, we need to create pairs of dimensions for plotting
      # If .dims is odd, we'll use the last dimension with the first
      if (.dims %% 2 == 1) {
        # For odd numbers, create pairs and handle the last dimension separately
        pairs = split(components[1:(.dims-1)], ceiling(seq_along(components[1:(.dims-1)])/2))
        # Add the last dimension paired with the first
        pairs[[length(pairs) + 1]] = c(components[.dims], components[1])
        components_list = pairs
      } else {
        # For even numbers, create pairs normally
        components_list = split(components, ceiling(seq_along(components)/2))
      }
    }
    
    # Loop over components list and calculate MDS. (I have to make this process more elegant)
    mds_object =
      components_list |>
      map(
        ~ .data |>
          limma::plotMDS(dim.plot = .x, plot = FALSE, top = top)
      )
    
    # Check if mds_object is NULL or contains NULL elements
    if (is.null(mds_object) || any(sapply(mds_object, is.null))) {
      stop("MDS calculation failed - check your data dimensions")
    }
    
    # Return
    list(
      raw_result = mds_object,
      result =
        map2_dfr(
          mds_object, components_list,
          ~ {
            
            # Use temporary variables for clarity
            mds_result = .x
            component_pair = .y
            
            # Extract sample names from the MDS result
            if ("distance.matrix.squared" %in% names(mds_result)) {
              my_rownames = rownames(mds_result$distance.matrix.squared)
            } else if ("distance.matrix" %in% names(mds_result)) {
              my_rownames = rownames(mds_result$distance.matrix)
            } else {
              # Fallback: use names from x component
              my_rownames = names(mds_result$x)
            }
            
            tibble(my_rownames, mds_result$x, mds_result$y) |>
              rename(
                sample := my_rownames,
                !!as.symbol(component_pair[1]) := `mds_result$x`,
                !!as.symbol(component_pair[2]) := `mds_result$y`
              ) |> gather(Component, `Component value`,-sample)
            
          }
          
          
        )  |>
        distinct() |>
        pivot_wider(names_from = Component, values_from = `Component value`) |>
        (\(.) {
          if (is.null(.)) {
            stop("MDS result is NULL - check data dimensions")
          }
          # Use temporary variables for clarity
          result_df = .
          col_names = colnames(result_df)
          print(sprintf("MDS result_df colnames: %s", paste(col_names, collapse=", ")))
          sample_col = col_names[1]
          dim_cols = col_names[-1]
          # Only keep the first .dims dimensions and ensure we have exactly .dims
          if (length(dim_cols) > .dims) {
            dim_cols = dim_cols[1:.dims]
          }
          new_names = paste0("Dim", 1:length(dim_cols))
          # Ensure we only return exactly .dims dimensions
          final_result = setNames(result_df[, c(sample_col, dim_cols)], c(sample_col, new_names))
          # Double-check that we have the right number of columns
          if (ncol(final_result) != .dims + 1) {
            stop(sprintf("MDS result has wrong number of columns: expected %d, got %d", .dims + 1, ncol(final_result)))
          }
          final_result
        })()
    )
    
    
  }





#' Get principal component information to a tibble using PCA
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#' @importFrom magrittr divide_by
#' @importFrom Matrix t
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean
#' @param ... Further parameters passed to the function prcomp
#'
#' @return A tibble with additional columns
#'
#'
get_reduced_dimensions_PCA_bulk_SE <-
  function(.data,
           .dims = 2,
           top = 500,
           of_samples = TRUE,
           transform = log1p,
           scale = FALSE,
           ...) {
    # Comply with CRAN NOTES
    . = NULL
    sdev = NULL
    x = NULL
    
    # Get components from dims
    components = 1:.dims
    
    
    
    # check that there are non-NA genes for enough samples
    if (nrow(.data) == 0) 
      stop(
        "tidybulk says: In calculating PCA there is no gene that have non NA values is all samples"
      )
    else if (nrow(.data) < 100) 
      warning(
        "
					tidybulk says: In PCA correlation there is < 100 genes that have non NA values is all samples.
The correlation calculation would not be reliable,
we suggest to partition the dataset for sample clusters.
					"
      )
    
    
    prcomp_obj =
      .data |>
      t() |>
      # Calculate principal components
      prcomp(scale = scale, ...)
    
    # Return
    list(
      raw_result = prcomp_obj,
      result =
        prcomp_obj |>
        # Anonymous function - Prints fraction of variance
        # input: PCA object
        # output: PCA object
        (\(.) {
          message("Fraction of variance explained by the selected principal components")
          (.) %$% sdev |> pow(2) |> # Eigen value
            unlist() |> divide_by(sum(unlist(.))) |>
            _[components] |>
            enframe() |>
            select(-name) |>
            rename(`Fraction of variance` = value) |>
            mutate(PC = components) |>
            capture.output() |> paste0(collapse = "\n") |> message()
          (.)
        })() %$%
        # Parse the PCA results to a tibble
        x |>
        as_tibble(rownames = "sample") |>
        select(sprintf("PC%s", components))
    )
    
    
  }

#' Get principal component information to a tibble using tSNE
#'
#' @keywords internal
#' @noRd
#'
#'
#'
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom Matrix t
#'
#' @param .data A tibble
#' @param .abundance A column symbol with the value the clustering is based on (e.g., `count`)
#' @param .dims A integer vector corresponding to principal components of interest (e.g., 1:6)
#' @param .feature A column symbol. The column that is represents entities to cluster (i.e., normally genes)
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param top An integer. How many top genes to select
#' @param of_samples A boolean
#' @param transform A function that will tranform the counts, by default it is log1p for RNA sequencing data, but for avoinding tranformation you can use identity
#' @param scale A boolean
#' @param ... Further parameters passed to the function Rtsne
#'
#' @return A tibble with additional columns
#'
get_reduced_dimensions_TSNE_bulk_SE <-
  function(.data,
           .dims = 2,
           top = 500,
           of_samples = TRUE,
           transform = log1p,
           scale = NULL, # This is only a dummy argument for making it compatibble with PCA
           ...) {
    # Comply with CRAN NOTES
    . = NULL
    Y = NULL
    
    # To avoid dplyr complications
    
    
    # Evaluate ...
    arguments <- list(...)
    if (!"check_duplicates" %in% names(arguments))
      arguments = arguments |> c(check_duplicates = FALSE)
    if (!"verbose" %in% names(arguments))
      arguments = arguments |> c(verbose = TRUE)
    if (!"dims" %in% names(arguments))
      arguments = arguments |> c(dims = .dims)
    
    
    # Check if package is installed, otherwise install
    check_and_install_packages("Rtsne")
    
    
    # Set perprexity to not be too high
    if (!"perplexity" %in% names(arguments))
      arguments = arguments |> c(perplexity = ((
        .data |> ncol() |> sum(-1)
      ) / 3 / 2) |> floor() |> min(30))
    
    # If not enough samples stop
    if (arguments$perplexity <= 2)
      stop("tidybulk says: You don't have enough samples to run tSNE")
    
    tsne_obj = .data |> t() |> as.matrix() |>  list() |> c(arguments) |> do.call(Rtsne::Rtsne, args = _)
    
    
    
    list(
      raw_result = tsne_obj,
      result = tsne_obj %$%
        Y |>
        as_tibble(.name_repair = "minimal") |>
        setNames(paste0("tSNE", seq_len(ncol(tsne_obj$Y)))) |>
        # add element name
        dplyr::mutate(sample = !!.data |> colnames()) |>
        select(-sample)
    )
    
  }

#' Get UMAP reduced dimensions
#'
#' @keywords internal
#' @noRd
#'
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#'
#' @param .data A SummarizedExperiment object containing the data to be reduced
#' @param .dims An integer specifying the number of UMAP dimensions to return (default = 2)
#' @param top An integer specifying how many top variable features to select
#' @param of_samples A boolean indicating whether to calculate for samples (TRUE) or features (FALSE)
#' @param transform A function to transform the abundance values (default = log1p)
#' @param scale A boolean indicating whether to scale the data before UMAP (default = NULL, only used for PCA compatibility)
#' @param calculate_for_pca_dimensions An integer specifying the number of PCA dimensions to use as UMAP input. 
#'                                    If NULL, variable features are used directly without PCA reduction.
#' @param ... Additional arguments passed to uwot::tumap()
#'
#' @return A list containing:
#'   \item{raw_result}{The raw UMAP object from uwot}
#'   \item{result}{A tibble with UMAP coordinates (UMAP1, UMAP2, etc.)}
#'
get_reduced_dimensions_UMAP_bulk_SE <-
  function(.data,
           .dims = 2,
           top = 500,
           of_samples = TRUE,
           transform = log1p,
           scale = NULL, # This is only a dummy argument for making it compatibble with PCA
           calculate_for_pca_dimensions = min(20, top),
           ...) {
    # Comply with CRAN NOTES
    . = NULL
    x = NULL
    
    # To avoid dplyr complications
    
    # Evaluate ...
    arguments <- list(...)
    # if (!"check_duplicates" %in% names(arguments))
    #   arguments = arguments |> c(check_duplicates = FALSE)
    if (!"dims" %in% names(arguments))
      arguments = arguments |> c(n_components = .dims)
    if (!"init" %in% names(arguments))
      arguments = arguments |> c(init = "spca")
    
    
    # Check if package is installed, otherwise install
    check_and_install_packages("uwot")
    
    
    # Calculate based on PCA
    if(!is.null(calculate_for_pca_dimensions))
      df_UMAP =
      .data |>
      t() |>
      # Calculate principal components
      prcomp(scale = scale) %$%
      # Parse the PCA results to a tibble
      x |>
      (\(.) .[,1:calculate_for_pca_dimensions])()
    
    # Calculate based on all features
    else
      df_UMAP = .data |> t() |>  as.matrix()
    
    umap_obj = do.call(uwot::tumap, c(list(df_UMAP), arguments))
    
    list(
      raw_result = umap_obj,
      result = umap_obj  |>
        as_tibble(.name_repair = "minimal") |>
        setNames(paste0("UMAP", seq_len(ncol(umap_obj)))) |>
        # add element name
        dplyr::mutate(sample = !!.data |> colnames()) |>
        select(-sample)
    )
    
  }
