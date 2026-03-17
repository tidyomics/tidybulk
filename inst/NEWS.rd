\name{NEWS}
\title{News for Package \pkg{tidybulk}}

\section{Changes in version 1.2.0, Bioconductor 3.12 Release}{
\itemize{
    \item Make gene filtering functionality `identify_abundance` explicit, a warning will be given if this has not been performed before the majority of workflow steps (e.g. `test_differential_abundance`).
    \item Add Automatic bibliography `get_bibliography`.
    \item Add DESeq2 and limma-voom to the methods for `test_differential_abundance` (method="DESeq2").
    \item Add prefix to test_differential_abundance for multi-methods analyses.
    \item Add other cell-type signature to `deconvolve_cellularity`.
    \item Add differential cellularity analyses `test_differential_cellularity`.
    \item Add gene descrption annotation utility `describe_transcript`.
    \item Add `nest` functionality for functional-coding transcriptomic analyses.
    \item Add gene overrepresentation functionality `test_gene_overrepresentation`.
    \item Add github website.
    \item Seep up data frame vadidation.
    \item Several bug fixes.
}}

\section{Changes in version 1.3.2, Bioconductor 3.13 Release}{
\itemize{
    \item Tidybulk now operates natively with SummarizedExperment data container, in a seamless way thanks to tidySummarisedExperiment 10.18129/B9.bioc.tidySummarizedExperiment
    \item Added robust edgeR as it outperforms many other methods as shown here doi.org/10.1093/nargab/lqab005
    \item Added test stratifiction cellularity, to easily calculate Kaplan-Meier curves
    \item Production of SummarizedExperiment from BAM or SAM files
    \item Added treat method to edgeR and voom differential transcription tests doi.org/10.1093/bioinformatics/btp053
    \item Added the method as_SummarizedExperiment
    \item Vastly improved test_gene_enrichment
    \item Added test_gene_rank, based on GSEA
    \item Several bug fixes.
}}

\section{Changes in version 1.5.5, Bioconductor 3.14 Release}{
\itemize{
    \item Added user-defined gene set for gene rank test
    \item Sped up aggregate_transcripts for large scale tibbles or SummarizedExperiment objects
    \item Allow passing additional arguments to DESeq2 method in test_differential_abundance
    \item Allow scale_abundance to run with a user-defined subset of genes (e.g. housekeeping genes)
    \item Add UMAP to reduce_dimensions()
    \item Several minor fixes, optimisations and documentation improvements
}}

\section{Changes in version 1.7.3, Bioconductor 3.15 Release}{
\itemize{
    \item Improve imputation and other features for sparse counts
    \item Cibersort deconvolution, check 0 counts
    \item Improve missing abundance with force scaling
    \item Other small fixes and messaging
}}

\section{Changes in version 1.7.4, Bioconductor 3.16 Dev}{
\itemize{
    \item Improved deconvolution robustness for SummarizedExperiment, edge cases
    \item Allow mapping of tidybulk_SAM_BAM to non-human genomes
    \item Adopt the vocabulary .feature, .sample, for conversion between SummarizedExperiment and tibble, similarly to tidySummarizedExperiment
    \item Deprecate .contrasts argument if favour of contrasts (with no dot)
    \item Make aggregate_duplicates more robust for tibble and SummarizedExperiment inputs
    \item Deprecate log_tranform argument for all methods for a more generic tranform argument that accepts arbitrary functions
}}

\section{Changes in version 1.9.2, Bioconductor 3.16 Dev}{
\itemize{
    \item Improve aggregate_duplicates for tibble and SummarizedExperiment
    \item Fix epic deconvolution when using DelayedMatrix
    \item Allow as_SummarizedExperiment with multiple columns identifiers for .sample and .feature
}}

\section{Changes in version 2.0.0, Bioconductor 3.22 Release}{
\itemize{
    \item Major refactoring to improve code maintainability and performance. This included the removal of all tbl methods in favor of SummarizedExperiment-based approaches.
    \item Replace deprecated pipe operator \%>\% with native |> operator for improved readability
    \item Refactor functions to use explicit function definitions instead of dplyr pipes
    \item Remove deprecated \code{action} parameter from multiple functions
    \item Remove deprecated \code{as_data_frame} function and clean up related code
    \item Remove all tbl methods in favor of SummarizedExperiment-based approaches
    \item Remove deprecated \code{.sample} and \code{.transcript} parameters from various methods
    \item Deprecate \code{.transcript} parameter in \code{aggregate_duplicates} and introduce \code{feature} parameter
    \item Add \code{feature_column} parameter to \code{deconvolve_cellularity} function
    \item Add \code{keep_integer} functionality to \code{aggregate_duplicates}
    \item Add \code{get_X_cibersort} function and update \code{deconvolve_cellularity} for improved reference handling
    \item Add \code{scale_x_log10_reverse} and \code{scale_y_log10_reverse} functions for enhanced axis transformations
    \item Add \code{formula_design} parameter to \code{identify_abundant} and \code{keep_abundant} methods
    \item Add \code{minimum_count_per_million} parameter to abundance-related methods
    \item Add \code{cores} parameter for dynamic specification of CPU cores during execution
    \item Add comprehensive tests for various functions in the analysis pipeline
    \item Add \code{filterByExpr} function and update methods to support \code{minimum_count_per_million}
    \item Enhance \code{as_SummarizedExperiment} function and deprecate \code{test_stratification_cellularity}
    \item Enhance parameter handling and deprecation warnings in abundance-related methods
    \item Enhance gene enrichment and testing functions for improved robustness and clarity
    \item Enhance \code{aggregate_duplicates} function to robustly handle numeric columns
    \item Update core functions to utilize dynamic core detection
    \item Update functions to use metadata instead of internals for better data handling
    \item Update differential analysis functions and documentation
    \item Update package dependencies and improve documentation
    \item Update examples and documentation to utilize the airway dataset
    \item Update vignettes to replace \code{.contrasts} with \code{contrasts} parameter
    \item Update vignettes to replace tidyverse with specific dplyr, tidyr, tibble, and purrr libraries
    \item Update README and documentation for clarity and completeness
    \item Fix CHECK issues and improve package validation
    \item Fix image paths in vignettes and README for better compatibility
    \item Remove logo image inclusion from vignettes for cleaner presentation
    \item Remove obsolete figures and streamline content
    \item Remove deprecated warnings and redundant messages
    \item Several bug fixes and optimizations
}}

\section{Changes in version 2.1.1, Bioconductor 3.23 Dev}{
\itemize{
    \item Maintenance release with minor bug fixes and performance improvements
    \item Bump minor version number following creation of RELEASE_3_22 branch
}}
