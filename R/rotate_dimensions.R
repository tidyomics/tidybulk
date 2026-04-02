#' Rotate two dimensions (e.g., principal components) of an arbitrary angle
#'
#' `r lifecycle::badge("maturing")`
#'
#' @description rotate_dimensions() takes as input a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo quo_name
#' @importFrom magrittr not
#' @importFrom dplyr between
#' @importFrom SummarizedExperiment colData rowData
#'
#'
#' @name rotate_dimensions
#'
#' @param .data A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))
#' @param .element The name of the element column (normally samples).
#'
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param of_samples A boolean. In case the input is a tidybulk object, it indicates Whether the element column will be sample or transcript column
#' @param dimension_1_column_rotated A character string. The column of the rotated dimension 1 (optional)
#' @param dimension_2_column_rotated A character string. The column of the rotated dimension 2 (optional)
#'
#' @details This function to rotate two dimensions such as the reduced dimensions.
#'
#' Underlying custom method:
#' rotation = function(m, d) {
#'   r = d * pi / 180
#'   ((bind_rows(
#'     c(`1` = cos(r), `2` = -sin(r)),
#'     c(`1` = sin(r), `2` = cos(r))
#'   ) |> as_matrix()) %*% m)
#' }
#'
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
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
#' counts.MDS =
#'  airway |>
#'  identify_abundant() |>
#'  reduce_dimensions(.abundance = counts, method="MDS", .dims = 3)
#'
#' counts.MDS.rotated =  rotate_dimensions(counts.MDS, `Dim1`, `Dim2`, rotation_degrees = 45, .element = sample)
#'
#' @references
#' Mangiola, S., Molania, R., Dong, R., Doyle, M. A., & Papenfuss, A. T. (2021). tidybulk: an R tidy framework for modular transcriptomic data analysis. Genome Biology, 22(1), 42. doi:10.1186/s13059-020-02233-7
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#' @export
#'
setGeneric("rotate_dimensions", function(.data,
                                         dimension_1_column,
                                         dimension_2_column,
                                         rotation_degrees,
                                         .element = NULL,
                                         of_samples = TRUE,
                                         dimension_1_column_rotated = NULL,
                                         dimension_2_column_rotated = NULL)
  standardGeneric("rotate_dimensions"))




#' Rotate two coordinate columns and append the rotated axes
#'
#' This internal helper applies a planar rotation to two numeric columns
#' that represent a low-dimensional embedding (for example PCA or UMAP
#' coordinates) stored in either `colData()` or `rowData()` of a
#' `SummarizedExperiment`.  It returns the original object with two
#' additional columns containing the rotated values.  The user specifies
#' the rotation angle in degrees and may provide custom names for the new
#' columns; otherwise sensible defaults are generated.
#'
#' @param .data  A `SummarizedExperiment` (or derivative) holding the
#'               coordinates to be rotated.
#' @param dimension_1_column Symbol or bare column name for the first axis
#'               (e.g. `UMAP_1`).
#' @param dimension_2_column Symbol or bare column name for the second axis
#'               (e.g. `UMAP_2`).
#' @param rotation_degrees   Numeric scalar in the closed interval
#'               \([-360, 360]\) indicating the anti-clockwise rotation
#'               angle.
#' @param .element Optional quoted column holding sample or feature labels
#'               (unused, retained for compatibility).
#' @param of_samples Logical.  If `TRUE` (default) the function rotates
#'               columns in `colData()`.  If `FALSE` it operates on
#'               `rowData()`.
#' @param dimension_1_column_rotated Optional symbol to name the new first
#'               rotated coordinate column.
#' @param dimension_2_column_rotated Optional symbol to name the new second
#'               rotated coordinate column.
#'
#' @return The input `SummarizedExperiment` with two extra metadata columns
#'         containing the rotated axes.
#'
#' @keywords internal
#'
#' @importFrom rlang enquo quo_name quo_is_null
#' @importFrom purrr when
#' @importFrom magrittr not
#' @importFrom dplyr between
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom stats setNames
.rotate_dimensions_se = function(.data,
                                 dimension_1_column,
                                 dimension_2_column,
                                 rotation_degrees,
                                 .element = NULL,
                                 
                                 of_samples = TRUE,
                                 dimension_1_column_rotated = NULL,
                                 dimension_2_column_rotated = NULL) {
  
  # Fix NOTEs
  . = NULL
  
  # Parse other colnames
  dimension_1_column = enquo(dimension_1_column)
  dimension_2_column = enquo(dimension_2_column)
  dimension_1_column_rotated = enquo(dimension_1_column_rotated)
  dimension_2_column_rotated = enquo(dimension_2_column_rotated)
  
  # Set default col names for rotated dimensions if not set
  if (quo_is_null(dimension_1_column_rotated))
    dimension_1_column_rotated = as.symbol(sprintf(
      "%s_rotated_%s",
      quo_name(dimension_1_column),
      rotation_degrees
    ))
  if (quo_is_null(dimension_2_column_rotated))
    dimension_2_column_rotated = as.symbol(sprintf(
      "%s_rotated_%s",
      quo_name(dimension_2_column),
      rotation_degrees
    ))
  
  # Sanity check of the angle selected
  if (rotation_degrees |> between(-360, 360) |> not())
    stop("tidybulk says: rotation_degrees must be between -360 and 360")
  
  
  # Return
  my_rotated_dimensions =
    .data |>
    
    # Select correct annotation
    (function(data_obj) {
      if (of_samples) {
        colData(data_obj)
      } else {
        rowData(data_obj)
      }
    })() |>
    
    # Select dimensions
    (\(.) .[,c(quo_name(dimension_1_column), quo_name(dimension_2_column))])() |>
    as.matrix() |>
    t() |>
    rotation(rotation_degrees) |>
    t() |>
    as.data.frame() |>
    setNames(c(
      quo_name(dimension_1_column_rotated),
      quo_name(dimension_2_column_rotated)
    ))
  
  
  .data |>
    
    # Add dimensions to metadata
    (function(data_obj) {
      if (of_samples) {
        colData(data_obj) <- colData(data_obj) |> cbind(my_rotated_dimensions)
        data_obj
      } else {
        rowData(data_obj) <- rowData(data_obj) |> cbind(my_rotated_dimensions)
        data_obj
      }
    })()
  
}

#' rotate_dimensions
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("rotate_dimensions",
          "SummarizedExperiment",
          .rotate_dimensions_se)

#' rotate_dimensions
#'
#' @docType methods
#' @rdname rotate_dimensions-methods
#'
#' @return A `SummarizedExperiment` object
#'
setMethod("rotate_dimensions",
          "RangedSummarizedExperiment",
          .rotate_dimensions_se)

