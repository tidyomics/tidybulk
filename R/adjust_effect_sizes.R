#' Adjust effect sizes of DESeq2 results
#'
#' @param .data A 'SummarizedExperiment'
#' @param DE_name the name of the DE results to shrink
#' @param type Character string, type of adjustment to perform
#' @param COEF a character string. Coefficient to adjust effect sizes for
#' @param DESIGN A formula or character string that can be converted to a formula for the DE calculations.
#' @param ... passed on to DESeq2::lfcShrink()
#'
#' @details
#' The function allows for effect size adjustments after DESeq2 DE calculations. By default, specified .data and DE_name is enough to perform the adjustment. In case the DE analysis does not exist, it will be generated using the specified coefficient.
#'
#' @return a SummarizedExperiment
#' @examples
#' T
#'
#' @export
#'
adjust_effect_sizes <- function(
  .data,
  DE_name,
  type = c("apeglm", "ashr", "normal"),
  COEF = NULL,
  DESIGN = NULL,
  ...
) {
  type <- match.arg(type, choices = c("apeglm", "ashr", "normal"))
  deseq_object <- metadata(.data)$tidybulk$DESeq2_object #legibility

  if (missing(DE_name)) {
    if (is.null(COEF)) {
      stop("Both DE_name and COEF are missing. Please specifiy at least one!")
    }
    # didn't know what the convention is here.
    # I just want the comparison to fail reliably below, forcing a generation of results
    incomingCoef <- "def_not_a_valid_coef"
    DE_name <- "comparison"
  } else {
    de_res <- metadata(.data)$tidybulk$DESeq2_Resultsobject[[DE_name]]
    if (!is(de_res, "DESeqResults")) {
      stop(
        "The specified DE result is not a valid DESeqResults object. Please respecify!"
      )
    }
    incomingCoef <- gsub(
      " ",
      "_",
      sub("log2 fold change \\(MLE\\): ", "", mcols(de_res)$description[2])
    )
  }

  incomingDesign <- BiocGenerics::design(deseq_object)

  if (!is.null(DESIGN) && as.formula(DESIGN) != incomingDesign) {
    design(deseq_object) <- as.formula(DESIGN)
    deseq_object <- DESeq2::DESeq(deseq_object)
    incomingDesign <- BiocGenerics::design(deseq_object)
  }

  coef_message <- paste0("The coefficient is assumed as ", incomingCoef, ".")
  if (!is.null(COEF) && COEF != incomingCoef) {
    deseq_object <- DESeq2::DESeq(deseq_object)
    de_res <- DESeq2::results(deseq_object, name = COEF)

    incomingCoef <- COEF
    coef_message <- paste0(
      "The coefficient was specified as",
      as.character(COEF)
    )
  }

  if (type == "ashr") {
    adjusted <- DESeq2::lfcShrink(
      deseq_object,
      res = de_res,
      type = "ashr",
      ...
    )
  } else if (type == "apeglm") {
    adjusted <- DESeq2::lfcShrink(
      deseq_object,
      coef = incomingCoef,
      type = "apeglm",
      ...
    )
  } else if (type == "normal") {
    adjusted <- DESeq2::lfcShrink(
      deseq_object,
      coef = incomingCoef,
      type = "normal",
      ...
    )
  } else {
    stop(
      "please specify a valid type. Valid options include apeglm, ahr and normal"
    )
  }

  metadata(.data)$tidybulk$DESeq2_object <- deseq_object

  message_text <- paste0(
    "Adjusting effects for the de results stored as ",
    DE_name,
    ". The design is ",
    stringr::str_flatten(incomingDesign),
    "."
  )

  res_df <- adjusted |>
    as_tibble(rownames = "transcript")
  cn <- colnames(res_df)
  cn_new <- c(
    cn[1],
    sprintf("%s%s", paste(DE_name, sep = "_"), cn[2:length(cn)])
  )
  colnames(res_df) <- cn_new

  stats_matrix <- res_df |>
    as_matrix(rownames = "transcript")
  stats_matrix <- stats_matrix[
    match(rownames(rowData(.data)), rownames(stats_matrix)),
    ,
    drop = FALSE
  ]
  rowData(.data) <- cbind(rowData(.data), stats_matrix)

  message(paste(message_text, coef_message))

  data_obj_intermediate <- attach_to_metadata(
    .data,
    adjusted,
    paste0(DE_name, "_adjusted")
  )

  rlang::inform(
    sprintf(
      "tidybulk says: to access the adjusted effects do `metadata(.)$tidybulk$%s_adjusted`",
      DE_name
    ),
    .frequency_id = sprintf("Access DE results %s", DE_name),
    .frequency = "always"
  )

  return(data_obj_intermediate)
}
