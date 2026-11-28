# Internal. edgeR tagwise (or trended if n >= 1000) written to rowData.
# dispersion_trend = "common" writes one shared dispersion instead.
# formula_abundance must be a fixed-effects formula; the user chooses it.

estimate_dispersion <- function(.data, formula_abundance, abundance = "counts",
                                dispersion_trend = c("tagwise", "common")) {
  dispersion_trend <- match.arg(dispersion_trend)
  check_and_install_packages("edgeR")
  check_formula(formula_abundance)
  if (grepl("|", paste(deparse(formula_abundance), collapse = ""), fixed = TRUE)) {
    stop(
      "tidybulk says: estimate_dispersion() uses edgeR, which has no random effects. ",
      "Pass a fixed-effects formula (e.g. ~ treatment + donor), not ",
      paste(deparse(formula_abundance), collapse = ""),
      ".",
      call. = FALSE
    )
  }

  n_sample <- ncol(.data)
  design <- stats::model.matrix(
    formula_abundance,
    data = droplevels(as.data.frame(colData(.data)))
  )
  if (n_sample <= ncol(design)) {
    stop(
      "tidybulk says: The dispersion design has ", ncol(design),
      " coefficients for ", n_sample,
      " samples, leaving no residual degrees of freedom.",
      call. = FALSE
    )
  }

  formula_txt <- paste(deparse(formula_abundance), collapse = "")
  if (identical(dispersion_trend, "common")) {
    message(
      "tidybulk says: calculating common dispersion ",
      "with edgeR::estimateDisp(trend.method = \"none\", tagwise = FALSE) using ",
      formula_txt, "."
    )
    fit <- edgeR::estimateDisp(
      assay(.data, abundance),
      design = design,
      trend.method = "none",
      tagwise = FALSE
    )
    disp <- trended <- fit$common.dispersion
    d_eff <- n_sample - ncol(design)
  } else if (n_sample < 1000L) {
    message(
      "tidybulk says: calculating tagwise (shrinked) and trended dispersion ",
      "with edgeR::estimateDisp() using ", formula_txt, "."
    )
    counts <- assay(.data, abundance)
    fit <- edgeR::estimateDisp(counts, design = design)
    disp <- fit$tagwise.dispersion
    trended <- fit$trended.dispersion
    d_eff <- (n_sample - ncol(design)) + fit$prior.df
  } else {
    sampled <- sample(seq_len(n_sample), size = min(n_sample, 2000L))
    se_sub <- .data[, sampled, drop = FALSE]
    design <- stats::model.matrix(
      formula_abundance,
      data = droplevels(as.data.frame(colData(se_sub)))
    )
    message(
      "tidybulk says: n >= 1000; calculating trended dispersion ",
      "with edgeR::estimateTrendedDisp()."
    )
    fit <- edgeR::estimateTrendedDisp(
      assay(se_sub, abundance),
      design = design,
      subset = 1000,
      rowsum.filter = 10
    )
    disp <- trended <- fit
    d_eff <- ncol(se_sub) - ncol(design)
  }

  rowData(.data)[["dispersion_shrinked"]] <- rep_len(as.numeric(disp), nrow(.data))
  rowData(.data)[["dispersion_trended"]] <- rep_len(as.numeric(trended), nrow(.data))
  rowData(.data)[["dispersion_degrees_freedom"]] <-
    rep_len(as.numeric(d_eff), nrow(.data))
  .data <- attach_to_metadata(.data, fit, "estimateDisp")
  memorise_methods_used(.data, "edger")
}
