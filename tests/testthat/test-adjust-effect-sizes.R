context('Effect Size Adjustment Functions')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]

library(dplyr)
library(SummarizedExperiment)
library(tidySummarizedExperiment)
library(tidybulk)

# DESeq2 differential abundance on airway_mini, stored under the prefix "de_dex"
de_mini <- function() {
  airway_mini |>
    identify_abundant() |>
    test_differential_abundance(~ dex, method = "DESeq2", prefix = "de_dex")
}

test_adjustment <- function(type) {
  skip_if_not_installed("DESeq2")
  # "normal" is built into DESeq2, the others need their own package
  if (type != "normal") skip_if_not_installed(type)

  se_de <- suppressWarnings(suppressMessages(de_mini()))
  original <- metadata(se_de)$tidybulk$DESeq2_Resultsobject$de_dex

  res <- suppressMessages(
    adjust_effect_sizes(
      se_de,
      DE_name = "de_dex",
      COEF = "dex_untrt_vs_trt",
      type = type
    )
  )

  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(dim(res), dim(se_de))

  # Adjusted results are stored in the metadata under <DE_name>_adjusted
  adjusted <- metadata(res)$tidybulk$de_dex_adjusted
  expect_s4_class(adjusted, "DESeqResults")
  expect_equal(nrow(adjusted), nrow(original))
  expect_true("log2FoldChange" %in% colnames(adjusted))

  # The DESeq2 object is kept in the metadata
  expect_s4_class(metadata(res)$tidybulk$DESeq2_object, "DESeqDataSet")

  # Adjusted statistics are appended to rowData with their own names
  expect_true("de_dex_adjusted_log2FoldChange" %in% colnames(rowData(res)))
  expect_true(all(colnames(rowData(se_de)) %in% colnames(rowData(res))))
  expect_false(any(duplicated(colnames(rowData(res)))))
  expect_equal(
    unname(rowData(res)$de_dex_adjusted_log2FoldChange),
    unname(adjusted$log2FoldChange)
  )

  # Shrinkage should not inflate the overall effect sizes
  expect_lte(
    sum(abs(adjusted$log2FoldChange), na.rm = TRUE),
    sum(abs(original$log2FoldChange), na.rm = TRUE)
  )

  invisible(res)
}

test_that("adjust_effect_sizes works with apeglm for DESeq2 results", {
  test_adjustment("apeglm")
})

test_that("adjust_effect_sizes works with ashr for DESeq2 results", {
  test_adjustment("ashr")
})

test_that("adjust_effect_sizes works with normal for DESeq2 results", {
  test_adjustment("normal")
})

test_that("adjust_effect_sizes reports the DE name and coefficient", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("apeglm")

  se_de <- suppressWarnings(suppressMessages(de_mini()))

  expect_message(
    adjust_effect_sizes(
      se_de,
      DE_name = "de_dex",
      COEF = "dex_untrt_vs_trt",
      type = "apeglm"
    ),
    "de_dex.*dex_untrt_vs_trt"
  )
})

test_that("adjust_effect_sizes reports a newly specified coefficient", {
  skip_if_not_installed("DESeq2")

  se_de <- suppressWarnings(suppressMessages(de_mini()))

  expect_message(
    adjust_effect_sizes(se_de, COEF = "dex_untrt_vs_trt", type = "normal"),
    "The coefficient was specified as dex_untrt_vs_trt"
  )
})

test_that("adjust_effect_sizes works when DESeq2 was called in lowercase", {
  skip_if_not_installed("DESeq2")

  se_de <- suppressWarnings(suppressMessages(
    airway_mini |>
      identify_abundant() |>
      test_differential_abundance(~ dex, method = "deseq2", prefix = "de_dex")
  ))

  res <- suppressMessages(
    adjust_effect_sizes(
      se_de,
      DE_name = "de_dex",
      COEF = "dex_untrt_vs_trt",
      type = "normal"
    )
  )

  expect_s4_class(metadata(res)$tidybulk$de_dex_adjusted, "DESeqResults")
  expect_s4_class(metadata(res)$tidybulk$deseq2_object, "DESeqDataSet")
  expect_null(metadata(res)$tidybulk$DESeq2_object)
})

test_that("adjusting twice replaces the adjusted columns instead of duplicating them", {
  skip_if_not_installed("DESeq2")

  se_de <- suppressWarnings(suppressMessages(de_mini()))

  once <- suppressMessages(
    adjust_effect_sizes(se_de, DE_name = "de_dex", type = "normal")
  )
  twice <- suppressMessages(
    adjust_effect_sizes(once, DE_name = "de_dex", type = "normal")
  )

  expect_equal(colnames(rowData(twice)), colnames(rowData(once)))
})

test_that("adjust_effect_sizes errors on missing or invalid inputs", {
  skip_if_not_installed("DESeq2")

  se_de <- suppressWarnings(suppressMessages(de_mini()))

  expect_error(
    adjust_effect_sizes(se_de),
    "Both DE_name and COEF are missing"
  )
  expect_error(
    adjust_effect_sizes(se_de, DE_name = "not_a_result"),
    "not a valid DESeqResults object"
  )
  expect_error(
    adjust_effect_sizes(se_de, DE_name = "de_dex", type = "not_a_type")
  )
  expect_error(
    adjust_effect_sizes(airway_mini, DE_name = "de_dex"),
    "No DESeq2 object found"
  )
})
