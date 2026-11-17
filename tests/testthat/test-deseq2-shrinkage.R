context('DESeq2 log2 fold change shrinkage')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]

library(dplyr)
library(SummarizedExperiment)
library(tidySummarizedExperiment)
library(tidybulk)

run_deseq2 <- function(method, ...) {
  suppressWarnings(suppressMessages(
    airway_mini |>
      identify_abundant() |>
      test_differential_abundance(~ dex, method = method, ...)
  ))
}

# One contrast with plain column names, as shrinkage requires explicit contrasts
run_one_contrast <- function(method, ...) {
  run_deseq2(method, contrasts = list(c("dex", "untrt", "trt")), omit_contrast_in_colnames = TRUE, ...)
}

test_shrinkage <- function(method, type) {
  skip_if_not_installed("DESeq2")
  # "normal" is built into DESeq2, the others need their own package
  if (type != "normal") skip_if_not_installed(type)

  unshrunk <- rowData(run_one_contrast("DESeq2"))
  shrunk <- rowData(run_one_contrast(method))

  expect_equal(colnames(shrunk), colnames(unshrunk))

  # p-values come from DESeq2::results(), only the effect sizes are shrunk
  expect_equal(shrunk$pvalue, unshrunk$pvalue)
  expect_equal(shrunk$padj, unshrunk$padj)
  expect_false(isTRUE(all.equal(shrunk$log2FoldChange, unshrunk$log2FoldChange)))
  expect_lte(
    sum(abs(shrunk$log2FoldChange), na.rm = TRUE),
    sum(abs(unshrunk$log2FoldChange), na.rm = TRUE)
  )
}

test_that("DESeq2_shrink_apeglm shrinks log2 fold changes", {
  test_shrinkage("DESeq2_shrink_apeglm", "apeglm")
})

test_that("DESeq2_shrink_ashr shrinks log2 fold changes", {
  test_shrinkage("DESeq2_shrink_ashr", "ashr")
})

test_that("DESeq2_shrink_normal shrinks log2 fold changes", {
  test_shrinkage("DESeq2_shrink_normal", "normal")
})

test_that("shrinkage works with explicit contrasts, including a non-reference denominator", {
  skip_if_not_installed("DESeq2")

  contrasts <- list(c("dex", "trt", "untrt"), c("dex", "untrt", "trt"))
  methods <- c(
    "DESeq2_shrink_normal",
    if (requireNamespace("apeglm", quietly = TRUE)) "DESeq2_shrink_apeglm",
    if (requireNamespace("ashr", quietly = TRUE)) "DESeq2_shrink_ashr"
  )
  for (method in methods) {
    res <- rowData(run_deseq2(method, contrasts = contrasts))
    expect_equal(
      res$`log2FoldChange___dex trt-untrt`,
      -res$`log2FoldChange___dex untrt-trt`,
      tolerance = 0.05
    )
  }
})

test_that("shrinkage handles multiple contrasts across different factors", {
  skip_if_not_installed("DESeq2")

  # The second contrast's denominator (N61311) is not the reference level of cell (N052611)
  contrasts <- list(c("dex", "trt", "untrt"), c("cell", "N061011", "N61311"))
  labels <- c("dex trt-untrt", "cell N061011-N61311")
  run_multi <- function(method) {
    suppressWarnings(suppressMessages(
      airway[1:100, ] |>
        identify_abundant() |>
        test_differential_abundance(~ dex + cell, method = method, contrasts = contrasts)
    ))
  }

  unshrunk <- rowData(run_multi("DESeq2"))
  methods <- c(
    "DESeq2_shrink_normal",
    if (requireNamespace("apeglm", quietly = TRUE)) "DESeq2_shrink_apeglm",
    if (requireNamespace("ashr", quietly = TRUE)) "DESeq2_shrink_ashr"
  )

  for (method in methods) {
    res <- run_multi(method)
    shrunk <- rowData(res)
    expect_equal(colnames(shrunk), colnames(unshrunk))

    stored <- metadata(res)$tidybulk[[paste0(method, "_shrunk")]]
    expect_equal(names(stored), labels)

    for (label in labels) {
      lfc <- sprintf("log2FoldChange___%s", label)
      pvalue <- sprintf("pvalue___%s", label)

      # Each contrast is shrunk on its own and keeps its own DESeq2 p-values
      expect_equal(shrunk[[pvalue]], unshrunk[[pvalue]])
      expect_false(isTRUE(all.equal(shrunk[[lfc]], unshrunk[[lfc]])))
      expect_lte(sum(abs(shrunk[[lfc]]), na.rm = TRUE), sum(abs(unshrunk[[lfc]]), na.rm = TRUE))

      expect_s4_class(stored[[label]], "DESeqResults")
      expect_equal(unname(stored[[label]][rownames(res), "log2FoldChange"]), unname(shrunk[[lfc]]))
    }

    expect_false(isTRUE(all.equal(
      shrunk[[sprintf("log2FoldChange___%s", labels[1])]],
      shrunk[[sprintf("log2FoldChange___%s", labels[2])]]
    )))
  }
})

test_that("shrunk DESeqResults are stored in the metadata", {
  skip_if_not_installed("DESeq2")

  expect_message(
    suppressWarnings(
      airway_mini |>
        identify_abundant() |>
        test_differential_abundance(~ dex, method = "DESeq2_shrink_normal", contrasts = list(c("dex", "untrt", "trt")))
    ),
    "metadata\\(\\.\\)\\$tidybulk\\$DESeq2_shrink_normal_shrunk"
  )

  res <- run_one_contrast("DESeq2_shrink_normal")
  shrunk <- metadata(res)$tidybulk$DESeq2_shrink_normal_shrunk
  expect_s4_class(shrunk, "DESeqResults")
  expect_equal(unname(shrunk[rownames(res), "log2FoldChange"]), unname(rowData(res)$log2FoldChange))

  contrasts <- list(c("dex", "trt", "untrt"), c("dex", "untrt", "trt"))
  res <- run_deseq2("DESeq2_shrink_normal", contrasts = contrasts)
  shrunk <- metadata(res)$tidybulk$DESeq2_shrink_normal_shrunk
  expect_equal(names(shrunk), c("dex trt-untrt", "dex untrt-trt"))
  expect_equal(
    unname(shrunk[["dex trt-untrt"]][rownames(res), "log2FoldChange"]),
    unname(rowData(res)$`log2FoldChange___dex trt-untrt`)
  )

  expect_null(metadata(run_deseq2("DESeq2"))$tidybulk$DESeq2_shrunk)
})

test_that("shrunk methods keep prefix and metadata conventions", {
  skip_if_not_installed("DESeq2")

  res <- run_one_contrast("DESeq2_shrink_normal", prefix = "de_")

  expect_true("de_log2FoldChange" %in% colnames(rowData(res)))
  expect_s4_class(metadata(res)$tidybulk$DESeq2_shrink_normal_object, "DESeqDataSet")
})

test_that("shrinkage refuses to run without explicit contrasts", {
  skip_if_not_installed("DESeq2")

  expect_error(run_deseq2("DESeq2_shrink_normal"), "requires `contrasts`")
  expect_error(run_deseq2("DESeq2_shrink_normal", contrasts = list("dex")), "requires `contrasts`")
})
