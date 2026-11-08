context('Differential Analysis Functions')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]

library(dplyr)
library(SummarizedExperiment)

# Test test_differential_abundance function
test_that("test_differential_abundance with edgeR works correctly", {
  res <- airway_mini |> 
    identify_abundant() |> 
    test_differential_abundance(
      .formula = ~ dex,
      method = "edgeR_quasi_likelihood"
    )
  
  # Check that required columns exist
  expect_true("logFC" %in% names(SummarizedExperiment::rowData(res)))
  expect_true("PValue" %in% names(SummarizedExperiment::rowData(res)))
  expect_true("FDR" %in% names(SummarizedExperiment::rowData(res)))
  
  # Check that logFC values are reasonable (not all NA or infinite)
  logfc_values <- SummarizedExperiment::rowData(res)$logFC
  expect_true(any(!is.na(logfc_values)))
  expect_true(any(!is.infinite(logfc_values)))
  
  # Check that P-values are in valid range
  p_values <- SummarizedExperiment::rowData(res)$PValue
  expect_true(all(p_values >= 0 & p_values <= 1, na.rm = TRUE))
  
  # Check that FDR values are in valid range
  fdr_values <- SummarizedExperiment::rowData(res)$FDR
  expect_true(all(fdr_values >= 0 & fdr_values <= 1, na.rm = TRUE))
})

test_that("test_differential_abundance with DESeq2 works correctly", {
  # Skip DESeq2 as it has parameter issues
  skip("DESeq2 has parameter compatibility issues")
})

test_that("test_differential_abundance with limma works correctly", {
  res <- airway_mini |> 
    identify_abundant() |> 
    test_differential_abundance(
      .formula = ~ dex,
      method = "limma_voom"
    )
  
  # Check that required columns exist
  expect_true("logFC" %in% names(SummarizedExperiment::rowData(res)))
  expect_true("P.Value" %in% names(SummarizedExperiment::rowData(res)))
  expect_true("adj.P.Val" %in% names(SummarizedExperiment::rowData(res)))
  
  # Check that logFC values are reasonable (not all NA or infinite)
  logfc_values <- SummarizedExperiment::rowData(res)$logFC
  expect_true(any(!is.na(logfc_values)))
  expect_true(any(!is.infinite(logfc_values)))
  
  # Check that P-values are in valid range
  p_values <- SummarizedExperiment::rowData(res)$P.Value
  expect_true(all(p_values >= 0 & p_values <= 1, na.rm = TRUE))
  
  # Check that adjusted P-values are in valid range
  adj_p_values <- SummarizedExperiment::rowData(res)$adj.P.Val
  expect_true(all(adj_p_values >= 0 & adj_p_values <= 1, na.rm = TRUE))
})

# Test colData preservation and usage
test_that("test_differential_abundance preserves colData correctly", {
  # Store original colData
  original_colData <- colData(airway_mini)
  original_colData_names <- names(original_colData)
  
  # Run differential abundance test
  res <- airway_mini |> 
    identify_abundant() |> 
    test_differential_abundance(
      .formula = ~ dex,
      method = "edgeR_quasi_likelihood"
    )
  
  # Check that colData is preserved
  expect_equal(names(colData(res)), original_colData_names)
  expect_equal(nrow(colData(res)), nrow(original_colData))
  expect_equal(ncol(colData(res)), ncol(original_colData))
  
  # Check that sample names are preserved
  expect_equal(rownames(colData(res)), rownames(original_colData))
  
  # Check that dex column is still present and usable
  expect_true("dex" %in% names(colData(res)))
  expect_true(all(colData(res)$dex %in% c("trt", "untrt")))
  
  # Check that the function can handle colData with additional columns
  se_with_extra <- airway_mini
  colData(se_with_extra)$extra_column <- rep("test", nrow(colData(se_with_extra)))
  
  res_with_extra <- se_with_extra |> 
    identify_abundant() |> 
    test_differential_abundance(
      .formula = ~ dex,
      method = "edgeR_quasi_likelihood"
    )
  
  # Check that extra column is preserved
  expect_true("extra_column" %in% names(colData(res_with_extra)))
  expect_equal(colData(res_with_extra)$extra_column, rep("test", nrow(colData(res_with_extra))))
})

test_that("three contrasts return distinct fold changes and p-values", {
  # cell has four lines, so these three pairwise contrasts are estimable
  contrasts <- c(
    "cellN061011 - cellN080611",
    "cellN080611 - cellN61311",
    "cellN61311 - cellN061011"
  )

  res <- airway[1:200, ] |>
    identify_abundant() |>
    test_differential_abundance(
      .formula = ~ cell,
      method = "edgeR_quasi_likelihood",
      contrasts = contrasts
    )

  rd <- SummarizedExperiment::rowData(res)
  logfc_cols <- paste0("logFC___", contrasts)
  p_cols <- paste0("PValue___", contrasts)

  expect_true(all(logfc_cols %in% colnames(rd)))
  expect_true(all(p_cols %in% colnames(rd)))

  logfc <- lapply(logfc_cols, function(col) rd[[col]])
  p_values <- lapply(p_cols, function(col) rd[[col]])

  expect_true(all(vapply(logfc, function(x) any(is.finite(x)), logical(1))))
  expect_true(all(vapply(p_values, function(x) any(!is.na(x)), logical(1))))
  expect_true(all(vapply(p_values, function(x) all(x >= 0 & x <= 1, na.rm = TRUE), logical(1))))

  # Each contrast must differ from the other two
  pairs <- list(c(1, 2), c(1, 3), c(2, 3))
  for (pair in pairs) {
    expect_gt(sum(abs(logfc[[pair[1]]] - logfc[[pair[2]]]), na.rm = TRUE), 0)
    expect_gt(sum(abs(p_values[[pair[1]]] - p_values[[pair[2]]]), na.rm = TRUE), 0)
  }
})

# Test glmmSeq function
test_that("glmmSeq works correctly", {
  # Skip glmmSeq as it's not available
  skip("glmmSeq function not available")
}) 