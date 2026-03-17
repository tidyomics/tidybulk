
library(SummarizedExperiment)
library(airway)

# Use testthat 3e style
test_that("filterByExpr filters by min.count and CPM.Cutoff (minimum_count_per_million)", {
  data(airway)
  se = airway
  mat = assay(se)
  filterByExpr_SE <- get("filterByExpr_SE", asNamespace("tidybulk"))
  filterByExpr <- get("filterByExpr", asNamespace("tidybulk"))

  # Test default min.count filtering on SummarizedExperiment
  keep_default_se <- filterByExpr_SE(se, min.count = 10)
  expect_type(keep_default_se, "logical")
  expect_length(keep_default_se, nrow(se))
  expect_true(any(!keep_default_se))

  # Test default min.count filtering on matrix
  keep_default <- filterByExpr(mat, min.count = 10)
  expect_type(keep_default, "logical")
  expect_length(keep_default, nrow(mat))
  expect_true(any(!keep_default))

  # Test CPM.Cutoff filtering on matrix: compare two CPM cutoffs
  keep_cpm10 <- filterByExpr(mat, CPM.Cutoff = 10)
  keep_cpm100 <- filterByExpr(mat, CPM.Cutoff = 100)
  expect_true(sum(keep_cpm100) <= sum(keep_cpm10))

  # If all counts are high, all should be kept
  mat_high <- matrix(1000, nrow = 3, ncol = 4)
  keep_all <- filterByExpr(mat_high, min.count = 10)
  expect_true(all(keep_all))
})

test_that("keep_abundant works with minimum_counts and minimum_count_per_million", {
  data(airway)
  se = airway
  colnames_cd <- colnames(colData(se))

  # Using minimum_counts
  se_abundant_10 <- tidybulk::keep_abundant(se, minimum_counts = 10)
  se_abundant_100 <- tidybulk::keep_abundant(se, minimum_counts = 100)
  expect_true(nrow(se_abundant_100) <= nrow(se_abundant_10))

  # Using minimum_count_per_million
  se_abundant_cpm10 <- tidybulk::keep_abundant(se, minimum_count_per_million = 10)
  se_abundant_cpm100 <- tidybulk::keep_abundant(se, minimum_count_per_million = 100)
  expect_true(nrow(se_abundant_cpm100) <= nrow(se_abundant_cpm10))

  # Using both together (minimum_count_per_million should override minimum_counts)
  se_abundant_both <- tidybulk::keep_abundant(se, minimum_counts = 100, minimum_count_per_million = 100)
  # When both are provided, minimum_count_per_million takes precedence
  # So this should behave the same as using only minimum_count_per_million = 100
  expect_equal(nrow(se_abundant_both), nrow(se_abundant_cpm100))
})

test_that("keep_abundant works with factor_of_interest", {
  data(airway)
  se = airway
  colnames_cd <- colnames(colData(se))
  
  # Use the first available column as factor_of_interest
  factor_col <- colnames_cd[1]
  
  # Test with factor_of_interest
  se_abundant_with_factor <- tidybulk::keep_abundant(se, minimum_counts = 10, factor_of_interest = !!rlang::sym(factor_col))
  se_abundant_without_factor <- tidybulk::keep_abundant(se, minimum_counts = 10)
  
  # Both should work without errors
  expect_true(nrow(se_abundant_with_factor) > 0)
  expect_true(nrow(se_abundant_without_factor) > 0)
})

test_that("identify_abundant works with minimum_counts and minimum_count_per_million", {
  data(airway)
  se = airway
  
  # Test with minimum_counts only
  se_abundant_10 <- tidybulk::identify_abundant(se, minimum_counts = 10)
  se_abundant_100 <- tidybulk::identify_abundant(se, minimum_counts = 100)
  expect_true(sum(rowData(se_abundant_100)$.abundant) <= sum(rowData(se_abundant_10)$.abundant))
  
  # Test with minimum_count_per_million only
  se_abundant_cpm10 <- tidybulk::identify_abundant(se, minimum_count_per_million = 10)
  se_abundant_cpm100 <- tidybulk::identify_abundant(se, minimum_count_per_million = 100)
  expect_true(sum(rowData(se_abundant_cpm100)$.abundant) <= sum(rowData(se_abundant_cpm10)$.abundant))
  
  # Test that minimum_count_per_million overrides minimum_counts
  se_abundant_both <- tidybulk::identify_abundant(se, minimum_counts = 100, minimum_count_per_million = 100)
  # When both are provided, minimum_count_per_million takes precedence
  # So this should behave the same as using only minimum_count_per_million = 100
  expect_equal(sum(rowData(se_abundant_both)$.abundant), sum(rowData(se_abundant_cpm100)$.abundant))
}) 

test_that("identify_abundant and keep_abundant work with design argument", {
  data(airway)
  se = airway
  # Use the first available column as a factor for the design
  design <- model.matrix(~ dex + cell, data = colData(se))

  # identify_abundant with design
  se_abundant_10 <- tidybulk::identify_abundant(se, minimum_counts = 10, design = design)
  se_abundant_100 <- tidybulk::identify_abundant(se, minimum_counts = 100, design = design)
  expect_true(sum(rowData(se_abundant_100)$.abundant) <= sum(rowData(se_abundant_10)$.abundant))

  se_abundant_cpm10 <- tidybulk::identify_abundant(se, minimum_count_per_million = 10, design = design)
  se_abundant_cpm100 <- tidybulk::identify_abundant(se, minimum_count_per_million = 100, design = design)
  expect_true(sum(rowData(se_abundant_cpm100)$.abundant) <= sum(rowData(se_abundant_cpm10)$.abundant))

  se_abundant_both <- tidybulk::identify_abundant(se, minimum_counts = 100, minimum_count_per_million = 100, design = design)
  expect_equal(sum(rowData(se_abundant_both)$.abundant), sum(rowData(se_abundant_cpm100)$.abundant))

  # keep_abundant with design
  se_keep_10 <- tidybulk::keep_abundant(se, minimum_counts = 10, design = design)
  se_keep_100 <- tidybulk::keep_abundant(se, minimum_counts = 100, design = design)
  expect_true(nrow(se_keep_100) <= nrow(se_keep_10))

  se_keep_cpm10 <- tidybulk::keep_abundant(se, minimum_count_per_million = 10, design = design)
  se_keep_cpm100 <- tidybulk::keep_abundant(se, minimum_count_per_million = 100, design = design)
  expect_true(nrow(se_keep_cpm100) <= nrow(se_keep_10))

  se_keep_both <- tidybulk::keep_abundant(se, minimum_counts = 100, minimum_count_per_million = 100, design = design)
  expect_equal(nrow(se_keep_both), nrow(se_keep_cpm100))
}) 

test_that("identify_abundant and keep_abundant work with formula_design argument", {
  data(airway)
  se = airway
  # Use a formula for design
  formula <- ~ dex + cell

  # identify_abundant with formula_design
  se_abundant_10 <- tidybulk::identify_abundant(se, minimum_counts = 10, formula_design = formula)
  se_abundant_100 <- tidybulk::identify_abundant(se, minimum_counts = 100, formula_design = formula)
  expect_true(sum(rowData(se_abundant_100)$.abundant) <= sum(rowData(se_abundant_10)$.abundant))

  se_abundant_cpm10 <- tidybulk::identify_abundant(se, minimum_count_per_million = 10, formula_design = formula)
  se_abundant_cpm100 <- tidybulk::identify_abundant(se, minimum_count_per_million = 100, formula_design = formula)
  expect_true(sum(rowData(se_abundant_cpm100)$.abundant) <= sum(rowData(se_abundant_cpm10)$.abundant))

  se_abundant_both <- tidybulk::identify_abundant(se, minimum_counts = 100, minimum_count_per_million = 100, formula_design = formula)
  expect_equal(sum(rowData(se_abundant_both)$.abundant), sum(rowData(se_abundant_cpm100)$.abundant))

  # keep_abundant with formula_design
  se_keep_10 <- tidybulk::keep_abundant(se, minimum_counts = 10, formula_design = formula)
  se_keep_100 <- tidybulk::keep_abundant(se, minimum_counts = 100, formula_design = formula)
  expect_true(nrow(se_keep_100) <= nrow(se_keep_10))

  se_keep_cpm10 <- tidybulk::keep_abundant(se, minimum_count_per_million = 10, formula_design = formula)
  se_keep_cpm100 <- tidybulk::keep_abundant(se, minimum_count_per_million = 100, formula_design = formula)
  expect_true(nrow(se_keep_cpm100) <= nrow(se_keep_10))

  se_keep_both <- tidybulk::keep_abundant(se, minimum_counts = 100, minimum_count_per_million = 100, formula_design = formula)
  expect_equal(nrow(se_keep_both), nrow(se_keep_cpm100))

  # Precedence: formula_design overrides design
  design <- model.matrix(~ dex, data = colData(se))
  expect_warning(
    se_abundant_precedence <- tidybulk::identify_abundant(se, minimum_counts = 10, design = design, formula_design = formula),
    "Both formula_design and design provided; formula_design will be used."
  )
  expect_equal(sum(rowData(se_abundant_precedence)$.abundant), sum(rowData(se_abundant_10)$.abundant))
})

test_that("identify_abundant and keep_abundant work with NULL rownames", {
  # Simulates get_pseudobulk() -> as(\"SummarizedExperiment\") where coercion drops feature rownames
  set.seed(1)
  mat <- matrix(rnbinom(100, mu = 10, size = 1), nrow = 10, ncol = 10)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    colData = data.frame(group = rep(c("A", "B"), 5))
  )
  rownames(se) <- NULL

  res <- tidybulk::identify_abundant(se, minimum_counts = 5)
  expect_true(".abundant" %in% colnames(SummarizedExperiment::rowData(res)))
  expect_length(SummarizedExperiment::rowData(res)$.abundant, nrow(se))
  expect_type(SummarizedExperiment::rowData(res)$.abundant, "logical")

  res_keep <- tidybulk::keep_abundant(se, minimum_counts = 5)
  expect_true(nrow(res_keep) <= nrow(se))
}) 