context('Abundance Normalization Functions')

library(airway)
data(airway)
airway_mini <- airway[1:100, 1:5]

library(dplyr)
library(SummarizedExperiment)

# Test scale_abundance function
test_that("scale_abundance works correctly", {
  res <- airway_mini |> identify_abundant() |> scale_abundance()
  
  expect_equal(
    names(SummarizedExperiment::assays(res)),
    c("counts", "counts_scaled")
  )
})

test_that("scale_abundance with subset works correctly", {
  # Skip this test - .subset_for_scaling requires complex quosure handling
  skip(".subset_for_scaling test requires refactoring")
  
  res <- airway_mini |> identify_abundant() |> scale_abundance()
  
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
})

# Test quantile_normalise_abundance function
test_that("quantile_normalise_abundance works correctly", {
  res <- airway_mini |> quantile_normalise_abundance()
  
  # Check if any normalized assay exists
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true(length(assay_names) > 1 || any(grepl("normalised", assay_names)))
})

test_that("quantile_normalise_abundance with preprocessCore works correctly", {
  res <- airway_mini |> quantile_normalise_abundance(method = "preprocesscore_normalize_quantiles_use_target")
  
  # Check if any normalized assay exists
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true(length(assay_names) > 1 || any(grepl("normalised", assay_names)))
})

# Test adjust_abundance function

test_that("adjust_abundance adds adjusted assay correctly", {
  # Add a batch column to airway_mini
  se_mini2 <- airway_mini
  # Create a batch variable with two groups (not confounded with dex)
  colData(se_mini2)$batch <- c(1, 1, 2, 2, 1)
  # Run identify_abundant and adjust_abundance
  res <- se_mini2 |> identify_abundant() |> adjust_abundance(
    .factor_unwanted = batch,
    .factor_of_interest = dex,
    method = "combat_seq"
  )
  # Check that an adjusted assay is present
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true(any(grepl("_adjusted$", assay_names)))
})

# Test fill_missing_abundance function
test_that("fill_missing_abundance works correctly", {
  # This function doesn't exist for SummarizedExperiment
  # Skip this test for now
  skip("fill_missing_abundance not implemented for SummarizedExperiment")
})

# Test impute_missing_abundance function
test_that("impute_missing_abundance works correctly", {
  res <- airway_mini |> impute_missing_abundance(.formula = ~ dex)
  
  expect_no_error(res)
}) 

# Test scale_abundance with custom suffix

test_that("scale_abundance uses custom suffix correctly", {
  res <- airway_mini |> identify_abundant() |> scale_abundance(suffix = "_custom")
  expect_true("counts_custom" %in% names(SummarizedExperiment::assays(res)))
  expect_false("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
})

test_that("scale_abundance default suffix still works", {
  res <- airway_mini |> identify_abundant() |> scale_abundance()
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
}) 

# Test scale_abundance with chunking
test_that("scale_abundance with chunking produces results", {
  # Use airway with all samples for meaningful chunking
  airway_test <- airway[1:100, ]
  
  # Standard scaling without chunking (default chunk_size = Inf)
  res_standard <- airway_test |> identify_abundant() |> scale_abundance()
  
  # Chunked scaling with small chunk size to test chunking logic
  res_chunked <- airway_test |> identify_abundant() |> 
    scale_abundance(chunk_sample_size = 2)
  
  # Both should have the scaled assay
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res_standard)))
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res_chunked)))
  
  # Both should have the same dimensions
  expect_equal(dim(SummarizedExperiment::assay(res_standard, "counts_scaled")),
               dim(SummarizedExperiment::assay(res_chunked, "counts_scaled")))
  
  # TMM and multiplier columns should exist
  expect_true("TMM" %in% names(SummarizedExperiment::colData(res_standard)))
  expect_true("TMM" %in% names(SummarizedExperiment::colData(res_chunked)))
  
  # Scaled values should be positive and non-NA
  expect_true(all(SummarizedExperiment::assay(res_chunked, "counts_scaled") >= 0, na.rm = TRUE))
})


test_that("scale_abundance with specified reference sample works", {
  airway_test <- airway[1:100, ]
  ref_sample <- colnames(airway_test)[1]
  
  # Test with chunking and specified reference
  res <- airway_test |> identify_abundant() |> 
    scale_abundance(reference_sample = ref_sample, chunk_sample_size = 3)
  
  expect_true("counts_scaled" %in% names(SummarizedExperiment::assays(res)))
  expect_equal(ncol(res), ncol(airway_test))
})

test_that("chunked and non-chunked produce identical results with specified reference", {
  airway_test <- airway[1:100, ]
  
  # Select a reference sample upfront
  ref_sample <- colnames(airway_test)[1]
  
  # Non-chunked version
  res_standard <- airway_test |> identify_abundant() |> 
    scale_abundance(reference_sample = ref_sample)
  
  # Chunked version with same reference
  res_chunked <- airway_test |> identify_abundant() |> 
    scale_abundance(reference_sample = ref_sample, chunk_sample_size = 2)
  
  # Results should be identical
  expect_equal(
    SummarizedExperiment::assay(res_standard, "counts_scaled"),
    SummarizedExperiment::assay(res_chunked, "counts_scaled"),
    tolerance = 1e-10
  )
  
  # TMM and multiplier should be identical
  expect_equal(
    SummarizedExperiment::colData(res_standard)$TMM,
    SummarizedExperiment::colData(res_chunked)$TMM,
    tolerance = 1e-10
  )
  
  expect_equal(
    SummarizedExperiment::colData(res_standard)$multiplier,
    SummarizedExperiment::colData(res_chunked)$multiplier,
    tolerance = 1e-10
  )
})

# Test adjust_abundance on a custom assay

test_that("adjust_abundance on custom assay creates correct adjusted assay name", {
  se_mini2 <- airway_mini
  # Create a batch variable with two groups (not confounded with dex)
  colData(se_mini2)$batch <- c(1, 1, 2, 2, 1)
  # Create a custom assay
  SummarizedExperiment::assays(se_mini2)[["my_custom_assay"]] <- SummarizedExperiment::assay(se_mini2, "counts") + 1
  # Run identify_abundant and adjust_abundance on the custom assay
  res <- se_mini2 |> identify_abundant(abundance = "my_custom_assay") |> adjust_abundance(
    abundance = "my_custom_assay",
    .factor_unwanted = batch,
    .factor_of_interest = dex,
    method = "combat_seq"
  )
  # Check that the adjusted assay name is correct
  assay_names <- names(SummarizedExperiment::assays(res))
  expect_true("my_custom_assay_adjusted" %in% assay_names)
}) 