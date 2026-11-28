estimate_dispersion <- function(...) tidybulk:::estimate_dispersion(...)

simulated_se <- function(n_genes = 100) {
  set.seed(1)
  n_donor <- 6L
  n_rep <- 2L
  n <- n_donor * 2L * n_rep
  donor <- factor(rep(paste0("d", seq_len(n_donor)), each = 2L * n_rep))
  treatment <- factor(rep(rep(c("ctrl", "trt"), each = n_rep), n_donor))
  true_dispersion <- exp(stats::rnorm(n_genes, mean = log(0.1), sd = 0.4))
  mu <- 80 * (1 + 0.4 * as.integer(treatment == "trt"))
  counts <- matrix(
    stats::rnbinom(n_genes * n, mu = rep(mu, each = n_genes), size = 1 / true_dispersion),
    nrow = n_genes
  )
  colnames(counts) <- paste0("s", seq_len(n))
  rownames(counts) <- paste0("g", seq_len(n_genes))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(treatment = treatment, donor = donor),
    rowData = S4Vectors::DataFrame(true_dispersion = true_dispersion)
  )
}

test_that("estimate_dispersion writes tagwise phi", {
  skip_if_not_installed("edgeR")
  se <- simulated_se()
  expect_message(
    out <- estimate_dispersion(se, ~ treatment + donor),
    "tagwise \\(shrinked\\) and trended dispersion"
  )
  rd <- SummarizedExperiment::rowData(out)
  expect_true(all(is.finite(rd$dispersion_shrinked) & rd$dispersion_shrinked > 0))
  expect_true(all(is.finite(rd$dispersion_trended) & rd$dispersion_trended > 0))
})

test_that("estimate_dispersion writes one shared phi when asked for common", {
  skip_if_not_installed("edgeR")
  se <- simulated_se()
  expect_message(
    out <- estimate_dispersion(se, ~ treatment + donor, dispersion_trend = "common"),
    "common dispersion"
  )
  rd <- SummarizedExperiment::rowData(out)
  tagwise <- SummarizedExperiment::rowData(
    estimate_dispersion(se, ~ treatment + donor)
  )$dispersion_shrinked
  expect_length(unique(rd$dispersion_shrinked), 1)
  expect_gt(length(unique(tagwise)), 1)
  expect_false(isTRUE(all.equal(tagwise, rd$dispersion_shrinked)))
  expect_true(all(is.finite(rd$dispersion_shrinked) & rd$dispersion_shrinked > 0))
  expect_equal(rd$dispersion_trended, rd$dispersion_shrinked)
  expect_length(rd$dispersion_shrinked, nrow(se))
})

test_that("a mixed-model formula is an error", {
  skip_if_not_installed("edgeR")
  se <- simulated_se()
  expect_error(
    estimate_dispersion(se, ~ treatment + (1 | donor)),
    "no random effects"
  )
})

test_that("a saturated design is an error", {
  skip_if_not_installed("edgeR")
  se <- simulated_se()
  se$sid <- factor(colnames(se))
  expect_error(
    estimate_dispersion(se, ~ treatment + sid),
    "no residual degrees of freedom"
  )
})
