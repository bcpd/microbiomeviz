testthat::test_that("valid inputs are cleaned, matched, and deterministically ordered", {
  data <- load_fixture_data()

  testthat::expect_identical(colnames(data$abundance), c("S3", "S1", "S6", "S2", "S5", "S4"))
  testthat::expect_identical(as.character(data$metadata$SampleID), colnames(data$abundance))
  testthat::expect_identical(data$metadata$Group, c("A", "A", "B", "A", "B", "B"))
  testthat::expect_identical(data$excluded$taxonomy_only, "Control")
  testthat::expect_identical(data$excluded$metadata_only, "S7")
  testthat::expect_true(any(grepl("Trimmed surrounding whitespace", data$notes)))
  testthat::expect_true(any(grepl("without metadata", data$notes)))
})

testthat::test_that("invalid abundance values produce feature- and sample-specific errors", {
  taxonomy <- utils::read.csv(fixture_path("taxonomy_data.csv"), check.names = FALSE)

  taxonomy[2L, "S1"] <- -1
  negative_path <- tempfile(fileext = ".csv")
  utils::write.csv(taxonomy, negative_path, row.names = FALSE, quote = FALSE)
  testthat::expect_error(
    read_taxonomy_file(negative_path),
    "Negative value.*g__Beta.*S1"
  )

  taxonomy[2L, "S1"] <- "not-a-number"
  nonnumeric_path <- tempfile(fileext = ".csv")
  utils::write.csv(taxonomy, nonnumeric_path, row.names = FALSE, quote = FALSE)
  testthat::expect_error(
    read_taxonomy_file(nonnumeric_path),
    "numeric and finite.*g__Beta.*S1"
  )
})

testthat::test_that("duplicate and unmatched sample IDs are rejected clearly", {
  metadata <- utils::read.csv(fixture_path("sample_metadata.csv"), check.names = FALSE)
  metadata$SampleID[2L] <- metadata$SampleID[1L]
  duplicate_path <- tempfile(fileext = ".csv")
  utils::write.csv(metadata, duplicate_path, row.names = FALSE, quote = FALSE)
  testthat::expect_error(read_metadata_file(duplicate_path), "sample IDs must be unique.*S1")

  unmatched <- metadata
  unmatched$SampleID <- paste0("X", seq_len(nrow(unmatched)))
  unmatched_path <- tempfile(fileext = ".csv")
  utils::write.csv(unmatched, unmatched_path, row.names = FALSE, quote = FALSE)
  testthat::expect_error(
    load_microbiome_data(unmatched_path, fixture_path("taxonomy_data.csv")),
    "only 0 matching sample ID.*at least 2"
  )
})
