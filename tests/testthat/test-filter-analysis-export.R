testthat::test_that("one filtered dataset drives metadata, abundance, and ordination", {
  data <- load_fixture_data()
  filtered <- filter_microbiome_data(
    data,
    group_variable = "Group",
    selected_groups = "A",
    min_library_size = 15,
    min_prevalence_percent = 1,
    min_total_abundance = 1
  )

  testthat::expect_identical(colnames(filtered$abundance), c("S3", "S1", "S2"))
  testthat::expect_identical(filtered$metadata$SampleID, colnames(filtered$abundance))
  testthat::expect_false(any(grepl("g__Zero$", rownames(filtered$abundance))))
  testthat::expect_identical(filtered$filter_summary$samples_after, 3L)

  ordination <- calculate_ordination(filtered)
  testthat::expect_identical(ordination$Sample, colnames(filtered$abundance))
  testthat::expect_equal(nrow(ordination), 3L)
  testthat::expect_length(attr(ordination, "explained"), 2L)
})

testthat::test_that("default filters preserve all matched samples and features", {
  data <- load_fixture_data()
  filtered <- filter_microbiome_data(data)

  testthat::expect_identical(filtered$abundance, data$abundance)
  testthat::expect_identical(filtered$metadata, data$metadata)
})

testthat::test_that("heatmap transformations are display-only and deterministic", {
  data <- filter_microbiome_data(load_fixture_data())
  original <- data$abundance
  raw <- prepare_heatmap_matrix(data, top_n = 3, transform = "raw")
  relative <- prepare_heatmap_matrix(data, top_n = 3, transform = "relative")

  testthat::expect_equal(nrow(raw), 3L)
  expected_relative <- sweep(
    original[rownames(relative), , drop = FALSE], 2L, colSums(original), "/"
  )
  testthat::expect_equal(relative, expected_relative, tolerance = 1e-12)
  testthat::expect_identical(data$abundance, original)
})

testthat::test_that("timeseries data retain selected microbial features without averaging", {
  data <- filter_microbiome_data(load_fixture_data())
  features <- rownames(data$abundance)[1:2]
  series <- prepare_timeseries_data(
    data,
    time_variable = "SamplingDate",
    features = features,
    transform = "relative",
    colour_variable = "Group"
  )

  testthat::expect_equal(nrow(series), 2L * ncol(data$abundance))
  testthat::expect_setequal(unique(series$Feature), features)
  testthat::expect_identical(series$Sample[1:6], data$metadata$SampleID)
  testthat::expect_true(all(series$Abundance >= 0 & series$Abundance <= 1))
  testthat::expect_identical(series$.colour[1:6], as.character(data$metadata$Group))

  same_time <- data
  same_time$metadata$SamplingDate <- "2024-01-01"
  testthat::expect_error(
    prepare_timeseries_data(same_time, "SamplingDate", features),
    "fewer than two distinct sampling times"
  )
})

testthat::test_that("download bundle contains matching filtered samples and notes", {
  filtered <- filter_microbiome_data(
    load_fixture_data(), group_variable = "Group", selected_groups = "B"
  )
  archive <- tempfile(fileext = ".zip")
  write_project_bundle(filtered, archive)
  testthat::expect_true(file.exists(archive))

  destination <- tempfile("unzip-")
  dir.create(destination)
  utils::unzip(archive, exdir = destination)
  testthat::expect_setequal(
    list.files(destination),
    c("analysis_notes.txt", "filtered_sample_metadata.csv", "filtered_taxonomy_data.csv")
  )
  metadata <- utils::read.csv(file.path(destination, "filtered_sample_metadata.csv"), check.names = FALSE)
  taxonomy <- utils::read.csv(file.path(destination, "filtered_taxonomy_data.csv"), check.names = FALSE)
  testthat::expect_identical(as.character(metadata$SampleID), names(taxonomy)[-1L])
  testthat::expect_match(paste(readLines(file.path(destination, "analysis_notes.txt")), collapse = "\n"),
                         "Selected groups: B")
})

testthat::test_that("analysis prerequisites fail with actionable messages", {
  data <- filter_microbiome_data(load_fixture_data(), group_variable = "Group", selected_groups = "A")
  data$abundance <- data$abundance[, 1:2, drop = FALSE]
  data$metadata <- data$metadata[1:2, , drop = FALSE]
  testthat::expect_error(calculate_ordination(data), "requires at least 3 samples.*2 remain")

  zero_data <- filter_microbiome_data(load_fixture_data())
  zero_data$abundance[, 1L] <- 0
  testthat::expect_error(calculate_ordination(zero_data), "cannot include zero-sum samples.*S3")

  if (!requireNamespace("maaslin3", quietly = TRUE)) {
    testthat::expect_error(
      run_maaslin3(
        filter_microbiome_data(load_fixture_data()), "Group", tempfile(),
        reference_levels = c(Group = "A")
      ),
      "Bioconductor package 'maaslin3'.*BiocManager::install"
    )
  }
})

testthat::test_that("MaAsLin 3 inputs are aligned and read depth remains pre-filter", {
  data <- load_fixture_data()
  filtered <- filter_microbiome_data(
    data,
    group_variable = "Group",
    selected_groups = "B",
    min_total_abundance = 5
  )
  inputs <- prepare_maaslin3_inputs(filtered, "Age", include_read_depth = TRUE)

  testthat::expect_identical(rownames(inputs$input_data), colnames(filtered$abundance))
  testthat::expect_identical(colnames(inputs$input_data), rownames(filtered$abundance))
  testthat::expect_identical(rownames(inputs$input_metadata), colnames(filtered$abundance))
  testthat::expect_false(filtered$sample_id_column %in% names(inputs$input_metadata))
  testthat::expect_identical(inputs$fixed_effects, c("Age", "sample_read_depth"))
  expected_depth <- unname(data$input_library_sizes[colnames(filtered$abundance)])
  testthat::expect_equal(inputs$input_metadata$sample_read_depth, expected_depth)

  testthat::expect_error(
    prepare_maaslin3_inputs(filtered, "Group"),
    "must have at least two distinct values"
  )

  all_groups <- filter_microbiome_data(data)
  all_groups$metadata$Group[all_groups$metadata$SampleID == "S6"] <- "C"
  testthat::expect_error(
    prepare_maaslin3_inputs(all_groups, "Group"),
    "Choose a reference level.*more than two levels"
  )
  referenced <- prepare_maaslin3_inputs(
    all_groups, "Group", reference_levels = c(Group = "B")
  )
  testthat::expect_identical(referenced$reference, "Group,B")
  testthat::expect_identical(levels(referenced$input_metadata$Group)[1L], "B")
})

testthat::test_that("MaAsLin 3 integration returns abundance and prevalence results", {
  testthat::skip_if_not_installed("maaslin3")
  data <- filter_microbiome_data(load_fixture_data(), min_prevalence_percent = 1)
  output_directory <- tempfile("maaslin3-test-")
  on.exit(unlink(output_directory, recursive = TRUE), add = TRUE)

  results <- run_maaslin3(
    data, "Group", output_directory,
    reference_levels = c(Group = "A"), include_read_depth = FALSE
  )
  expected_columns <- c(
    "feature", "name", "coef", "qval_individual", "qval_joint", "model"
  )
  testthat::expect_s3_class(results$all_results, "data.frame")
  testthat::expect_true(all(expected_columns %in% names(results$all_results)))
  testthat::expect_setequal(unique(results$all_results$model), c("abundance", "prevalence"))
  testthat::expect_true(file.exists(file.path(output_directory, "significant_results.tsv")))
})
