testthat::test_that("Shiny filter selections propagate to the shared reactive dataset", {
  app_environment <- new.env(parent = globalenv())
  sys.source(project_path("app.R"), envir = app_environment)

  server_under_test <- function(input, output, session) {
    app_environment$app_server(
      input, output, session,
      default_metadata = fixture_path("sample_metadata.csv"),
      default_taxonomy = fixture_path("taxonomy_data.csv")
    )
  }

  shiny::testServer(server_under_test, {
      session$setInputs(
        group_variable = "",
        min_library_size = 0,
        min_prevalence = 0,
        min_total_abundance = 0,
        run_ordination = 0,
        run_differential = 0
      )
      state <- session$userData$microbiomeviz
      testthat::expect_true(state$data_result()$ok)
      testthat::expect_equal(ncol(state$filtered_data()$abundance), 6L)

      session$setInputs(group_variable = "Group", selected_groups = "B")
      filtered <- state$filtered_data()
      testthat::expect_identical(colnames(filtered$abundance), c("S6", "S5", "S4"))
      testthat::expect_identical(filtered$metadata$SampleID, colnames(filtered$abundance))

      invisible(try(state$ordination_result(), silent = TRUE))
      session$setInputs(run_ordination = 1)
      session$flushReact()
      result <- state$ordination_result()
      testthat::expect_identical(result$value$Sample, c("S6", "S5", "S4"))
      testthat::expect_identical(result$key, state$analysis_key())

      session$setInputs(selected_groups = "A")
      testthat::expect_false(identical(result$key, state$analysis_key()))

      invalid_taxonomy <- tempfile(fileext = ".csv")
      writeLines(c("FeatureID,S1,S2", "F1,1,not-numeric"), invalid_taxonomy)
      session$setInputs(taxonomy_file = list(
        name = "invalid.csv", size = file.info(invalid_taxonomy)$size,
        type = "text/csv", datapath = invalid_taxonomy
      ))
      invalid_result <- state$data_result()
      testthat::expect_false(invalid_result$ok)
      testthat::expect_match(invalid_result$message, "numeric and finite.*F1.*S2")
  })
})
