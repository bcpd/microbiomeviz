if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Tests require the 'testthat' package.", call. = FALSE)
}

testthat::test_dir("tests/testthat", reporter = "summary")
