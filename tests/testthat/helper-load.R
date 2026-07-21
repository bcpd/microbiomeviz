project_path <- function(...) {
  testthat::test_path("..", "..", ...)
}

source(project_path("R", "data.R"))
source(project_path("R", "analysis.R"))

fixture_path <- function(file) {
  testthat::test_path("fixtures", file)
}

load_fixture_data <- function() {
  load_microbiome_data(
    fixture_path("sample_metadata.csv"),
    fixture_path("taxonomy_data.csv")
  )
}
