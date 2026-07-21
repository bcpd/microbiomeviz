is_blank <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

read_metadata_file <- function(path) {
  if (!file.exists(path)) {
    stop("Metadata file does not exist: ", path, call. = FALSE)
  }

  metadata <- tryCatch(
    utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA")
    ),
    error = function(e) stop("Could not read metadata CSV: ", conditionMessage(e), call. = FALSE)
  )

  if (ncol(metadata) < 1L) {
    stop("Metadata must contain a sample-ID column.", call. = FALSE)
  }

  notes <- character()
  blank_rows <- apply(metadata, 1L, function(row) all(is_blank(row)))
  if (any(blank_rows)) {
    metadata <- metadata[!blank_rows, , drop = FALSE]
    notes <- c(notes, sprintf("Removed %d completely blank metadata row(s).", sum(blank_rows)))
  }

  blank_unnamed_columns <- vapply(
    seq_len(ncol(metadata)),
    function(i) is_blank(names(metadata)[i]) && all(is_blank(metadata[[i]])),
    logical(1)
  )
  if (any(blank_unnamed_columns)) {
    metadata <- metadata[, !blank_unnamed_columns, drop = FALSE]
    notes <- c(
      notes,
      sprintf("Removed %d unnamed, completely blank metadata column(s).", sum(blank_unnamed_columns))
    )
  }

  names(metadata) <- trimws(names(metadata))
  if (any(!nzchar(names(metadata)))) {
    stop("Every non-empty metadata column must have a name.", call. = FALSE)
  }
  if (anyDuplicated(names(metadata))) {
    duplicate_names <- unique(names(metadata)[duplicated(names(metadata))])
    stop("Metadata column names must be unique. Duplicated: ",
         paste(duplicate_names, collapse = ", "), call. = FALSE)
  }

  character_columns <- vapply(metadata, is.character, logical(1))
  whitespace_changed <- 0L
  for (column in names(metadata)[character_columns]) {
    original <- metadata[[column]]
    trimmed <- trimws(original)
    whitespace_changed <- whitespace_changed + sum(original != trimmed, na.rm = TRUE)
    metadata[[column]] <- trimmed
  }
  if (whitespace_changed > 0L) {
    notes <- c(notes, sprintf("Trimmed surrounding whitespace in %d metadata value(s).", whitespace_changed))
  }

  sample_id_column <- names(metadata)[1L]
  sample_ids <- metadata[[sample_id_column]]
  if (nrow(metadata) < 1L) {
    stop("Metadata contains no non-empty sample rows.", call. = FALSE)
  }
  if (any(is_blank(sample_ids))) {
    stop("Metadata sample IDs (first column, '", sample_id_column,
         "') must not be blank.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    duplicated_ids <- unique(sample_ids[duplicated(sample_ids)])
    stop("Metadata sample IDs must be unique. Duplicated: ",
         paste(duplicated_ids, collapse = ", "), call. = FALSE)
  }

  list(data = metadata, sample_id_column = sample_id_column, notes = notes)
}

read_taxonomy_file <- function(path) {
  if (!file.exists(path)) {
    stop("Taxonomy/abundance file does not exist: ", path, call. = FALSE)
  }

  taxonomy <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) stop("Could not read taxonomy/abundance CSV: ", conditionMessage(e), call. = FALSE)
  )
  if (ncol(taxonomy) < 2L) {
    stop("Taxonomy/abundance data must contain a feature-ID column followed by at least one sample column.",
         call. = FALSE)
  }
  if (nrow(taxonomy) < 1L) {
    stop("Taxonomy/abundance data contains no features.", call. = FALSE)
  }

  feature_column <- names(taxonomy)[1L]
  if (is_blank(feature_column)) {
    stop("The first taxonomy/abundance column (feature IDs) must have a name.", call. = FALSE)
  }
  feature_ids <- trimws(as.character(taxonomy[[1L]]))
  sample_ids <- trimws(names(taxonomy)[-1L])
  if (any(!nzchar(feature_ids))) {
    stop("Feature IDs (first taxonomy/abundance column) must not be blank.", call. = FALSE)
  }
  if (anyDuplicated(feature_ids)) {
    duplicated_ids <- unique(feature_ids[duplicated(feature_ids)])
    stop("Feature IDs must be unique. Duplicated: ",
         paste(utils::head(duplicated_ids, 5L), collapse = ", "), call. = FALSE)
  }
  if (any(!nzchar(sample_ids))) {
    stop("Taxonomy/abundance sample column names must not be blank.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    duplicated_ids <- unique(sample_ids[duplicated(sample_ids)])
    stop("Taxonomy/abundance sample names must be unique. Duplicated: ",
         paste(duplicated_ids, collapse = ", "), call. = FALSE)
  }

  abundance_columns <- taxonomy[-1L]
  abundance <- vapply(
    abundance_columns,
    function(column) suppressWarnings(as.numeric(column)),
    numeric(nrow(taxonomy))
  )
  if (is.null(dim(abundance))) {
    abundance <- matrix(
      abundance,
      nrow = nrow(taxonomy),
      ncol = length(abundance_columns)
    )
  }
  rownames(abundance) <- feature_ids
  colnames(abundance) <- sample_ids

  invalid <- !is.finite(abundance)
  if (any(invalid)) {
    first_invalid <- which(invalid, arr.ind = TRUE)[1L, ]
    stop(
      sprintf(
        "Abundances must be numeric and finite. Invalid value for feature '%s', sample '%s'.",
        rownames(abundance)[first_invalid[1L]], colnames(abundance)[first_invalid[2L]]
      ),
      call. = FALSE
    )
  }
  if (any(abundance < 0)) {
    first_negative <- which(abundance < 0, arr.ind = TRUE)[1L, ]
    stop(
      sprintf(
        "Abundances must be non-negative. Negative value for feature '%s', sample '%s'.",
        rownames(abundance)[first_negative[1L]], colnames(abundance)[first_negative[2L]]
      ),
      call. = FALSE
    )
  }

  storage.mode(abundance) <- "double"
  list(data = abundance, feature_column = feature_column, notes = character())
}

load_microbiome_data <- function(metadata_path, taxonomy_path) {
  metadata_input <- read_metadata_file(metadata_path)
  taxonomy_input <- read_taxonomy_file(taxonomy_path)

  metadata <- metadata_input$data
  abundance <- taxonomy_input$data
  sample_id_column <- metadata_input$sample_id_column
  metadata_ids <- as.character(metadata[[sample_id_column]])
  matched_ids <- colnames(abundance)[colnames(abundance) %in% metadata_ids]

  if (length(matched_ids) < 2L) {
    stop(
      sprintf(
        paste0(
          "Metadata and taxonomy/abundance data have only %d matching sample ID(s). ",
          "IDs are matched between the first metadata column ('%s') and abundance column names; at least 2 are required."
        ),
        length(matched_ids), sample_id_column
      ),
      call. = FALSE
    )
  }

  taxonomy_only <- setdiff(colnames(abundance), metadata_ids)
  metadata_only <- setdiff(metadata_ids, colnames(abundance))
  notes <- c(metadata_input$notes, taxonomy_input$notes)
  if (length(taxonomy_only)) {
    notes <- c(notes, sprintf(
      "Excluded %d abundance sample(s) without metadata: %s.",
      length(taxonomy_only), paste(utils::head(taxonomy_only, 5L), collapse = ", ")
    ))
  }
  if (length(metadata_only)) {
    notes <- c(notes, sprintf(
      "Excluded %d metadata sample(s) without abundance data: %s.",
      length(metadata_only), paste(utils::head(metadata_only, 5L), collapse = ", ")
    ))
  }

  abundance <- abundance[, matched_ids, drop = FALSE]
  metadata <- metadata[match(matched_ids, metadata_ids), , drop = FALSE]
  rownames(metadata) <- matched_ids

  list(
    metadata = metadata,
    abundance = abundance,
    input_library_sizes = colSums(abundance),
    sample_id_column = sample_id_column,
    feature_column = taxonomy_input$feature_column,
    notes = notes,
    excluded = list(taxonomy_only = taxonomy_only, metadata_only = metadata_only)
  )
}

filter_microbiome_data <- function(data, group_variable = NULL, selected_groups = NULL,
                                   min_library_size = 0, min_prevalence_percent = 0,
                                   min_total_abundance = 0) {
  stopifnot(is.list(data), is.matrix(data$abundance), is.data.frame(data$metadata))
  if (!is.numeric(min_library_size) || length(min_library_size) != 1L ||
      !is.finite(min_library_size) || min_library_size < 0) {
    stop("Minimum library size must be one non-negative number.", call. = FALSE)
  }
  if (!is.numeric(min_prevalence_percent) || length(min_prevalence_percent) != 1L ||
      !is.finite(min_prevalence_percent) || min_prevalence_percent < 0 ||
      min_prevalence_percent > 100) {
    stop("Minimum prevalence must be between 0 and 100 percent.", call. = FALSE)
  }
  if (!is.numeric(min_total_abundance) || length(min_total_abundance) != 1L ||
      !is.finite(min_total_abundance) || min_total_abundance < 0) {
    stop("Minimum total abundance must be one non-negative number.", call. = FALSE)
  }

  keep_samples <- rep(TRUE, ncol(data$abundance))
  if (!is.null(group_variable) && nzchar(group_variable)) {
    if (!group_variable %in% names(data$metadata)) {
      stop("Selected grouping variable is not present in metadata: ", group_variable, call. = FALSE)
    }
    group_values <- as.character(data$metadata[[group_variable]])
    group_values[is.na(group_values) | !nzchar(group_values)] <- "__MISSING__"
    if (!is.null(selected_groups)) {
      keep_samples <- keep_samples & group_values %in% selected_groups
    }
  }

  library_sizes <- colSums(data$abundance)
  keep_samples <- keep_samples & library_sizes >= min_library_size
  abundance <- data$abundance[, keep_samples, drop = FALSE]
  metadata <- data$metadata[keep_samples, , drop = FALSE]

  if (ncol(abundance) == 0L) {
    keep_features <- rep(FALSE, nrow(abundance))
  } else {
    prevalence <- rowMeans(abundance > 0) * 100
    keep_features <- prevalence >= min_prevalence_percent &
      rowSums(abundance) >= min_total_abundance
  }
  abundance <- abundance[keep_features, , drop = FALSE]

  data$metadata <- metadata
  data$abundance <- abundance
  data$filter_summary <- list(
    samples_before = length(keep_samples),
    samples_after = sum(keep_samples),
    features_before = length(keep_features),
    features_after = sum(keep_features),
    group_variable = group_variable,
    selected_groups = selected_groups,
    min_library_size = min_library_size,
    min_prevalence_percent = min_prevalence_percent,
    min_total_abundance = min_total_abundance
  )
  data
}

write_project_bundle <- function(data, file) {
  if (ncol(data$abundance) < 1L || nrow(data$abundance) < 1L) {
    stop("Cannot download an empty filtered dataset.", call. = FALSE)
  }

  bundle_dir <- tempfile("microbiomeviz-bundle-")
  dir.create(bundle_dir)
  on.exit(unlink(bundle_dir, recursive = TRUE), add = TRUE)

  metadata_path <- file.path(bundle_dir, "filtered_sample_metadata.csv")
  taxonomy_path <- file.path(bundle_dir, "filtered_taxonomy_data.csv")
  notes_path <- file.path(bundle_dir, "analysis_notes.txt")
  utils::write.csv(data$metadata, metadata_path, row.names = FALSE, na = "")
  taxonomy_export <- data.frame(
    feature = rownames(data$abundance), data$abundance,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  names(taxonomy_export)[1L] <- data$feature_column
  utils::write.csv(taxonomy_export, taxonomy_path, row.names = FALSE)

  filter_summary <- data$filter_summary
  notes <- c(
    "Microbiome Visualizer filtered project export",
    sprintf("Created: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "",
    "Data alignment and cleaning messages:",
    if (length(data$notes)) paste0("- ", data$notes) else "- None",
    "",
    "Active filters:",
    sprintf("- Samples: %d of %d", filter_summary$samples_after, filter_summary$samples_before),
    sprintf("- Features: %d of %d", filter_summary$features_after, filter_summary$features_before),
    sprintf("- Grouping variable: %s", filter_summary$group_variable %||% "None"),
    sprintf("- Selected groups: %s", paste(filter_summary$selected_groups %||% "All", collapse = ", ")),
    sprintf("- Minimum library size: %s", filter_summary$min_library_size),
    sprintf("- Minimum prevalence: %s%% (abundance > 0)", filter_summary$min_prevalence_percent),
    sprintf("- Minimum total abundance: %s", filter_summary$min_total_abundance),
    "",
    "Scientific defaults:",
    "- Filtering uses raw, non-negative abundances.",
    "- No normalization or transformation is applied to exported values."
  )
  writeLines(notes, notes_path)

  files <- c(metadata_path, taxonomy_path, notes_path)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zipr(file, files = files, root = bundle_dir)
  } else {
    old_wd <- setwd(bundle_dir)
    on.exit(setwd(old_wd), add = TRUE)
    status <- utils::zip(file, basename(files), flags = "-q")
    if (!identical(status, 0L)) {
      stop("Could not create ZIP archive. Install the 'zip' R package and retry.", call. = FALSE)
    }
  }
  invisible(file)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}
