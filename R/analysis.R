validate_analysis_data <- function(data, minimum_samples = 2L, minimum_features = 1L,
                                   analysis_name = "Analysis") {
  if (ncol(data$abundance) < minimum_samples) {
    stop(sprintf("%s requires at least %d samples after filtering; %d remain.",
                 analysis_name, minimum_samples, ncol(data$abundance)), call. = FALSE)
  }
  if (nrow(data$abundance) < minimum_features) {
    stop(sprintf("%s requires at least %d features after filtering; %d remain.",
                 analysis_name, minimum_features, nrow(data$abundance)), call. = FALSE)
  }
  invisible(TRUE)
}

calculate_ordination <- function(data) {
  validate_analysis_data(data, minimum_samples = 3L, analysis_name = "PCoA")
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("PCoA requires the 'vegan' R package. Install it with install.packages('vegan').",
         call. = FALSE)
  }
  if (any(colSums(data$abundance) == 0)) {
    zero_samples <- colnames(data$abundance)[colSums(data$abundance) == 0]
    stop("Bray-Curtis PCoA cannot include zero-sum samples. Adjust filters or remove: ",
         paste(zero_samples, collapse = ", "), call. = FALSE)
  }
  if (all(apply(data$abundance, 1L, function(x) length(unique(x)) == 1L))) {
    stop("PCoA requires abundance variation among samples; all features are constant.", call. = FALSE)
  }

  distance <- vegan::vegdist(t(data$abundance), method = "bray")
  if (any(!is.finite(distance)) || all(distance == 0)) {
    stop("Bray-Curtis distances are undefined or all zero for the filtered data.", call. = FALSE)
  }
  pcoa <- stats::cmdscale(distance, eig = TRUE, k = 2, add = FALSE)
  positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
  explained <- if (length(positive_eigenvalues)) {
    100 * pcoa$eig[seq_len(2L)] / sum(positive_eigenvalues)
  } else {
    c(NA_real_, NA_real_)
  }
  points <- data.frame(
    Sample = rownames(pcoa$points),
    Axis1 = pcoa$points[, 1L],
    Axis2 = pcoa$points[, 2L],
    stringsAsFactors = FALSE
  )
  metadata <- data$metadata[match(points$Sample, data$metadata[[data$sample_id_column]]), , drop = FALSE]
  cbind(points, metadata[, setdiff(names(metadata), data$sample_id_column), drop = FALSE]) |>
    structure(explained = explained, eigenvalues = pcoa$eig)
}

prepare_heatmap_matrix <- function(data, top_n = 50L, transform = c("raw", "relative")) {
  transform <- match.arg(transform)
  validate_analysis_data(data, minimum_samples = 1L, analysis_name = "Heatmap")
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
    stop("Heatmap feature limit must be a positive number.", call. = FALSE)
  }

  abundance <- data$abundance
  if (transform == "relative") {
    sample_totals <- colSums(abundance)
    if (any(sample_totals == 0)) {
      stop("Relative-abundance heatmap cannot include zero-sum samples.", call. = FALSE)
    }
    abundance <- sweep(abundance, 2L, sample_totals, "/")
  }
  feature_order <- order(rowSums(abundance), decreasing = TRUE)
  abundance[utils::head(feature_order, min(as.integer(top_n), length(feature_order))), , drop = FALSE]
}

parse_time_variable <- function(values, variable_name) {
  parsed_time <- tryCatch(
    suppressWarnings(as.POSIXct(values, tz = "UTC")),
    error = function(e) rep(as.POSIXct(NA), length(values))
  )
  if (all(is.na(parsed_time))) {
    parsed_time <- tryCatch(
      suppressWarnings(as.Date(values)),
      error = function(e) rep(as.Date(NA), length(values))
    )
  }
  if (any(is.na(parsed_time))) {
    stop(sprintf("'%s' contains values that cannot be parsed consistently as dates/times.",
                 variable_name), call. = FALSE)
  }
  parsed_time
}

prepare_timeseries_data <- function(data, time_variable, features,
                                    transform = c("raw", "relative"),
                                    colour_variable = NULL, trajectory_variable = NULL) {
  transform <- match.arg(transform)
  validate_analysis_data(data, minimum_samples = 2L, analysis_name = "Timeseries")
  if (length(time_variable) != 1L || !time_variable %in% names(data$metadata)) {
    stop("Choose one valid date/time metadata column.", call. = FALSE)
  }
  if (!length(features) || any(!features %in% rownames(data$abundance))) {
    stop("Choose one or more features present in the filtered abundance table.", call. = FALSE)
  }
  if (!is.null(colour_variable) && nzchar(colour_variable) &&
      !colour_variable %in% names(data$metadata)) {
    stop("Selected colour variable is not present in metadata.", call. = FALSE)
  }
  if (!is.null(trajectory_variable) && nzchar(trajectory_variable) &&
      !trajectory_variable %in% names(data$metadata)) {
    stop("Selected trajectory ID is not present in metadata.", call. = FALSE)
  }

  parsed_time <- parse_time_variable(data$metadata[[time_variable]], time_variable)
  if (length(unique(parsed_time)) < 2L) {
    stop(sprintf("'%s' has fewer than two distinct sampling times after filtering; a temporal trend cannot be shown.",
                 time_variable), call. = FALSE)
  }

  abundance <- data$abundance
  if (transform == "relative") {
    sample_totals <- colSums(abundance)
    if (any(sample_totals == 0)) {
      stop("Relative-abundance timeseries cannot include zero-sum samples.", call. = FALSE)
    }
    abundance <- sweep(abundance, 2L, sample_totals, "/")
  }

  selected_abundance <- abundance[features, , drop = FALSE]
  sample_count <- ncol(selected_abundance)
  feature_count <- nrow(selected_abundance)
  plot_data <- data.frame(
    Time = rep(parsed_time, times = feature_count),
    Feature = rep(rownames(selected_abundance), each = sample_count),
    Abundance = as.vector(t(selected_abundance)),
    Sample = rep(data$metadata[[data$sample_id_column]], times = feature_count),
    stringsAsFactors = FALSE
  )

  if (!is.null(colour_variable) && nzchar(colour_variable)) {
    plot_data$.colour <- rep(as.character(data$metadata[[colour_variable]]), times = feature_count)
  }
  if (!is.null(trajectory_variable) && nzchar(trajectory_variable)) {
    trajectory <- as.character(data$metadata[[trajectory_variable]])
    missing_trajectory <- is.na(trajectory) | !nzchar(trajectory)
    trajectory[missing_trajectory] <- paste0("sample:", data$metadata[[data$sample_id_column]][missing_trajectory])
    plot_data$.trajectory <- rep(trajectory, times = feature_count)
  }

  plot_data
}

prepare_maaslin3_inputs <- function(data, fixed_effects, include_read_depth = FALSE) {
  validate_analysis_data(data, minimum_samples = 3L, analysis_name = "MaAsLin 3")
  if (length(fixed_effects) < 1L || any(!fixed_effects %in% names(data$metadata))) {
    stop("Choose one or more valid metadata variables as MaAsLin 3 fixed effects.", call. = FALSE)
  }
  invalid_names <- fixed_effects[make.names(fixed_effects) != fixed_effects]
  if (length(invalid_names)) {
    stop(
      "MaAsLin 3 model-variable names must not contain spaces or special characters. Rename: ",
      paste(invalid_names, collapse = ", "),
      call. = FALSE
    )
  }
  for (variable in fixed_effects) {
    effect <- data$metadata[[variable]]
    if (any(is.na(effect) | (is.character(effect) & !nzchar(effect)))) {
      stop("MaAsLin 3 fixed effect '", variable,
           "' contains missing values in filtered samples.", call. = FALSE)
    }
    if (length(unique(effect)) < 2L) {
      stop("MaAsLin 3 fixed effect '", variable,
           "' must have at least two distinct values after filtering.", call. = FALSE)
    }
  }
  if (nrow(data$abundance) < 2L) {
    stop("MaAsLin 3 requires at least two features after filtering.", call. = FALSE)
  }

  model_metadata <- data$metadata[, setdiff(names(data$metadata), data$sample_id_column), drop = FALSE]
  rownames(model_metadata) <- data$metadata[[data$sample_id_column]]
  model_effects <- fixed_effects
  if (isTRUE(include_read_depth)) {
    read_depth_name <- "sample_read_depth"
    if (read_depth_name %in% names(model_metadata)) {
      stop("Metadata already contains the reserved column 'sample_read_depth'. Rename it or disable the derived read-depth covariate.",
           call. = FALSE)
    }
    if (is.null(data$input_library_sizes)) {
      stop("Input library sizes are unavailable; the derived read-depth covariate cannot be used.",
           call. = FALSE)
    }
    sample_ids <- colnames(data$abundance)
    library_sizes <- data$input_library_sizes[match(sample_ids, names(data$input_library_sizes))]
    if (any(!is.finite(library_sizes))) {
      stop("Input library sizes could not be matched to all filtered samples.", call. = FALSE)
    }
    model_metadata[[read_depth_name]] <- unname(library_sizes)
    model_effects <- c(model_effects, read_depth_name)
  }

  list(
    input_data = as.data.frame(t(data$abundance), check.names = FALSE),
    input_metadata = model_metadata,
    fixed_effects = model_effects
  )
}

run_maaslin3 <- function(data, fixed_effects, output_directory, include_read_depth = FALSE) {
  inputs <- prepare_maaslin3_inputs(data, fixed_effects, include_read_depth)
  if (!requireNamespace("maaslin3", quietly = TRUE)) {
    stop(
      paste0(
        "Differential abundance requires the Bioconductor package 'maaslin3'. ",
        "Install it with BiocManager::install('maaslin3'), then restart the app."
      ),
      call. = FALSE
    )
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  on.exit(
    if ("maaslin_log_reset" %in% getNamespaceExports("maaslin3")) maaslin3::maaslin_log_reset(),
    add = TRUE
  )
  maaslin3::maaslin3(
    input_data = inputs$input_data,
    input_metadata = inputs$input_metadata,
    output = output_directory,
    fixed_effects = inputs$fixed_effects,
    min_abundance = 0,
    min_prevalence = 0,
    max_prevalence = 1.01,
    zero_threshold = 0,
    min_variance = 0,
    max_significance = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = TRUE,
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    warn_prevalence = TRUE,
    augment = TRUE,
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    verbosity = "WARN"
  )
  all_results_path <- file.path(output_directory, "all_results.tsv")
  significant_results_path <- file.path(output_directory, "significant_results.tsv")
  if (!file.exists(all_results_path) || !file.exists(significant_results_path)) {
    stop("MaAsLin 3 completed without producing its expected all_results.tsv and significant_results.tsv files.",
         call. = FALSE)
  }
  list(
    all_results = utils::read.delim(all_results_path, check.names = FALSE, stringsAsFactors = FALSE),
    significant_results = utils::read.delim(
      significant_results_path, check.names = FALSE, stringsAsFactors = FALSE
    )
  )
}
