required_packages <- c("shiny", "DT", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them before starting the app.",
    call. = FALSE
  )
}

find_app_dir <- function(start = getwd()) {
  candidate <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(candidate, "R", "data.R")) &&
        file.exists(file.path(candidate, "R", "analysis.R"))) {
      return(candidate)
    }
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop("Could not locate the application root containing R/data.R and R/analysis.R.", call. = FALSE)
}
app_dir <- find_app_dir()
source(file.path(app_dir, "R", "data.R"), local = TRUE)
source(file.path(app_dir, "R", "analysis.R"), local = TRUE)

app_ui <- shiny::fluidPage(
  shiny::titlePanel("Microbiome Data Analysis"),
  shiny::fluidRow(
    shiny::column(
      width = 3,
      shiny::wellPanel(
        shiny::h4("Data and filters"),
        shiny::fileInput("metadata_file", "Sample metadata (CSV)", accept = ".csv"),
        shiny::fileInput("taxonomy_file", "Taxonomy/abundance data (CSV)", accept = ".csv"),
        shiny::helpText("Leave both inputs empty to use the bundled example files. Sample IDs are matched exactly."),
        shiny::uiOutput("data_status"),
        shiny::hr(),
        shiny::selectInput("group_variable", "Filter samples by", choices = c("None" = "")),
        shiny::uiOutput("group_filter_ui"),
        shiny::numericInput("min_library_size", "Minimum sample library size", value = 0, min = 0),
        shiny::sliderInput(
          "min_prevalence", "Minimum feature prevalence (%)",
          min = 0, max = 100, value = 0, step = 1
        ),
        shiny::numericInput("min_total_abundance", "Minimum feature total abundance", value = 0, min = 0),
        shiny::helpText("Prevalence means the percentage of retained samples with abundance > 0. Defaults retain all samples and features."),
        shiny::uiOutput("filter_status")
      )
    ),
    shiny::column(
      width = 9,
      shiny::tabsetPanel(
        id = "main_tab",
        shiny::tabPanel(
          "Data summary", value = "summary",
          shiny::h3("Data summary"),
          shiny::uiOutput("summary_content")
        ),
        shiny::tabPanel(
          "Sample metadata", value = "metadata",
          shiny::h3("Filtered sample metadata"),
          DT::DTOutput("metadata_table")
        ),
        shiny::tabPanel(
          "Taxonomy", value = "taxonomy",
          shiny::h3("Filtered taxonomy/abundance table"),
          DT::DTOutput("taxonomy_table")
        ),
        shiny::tabPanel(
          "Heatmap", value = "heatmap",
          shiny::h3("Abundance heatmap"),
          shiny::fluidRow(
            shiny::column(4, shiny::numericInput("heatmap_top", "Most abundant features shown", 50, min = 1, max = 500)),
            shiny::column(4, shiny::selectInput(
              "heatmap_transform", "Display transformation",
              choices = c("Raw abundance" = "raw", "Relative abundance per sample" = "relative")
            )),
            shiny::column(4, shiny::checkboxInput(
              "heatmap_all_features", "Show all filtered features (may be slow)", value = FALSE
            ))
          ),
          shiny::helpText("Feature ranking and the optional transformation affect this display only. No scaling is applied by pheatmap."),
          shiny::plotOutput("heatmap_plot", height = "700px")
        ),
        shiny::tabPanel(
          "Ordination", value = "ordination",
          shiny::h3("Bray-Curtis PCoA"),
          shiny::selectInput("ordination_colour", "Colour points by", choices = c("None" = "")),
          shiny::actionButton("run_ordination", "Run PCoA", class = "btn-primary"),
          shiny::helpText("Uses raw filtered abundances, Bray-Curtis dissimilarity, and classical multidimensional scaling (cmdscale)."),
          shiny::plotOutput("ordination_plot", height = "550px"),
          shiny::verbatimTextOutput("ordination_message")
        ),
        shiny::tabPanel(
          "Timeseries", value = "timeseries",
          shiny::h3("Microbial features over time"),
          shiny::fluidRow(
            shiny::column(4, shiny::selectInput("time_variable", "Date/time metadata column", choices = c("None" = ""))),
            shiny::column(4, shiny::selectizeInput(
              "timeseries_features", "Features to display", choices = character(), multiple = TRUE,
              options = list(placeholder = "Choose one or more features")
            )),
            shiny::column(4, shiny::selectInput(
              "timeseries_transform", "Abundance scale",
              choices = c("Raw abundance" = "raw", "Relative abundance per sample" = "relative")
            ))
          ),
          shiny::fluidRow(
            shiny::column(4, shiny::selectInput("time_colour", "Colour points by", choices = c("None" = ""))),
            shiny::column(4, shiny::selectInput(
              "time_trajectory", "Repeated-measures ID for lines (optional)", choices = c("None" = "")
            ))
          ),
          shiny::helpText(
            "Each point is one retained sample's selected feature abundance. Lines are drawn only when a repeated-measures ID is selected; no samples are averaged."
          ),
          shiny::plotOutput("timeseries_plot", height = "500px")
        ),
        shiny::tabPanel(
          "Differential abundance", value = "differential",
          shiny::h3("Differential abundance with MaAsLin 3"),
          shiny::selectizeInput(
            "fixed_effect", "Fixed effect(s)", choices = character(),
            multiple = TRUE, options = list(placeholder = "Choose one or more metadata variables")
          ),
          shiny::checkboxInput(
            "include_read_depth",
            "Use input column sum as read-depth covariate (raw counts only)",
            value = FALSE
          ),
          shiny::actionButton("run_differential", "Run MaAsLin 3", class = "btn-primary"),
          shiny::helpText(
            paste(
              "The app passes the raw filtered feature table to MaAsLin 3, which models both abundance and prevalence.",
              "Settings are explicit: TSS normalization, log2 transformation, BH correction, metadata standardization,",
              "abundance median comparison on, prevalence median comparison off, and q-value threshold 0.1."
            )
          ),
          shiny::verbatimTextOutput("differential_message"),
          shiny::plotOutput("differential_plot", height = "550px"),
          DT::DTOutput("maaslin3_results")
        ),
        shiny::tabPanel(
          "Download project data", value = "download",
          shiny::h3("Download filtered project data"),
          shiny::p("The ZIP contains aligned filtered metadata, raw abundance data, active filters, and data-cleaning notes."),
          shiny::downloadButton("download_data", "Download ZIP")
        )
      )
    )
  )
)

app_server <- function(input, output, session,
                       default_metadata = file.path(app_dir, "sample_metadata.csv"),
                       default_taxonomy = file.path(app_dir, "taxonomy_data.csv")) {
  data_result <- shiny::reactive({
    metadata_path <- if (is.null(input$metadata_file)) default_metadata else input$metadata_file$datapath
    taxonomy_path <- if (is.null(input$taxonomy_file)) default_taxonomy else input$taxonomy_file$datapath
    tryCatch(
      list(ok = TRUE, data = load_microbiome_data(metadata_path, taxonomy_path), message = NULL),
      error = function(e) list(ok = FALSE, data = NULL, message = conditionMessage(e))
    )
  })

  current_data <- shiny::reactive({
    result <- data_result()
    shiny::validate(shiny::need(result$ok, result$message))
    result$data
  })

  output$data_status <- shiny::renderUI({
    result <- data_result()
    if (!result$ok) {
      return(shiny::div(class = "alert alert-danger", shiny::strong("Input error: "), result$message))
    }
    data <- result$data
    shiny::tagList(
      shiny::div(
        class = "alert alert-success",
        sprintf("Loaded %d matched samples and %s features.",
                ncol(data$abundance), format(nrow(data$abundance), big.mark = ","))
      ),
      if (length(data$notes)) {
        shiny::div(
          class = "alert alert-warning",
          shiny::strong("Input notes"),
          shiny::tags$ul(lapply(data$notes, shiny::tags$li))
        )
      }
    )
  })

  shiny::observeEvent(current_data(), {
    data <- current_data()
    metadata_variables <- setdiff(names(data$metadata), data$sample_id_column)
    selected <- if ("Group" %in% metadata_variables) "Group" else ""
    shiny::updateSelectInput(
      session, "group_variable",
      choices = c("None" = "", metadata_variables), selected = selected
    )
    shiny::updateSelectInput(
      session, "ordination_colour",
      choices = c("None" = "", metadata_variables), selected = selected
    )
    shiny::updateSelectInput(
      session, "time_colour",
      choices = c("None" = "", metadata_variables), selected = selected
    )
    shiny::updateSelectInput(
      session, "time_trajectory",
      choices = c("None" = "", metadata_variables), selected = ""
    )
    date_guess <- metadata_variables[grepl("date|time", metadata_variables, ignore.case = TRUE)]
    shiny::updateSelectInput(
      session, "time_variable",
      choices = c("None" = "", metadata_variables),
      selected = if (length(date_guess)) date_guess[1L] else ""
    )
    top_features <- rownames(data$abundance)[order(rowSums(data$abundance), decreasing = TRUE)]
    shiny::updateSelectizeInput(
      session, "timeseries_features",
      choices = rownames(data$abundance), selected = utils::head(top_features, 5L), server = TRUE
    )
    shiny::updateSelectInput(
      session, "fixed_effect",
      choices = metadata_variables, selected = selected
    )
  }, ignoreInit = FALSE)

  output$group_filter_ui <- shiny::renderUI({
    data <- current_data()
    group_variable <- input$group_variable %||% ""
    if (!nzchar(group_variable) || !group_variable %in% names(data$metadata)) {
      return(NULL)
    }
    values <- as.character(data$metadata[[group_variable]])
    values[is.na(values) | !nzchar(values)] <- "__MISSING__"
    choices <- unique(values)
    labels <- ifelse(choices == "__MISSING__", "(missing)", choices)
    shiny::checkboxGroupInput(
      "selected_groups", "Retain values",
      choices = stats::setNames(choices, labels), selected = choices
    )
  })

  filtered_data <- shiny::reactive({
    data <- current_data()
    filter_microbiome_data(
      data,
      group_variable = input$group_variable %||% "",
      selected_groups = input$selected_groups,
      min_library_size = input$min_library_size %||% 0,
      min_prevalence_percent = input$min_prevalence %||% 0,
      min_total_abundance = input$min_total_abundance %||% 0
    )
  })

  analysis_key <- shiny::reactive({
    metadata_source <- if (is.null(input$metadata_file)) default_metadata else input$metadata_file$datapath
    taxonomy_source <- if (is.null(input$taxonomy_file)) default_taxonomy else input$taxonomy_file$datapath
    paste(
      metadata_source,
      taxonomy_source,
      input$group_variable %||% "",
      paste(sort(input$selected_groups %||% character()), collapse = ","),
      input$min_library_size %||% 0,
      input$min_prevalence %||% 0,
      input$min_total_abundance %||% 0,
      isTRUE(input$include_read_depth),
      sep = "|"
    )
  })

  output$filter_status <- shiny::renderUI({
    data <- filtered_data()
    summary <- data$filter_summary
    class <- if (summary$samples_after > 0 && summary$features_after > 0) "alert alert-info" else "alert alert-danger"
    shiny::div(
      class = class,
      sprintf("Retained %d/%d samples and %s/%s features.",
              summary$samples_after, summary$samples_before,
              format(summary$features_after, big.mark = ","),
              format(summary$features_before, big.mark = ","))
    )
  })

  output$summary_content <- shiny::renderUI({
    data <- filtered_data()
    shiny::tagList(
      shiny::tableOutput("summary_table"),
      shiny::h4("Scientific and data-handling assumptions"),
      shiny::tags$ul(
        shiny::tags$li("The first metadata column is the sample ID; the first abundance column is the feature ID."),
        shiny::tags$li("Only exact sample-ID matches are analyzed; both tables are reordered to the abundance-file sample order."),
        shiny::tags$li("Filters operate on raw, non-negative abundances. A feature is present when abundance > 0."),
        shiny::tags$li("Heatmap and timeseries display transformations do not affect ordination, tables, differential abundance, or downloads."),
        shiny::tags$li("Bray-Curtis PCoA and MaAsLin 3 run only when requested and use the current filtered dataset."),
        shiny::tags$li("MaAsLin 3 models abundance and prevalence. The optional read-depth covariate uses pre-filter input column sums and is valid only when abundances are raw counts.")
      ),
      if (length(data$notes)) shiny::tagList(shiny::h4("Input notes"), shiny::tags$ul(lapply(data$notes, shiny::tags$li)))
    )
  })

  output$summary_table <- shiny::renderTable({
    data <- filtered_data()
    library_sizes <- colSums(data$abundance)
    data.frame(
      Metric = c("Samples", "Features", "Minimum library size", "Median library size", "Maximum library size", "Total abundance"),
      Value = c(
        ncol(data$abundance), nrow(data$abundance),
        if (length(library_sizes)) min(library_sizes) else NA,
        if (length(library_sizes)) stats::median(library_sizes) else NA,
        if (length(library_sizes)) max(library_sizes) else NA,
        sum(data$abundance)
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  output$metadata_table <- DT::renderDT({
    DT::datatable(filtered_data()$metadata, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)

  output$taxonomy_table <- DT::renderDT({
    data <- filtered_data()
    taxonomy <- data.frame(feature = rownames(data$abundance), data$abundance, check.names = FALSE)
    names(taxonomy)[1L] <- data$feature_column
    DT::datatable(taxonomy, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)

  output$heatmap_plot <- shiny::renderPlot({
    data <- filtered_data()
    shiny::validate(shiny::need(nrow(data$abundance) > 0, "No features remain. Relax the feature filters."))
    shiny::validate(shiny::need(ncol(data$abundance) > 0, "No samples remain. Relax the sample filters."))
    shiny::validate(shiny::need(
      requireNamespace("pheatmap", quietly = TRUE),
      "Heatmaps require the 'pheatmap' package. Install it with install.packages('pheatmap')."
    ))
    matrix <- prepare_heatmap_matrix(
      data,
      top_n = if (isTRUE(input$heatmap_all_features)) nrow(data$abundance) else input$heatmap_top %||% 50,
      transform = input$heatmap_transform %||% "raw"
    )
    result <- pheatmap::pheatmap(
      matrix,
      cluster_rows = nrow(matrix) > 1L,
      cluster_cols = ncol(matrix) > 1L,
      show_rownames = nrow(matrix) <= 75L,
      silent = TRUE,
      main = sprintf("Top %d features (%s)", nrow(matrix),
                     if ((input$heatmap_transform %||% "raw") == "raw") "raw abundance" else "relative abundance")
    )
    grid::grid.newpage()
    grid::grid.draw(result$gtable)
  })

  ordination_result <- shiny::eventReactive(input$run_ordination, {
    shiny::withProgress(message = "Calculating Bray-Curtis PCoA", value = 0.5, {
      list(value = calculate_ordination(filtered_data()), key = analysis_key())
    })
  }, ignoreInit = TRUE)

  output$ordination_plot <- shiny::renderPlot({
    shiny::validate(shiny::need(input$run_ordination > 0, "Select filters, then click 'Run PCoA'."))
    result <- ordination_result()
    shiny::validate(shiny::need(
      identical(result$key, analysis_key()),
      "Filters or input data changed. Click 'Run PCoA' to calculate current results."
    ))
    points <- result$value
    explained <- attr(points, "explained")
    colour_variable <- input$ordination_colour %||% ""
    if (nzchar(colour_variable) && colour_variable %in% names(points)) {
      points$.colour <- as.factor(points[[colour_variable]])
      plot <- ggplot2::ggplot(points, ggplot2::aes(x = Axis1, y = Axis2, colour = .colour)) +
        ggplot2::labs(colour = colour_variable)
    } else {
      plot <- ggplot2::ggplot(points, ggplot2::aes(x = Axis1, y = Axis2))
    }
    plot +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = "PCoA of Bray-Curtis dissimilarity",
        x = sprintf("PCoA1 (%.1f%% of positive eigenvalues)", explained[1L]),
        y = sprintf("PCoA2 (%.1f%% of positive eigenvalues)", explained[2L])
      )
  })

  output$ordination_message <- shiny::renderText({
    if (input$run_ordination == 0) {
      return("PCoA has not been run. Results are recalculated only when you click the button.")
    }
    result <- ordination_result()
    if (!identical(result$key, analysis_key())) {
      return("Filters or input data changed. Click 'Run PCoA' to calculate current results.")
    }
    points <- result$value
    eigenvalues <- attr(points, "eigenvalues")
    negative <- sum(eigenvalues < -sqrt(.Machine$double.eps))
    sprintf("PCoA used %d samples and %d features. Negative eigenvalues: %d.",
            nrow(points), nrow(filtered_data()$abundance), negative)
  })

  output$timeseries_plot <- shiny::renderPlot({
    data <- filtered_data()
    time_variable <- input$time_variable %||% ""
    features <- input$timeseries_features %||% character()
    shiny::validate(shiny::need(nzchar(time_variable), "Choose a date/time metadata column."))
    shiny::validate(shiny::need(length(features), "Choose one or more features to display."))
    colour_variable <- input$time_colour %||% ""
    trajectory_variable <- input$time_trajectory %||% ""
    plot_data <- tryCatch(
      prepare_timeseries_data(
        data,
        time_variable = time_variable,
        features = features,
        transform = input$timeseries_transform %||% "raw",
        colour_variable = colour_variable,
        trajectory_variable = trajectory_variable
      ),
      error = function(e) shiny::validate(shiny::need(FALSE, conditionMessage(e)))
    )

    if (nzchar(colour_variable)) {
      plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Abundance, colour = .colour)) +
        ggplot2::labs(colour = colour_variable)
    } else {
      plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Abundance))
    }

    if (nzchar(trajectory_variable)) {
      plot <- plot + ggplot2::geom_line(
        ggplot2::aes(group = interaction(Feature, .trajectory)), alpha = 0.55
      )
    }

    plot +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::facet_wrap(~Feature, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = "Selected microbial features over time",
        x = time_variable,
        y = if ((input$timeseries_transform %||% "raw") == "raw") "Raw abundance" else "Relative abundance"
      )
  })

  differential_result <- shiny::eventReactive(input$run_differential, {
    fixed_effect <- input$fixed_effect %||% character()
    output_directory <- file.path(tempdir(), paste0("microbiomeviz-maaslin3-", session$token))
    if (dir.exists(output_directory)) unlink(output_directory, recursive = TRUE)
    shiny::withProgress(message = "Running MaAsLin 3", value = 0.5, {
      list(
        value = run_maaslin3(
          filtered_data(), fixed_effect, output_directory,
          include_read_depth = isTRUE(input$include_read_depth)
        ),
        key = analysis_key()
      )
    })
  }, ignoreInit = TRUE)

  session$onSessionEnded(function() {
    output_directory <- file.path(tempdir(), paste0("microbiomeviz-maaslin3-", session$token))
    if (dir.exists(output_directory)) unlink(output_directory, recursive = TRUE)
  })

  output$differential_message <- shiny::renderText({
    if (input$run_differential == 0) {
      if (!requireNamespace("maaslin3", quietly = TRUE)) {
        return("MaAsLin 3 is not installed. Install it with BiocManager::install('maaslin3') to enable this workflow.")
      }
      return("Choose one or more fixed effects and click 'Run MaAsLin 3'.")
    }
    result <- differential_result()
    if (!identical(result$key, analysis_key())) {
      return("Filters or input data changed. Click 'Run MaAsLin 3' to calculate current results.")
    }
    sprintf("MaAsLin 3 returned %d abundance/prevalence coefficients; %d are in significant_results.tsv.",
            nrow(result$value$all_results), nrow(result$value$significant_results))
  })

  output$differential_plot <- shiny::renderPlot({
    shiny::validate(shiny::need(input$run_differential > 0, "Run MaAsLin 3 to view results."))
    result <- differential_result()
    shiny::validate(shiny::need(
      identical(result$key, analysis_key()),
      "Filters or input data changed. Run MaAsLin 3 again to view current results."
    ))
    results <- result$value$significant_results
    required <- c("feature", "name", "coef", "qval_individual", "qval_joint", "model")
    shiny::validate(shiny::need(all(required %in% names(results)),
      "MaAsLin 3 results do not contain the expected association and q-value columns."))
    significant <- results
    significant$.min_q <- pmin(significant$qval_individual, significant$qval_joint, na.rm = TRUE)
    significant$.min_q[!is.finite(significant$.min_q)] <- NA_real_
    shiny::validate(shiny::need(nrow(significant) > 0, "No associations meet the MaAsLin 3 q-value threshold of 0.1."))
    significant <- significant[order(significant$.min_q, -abs(significant$coef)), , drop = FALSE]
    significant <- utils::head(significant, 30L)
    significant$association <- paste(significant$feature, significant$name, sep = " — ")
    significant$association <- factor(significant$association, levels = rev(unique(significant$association)))
    ggplot2::ggplot(significant, ggplot2::aes(x = association, y = coef, fill = model)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~model, scales = "free") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = "MaAsLin 3 significant associations (q-value <= 0.1)",
        x = "Feature — model term", y = "Coefficient (model-specific scale)"
      )
  })

  output$maaslin3_results <- DT::renderDT({
    shiny::validate(shiny::need(input$run_differential > 0, "Run MaAsLin 3 to view the result table."))
    result <- differential_result()
    shiny::validate(shiny::need(
      identical(result$key, analysis_key()),
      "Filters or input data changed. Run MaAsLin 3 again to view current results."
    ))
    DT::datatable(result$value$significant_results, rownames = FALSE, filter = "top",
                  options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)

  output$download_data <- shiny::downloadHandler(
    filename = function() sprintf("project_data-%s.zip", Sys.Date()),
    content = function(file) {
      data <- filtered_data()
      if (ncol(data$abundance) == 0L || nrow(data$abundance) == 0L) {
        stop("No filtered data are available to download. Relax the filters.", call. = FALSE)
      }
      write_project_bundle(data, file)
    },
    contentType = "application/zip"
  )

  server_state <- list(
    data_result = data_result,
    current_data = current_data,
    filtered_data = filtered_data,
    analysis_key = analysis_key,
    ordination_result = ordination_result,
    differential_result = differential_result
  )
  session$userData$microbiomeviz <- server_state
  server_state
}

app <- shiny::shinyApp(ui = app_ui, server = app_server)
app
