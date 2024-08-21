library(shiny)
library(shinydashboard)
library(DT)
library(pheatmap)
library(vegan)
library(ggplot2)

# Define the UI
ui <- dashboardPage(
  dashboardHeader(title = "Microbiome Data Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data summary", tabName = "summary", icon = icon("info-circle")),
      menuItem("Filter data", tabName = "filter", icon = icon("filter")),
      menuSubItem("Sample metadata", tabName = "metadata", icon = icon("table")),
      menuSubItem("Taxonomy", tabName = "taxonomy", icon = icon("table")),
      menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
      menuItem("Ordination", tabName = "ordination", icon = icon("chart-pie")),
      menuItem("Timeseries", tabName = "timeseries", icon = icon("line-chart")),
      menuItem("Differential Abundance", tabName = "differential", icon = icon("chart-bar")),
      menuItem("Download project data", tabName = "download", icon = icon("download"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "summary",
              h2("Data Summary"),
              # Add UI elements to summarize your data
      ),
      tabItem(tabName = "filter",
              h2("Filter Data"),
              # Add UI elements for data filtering
      ),
      tabItem(tabName = "metadata",
              h2("Sample Metadata"),
              DTOutput("metadata_table")
      ),
      tabItem(tabName = "taxonomy",
              h2("Taxonomy"),
              DTOutput("taxonomy_table")
      ),
      tabItem(tabName = "heatmap",
              h2("Heatmap"),
              plotOutput("heatmap_plot")
      ),
      tabItem(tabName = "ordination",
              h2("Ordination"),
              plotOutput("ordination_plot")
      ),
      tabItem(tabName = "timeseries",
              h2("Timeseries"),
              plotOutput("timeseries_plot")
      ),
      tabItem(tabName = "differential",
              h2("Differential Abundance"),
              plotOutput("differential_plot"),
              DTOutput("maaslin2_results")
      ),
      tabItem(tabName = "download",
              h2("Download Project Data"),
              downloadButton("download_data", "Download")
      )
    )
  )
)


####
# Define server logic
server <- function(input, output) {
  
  # Load and process data
  sample_data <- reactive({
    # Load your sample metadata here
    read.csv("sample_metadata.csv")
  })
  
  taxonomy_data <- reactive({
    # Load your taxonomy data here
    read.csv("taxonomy_data.csv", row.names = 1) # Assume the first column has row names (taxa/features)
  })
  
  output$metadata_table <- renderDT({
    sample_data()
  })
  
  output$taxonomy_table <- renderDT({
    taxonomy_data()
  })
  
  # Heatmap generation using pheatmap
  output$heatmap_plot <- renderPlot({
    pheatmap(taxonomy_data(), cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
  })
  
  # Ordination plot (example using PCoA from the vegan package)
  output$ordination_plot <- renderPlot({
    dist_matrix <- vegdist(t(taxonomy_data()), method = "bray") # Bray-Curtis dissimilarity
    pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
    pcoa_data <- data.frame(Sample = rownames(pcoa_result$points),
                            Axis1 = pcoa_result$points[, 1],
                            Axis2 = pcoa_result$points[, 2])
    ggplot(pcoa_data, aes(x = Axis1, y = Axis2)) +
      geom_point() +
      theme_minimal() +
      labs(title = "PCoA of Microbiome Data", x = "PCoA1", y = "PCoA2")
  })
  
  # Differential abundance using Maaslin2
  output$differential_plot <- renderPlot({
    library(Maaslin2)
    maaslin2_results <- Maaslin2(
      input_data = taxonomy_data(),
      input_metadata = sample_data(),
      output = "maaslin2_output",
      fixed_effects = colnames(sample_data())[2:ncol(sample_data())] # Assuming first column is Sample ID
    )
    
    # Example of visualizing one of the significant results
    significant_results <- read.table("maaslin2_output/significant_results.tsv", header = TRUE, sep = "\t")
    ggplot(significant_results, aes(x = feature, y = coef)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      coord_flip() +
      labs(title = "Differential Abundance Analysis", x = "Feature", y = "Effect Size")
  })
  
  output$maaslin2_results <- renderDT({
    read.table("maaslin2_output/significant_results.tsv", header = TRUE, sep = "\t")
  })
  
  # Download data
  output$download_data <- downloadHandler(
    filename = function() {
      paste("project_data-", Sys.Date(), ".zip", sep="")
    },
    content = function(file) {
      # Bundle and zip project data for download
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
