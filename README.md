# Microbiome Visualizer 

MicrobiomeVisualizer2 is an interactive R Shiny application designed to facilitate comprehensive microbiome data analysis. The app provides a user-friendly dashboard interface that allows researchers to explore, visualize, and analyze their microbiome data through a series of interactive tabs.

## Key Features

- **Data Summary:**  
  Quickly review an overview of your microbiome data, including key statistics and summaries.

- **Data Filtering:**  
  Apply filters to your data based on sample metadata or taxonomy, helping to narrow down the datasets for specific analyses.

- **Interactive Tables:**  
  View sample metadata and taxonomy information in interactive tables powered by the `DT` package.

- **Heatmap Visualization:**  
  Generate clustered heatmaps of your taxonomy data using the `pheatmap` package, which helps visualize patterns and groupings among taxa.

- **Ordination Analysis (PCoA):**  
  Perform Principal Coordinates Analysis (PCoA) using Bray-Curtis dissimilarity metrics (via the `vegan` package) and visualize the results with `ggplot2`.

- **Timeseries Visualization:**  
  (Placeholder) Designed to plot timeseries data, enabling the analysis of temporal trends in microbiome profiles.

- **Differential Abundance Analysis:**  
  Run differential abundance analysis using the Maaslin2 package. The app generates both a visualization (bar plot) of significant features and a detailed results table.

- **Download Project Data:**  
  Bundle and download project data for further analysis or sharing.


