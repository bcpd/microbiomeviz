# Microbiome Visualizer

Microbiome Visualizer is an R Shiny application for validating, filtering, exploring, and exporting feature-by-sample microbiome abundance data with sample metadata.

## Run the app

Required packages:

```r
install.packages(c("shiny", "DT", "ggplot2", "pheatmap", "vegan", "testthat"))
```

MaAsLin 3 is optional and is needed only for differential abundance:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("maaslin3")
```

Start from the repository root:

```r
shiny::runApp()
```

The bundled `sample_metadata.csv` and `taxonomy_data.csv` load by default. CSV files can also be uploaded in the app.

## Input contract

- Metadata: the first column contains unique, non-empty sample IDs.
- Taxonomy/abundance: the first column contains unique, non-empty feature IDs; remaining column names are sample IDs.
- Abundances must be numeric, finite, and non-negative.
- At least two sample IDs must occur in both files. Exact matching is used.
- Wholly blank metadata rows and unnamed, wholly blank metadata columns are removed and reported. Surrounding metadata whitespace is trimmed and reported.
- Samples present in only one file are excluded and listed. Matched metadata is reordered to the abundance-file sample order.

## Analysis behaviour

Filter defaults preserve all matched samples and features. Filters operate on raw values. Prevalence is the percentage of retained samples where abundance is greater than zero.

- Heatmap: shows the most abundant retained features. Raw abundance is the default; relative abundance is an explicit display-only option.
- Ordination: uses raw filtered abundances, Bray-Curtis dissimilarity (`vegan::vegdist`), and two-dimensional classical multidimensional scaling (`stats::cmdscale`). It runs only when requested.
- Timeseries: plots selected microbial features at every sampling time. Raw abundance is the default; relative abundance is an explicit display option. It does not average samples. Lines are drawn only after selecting a repeated-measures ID, so unrelated samples are never connected as a trajectory.
- Differential abundance: passes samples-by-features raw filtered data and aligned metadata to MaAsLin 3, which fits abundance and prevalence models. The app explicitly uses TSS normalization, log2 transformation, BH correction, metadata standardization, abundance median comparison, no prevalence median comparison, and a 0.1 q-value threshold. An opt-in option uses pre-filter input column sums as a read-depth covariate; enable it only for raw-count inputs. It runs only when requested in a session-specific temporary directory.
- Download: exports the same aligned filtered metadata and raw abundance data used by tables and analyses, plus filter and cleaning notes.

The bundled data currently produce visible notes: one blank metadata row, five blank unnamed columns, trailing whitespace in `Group`, and two abundance-only control samples. These are not silently analyzed.

## Tests

A small deterministic dataset under `tests/testthat/fixtures/` covers validation, matching, filtering, ordination, heatmap and timeseries preparation, Shiny reactive propagation, and download consistency.

```r
Rscript tests/testthat.R
```

MaAsLin 3 itself is not exercised in the default test suite because it is an optional Bioconductor dependency. Its input prerequisites and missing-package message are validated by the application layer.
