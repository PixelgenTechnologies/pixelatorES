# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `get_hash_stats()` to summarise sample hashing from per-cell `hash_counts` metadata (purity, heatmap inputs, sample-level stats).
- `component_hashing()` to render hash purity violin, summary table, and hash purity / hash fraction heatmaps.
- `sample_calling` as a recognised `pixelator` pipeline stage in `find_stage()` / `get_file_paths()` file resolution.
- Quarto parameters on `inst/quarto/pixelatorES.qmd`: `include_pool_hash_metrics`, `include_sample_hashing_section`, and `key_metrics_detailed` to control key-metrics and hashing sections in the bundled report.
- `inst/quarto/quality_metrics.qmd` now supports hashed experiments: optional pool vs sample key-metric layouts, conditional sample hashing subsection, pool-aware sequencing and QC components, expanded graph metrics (including crossing edges and k-core heatmaps), bleedover noise, and control markers with both violin plots where applicable.
- `pipeline_pool_stages` containing the stages of the nf-core/pixelator pipeline that produce pool-level QC files.
- Added a Samples tab `inst/quarto/samples.qmd`.
- Added render utilities `section_table` and `section_intro`.

### Changed

- `read_qc_files()` returns a nested list `list(qc_files = <named list per sample>, pool_qc_files = NULL | <named list per pool>)`. Passing a data frame of QC paths (legacy) is still supported; the return value is wrapped in the same structure.
- `key_metric_table()` returns `list(sample = <datatables>, pool = <datatables or NULL>)` for HTML output, and `return_data = TRUE` returns the corresponding pre-styled wide tibbles in the same shape. Hashing columns are included in definitions when `sample_hash_stats` is present.
- `get_file_paths()` accepts either `sample_aliases` (previous behaviour) or `sample_sheet` (maps files via `sample`, detects pool-level JSON/PXL, returns `pool_qc_files`). The return value always includes `pool_qc_files` (set to `NULL` when using `sample_aliases` only).
- `read_samplesheet()` adds a `pool` column (all `NA`) when the CSV has no pool column, so downstream code can rely on a stable schema.
- `merge_data()` joins optional `pool` from the sample sheet into object metadata when that column exists.
- `get_qc_metrics()` adds optional `sample_hash_stats`, normalises nested QC input, and orders `sample_alias` / `pool` factors using levels that actually appear in each derived table (avoids errors when the Seurat object does not contain every sample in the sheet).
- `get_seq_saturation()` uses an inner join to object-derived graph totals so QC rows for samples absent from the object do not produce undefined saturations.
- `get_crossing_edges()` and `get_denoising_data()` accept nested QC lists from `read_qc_files()`; legacy flat per-sample lists are still accepted via normalisation.
- `extract_sample_qc_metrics()` automatically uses `qc_input$qc_files` when given the full return value of `read_qc_files()`.
- `component_sequencing_reads_and_molecules()`, `component_cell_recovery()`, and `component_qc_molecule_rank_plot()` select pool-level vs sample-level QC JSON lists when `pool_qc_files` is present.
- `component_crossing_edges()` and `component_sequencing_saturation()` rename the first grouping column to `sample_alias` for plotting so pool-level metrics render on the same code paths.
- `inst/quarto/preprocessing.qmd` uses `get_file_paths(data_folder, sample_sheet = sample_sheet)`, extends `sample_aliases` with pool IDs when a `pool` column exists, and calls `read_qc_files(file_paths, sample_sheet)`.
- `inst/quarto/quality_metrics.qmd` still omits dispersion and “percent nodes in largest component” subsections; the corresponding helpers are not part of this package yet.
- `test_samplesheet()`, `test_data_folder()`, and `get_test_qc_metrics()`, now return two types of content, controlled by `type = c("default", "hashing")`.

### Fixed

- Removed a duplicate `read_qc_files()` invocation in `get_test_qc_metrics()`.
- `downsample_data` will now give an informative error if `control_markers` are missing from data.
- Bug in `order_cd_markers` which would return duplicate markers if a control marker string started with "CD".

## [0.8.6] 2026-03-23

### Updated

- `pixelatorES` now depends on `pixelatorR` >= v0.17.0.

## [0.8.5] 2026-03-16

### Updated

- Added a `ci` quarto profile that is used for CI smoketests, disabling chunk error tolerance.
- Added global option `pixelatorES.dev_mode = isTRUE(params$test_mode) || isTRUE(params$debug_mode)` to control whether the ES should be rendered in "dev mode", which is a mode where plot build errors in `title_plotlist()` are not tolerated and will cause the render to fail (overall chunk error tolerance is still controlled by Quarto's `execute: error:` setting).

## [0.8.4] 2026-03-09

### Fixed

- Fixed another issue in per-cell type abundance violin plots where the plots would fail to render due to that the `condition` column is missing.

## [0.8.3] 2026-03-06

### Fixed

- Fixed an issue where per-cell type abundance violin plots would fail to render with more than 20 samples.

## [0.8.2] 2026-03-05

### Fixed

- Fixed an issue where cells annotated as "DC" would fail to get aggregated to "Mono & DC" in the cell type aggregation step.

## [0.8.1] 2026-03-05

### Fixed

- Fixed an issue where spatial score clustering violin plots would fail to render with more than 20 samples.

## [0.8.0] 2026-03-05

### Added

- Added Q30 Phred score to key metrics table. 

### Fixed

- Removed unused imports from NAMESPACE

### Updated

- Cell types are now aggregated to "T", "B", "Mono & DC", "NK", and "Platelets".
- Cell annotation reference is now pulled from the `pixelatorR` package instead of being stored in the ES.
- The `annotate_cells` function has been removed from the `pixelatorES` package and is now contained in the `pixelatorR` package.


## [0.7.0] 2026-01-19

### Added
- Added `print_pixelator_version` to print a table with the Pixelator version and panel file for each sample.

### Updated
- The key metrics table is now in long format to avoid horizontal scrolling with high number of samples.
- The `Run settings` tab is now its own module.

### Fixed
- The ES now tolerates chunk errors and will in most cases render a report even if some chunks fail.
- Fixed issue where the ES would fail if `pg_data` contains fewer cells than the number of PCs to compute.
- Fixed DuckDB issue with temporary directory
- Fixed issue with violin median lines displaying per x coordinate rather than per group when also using the `fill` argument. 
- Fixed issue where DT tables would not render in full until switching back and fourth between tabs.

## [0.6.0] 2025-11-28

### Updated

- `order_sample_alias_factors` can now be used to order factors of any column. 
- The protein selection for colocalization heatmaps is now done by abundance instead of colocalization variance. 
- Switched to NMF algorithm for cell annotation. 
- The cell annotation section has been reorganized. 

### Fixed
- Fixed an issue where the plot legend of dimensionality reduction plots would overlap with the plot x-axis title. 
- Fixed an issue where the screen would "jump" to the top of the page when changing tabs.
- Removed `draw_quantiles` in favor of `ggplot2` functionality doing the same thing but less error prone. 
- Fixed a bug in `component_proximity_per_marker` that would cause samples order alphabetically. 

## [0.5.0] 2025-09-26

### Updated

- Proximity Score filtering is now harsher, requiring markers to be detected higher than the 90th percentile of `isotype_fraction` and 0.25% of counts. 
- The protein selection for colocalization heatmaps are now done per each cell type. 
- The colocalization heatmap color legend range has been decreased to range -0.8 to 0.8 to better visualize the differences in colocalization.

## [0.4.5] 2025-09-19

### Updated

- The ES is now more lenient with the file input, allowing some QC data to be missing, and allowing `sample_alias` and `condition` columns to be missing or empty in the sample sheet. 

## [0.4.4] 2025-09-05

### Added

- Added Samplewise Dimensionality Reduction plots and restyled the Cell annotation tab somewhat. 

### Fixed

- Fixed an issue in ANOVA components that would throw an error when only a single sample was run, as ANOVA requires at least two groups to compare. These functions will now default to using only the `seurat_clusters` variable for ANOVA if only a single sample is detected in the data. 

## [0.4.3] 2025-08-26

### Updated

- Renamed report to "Proxiome Experiment Summary". 

## [0.4.2] 2025-08-21

### Fixed

- Fixed an issue in colocalization heatmaps that would throw an error when running tests. 

## [0.4.1] 2025-08-21

### Fixed

- Fixed a bug in `plot_violin` that would throw an error if there is only a single point to be plotted.
- Removed equations from key metrics table tooltips, as these could not be rendered properly.

### Updated

- Images are now converted to WebP to limit file size somewhat.

### Added

- Workflow to label images with the version upon release.
- Push containers to quay.io

## [0.4.0] 2025-08-18

### Updates

- More extensive document styling has been added.
- Default color palettes for cell types and samples have been modified.
- `pixelatorES.qmd` now calls upon modular qmd files for each tab.
- Slimmed down content in "Quality metrics".
- `key_metric_table` now has an argument `detailed` to control whether all key metrics should be shown.
- Total and sample-wise proximity score violin plots are now plotted together in a tabset to avoid "nested" structures in the ES.
- The Proximity tab is now called "Spatial metrics" and "Self-proximity" is now called "Clustering".
- The proximity score content has had a major overhaul to make calculations within the functions that generate the plots, rather than in previous steps. This allows for more flexibility in the plots and avoids unexpected dependencies of different functions.
- Colocalization heatmaps are now showing the 40 markers with the highest F statistic in proximity ANOVA.

### Added

- `component_sequencing_saturation_curve` and sequencing saturation curves.
- Pop up tooltips for all key metrics.
- Horizontal line (`hline`) at 0 for poximity score violins.
- Functions for performing ANOVA on abundance and proximity metrics.
- `complete_proximity_scores` function to add 0's to proximity scores for markers that are missing for each cell.

## [0.3.2] 2025-08-11

### Added

- `cluster_palette` for coloring clusters.
- `plot_void` function to plot empty `ggplot2` plots.
- `extract_legend` function to extract legends from `ggplot2` plots.

### Updates

- Lowered resolution to 96 dpi.
- `tabset_nested_plotlist` and `tabset_plotlist` now have an argument `close` to control whether the tabset should be closed or not.

### Fixed

- `draw_quantiles` now works when there are missing values or few data points.
- Overflowing legends should happen less often in Dimensionality Reduction plots.
- Bug causing tabs in the "Selected contrasts" section to lack names.
- `annotation_dimred_heatmap` is now given more space and is tabset together with "Automated annotation".

## [0.3.1] 2025-08-07

### Added

- `displayed_cell_types` vector to hold cell type levels to plot, when plotted cell types should be restricted to a smaller selection.

### Updates

- `component_proximity_selected` and `component_abundance_per_celltype` now only plot the displayed cell types.

### Fixed

- `draw_quantiles` now works with `facet_var` in `plot_violin`.
- `read_samplesheet` will now work when a sample consists of concatenated fastq files.
- `draw_quantiles` now works with log scale.
- `draw_quantiles` now plots the quantile at the quantile of the data rather than of the density function (note: which is default of ggplot2).
- `pixelatorES.qmd` file is now runnable with the output of `nf-core/pixelator`.

## [0.3.0] 2025-08-06

### Added

- `pixelatorES.qmd` file added.

### Updated

- Removed the `metadata` variable.

## [0.2.0] 2025-08-04

### Added

- nav-pills to style tabsets
- `get_test_qc_metrics` to get qc metrics quickly for testing.
- `component_qc_molecule_rank_plot` to make the molecule rank plot from qc data.
- `get_test_data` to make a small Seurat object for testing components.
- `Graph Nodes [M]` and `Graph Edges [M]` to the key metrics table.

### Updated

- `process_data` now uses the public `AnnotateCells` function from `pixelatorR`. `process_data` can now switch annotation method to "nmf" by specifying `params$annotation_method = "nmf"`.
- `component_cell_recovery` now returns a tabsetted Molecule Rank Plot.
- `violin_plot` is now borderless by default.
- `component_crossing_edges` now returns a plot without absolute numbers and with fewer decimals to avoid overplotting.
- `component_annotation` now outputs an additional `celltype_numbers_table` with a table of cell type numbers per sample.

### Fixed

- Disabled the `draw_quantile` option when a `facet_var` is provided to `plot_violin`. Otherwise, the medians are placed incorrectly.
- `title_plotlist` now returns an empty plot if the render fails.
- Overlapping text in `component_cell_recovery` plots.
- All violin plots are now using `scale = "width"`.

## [0.1.0] 2025-07-22

### Added

- Initiated repository.
