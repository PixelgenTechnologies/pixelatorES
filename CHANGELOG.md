# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Changed

- Files whose pipeline stage cannot be determined no longer abort Experiment Summary ingestion. They are skipped and recorded as `file_discovery` diagnostics.
- Colocalization heatmap text now explains that proteins are selected by mean abundance (up to 40 markers), globally across samples or within each cell type.
- The Cell recovery section now explains how to read the molecule rank plot.
- Experiment Summary labels now use "isotype control markers" consistently.
- Molecule rank plots mark components excluded by the component size thresholds with a lighter shade of their sample's colour, so it is visible which components the report drops.

## [0.14.0] 2026-08-25

### Changed

- Breaking: stage vocabularies are workflow-owned — `register_es_data_workflow()` and discovery helpers require `stages`, read via `get_es_workflow_stages()`.

## [0.13.0] 2026-08-10

### Added

- Workflows register a Quarto report recipe (`preamble` + `sections`) via `register_es_data_workflow()`, retrieved with `get_es_workflow_report()`. Child paths are validated for existence at registration.
- `test_es_data()` helper for building lightweight `es_data` fixtures in tests and downstream packages.
- Relative QC stage completeness diagnostics: samples or pools missing stages that peers have receive a `qc_load` diagnostic while keeping their partial QC data.
- Samples page report-data callout listing loading and analysis-step diagnostics, with warning markers on samples affected by loading issues.
- Run info tab (renamed from Run settings) with a conditional Diagnostics section.

### Changed

- Breaking: `register_es_data_workflow()` now requires a `report` factory. Callers that only registered extractors in 0.12.0 must supply a valid report recipe (`preamble` + `sections`).
- The Quarto shell dispatches report sections from the registered workflow recipe. Shared children live under `inst/quarto/shared/`; workflow-owned children live under `inst/quarto/workflows/<id>/`.

### Fixed

- `print_metadata_table()` no longer requires `pxl_data_processed` for hashed experiments; `% of pool` is omitted when processed data is unavailable.

## [0.12.0] 2026-08-06

### Added

- `build_es_data()` and an `es_data` object as the Experiment Summary ingestion entry point. Extension packages can register additional workflows with `register_es_data_workflow()`.

### Changed

- Ingestion is resilient: failed sample or QC loads no longer abort the whole report. Only an unreadable samplesheet stops the build.
- Breaking: `component_*()`, `key_metric_table()`, and `print_*()` helpers now take `es_data` as their first argument.

## [0.11.5] 2026-07-29

### Fixed

- `run_ES_docker.sh` now supports relative/nested output paths (`-n`), optional samplesheet in debug mode (`-D`), and correctly forwards Quarto `-P` parameters and output file ownership.
- `experiment-summary` exits with an error when Quarto rendering fails.

## [0.11.4] 2026-07-21

### Fixed

- Bug which would cause the sample confidence / hashing enrichment ratio of undetermined to not be displayed for undetermined cells in the sample confidence plot. 

## [0.11.3] 2026-06-23

### Added

- In-report deep-link anchor IDs for plot tabsets (`title_plotlist`, `tabset_plotlist`, `tabset_nested_plotlist`). Anchor IDs are derived from the knitr chunk label and tab name (e.g. `#abundance-per-marker-cd3`, `#qc-metrics-sequencing-saturation-s11`).
- ES report hash navigation in `custom.html`: opens nested tab panes for `#` links, scrolls to the target, refreshes HTMLWidgets and DataTables in newly shown tabs, and restores tab/scroll state on browser back.

### Fixed
 
- `downsample_data` will now default to the number of markers or cells in the data if `n_markers` or `n_cells` are higher than available in the data.
- Bug in `component_denoising` in the isotype reduction plot, where components with zero isotype counts would cause the summary median to become `NAs`.

## [0.11.2] 2026-06-15

### Changed

- `component_hashing` now returns a sample confidence plot with either hash purity or hash enrichment factor, depending on which metric is present in the data.
 
### Removed
 
- `harmony` has been removed, and `do_harmonize` is no longer an option. 
 
## [0.11.1] 2026-06-09

### Changed

- `read_samplesheet` now also returns "lot_role", "kit_lot_id" columns if they exist.

## [0.11.0] 2026-06-08

### Added

- `component_denoising` with plots and tables for UMIs removed, denoising by method, and isotype reduction.
- `get_denoising_detail_data()` for per-component denoising metadata from Seurat objects.
- `component_abundance_per_marker` combining sample-wise and per-celltype abundance plots.
- `DEVELOPERS.md` with linting, styling, and test data update instructions.

### Removed

- `component_abundance_per_sample` and `component_abundance_per_celltype` (replaced by `component_abundance_per_marker`).

### Changed

- Renamed `component_bleedover_noise` to `component_denoising`.
- Test layout PXL files updated with denoising metadata columns in `__adata__obs`.
- Replaced `component_abundance_per_sample` and `component_abundance_per_celltype` with `component_abundance_per_marker` (mirrors `component_proximity_per_marker`).
- `plot_violin()` now uses `ggbeeswarm::geom_quasirandom` with dodge alignment and a fixed contrasting point color instead of opaque `geom_jitter`.
- The package now depends on `pixelatorR` >= v0.18.0.

### Fixed

- Bug in `get_hash_stats` which would cause an error to be thrown for hashed samples if the samplesheet didn't contain all hash groups that exist in the QC files.
- `get_qc_metrics` now reports which data causes an error when data fails in formatting.
- Changed `run_proximity_anova` to use ":" instead of "/" as a marker pair separator.

## [0.10.5] 2026-05-21

### Changed

- `pool` IDs with duplicates in `sample` or `sample_alias` are no longer tolerated. 
- The chunk `read_qc_metrics` will now always throw a hard error if there is an error, regardless of which Quarto profile is used. This is to avoid silent failures in this crucial step of the data processing.
- The chunk `read_qc_metrics` will now always output messages and warnings.
- `get_qc_metrics` is now more lenient and will return `NULL` for extracted tables with zero rows rather than throwing an error. 

### Fixed

- Added an explicit title level for the tab set `component$sample_confidence_plots`.

## [0.10.4] 2026-05-19

### Fixed

- Bug in `component_bleedover_noise` that would cause a crash if the "pool" column exists in the indata instead of "sample_alias"

## [0.10.3] 2026-05-18

### Fixed

- Bug in `component_bleedover_noise` related to name change from `ratio` to `percent_umis_denoised`.

## [0.10.2] 2026-05-13

### Changed

- Added argument `package` to test data pointer functions `get_test_data`, `test_samplesheet` and `test_data_folder`, to allow loading test data from other packages if needed. 
- Removed redundant function `get_test_qc_metrics`.

## [0.10.1] 2026-05-13

### Fixed

- Fixed an issue that would throw an error if denoise data belonged to a pool rather than a sample.
- Bug where the component sample confidence plot would fail to display if there are no undetermined cells.
- Bug in retrieving sample calling data that would throw an error if samples are different from sample aliases.

## [0.10.0] 2026-05-12

### Added

- Molecule Rank Plots for individual samples for hashed experiments.
- `tabset_nested` mode to `tabset_figure_table` to enable nested plot lists to be displayed together with a table.

### Changed

- Switched from `SequenceSaturationCurve` to `approximate_saturation_curve` for better performance
- Removed `tests/testdata/es_test_data/` in favor of using `inst/extdata/` for test data.
- Renamed the metric `ratio` to a more descriptive name: `percent_umis_denoised`.
- `fraction_graph_reads` and `percent_umis_denoised` are now always displayed in the Key Metrics Table.
- Changed Key Metrics Table titles of "Valid reads fraction [%]" and "Graph reads fraction [%]" to "% Valid reads" and "% Graph reads", respectively.
- Removed some superfluous sequencing saturation barplot content from Quality Metrics.
- The Key Metrics Table is now horizontal again.
- Removed unused function `component_molecule_rank_plot`.
- Moved function `component_qc_molecule_rank_plot` to inside of `component_cell_recovery`.


## [0.9.2] 2026-05-08

### Added 

- Added `get_degree_distribution` to load degree distribution from QC files.
- `get_hash_stats` now collects sample confidence data.
- Metrics `Number of cells` and `% Sample called cells` have been added to the Key Metrics Table.
- Component "Component sample confidence" added to the Quality Metrics tab.

### Changed

- Added new slots to example data jsons with degree distribution. 
- .pxl files are now tracked with Git LFS.
- `get_test_data` now loads data from `extdata` to allow easier update of test data and more representative test data.
- `get_test_data` now creates a temporary folder where it puts the .pxl files to be read. This avoids `Duckdb` file accession collisions in CI.
- The Hash Purity plot's y-axis now ranges 90-100% by default.

## [0.9.1] 2026-04-30

### Added

- Expanded component and data reading tests to include hashing data.

### Fixed

- Fixed a bug in `component_crossing_edges` that would throw an error for hashed data.

## [0.9.0] 2026-04-27

### Added

- `get_hash_stats()` to summarise sample hashing from per-cell `hash_counts` metadata (purity, heatmap inputs, sample-level stats).
- `component_hashing()` to render hash purity violin, summary table, and hash purity / hash fraction heatmaps.
- `sample_calling` as a recognised `pixelator` pipeline stage in `find_stage()` / `get_file_paths()` file resolution.
- `inst/quarto/quality_metrics.qmd` now supports hashed experiments: optional pool vs sample key-metric layouts, conditional sample hashing subsections.
- `pipeline_pool_stages` containing the stages of the nf-core/pixelator pipeline that produce pool-level QC files.
- Added a Samples tab `inst/quarto/samples.qmd`.
- Added render utilities `section_table` and `section_intro`.

### Changed

- `read_qc_files()` returns a nested list `list(qc_files = <named list per sample>, pool_qc_files = NULL | <named list per pool>)`. Passing a data frame of QC paths (legacy) is still supported; the return value is wrapped in the same structure.
- `key_metric_table()` returns `list(sample = <datatables>, pool = <datatables or NULL>)` for HTML output, and `return_data = TRUE` returns the corresponding pre-styled wide tibbles in the same shape. Hashing columns are included in definitions when `sample_hash_stats` is present.
- `get_file_paths()` now accepts `sample_sheet` instead of `sample_aliases` (previous behaviour).
- `read_samplesheet()` adds a `pool` column when the CSV has a pool column.
- `merge_data()` joins optional `pool` from the sample sheet into object metadata when that column exists.
- `get_qc_metrics()` adds optional `sample_hash_stats`.
- `get_crossing_edges()` and `get_denoising_data()` accept nested QC lists from `read_qc_files()`.
- `extract_sample_qc_metrics()` automatically extracts content from the appropriate slot `qc_files` or `pool_qc_files` when given the full return value of `read_qc_files()`.
- `component_sequencing_reads_and_molecules()`, `component_cell_recovery()`, and `component_qc_molecule_rank_plot()` select pool-level vs sample-level QC JSON lists when `pool_qc_files` is present.
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
