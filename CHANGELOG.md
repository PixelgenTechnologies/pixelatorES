# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- Fixed a bug in `plot_violin` that would throw an error if there is only a single point to be plotted. 
- Removed equations from key metrics table tooltips, as these could not be rendered properly.

### Added
- Workflow to label images with the version upon release. 

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
