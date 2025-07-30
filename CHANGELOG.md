# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- nav-pills to style tabsets
- `get_test_qc_metrics` to get qc metrics quickly for testing.
- `component_qc_molecule_rank_plot` to make the molecule rank plot from qc data.
- `get_test_data` to make a small Seurat object for testing components.

### Updated
- `process_data` now uses the public `AnnotateCells` function from `pixelatorR`. `process_data` can now switch annotation method to "nmf" by specifying `params$annotation_method = "nmf"`.
- `component_cell_recovery` now returns a tabsetted Molecule Rank Plot.
- `violin_plot` is now borderless by default.

### Fixed
- Disabled the `draw_quantile` option when a `facet_var` is provided to `plot_violin`. Otherwise, the medians are placed incorrectly.
- `title_plotlist` now returns an empty plot if the render fails.
- Overlapping text in `component_cell_recovery` plots.
- All violin plots are now using `scale = "width"`.

## [0.1.0] 2025-07-22

### Added
- Initiated repository.
