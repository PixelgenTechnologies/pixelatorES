# Declarations used in package check
globalVariables(
  names = c(
    ".", "A", "After filtering", "anno", "attached", "B", "Before filtering",
    "cell_id", "cherry_gradient", "cluster", "COa", "COb", "column",
    "comp_id", "component", "component_n_post_filtering", "component_n_pre_filtering",
    "condition", "contrast", "control_outlier", "control_perc", "control_score",
    "coreness", "deduped_valid_reads", "Default", "degree", "dispersion_gini",
    "dispersion_outlier", "dispersion_tau", "Edge count after recovery",
    "Edge count before recovery", "edge_count_post_recovery", "edge_count_pre_recovery",
    "edges", "file_basename", "file_ext", "filename", "filtered",
    "fload", "fraction_graph_reads", "fraction_nodes_in_largest_component_post_recovery",
    "fraction_nodes_in_largest_component_pre_recovery", "fraction_valid_reads",
    "graph_edge_saturation", "graph_edges", "graph_node_saturation",
    "graph_proteins", "graph_reads", "group_id", "i", "id", "input",
    "intracellular_fraction", "is_default", "isotype_counts", "isotype_fraction",
    "join_count_expected_mean_list", "join_count_list", "k_core",
    "l1_annotation", "l1_annotation_summary", "l2_annotation", "l2_annotation_summary",
    "label", "loadedversion", "M_reads", "marker", "marker_1", "marker_2",
    "marker_i", "max_size_theshold", "mean_coreness", "mean_join_count_z",
    "mean_log2_ratio", "Median dangling nodes [%]", "Median mean coreness",
    "Median percent control markers", "Median reads per cell [thousands]",
    "Median well connected nodes [%]", "median_abs_per_cell", "median_intracellular_count_pct",
    "median_isotype_count_pct", "median_isotype_counts", "median_isotype_percent",
    "median_mean_coreness", "median_percent_dangling_nodes", "median_percent_well_connected_nodes",
    "median_reads_in_component", "median_reads_per_cell", "meta",
    "metadata", "min_size_theshold", "molecules", "n_cells", "n_cells_detected",
    "n_cells_over10k", "n_edges", "n_true", "n_umi", "name", "Node count after recovery",
    "Node count before recovery", "node_count_post_recovery", "node_count_pre_recovery",
    "nodes", "number_of_umis_removed", "output", "package", "percent",
    "percent_dangling_nodes", "percent_nodes", "percent_nodes_in_largest_post",
    "percent_nodes_in_largest_pre", "percent_well_connected_nodes",
    "plot_data", "QC Cell\nfiltering", "QC Cell filtering", "ratio",
    "reads_in_component", "reads_per_component", "removed_total",
    "Sample ID", "sample_alias", "sample_component", "sample_frac",
    "saturation", "selected_sample", "size_outlier", "stage", "stage_i",
    "top_markers", "top3_fraction", "top5_fraction", "total_edges_in",
    "total_reads", "type", "V1", "V2", "valid_reads", "valid_reads_saturation",
    "value", "z", "Cell annotation", "cell_annotation",
    "doublet_prediction", "frac", "normcount", "seurat_clusters"
  ),
  package = "pixelatorES",
  add = TRUE
)

expect_quarto <- function(...) {
  rlang::check_installed("quarto", ...)
}

expect_rstudioapi <- function(...) {
  rlang::check_installed("rstudioapi", ...)
}

expect_DBI <- function(...) {
  rlang::check_installed("DBI", ...)
}

expect_fs <- function(...) {
  rlang::check_installed("fs", ...)
}

expect_duckdb <- function(...) {
  rlang::check_installed("duckdb", ...)
}
