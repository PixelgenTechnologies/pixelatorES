# Tests mutate the run, adding reports and rewriting sample_id, so they work on
# a temp copy rather than the installed fixture. The fixture ships the layout
# PXL a real run hands to the report; the upstream stages that also emit a PXL
# are added here as empty placeholders, since discovery never opens a PXL and
# their presence is what gives `pxl_preference` something to rank.
materialise_edgelist_aggregate_run <- function() {
  run_dir <- tempfile("edgelist_aggregate_run_")
  dir.create(run_dir)

  file.copy(
    list.files(test_data_folder("edgelist_aggregate"), full.names = TRUE),
    run_dir,
    recursive = TRUE
  )

  placeholders <- c(
    file.path("graph", paste0(c("pool1", "pool2"), ".graph.pxl")),
    file.path("analysis", paste0("sample", 1:4, ".analysis.pxl")),
    file.path("post_analysis", paste0("sample", 1:4, ".post_analysis.pxl")),
    file.path("sample_calling", paste0("sample", 1:4, ".dehashed.pxl"))
  )

  for (rel in file.path(run_dir, placeholders)) {
    dir.create(dirname(rel), recursive = TRUE, showWarnings = FALSE)
    file.create(rel)
  }

  return(run_dir)
}

edgelist_aggregate_file_paths <- function(run_dir) {
  return(get_file_paths(
    data_folder = run_dir,
    sample_sheet = read_samplesheet(test_samplesheet("edgelist_aggregate")),
    stages = get_es_workflow_stages("edgelist_aggregate")
  ))
}

edgelist_aggregate_es_data_stub <- function(run_dir) {
  return(list(
    params = list(data_folder = run_dir),
    samplesheet = read_samplesheet(test_samplesheet("edgelist_aggregate")),
    stages = get_es_workflow_stages("edgelist_aggregate"),
    file_paths = edgelist_aggregate_file_paths(run_dir)
  ))
}

edgelist_aggregate_qc_raw <- function(run_dir) {
  return(pixelatorES:::.extract_qc_raw(
    edgelist_aggregate_es_data_stub(run_dir)
  )$value)
}
