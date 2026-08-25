test_that("Sample diagnostics work as expected", {
  samplesheet <- tibble(
    pool = c("pool1", "pool1"),
    sample = c("one", "two"),
    sample_alias = c("S1", "S2"),
    condition = c("A", "B")
  )

  clean <- test_es_data(samplesheet = samplesheet)
  expect_equal(
    list(
      tibble = diagnostics_to_tibble(clean),
      has_sample = has_sample_diagnostics(clean),
      targets = pixelatorES:::sample_diagnostic_targets(clean),
      summary = format_sample_diagnostics_summary(clean),
      callout = format_sample_diagnostics_callout(clean)
    ),
    list(
      tibble = tibble(
        type = character(),
        target = character(),
        message = character()
      ),
      has_sample = FALSE,
      targets = character(),
      summary = NULL,
      callout = NULL
    )
  )

  loading <- test_es_data(
    samplesheet = samplesheet,
    diagnostics = list(
      list(
        type = "pxl_load",
        target = "S1",
        message = "No PXL file was found."
      ),
      list(
        type = "qc_load",
        target = "pool1",
        message = "No QC files were found."
      )
    )
  )
  expect_equal(
    list(
      has_sample = has_sample_diagnostics(loading),
      targets = pixelatorES:::sample_diagnostic_targets(loading),
      summary = format_sample_diagnostics_summary(loading),
      callout = format_sample_diagnostics_callout(loading)
    ),
    list(
      has_sample = TRUE,
      targets = c("S1", "S2"),
      summary = paste(
        "- **S1** (PXL loading): No PXL file was found.",
        "- **pool1** (QC loading): No QC files were found.",
        sep = "\n"
      ),
      callout = paste0(
        '::: {.callout-important title="Report data issues"}\n',
        "Some input data could not be loaded or some analyses could not be ",
        "completed, and the metrics in this report are therefore incomplete. ",
        "Affected samples are marked with a warning symbol in the table below.",
        "\n\n",
        "- **S1** (PXL loading): No PXL file was found.\n",
        "- **pool1** (QC loading): No QC files were found.",
        "\n\nSee the Diagnostics section under Run info for the complete list.\n",
        ":::\n"
      )
    )
  )

  extractor_only <- test_es_data(
    samplesheet = samplesheet,
    diagnostics = list(list(
      type = "extractor",
      target = "pxl_data_processed",
      message = "You've supplied a <NULL> object."
    ))
  )
  expect_equal(
    list(
      has_sample = has_sample_diagnostics(extractor_only),
      targets = pixelatorES:::sample_diagnostic_targets(extractor_only),
      summary = format_sample_diagnostics_summary(extractor_only),
      callout = format_sample_diagnostics_callout(extractor_only)
    ),
    list(
      has_sample = FALSE,
      targets = character(),
      summary = paste0(
        "- **pxl_data_processed** (Analysis step): ",
        "You've supplied a &lt;NULL&gt; object."
      ),
      callout = paste0(
        '::: {.callout-important title="Report data issues"}\n',
        "Some input data could not be loaded or some analyses could not be ",
        "completed, and the metrics in this report are therefore incomplete.",
        "\n\n",
        "- **pxl_data_processed** (Analysis step): ",
        "You've supplied a &lt;NULL&gt; object.",
        "\n\nSee the Diagnostics section under Run info for the complete list.\n",
        ":::\n"
      )
    )
  )

  mixed <- test_es_data(
    samplesheet = samplesheet,
    diagnostics = list(
      list(
        type = "pxl_load",
        target = "S1",
        message = "No PXL file was found."
      ),
      list(
        type = "extractor",
        target = "proximity",
        message = "Unavailable."
      )
    )
  )
  expect_equal(
    list(
      has_sample = has_sample_diagnostics(mixed),
      targets = pixelatorES:::sample_diagnostic_targets(mixed),
      summary = format_sample_diagnostics_summary(mixed)
    ),
    list(
      has_sample = TRUE,
      targets = "S1",
      summary = paste(
        "- **S1** (PXL loading): No PXL file was found.",
        "- **proximity** (Analysis step): Unavailable.",
        sep = "\n"
      )
    )
  )
})

test_that("Relative QC stage completeness diagnostics work as expected", {
  equal_values <- list(
    S1 = list(analysis = list(), denoise = list()),
    S2 = list(analysis = list(), denoise = list())
  )
  expect_equal(
    pixelatorES:::.qc_stage_completeness_diagnostics(equal_values, "sample"),
    list()
  )

  expect_warning(
    incomplete <- pixelatorES:::.qc_stage_completeness_diagnostics(
      list(
        S1 = list(analysis = list(), denoise = list(), graph = list()),
        S2 = list(analysis = list())
      ),
      "sample"
    ),
    "Incomplete QC data for sample"
  )
  expect_equal(
    incomplete,
    list(list(
      type = "qc_load",
      target = "S2",
      message = "Missing QC stages: denoise, graph."
    ))
  )

  sample_sheet <- read_samplesheet(test_samplesheet(type = "hashing"))
  data_folder <- tempfile("pixelatorES_test_hashing_")
  dir.create(data_folder)
  stopifnot(all(file.copy(
    list.files(test_data_folder(type = "hashing"), full.names = TRUE),
    data_folder,
    recursive = TRUE
  )))
  file.remove(file.path(data_folder, "denoise", "S2.report.json"))
  file.remove(file.path(data_folder, "graph", "pool2.report.json"))
  file_paths <- get_file_paths(
    data_folder = data_folder,
    sample_sheet = sample_sheet,
    stages = get_es_workflow_stages("amplicon_demux")
  )

  expect_warning(
    sample_result <- pixelatorES:::.read_qc_groups_soft(
      files = file_paths$qc_files,
      aliases = sort(sample_sheet$sample_alias),
      alias_column = "sample_alias",
      target_label = "sample"
    ),
    "Incomplete QC data for sample"
  )
  expect_equal(
    lapply(
      Filter(
        function(diagnostic) {
          return(grepl("^Missing QC stages", diagnostic$message))
        },
        sample_result$diagnostics
      ),
      function(diagnostic) {
        return(diagnostic[c("type", "target", "message")])
      }
    ),
    list(list(
      type = "qc_load",
      target = "S2",
      message = "Missing QC stages: denoise."
    ))
  )

  expect_warning(
    pool_result <- pixelatorES:::.read_qc_groups_soft(
      files = file_paths$pool_qc_files,
      aliases = sort(unique(sample_sheet$pool)),
      alias_column = "pool_alias",
      target_label = "pool"
    ),
    "Incomplete QC data for pool"
  )
  expect_equal(
    lapply(
      Filter(
        function(diagnostic) {
          return(grepl("^Missing QC stages", diagnostic$message))
        },
        pool_result$diagnostics
      ),
      function(diagnostic) {
        return(diagnostic[c("type", "target", "message")])
      }
    ),
    list(list(
      type = "qc_load",
      target = "pool2",
      message = "Missing QC stages: graph."
    ))
  )

  expect_warning(
    zero_qc_result <- pixelatorES:::.read_qc_groups_soft(
      files = file_paths$qc_files[
        file_paths$qc_files$sample_alias == "S1",
        ,
        drop = FALSE
      ],
      aliases = c("S1", "S2"),
      alias_column = "sample_alias",
      target_label = "sample"
    ),
    "No QC files were found"
  )
  expect_equal(
    lapply(
      zero_qc_result$diagnostics,
      function(diagnostic) {
        return(diagnostic[c("type", "target", "message")])
      }
    ),
    list(list(
      type = "qc_load",
      target = "S2",
      message = "No QC files were found."
    ))
  )

  default_sheet <- read_samplesheet(test_samplesheet(type = "default"))
  default_paths <- get_file_paths(
    data_folder = test_data_folder(type = "default"),
    sample_sheet = default_sheet,
    stages = get_es_workflow_stages("amplicon_demux")
  )
  default_result <- pixelatorES:::.read_qc_groups_soft(
    files = default_paths$qc_files,
    aliases = sort(default_sheet$sample_alias),
    alias_column = "sample_alias",
    target_label = "sample"
  )
  expect_equal(
    Filter(
      function(diagnostic) {
        return(grepl("^Missing QC stages", diagnostic$message))
      },
      default_result$diagnostics
    ),
    list()
  )
})
