test_that("Parameter tables work as expected", {
  expect_no_error(
    print_params(structure(list(
      sample_sheet = "../../../ES_test_data/test_samplesheet.csv",
      data_folder = "../../../ES_test_data",
      control_markers = "c('mIgG1', 'mIgG2a', 'mIgG2b')",
      molecule_rank_cutoff = 10000L,
      gini_dispersion_cutoff = 0.7,
      control_score_cutoff = 1.5,
      do_harmonize = FALSE,
      harmonization_vars = "'condition'",
      norm_method = "CLR",
      clustering_resolution = 1L,
      run_proximity_analysis = TRUE,
      proximity_count_cutoff = 25L,
      run_qc_metrics = TRUE
    ), class = "knit_param_list"))
  )

  expect_no_error(print_metadata_table(structure(list(
    sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2"), design = c("pna-2", "pna-2"), panel = c(
      "pna-rnd-158plex-final",
      "pna-rnd-158plex-final"
    ), fastq_1 = c(
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample1_S1_R1_001.fastq.gz",
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample2_S2_R1_001.fastq.gz"
    ), fastq_2 = c(
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample1_S1_R2_001.fastq.gz",
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample2_S2_R2_001.fastq.gz"
    ), sample_alias = c("S1", "S2"), condition = c("PBMC", "PHA"),
    reference_condition = c(TRUE, FALSE), alternative_condition = c(
      FALSE,
      TRUE
    )
  ), row.names = c(NA, -2L), class = c(
    "tbl_df", "tbl",
    "data.frame"
  ))))

  expect_no_error(print_session_info())
})
