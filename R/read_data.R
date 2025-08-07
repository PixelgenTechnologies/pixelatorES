#' nf-core/pixelator stages
#'
pipeline_stages <- c(
  "amplicon", "collapse",
  "demux", "denoise",
  "graph", "analysis",
  "post_analysis", "layout"
)


#' Find pixelator stage of file
#'
#' Determine which pixelator stage a file comes from
#'
#' @param filepath A file path
#'
#' @return The name of the stage it comes from
#'
#' @export
#'
find_stage <-
  function(filepath) {
    pixelatorR:::assert_single_value(filepath, type = "string")

    suffix <-
      basename(filepath) %>%
      str_extract("\\..*\\.") %>%
      str_remove("\\.dataset") %>%
      str_remove_all("\\.")

    filedir <-
      dirname(filepath) %>%
      str_remove(".*/")

    if (suffix %in% pipeline_stages) {
      return(suffix)
    } else if (filedir %in% pipeline_stages) {
      return(filedir)
    } else {
      cli_abort(
        c(
          "x" = "Could not determine stage for file: {.val {filepath}}"
        )
      )
    }
  }

#' Get file paths for data and QC files
#'
#' Retrieve file paths for data and QC files from a specified data folder.
#'
#' @param data_folder A character string specifying the path to the data folder.
#' @param file_paths A character vector with file paths including all data files.
#' @param sample_aliases A named character vector mapping sample aliases to their actual names.
#'
#' @return A list containing two data frames: `data_files` and `qc_files`.
#'
#' @export
#'
get_file_paths <-
  function(data_folder = NULL, file_paths = NULL, sample_aliases) {
    pixelatorR:::assert_single_value(data_folder, type = "string", allow_null = TRUE)
    pixelatorR:::assert_vector(file_paths, "character", allow_null = TRUE)
    pixelatorR:::assert_vector(sample_aliases, "character", allow_null = TRUE)
    pixelatorR:::assert_vector(names(sample_aliases), "character", allow_null = FALSE)

    if (is.null(file_paths)) {
      file_paths <- list.files(data_folder, recursive = TRUE, full.names = TRUE)
    }

    # Look for files in data folder and filter for sample files
    all_files <-
      file_paths %>%
      enframe("i", "filename") %>%
      mutate(file_ext = str_remove(filename, ".*\\.")) %>%
      filter(file_ext %in% c("json", "pxl")) %>%
      filter(!str_detect(filename, "pipeline_info")) %>%
      mutate(
        file_basename = basename(filename),
        stage = sapply(filename, find_stage),
        sample_alias = str_remove(file_basename, "\\..*")
      )

    if (!is.null(sample_aliases)) {
      all_files <-
        all_files %>%
        mutate(sample_alias = sample_aliases[sample_alias]) %>%
        filter(!is.na(sample_alias))
    }

    # Filter for data and QC files
    data_files <-
      all_files %>%
      filter(
        str_detect(file_basename, "\\.pxl$"),
        stage %in% c("graph", "analysis", "post_analysis", "layout")
      ) %>%
      mutate(stage_i = unclass(factor(stage, c("graph", "analysis", "post_analysis", "layout")))) %>%
      group_by(sample_alias) %>%
      top_n(1, stage_i) %>%
      ungroup() %>%
      select(sample_alias, filename)

    if (nrow(data_files) > length(unique(data_files$sample_alias))) {
      dup_files <- data_files$sample_alias[duplicated(data_files$sample_alias)]
      cli_abort(
        c(
          "x" = "Some samples match multiple data files: {.val {dup_files}}"
        )
      )
    }

    qc_files <-
      all_files %>%
      filter(str_detect(file_basename, "\\.report.json$")) %>%
      filter(!(str_detect(filename, "\\.part_\\d{3}\\.") & stage == "collapse")) %>%
      select(sample_alias, filename, stage)

    return(list(data_files = data_files, qc_files = qc_files))
  }


#' Load pxl files into a list
#'
#' Load pxl files into a list as specified by identified data paths and sample sheet.
#'
#' @param data_folder Path to the folder containing the data files.
#' @param data_files Paths to data files.
#' @param sample_sheet Sample sheet.
#'
#' @return A list of Seurat objects
#'
#' @export
#'
load_pxl_data_list <-
  function(data_folder, data_files, sample_sheet) {
    # Check that all expected data files exist
    missing_samples <-
      sample_sheet %>%
      filter(!sample_alias %in% data_files$sample_alias) %>%
      pull(sample)

    if (length(missing_samples) > 0) {
      found_samples <-
        list.files(data_folder, recursive = TRUE, pattern = ".pxl") %>%
        basename() %>%
        unique() %>%
        str_remove("\\..*")

      cli_abort("The following samples are missing data files: {missing_samples}\n
                The following samples were found: {found_samples}")
    }


    # Create a list of Seurat objects
    pg_data <-
      data_files %>%
      select(sample_alias, filename) %>%
      deframe() %>%
      lapply(ReadPNA_Seurat,
        load_proximity_scores = FALSE
      )

    if (length(pg_data) != nrow(sample_sheet)) {
      cli_abort("The number of loaded samples do not match the expected number.\n
                {length(pg_data)} samples are loaded: {names(pg_data)}.\n
                The samplesheet contains {nrow(sample_sheet)} samples: {sample_sheet$sample_alias}")
    }

    return(pg_data)
  }


#' Merge data from multiple samples into a single Seurat object
#'
#' Merges data from list of multiple samples into a single Seurat object, adding metadata and ranks based on UMI counts.
#'
#' @param pg_data A list of Seurat objects, each representing a sample.
#' @param sample_sheet A data frame containing sample metadata, including `sample_alias` and `condition`.
#'
#' @return A merged Seurat object with added metadata.
#'
#' @export
#'
merge_data <-
  function(pg_data,
             sample_sheet) {
    pg_data <-
      pg_data %>%
      map(. %>%
        AddMetaData(
          FetchData(., "n_umi") %>%
            mutate(rank = rank(-n_umi, ties.method = "random"))
        ))


    pg_data <-
      merge(
        pg_data[[1]],
        y = pg_data[-1],
        add.cell.ids = names(pg_data)
      ) %>%
      JoinLayers()

    metadata <-
      sample_sheet %>%
      select(sample_alias, condition)

    pg_data <-
      pg_data %>%
      AddMetaData(colnames(pg_data) %>%
        enframe("i", "comp_id") %>%
        mutate(
          sample_alias = str_extract(comp_id, ".*_") %>%
            str_remove("_$")
        ) %>%
        left_join(
          metadata,
          by = "sample_alias"
        ) %>%
        mutate(
          sample_alias = factor(sample_alias, metadata$sample_alias),
          condition = factor(condition, unique(metadata$condition))
        ) %>%
        select(-i) %>%
        column_to_rownames("comp_id"))

    return(pg_data)
  }

#' Downsample data to a specified number of cells and markers
#'
#' Downsamples the data to a specified number of cells and markers, ensuring that control markers are always included.
#'
#' @param pg_data A Seurat object containing the data to be downsampled.
#' @param control_markers A character vector of control markers to always include in the downsampled data.
#' @param n_cells An integer specifying the number of cells to keep in each sample.
#' @param n_markers An integer specifying the total number of markers to keep in the downsampled data.
#'
#' @return A downsampled Seurat object with the specified number of cells and markers.
#'
#' @export
#'
downsample_data <-
  function(pg_data,
             control_markers = NULL,
             n_cells = 50,
             n_markers = 20) {
    set.seed(37)

    keep_cells <-
      FetchData(pg_data, "sample_alias") %>%
      as_tibble(rownames = "cell_id") %>%
      group_by(sample_alias) %>%
      slice_sample(n = n_cells) %>%
      pull(cell_id)


    keep_markers <-
      rownames(pg_data) %>%
      {
        .[!. %in% control_markers]
      } %>%
      {
        .[sample(seq_len(length(.)), size = n_markers - length(control_markers), replace = FALSE)]
      } %>%
      union(control_markers)

    pg_data <-
      pg_data[keep_markers, keep_cells]

    return(pg_data)
  }

#' Read QC files and return metrics
#'
#' Reads QC files and returns a list of metrics for each sample.
#'
#' @param qc_files A data frame containing the file paths of json QC files.
#' @param sample_sheet A data frame containing sample metadata, including `sample_alias`.
#'
#' @return A list of QC metrics for each sample, where each element is a named list of metrics.
#'
#' @export
read_qc_files <-
  function(qc_files,
             sample_sheet) {
    pixelatorR:::assert_within_limits(nrow(qc_files), c(1, Inf))
    pixelatorR:::assert_class(qc_files, "data.frame")
    for (f in qc_files$filename) pixelatorR:::assert_file_exists(f)

    temp_data <-
      qc_files %>%
      group_by(sample_alias)

    sample_qc_metrics <-
      temp_data %>%
      group_split() %>%
      set_names(group_keys(temp_data)$sample_alias) %>%
      lapply(function(sample_files) {
        sample_files %>%
          select(stage, filename) %>%
          deframe() %>%
          lapply(fload)
      })

    # If some qc metrics are missing
    if (!all(sample_sheet$sample_alias %in% qc_files$sample_alias)) {
      cli_abort(
        "Some samples are missing .json files containing QC metrics.
        Please check the following samples:
        {sample_sheet$sample_alias[!sample_sheet$sample_alias %in%
        qc_files$sample_alias]}"
      )
    }

    return(sample_qc_metrics)
  }

#' Extract sample QC metrics from a list of sample QC data
#'
#' This function extracts specific quality control metrics from a list of sample
#' QC data (as generated by `read_qc_files`). It returns a tibble with the
#' sample alias and the specified metrics.
#'
#' @param sample_qc_metrics A list of sample QC data, where each element is a named list
#' containing QC metrics for a sample.
#' @param vars A character vector of variable names to extract from the sample QC data.
#' @param stage A character string specifying the stage from which to extract the metrics.
#'
#' @return A tibble with the sample alias and the specified metrics.
#'
#' @export
#'
extract_sample_qc_metrics <-
  function(sample_qc_metrics, vars, stage = "graph") {
    pixelatorR:::assert_vector(vars, "character", n = 1)
    pixelatorR:::assert_single_value(stage, type = "string")
    pixelatorR:::assert_class(sample_qc_metrics, "list")
    stage <- match.arg(stage, pipeline_stages)

    extracted_data <-
      sample_qc_metrics %>%
      lapply(function(sample_qc_data) {
        if (any(sapply(sample_qc_data[[stage]][vars], length) > 1)) {
          cli_abort("The following variables contain multiple values for a
                    single sample:
                    {.val {vars[sapply(sample_qc_data[[stage]][vars],
                    length) > 1]}}.")
        }

        sample_qc_data[[stage]][vars] %>%
          map(. %>%
            {
              ifelse(is.null(.), NA, .)
            }) %>%
          enframe() %>%
          unnest(value) %>%
          pivot_wider(names_from = name, values_from = value)
      }) %>%
      bind_rows(.id = "sample_alias")

    if (!is.null(names(vars))) {
      extracted_data <-
        extracted_data %>%
        rename(!!!vars)
    }

    return(extracted_data)
  }

#' Read a sample sheet from a CSV file
#'
#' Reads a sample sheet from a CSV file and returns it as a tibble.
#'
#' @param filepath A character string specifying the path to the sample sheet CSV file.
#'
#' @return A tibble containing the sample sheet data with columns `sample`, `sample_alias`, and `condition`.
#'
#' @export
#'
read_samplesheet <-
  function(filepath) {
    pixelatorR:::assert_single_value(filepath, type = "string")
    pixelatorR:::assert_file_exists(filepath)

    read_csv(filepath,
      col_types = cols(
        sample = col_character(),
        sample_alias = col_character(),
        condition = col_character()
      )
    ) %>%
      select(sample, sample_alias, condition) %>%
      mutate(sample_alias = ifelse(is.na(sample_alias), sample, sample_alias)) %>%
      distinct()
  }


#' Get the path to the test samplesheet
#'
#' Returns the path to the test samplesheet CSV file included in the package.
#'
#' @return A character string with the path to the test samplesheet CSV file.
#'
#' @export
#'
test_samplesheet <-
  function() {
    system.file(
      "extdata", "test_samplesheet.csv",
      package = "pixelatorES"
    )
  }

#' Get the path to the test data folder
#'
#' Returns the path to the test data folder containing QC JSON files included in the package.
#'
#' @return A character string with the path to the test data folder.
#'
#' @export
#'
test_data_folder <-
  function() {
    system.file(
      "extdata", "qc_jsons",
      package = "pixelatorES"
    )
  }

#' Get test QC metrics
#'
#' Reads the test samplesheet and retrieves QC metrics from the test data folder.
#'
#' @return A list of QC metrics for each sample, where each element is a named list of metrics.
#'
#' @export
#'
get_test_qc_metrics <-
  function() {
    sample_sheet <- read_samplesheet(test_samplesheet())

    data_paths <-
      get_file_paths(
        data_folder = test_data_folder(),
        sample_aliases = sample_sheet %>%
          select(sample, sample_alias) %>%
          deframe()
      )

    sample_qc_metrics <-
      read_qc_files(data_paths$qc_files, sample_sheet)

    sample_qc_metrics <- read_qc_files(data_paths$qc_files, sample_sheet)

    return(sample_qc_metrics)
  }



#' Get test data
#'
#' Generates a minimal Seurat object for testing purposes.
#'
#' @param concatenate A logical indicating whether to concatenate the data 6 times (default is TRUE).
#'
#' @return A Seurat object containing test data with normalized and scaled data, PCA results, and merged layers.
#'
#' @export
#'
get_test_data <-
  function(
    concatenate = TRUE
  ) {
    seur <-
      minimal_pna_pxl_file() %>%
      ReadPNA_Seurat(load_proximity_scores = FALSE)

    if (concatenate) {
      seur <-
        merge(
          seur,
          y = list(seur, seur, seur, seur, seur)
        ) %>%
        JoinLayers(verbose = FALSE)
    }

    seur[[]]$sample_alias <-
      "S1"

    seur <-
      seur %>%
      NormalizeData(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(
        features = rownames(.)[1:10],
        npcs = 2,
        verbose = FALSE,
        approx = TRUE
      )

    return(seur)
  }
