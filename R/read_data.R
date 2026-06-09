#' nf-core/pixelator stages
#'
pipeline_stages <- c(
  "amplicon", "collapse",
  "demux", "denoise",
  "graph", "analysis",
  "sample_calling",
  "post_analysis", "layout"
)

#' nf-core/pixelator pool stages
#'
pipeline_pool_stages <- c(
  "amplicon", "collapse",
  "demux", "graph"
)

#' Find pixelator stage of file
#'
#' Determine which pixelator stage a file comes from
#'
#' @param filepath A file path
#' @param allow_unknown Whether to allow unknown stages (default: FALSE)
#'
#' @return The name of the stage it comes from
#'
#' @export
#'
find_stage <-
  function(filepath, allow_unknown = FALSE) {
    pixelatorR:::assert_single_value(filepath, type = "string")
    pixelatorR:::assert_single_value(allow_unknown, type = "bool")

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
    } else if (allow_unknown) {
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
#' @param sample_sheet A sample sheet [tibble::tibble()] including `sample`, `sample_alias`, and
#'   optionally `pool` for hashed experiments. When supplied, pool-level QC files are detected.
#' @param allow_unknown A logical indicating whether to allow files from unknown stages (default: FALSE).
#'
#' @return A list with `data_files`, `qc_files`, and `pool_qc_files` (the last is `NULL` when absent).
#'
#' @export
#'
get_file_paths <-
  function(data_folder = NULL, file_paths = NULL, sample_sheet = NULL, allow_unknown = FALSE) {
    pixelatorR:::assert_single_value(data_folder, type = "string", allow_null = TRUE)
    pixelatorR:::assert_vector(file_paths, "character", allow_null = TRUE)
    pixelatorR:::assert_single_value(allow_unknown, type = "bool")

    if (is.null(file_paths)) {
      file_paths <- list.files(data_folder, recursive = TRUE, full.names = TRUE)
    }

    pixelatorR:::assert_class(sample_sheet, "tbl_df", allow_null = TRUE)

    all_files <-
      file_paths %>%
      enframe("i", "filename") %>%
      mutate(file_ext = str_remove(filename, ".*\\.")) %>%
      filter(file_ext %in% c("json", "pxl")) %>%
      filter(!str_detect(filename, "pipeline_info")) %>%
      mutate(
        file_basename = basename(filename),
        stage = sapply(filename, find_stage, allow_unknown),
        file_alias = str_remove(file_basename, "\\..*")
      )

    if ("pool" %in% names(sample_sheet)) {
      pool_ids <- unique(sample_sheet$pool)
    } else {
      pool_ids <- c()
    }

    all_files <-
      all_files %>%
      left_join(select(sample_sheet, sample, sample_alias),
        by = c("file_alias" = "sample")
      ) %>%
      mutate(is_pool = ifelse(file_alias %in% pool_ids, TRUE, FALSE)) %>%
      filter(is_pool | !is.na(sample_alias))

    data_files <-
      all_files %>%
      filter(!is_pool) %>%
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

    all_qc_files <-
      all_files %>%
      filter(str_detect(file_basename, "\\.report.json$")) %>%
      filter(!(str_detect(filename, "\\.part_\\d{3}\\.") & stage == "collapse"))

    qc_files <-
      all_qc_files %>%
      filter(!is_pool) %>%
      select(sample_alias, filename, stage)

    if (any(all_files$is_pool)) {
      pool_qc_files <-
        all_qc_files %>%
        filter(is_pool) %>%
        select(pool_alias = file_alias, filename, stage)
    } else {
      pool_qc_files <- NULL
    }

    return(
      list(
        data_files = data_files,
        qc_files = qc_files,
        pool_qc_files = pool_qc_files
      )
    )
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
#' @param sample_sheet A data frame containing sample metadata, including `sample_alias` and
#'   `condition`, and optionally `pool` for hashed experiments.
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
      select(sample_alias, any_of("pool"), condition)

    metadata_to_join <-
      tibble(comp_id = colnames(pg_data)) %>%
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
      )

    if ("pool" %in% names(metadata)) {
      metadata_to_join <-
        metadata_to_join %>%
        mutate(pool = factor(pool, unique(metadata$pool)))
    }


    pg_data <-
      pg_data %>%
      AddMetaData(column_to_rownames(metadata_to_join, "comp_id"))

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

    pixelatorR:::assert_x_in_y(control_markers, rownames(pg_data))

    keep_markers <-
      rownames(pg_data) %>%
      {
        .[!. %in% control_markers]
      } %>%
      {
        .[sample(seq_along(.), size = n_markers - length(control_markers), replace = FALSE)]
      } %>%
      union(control_markers)

    pg_data <-
      pg_data[keep_markers, keep_cells]

    return(pg_data)
  }


#' Utility function to add percentage of pool to the sample sheet
#'
#' Adds a column to the sample sheet indicating how much of the
#' pool each sample makes up.
#'
#' @param sample_sheet A data frame containing sample metadata,
#' including `sample_alias` and `pool` for a hashed sample.
#' @param object A Seurat object containing the sample metadata
#' in its `meta.data` slot.
#'
#' @return A modified sample sheet with an additional column `pool_fraction`.
#'
#' @export
add_pct_of_pool_to_samplesheet <-
  function(sample_sheet, object) {
    sample_sheet <- sample_sheet %>%
      left_join(
        object[[]] %>%
          select(sample_alias) %>%
          group_by(sample_alias) %>%
          count(),
        by = "sample_alias"
      ) %>%
      group_by(pool) %>%
      mutate(
        pool_fraction = round((n / sum(n)) * 100, 2)
      ) %>%
      select(-n) %>%
      ungroup()
    return(sample_sheet)
  }

#' Read QC files and return metrics
#'
#' Reads QC files and returns a nested list suitable for [get_qc_metrics()].
#'
#' @param qc_input Either a data frame of per-sample QC paths (legacy), or the list returned by
#'   [get_file_paths()] (with `qc_files` and optionally `pool_qc_files`).
#' @param sample_sheet A data frame containing sample metadata, including `sample_alias`.
#'
#' @return A list with `qc_files` (named list per sample) and `pool_qc_files` (`NULL` or named list per pool).
#'
#' @export
read_qc_files <-
  function(qc_input,
           sample_sheet) {
    if (is.data.frame(qc_input)) {
      qc_files <- qc_input
      pixelatorR:::assert_within_limits(nrow(qc_files), c(1, Inf))
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

      if (!all(sample_sheet$sample_alias %in% qc_files$sample_alias)) {
        cli_abort(
          "Some samples are missing .json files containing QC metrics.
        Please check the following samples:
        {sample_sheet$sample_alias[!sample_sheet$sample_alias %in%
        qc_files$sample_alias]}"
        )
      }

      return(list(qc_files = sample_qc_metrics, pool_qc_files = NULL))
    }

    pixelatorR:::assert_class(qc_input, "list")
    pixelatorR:::assert_x_in_y(
      c("qc_files", "pool_qc_files"),
      names(qc_input)
    )

    files_to_read <-
      qc_input[c("qc_files", "pool_qc_files")]

    for (f in files_to_read$qc_files$filename) {
      pixelatorR:::assert_file_exists(f)
    }
    if (!is.null(files_to_read$pool_qc_files)) {
      for (f in files_to_read$pool_qc_files$filename) {
        pixelatorR:::assert_file_exists(f)
      }
    }

    sample_qc_metrics <-
      lapply(
        files_to_read,
        function(files) {
          if (is.null(files)) {
            return(NULL)
          }

          temp_data <-
            files %>%
            rename(alias = any_of(c("sample_alias", "pool_alias"))) %>%
            group_by(alias)

          temp_data %>%
            group_split() %>%
            set_names(group_keys(temp_data)$alias) %>%
            lapply(function(sample_files) {
              sample_files %>%
                select(stage, filename) %>%
                deframe() %>%
                lapply(fload)
            })
        }
      )


    samples_in_data <- c(names(sample_qc_metrics$qc_files), names(sample_qc_metrics$pool_qc_files))
    if (!all(sample_sheet$sample_alias %in% samples_in_data)) {
      cli_abort(
        "Some samples are missing .json files containing QC metrics.
        Please check the following samples:
        {sample_sheet$sample_alias[!sample_sheet$sample_alias %in%
        samples_in_data]}"
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

    # Check if the sample_qc_metrics a sample qc list, if so extract the correct sample/pool type
    if ("qc_files" %in% names(sample_qc_metrics)) {
      if (stage %in% pipeline_pool_stages && !is.null(sample_qc_metrics$pool_qc_files)) {
        sample_qc_metrics <- sample_qc_metrics$pool_qc_files
      } else {
        sample_qc_metrics <- sample_qc_metrics$qc_files
      }
    }

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
#' @return A tibble with columns `sample`, `sample_alias`, `condition`, and `pool` (the latter all
#'   `NA` when the file has no pool column).
#'
#' @export
#'
read_samplesheet <-
  function(filepath, additional_columns = NULL) {
    pixelatorR:::assert_single_value(filepath, type = "string")
    pixelatorR:::assert_file_exists(filepath)

    sample_sheet <-
      read_csv(filepath,
        col_types = cols(
          sample = col_character(),
          sample_alias = col_character(),
          condition = col_character()
        )
      )

    if (!"sample_alias" %in% names(sample_sheet)) {
      sample_sheet$sample_alias <- NA
    }
    if (!"condition" %in% names(sample_sheet)) {
      sample_sheet$condition <- NA
    }

    sample_sheet <-
      sample_sheet %>%
      select(sample, sample_alias, condition, any_of(c("pool", "lot_role", "kit_lot_id"))) %>%
      mutate(
        sample_alias = ifelse(is.na(sample_alias), sample, sample_alias),
        condition = ifelse(is.na(condition), sample_alias, condition)
      ) %>%
      distinct()

    if ("pool" %in% names(sample_sheet)) {
      if (any(c(sample_sheet$sample, sample_sheet$sample_alias) %in% sample_sheet$pool)) {
        cli_abort("The `pool` column in the sample sheet cannot contain values that are also present in the `sample`
                  or `sample_alias` columns.")
      }
    }

    return(sample_sheet)
  }


#' Get the path to the test samplesheet
#'
#' Returns the path to the test samplesheet CSV file included in the package.
#'
#' @param type A character string specifying the type of test samplesheet to return. Options are "default"
#' (the standard test samplesheet) and "hashing" (a test samplesheet for hashing experiments). Default is "default".
#' @param package A character string specifying the name of the package where the test samplesheet is located.
#' Default is "pixelatorES".
#' @return A character string with the path to the test samplesheet CSV file.
#'
#' @export
#'
test_samplesheet <-
  function(type = c("default", "hashing"),
           package = "pixelatorES") {
    pixelatorR:::assert_vector(type, type = "character", n = 1)
    pixelatorR:::assert_single_value(package, type = "string")
    type <- match.arg(type, c("default", "hashing"))

    csv_file <- switch(type,
      "default" = "test_samplesheet.csv",
      "hashing" = "test_samplesheet_hashing.csv"
    )

    system.file(
      "extdata", csv_file,
      package = package
    )
  }

#' Get the path to the test data folder
#'
#' Returns the path to the test data folder containing QC JSON files included in the package.
#'
#' @param type A character string specifying the type of test data folder to return. Options are "default"
#' (the standard test data folder) and "hashing" (a test data folder for hashing experiments). Default is "default".
#' @param package A character string specifying the name of the package where the test data folder is located.
#' Default is "pixelatorES".
#' @return A character string with the path to the test data folder.
#'
#' @export
#'
test_data_folder <-
  function(type = c("default", "hashing"),
           package = "pixelatorES") {
    pixelatorR:::assert_vector(type, type = "character", n = 1)
    pixelatorR:::assert_single_value(package, type = "string")
    type <- match.arg(type, c("default", "hashing"))

    foldr <- switch(type,
      "default" = "qc_jsons",
      "hashing" = "qc_jsons_hashing"
    )

    system.file(
      "extdata", foldr,
      package = package
    )
  }

#' Get test data
#'
#' Generates a minimal Seurat object for testing purposes.
#'
#' @param type A character string specifying the type of test data to generate. Options are "default"
#' (the standard test data) and "hashing" (test data for hashing experiments).
#' @param package A character string specifying the name of the package where the test data is located.
#' Default is "pixelatorES".
#'
#' @return A Seurat object containing test data with normalized and scaled data, PCA results, and merged layers.
#'
#' @export
#'
get_test_data <-
  function(
    type = c("default", "hashing"),
    package = "pixelatorES"
  ) {
    pixelatorR:::assert_vector(type, type = "character", n = 1)
    type <- match.arg(type, c("default", "hashing"))

    samplesheet <-
      read_samplesheet(test_samplesheet(type = type, package = package))

    data_folder <- test_data_folder(type = type, package = package)

    data_files <- get_file_paths(data_folder, sample_sheet = samplesheet)$data_files

    # Copy PXL files to a unique tempdir to avoid duckdb
    # "blocked by another connection" errors when get_test_data() is called
    # multiple times (each call gets fresh paths so connections never collide)
    # and to ensure duckdb can write its sidecar files in a writable location.
    tmp_dir <- tempfile("pixelatorES_test_data_")
    dir.create(tmp_dir)
    data_files <- data_files %>%
      mutate(
        filename = vapply(filename, function(f) {
          dest <- file.path(tmp_dir, basename(f))
          copied <- file.copy(f, dest, overwrite = TRUE)
          if (!isTRUE(copied)) {
            stop(
              sprintf(
                "Failed to copy PXL file from '%s' to temporary location '%s'.",
                f, dest
              ),
              call. = FALSE
            )
          }
          dest
        }, character(1))
      )

    seur <-
      load_pxl_data_list(
        data_folder = tmp_dir,
        data_files = data_files,
        sample_sheet = samplesheet
      ) %>%
      merge_data(sample_sheet = samplesheet)

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

    set.seed(37)
    seur[["seurat_clusters"]] <-
      sample(1:3, ncol(seur), replace = TRUE) %>%
      factor()
    seur[["condition"]] <-
      sample(paste("Cond", 1:3), ncol(seur), replace = TRUE) %>%
      factor()

    return(seur)
  }
