#' Runs a ridiculously minimal test of the pixelatorESinternal
#'
#' This function is intended for development to quickly test if
#' new features break. It uses a combinations of json files from
#' `pixelatorESinternal` and PXL files from `pixelatorR`. The PXL
#' files only contain 5 cells and the data is therefore not useful
#' to evaluate most plots.
#'
#' It is possible to test with a larger dataset by providing a custom
#' `params` list where the `sample_sheet` and `data_folder` are changed.
#' Note that you should also set `use_tiny_test = FALSE` to use more
#' sutiable defaults for the test in this case.
#'
#' @section Quarto path:
#' It's important to make sure that quarto is isntalled and that `quarto_render`
#' uses the correct quarto path. The easiest way to set this is by using:
#' `Sys.setenv(QUARTO_PATH = "<path to quarto binary>")`
#'
#' @param output_file An optional specifying the output file where the
#' rendered document will be saved.
#' @param use_tiny_test A boolean indicating whether to use a tiny test dataset
#' with customized parameters (default is TRUE).
#' @param params A list of parameters to override the default settings.
#' @param view_ES A boolean indicating whether to open the rendered HTML file
#' @param ... Additional arguments passed to the `quarto_render` function.
#'
#' @return Nothing. Only used for its side effects
#'
#' @export
#'
test_ES <- function(
  output_file = NULL,
  use_tiny_test = TRUE,
  params = NULL,
  view_ES = FALSE,
  ...
) {
  expect_quarto()
  expect_DBI()
  expect_duckdb()
  expect_fs()

  pixelatorR:::assert_single_value(output_file, "string", allow_null = TRUE)
  pixelatorR:::assert_file_ext(output_file, "html", allow_null = TRUE)
  pixelatorR:::assert_single_value(use_tiny_test, "bool")
  pixelatorR:::assert_class(params, "list", allow_null = TRUE)

  if (view_ES && is.null(output_file)) {
    cli::cli_abort(
      c("x" = "To preview the ES, you must provide an `output_file` path.")
    )
  }

  qmd_path <- system.file("quarto/pixelatorESinternal.qmd", package = "pixelatorESinternal")
  samplesheet_path <- system.file("extdata/test_samplesheet.csv", package = "pixelatorESinternal")
  data_path <- system.file("extdata/qc_jsons", package = "pixelatorESinternal")

  # Move data to temp
  tmp_dir <- fs::file_temp() %>% stringr::str_replace_all("file", "dir")
  fs::dir_create(tmp_dir)
  qmd_temp_path <- file.path(tmp_dir, "pixelatorESinternal.qmd")
  fs::file_copy(qmd_path, qmd_temp_path)
  fs::dir_copy(data_path, tmp_dir, overwrite = TRUE)
  pxl1 <- file.path(tmp_dir, "S01_PBMC_unstimulated_S1.layout.pxl")
  pxl2 <- file.path(tmp_dir, "S02_PHA_S2.layout.pxl")
  pxl_files <- c(pxl1, pxl2)
  fs::file_copy(minimal_pna_pxl_file(), pxl1, overwrite = TRUE)
  fs::file_copy(minimal_pna_pxl_file(), pxl2, overwrite = TRUE)

  # We need to update the edgelist for compatibility with the
  # seuquencing saturation calculation
  # We also need to add annotations to the obs data
  ann <- data.frame(
    index = c("0a45497c6bfbfb22", "2708240b908e2eba", "c3c393e9a17c1981", "d4074c845bb62800", "efe0ed189cb499fc"),
    l1_annotation = rep(c("Mono", "DC", "CD4 T", "CD4 T", "CD4 T")),
    l2_annotation = rep(c("CD16 Mono", "pDC", "CD4 T", "CD4 T", "CD4 T"))
  ) %>%
    mutate(
      l1_annotation_summary = l1_annotation,
      l2_annotation_summary = l2_annotation
    )
  join_query <- "
    CREATE OR REPLACE TABLE __adata__obs AS
    SELECT
      o.*,
      a.*
    FROM
      __adata__obs AS o
    LEFT JOIN
      ann AS a ON o.index = a.index;
  "
  for (f in pxl_files) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = f, read_only = FALSE)
    DBI::dbWriteTable(con, "ann", ann, temporary = TRUE)
    DBI::dbExecute(con, join_query)
    DBI::dbExecute(con, "ALTER TABLE edgelist ALTER COLUMN read_count TYPE SMALLINT")
    DBI::dbDisconnect(con)
  }

  if (is.null(params)) {
    params <- list(
      sample_sheet = samplesheet_path,
      data_folder = tmp_dir,
      test_mode = TRUE
    )
  }

  if (!is.null(params) && use_tiny_test) {
    params[["sample_sheet"]] <- samplesheet_path
    params[["data_folder"]] <- tmp_dir
    params[["test_mode"]] <- TRUE
  }

  quarto::quarto_render(
    input = qmd_temp_path,
    execute_params = params
  )

  html_path <- file.path(tmp_dir, "pixelatorESinternal.html")

  if (!is.null(output_file)) {
    fs::file_copy(html_path, output_file, overwrite = TRUE)
  }

  if (view_ES) {
    expect_rstudioapi()
    rstudioapi::viewer(output_file)
  }

  return(invisible(NULL))
}
