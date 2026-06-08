# Developer Guide

This document contains instructions for developers working on the Proxiome Experiment Summary (`pixelatorES`) package.

## Table of Contents

- [Code Style](#code-style)
- [Linting](#linting)
- [Styler](#styler)
- [Type assertions and error messages](#type-assertions-and-error-messages)
- [Package dev tasks](#package-dev-tasks)
- [Updating test data](#updating-test-data)

---

## Code Style

This project uses `lintr` for linting and `styler` for code formatting. The configuration is compatible with `pixelatorR` style conventions.

---

## Linting

To run the linter, you need to install `lintr`. Then you can use one of the the following commands:

```r
# Lint entire package
lintr::lint_package()

# Lint single file
lintr::lint("path/to/file.R")
```

Alternatively, you can run the linter from RStudio through Addins -> Lint current file or Addins -> Lint current package.

The configuration file `.lintr` is used to specify the rules that the linter should follow. For compatibility with styler, some linting rules have been disabled.

CI runs `lintr::lint_package()` with `LINTR_ERROR_ON_LINT: true` on pull requests (see `.github/workflows/lint.yaml`).

---

## Styler

To style the code, you need to install `styler`. Then you can use one of the the following commands:

```r
# Style entire package
styler::style_pkg(transformers = pixelatorR::pixelatorR_style())

# Style single file
styler::style_file("R/components.R", transformers = pixelatorR::pixelatorR_style())
```

Alternatively, you can run the styler from RStudio. Configure `styler` to use the style guide provided in `pixelatorR`. Go to Addins -> Styler -> Set style and set `pixelatorR::pixelatorR_style()` as the style guide. Then you can use Addins -> Styler -> Style active file or Addins -> Styler -> Style active package.

---

## Type assertions and error messages

`pixelatorR` provides a set of type assertions functions to check function arguments. See `?type-checks` for a full list.

These functions are designed to cover the most common type checks that are needed in `pixelatorR` to ensure consistent error messages. They are also designed according to the tidyverse style guide (see [tidyverse style guide](https://style.tidyverse.org/errors.html) and [Including function calls in error messages](https://rlang.r-lib.org/reference/topic-error-call.html)).

### Example

Assert that the input is a single character string:

```
my_var <- 1
assert_single_value(my_var, type = "string")
```

Output:

```
Error:
ℹ `my_var` must be a single string.
✖ You've supplied a <numeric>.
```

Note that the error message includes the name of the variable, i.e. we don't need to hard code the name in the error message.

Another useful feature is that the assertion function can pass the caller environment to the error message to signal where the error happened. This is useful when the assertion is used inside a function.

```
my_var <- 1
my_func <- function(x) {
  assert_single_value(x, type = "string")
}
my_func(my_var)
```

Output:

```
Error in `my_func()`:
ℹ `x` must be a single string.
✖ You've supplied a <numeric>.
```

### Writing error messages

The assertion functions do not cover all possible type checks and are not necessarily useful for error messages that require more verbosity. If you need to write a custom error message, we suggest using `cli::cli_abort()`.

```
my_func <- function(x) {
  if (!is.finite(x)) {
    cli::cli_abort(
      c("i" = "{.arg x} must be a finite number.",
        "x" = "x = {.val {x}}")
    )
  }
}
my_func(Inf)
```

Output:

```
Error in `my_func()`:
ℹ `x` must be a finite number.
✖ x = Inf
```

---

## Package dev tasks

Generate roxygen2 documentation:

```r
devtools::document()
```

Run all package tests:

```r
devtools::test()
```

Run a subset of tests:

```r
devtools::test(filter = "component|key_table")
```

---

## Updating test data

Test fixtures live in `inst/extdata/`:

- `qc_jsons/` and `qc_jsons_hashing/` — QC JSON files and layout `.pxl` files used by `get_test_data()`, `test_data_folder()`, and component tests
- `test_samplesheet.csv` and `test_samplesheet_hashing.csv` — sample sheets for tests

`.pxl` files are tracked with **Git LFS**. Ensure Git LFS is configured before committing modified `.pxl` files.

### Adding obs metadata columns to PXL files

PXL files store component metadata in the DuckDB table `__adata__obs`. To add or update columns, copy the file to a temporary location, modify via DuckDB, then copy back:

```r
library(DBI)
library(duckdb)
library(fs)
library(dplyr)

pxl_path <- "inst/extdata/qc_jsons/S01_PBMC_unstimulated_S1.layout.pxl"

tmp_pxl <- fs::file_temp(ext = ".pxl")
fs::file_copy(pxl_path, tmp_pxl, overwrite = TRUE)

con <- DBI::dbConnect(
  duckdb::duckdb(),
  bigint = "integer64",
  dbdir = tmp_pxl,
  read_only = FALSE
)

DBI::dbGetQuery(con, "SHOW ALL TABLES")
obs <- DBI::dbGetQuery(con, "SELECT * FROM __adata__obs")

# Example: add new columns
obs <- obs %>%
  mutate(
    a_new_component_metric = pmax(0L, round(n_umi * 0.003)),
    another_component_metric = pmax(0L, round(n_umi * 0.002)),
  )

DBI::dbWriteTable(con, "__adata__obs", obs, overwrite = TRUE)
DBI::dbDisconnect(con, shutdown = TRUE)

fs::file_copy(tmp_pxl, pxl_path, overwrite = TRUE)

# Verify
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = pxl_path, read_only = TRUE)
head(DBI::dbGetQuery(con, "SELECT * FROM __adata__obs"))
DBI::dbDisconnect(con, shutdown = TRUE)
```

Commit the updated LFS-tracked `.pxl` files together with any related test expectations.
