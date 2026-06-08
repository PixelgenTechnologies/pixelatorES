# Developer Guide

This document contains instructions for developers working on the Proxiome Experiment Summary (`pixelatorES`) package.

## Table of Contents

- [Code Style](#code-style)
- [Linting](#linting)
- [Styler](#styler)
- [Type assertions and error messages](#type-assertions-and-error-messages)
- [Package dev tasks](#package-dev-tasks)
- [Quarto report rendering](#quarto-report-rendering)
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

## Quarto report rendering

The ES report is built with Quarto (see `inst/quarto/`). Components return `ggplot` objects and/or `DT::datatables` objects from [`style_table()`](R/tables.R). To place them in the HTML report, `.qmd` chunks call helper functions that emit Quarto markdown via `cat()`.

**Requirement:** any chunk that uses these helpers must set `#| results: 'asis'` so printed output is passed through as raw markdown/HTML rather than wrapped in a code block.

### How tabsets work

Tab UI is created with Quarto fenced divs:

```markdown
::: {.panel-tabset .nav-pills}
### Tab one
(content)
### Tab two
(content)
:::
```

The `tabset_*` helpers write these divs for you. Tab titles come from markdown headings (`#` level set by the `level` argument). Use `level` consistently within a section so tab nesting matches the surrounding heading hierarchy (e.g. `level = 5` under a `####` markdown heading).

### `tabset_plotlist()`

Renders a **named list of plots** as a tabset — one tab per plot.

```r
#| results: 'asis'
plots <- component_sequencing_saturation_curve(pg_data, data_files, sample_levels = sample_aliases)
tabset_plotlist(plots, level = 5)
```

- **Input:** named `list` of `ggplot` or `datatables` objects.
- **`close`:** defaults to `TRUE`; the function opens and closes the tabset div.
- **Use when:** a component returns multiple plots with no paired summary table, or tables are rendered separately.

See `inst/quarto/quality_metrics.qmd` (sequencing saturation curves) and `inst/quarto/abundance.qmd`.

### `tabset_nested_plotlist()`

Like `tabset_plotlist()`, but each list element may itself be a `list` of plots, producing nested tabsets.

```r
#| results: 'asis'
tabset_figure_table(
  component$plots,   # named list of ggplot objects
  component$table,
  level = 5,
  mode = "tabset_nested"
)
close_tabset()
```

- **Use when:** one table accompanies several plots that should be grouped in sub-tabs (e.g. cell recovery molecule-rank plots).

See `inst/quarto/quality_metrics.qmd` (`qc_metrics_molrank_plot`) and `inst/quarto/spatial.qmd`.

### `tabset_figure_table()`

Renders a **figure + table** pair as a two-tab set ("Figure" and "Table").

```r
#| results: 'asis'
component <- component_control_markers(pg_data)

tabset_figure_table(
  component$p1,
  component$tabl,
  level = 5,
  mode = "title"
)
close_tabset()
```

| Argument | Description |
|---|---|
| `figure` | A single `ggplot`, or a `list` of ggplots when `mode` is `"tabset"` or `"tabset_nested"`. |
| `table` | A `datatables` object from `style_table()`. |
| `level` | Markdown heading level for the "Figure" / "Table" tab labels. |
| `mode` | `"tabset"` — sub-tab per plot in a list; `"tabset_nested"` — nested sub-tabs; `"title"` — print a single figure directly (no sub-tabset). |

**Important:** `tabset_figure_table()` opens a `::: {.panel-tabset}` div but does **not** close it. Always call `close_tabset()` immediately after to emit the closing `:::`.

The function **returns** the `table` object; in an `results: 'asis'` chunk the return value is printed automatically, which renders the table in the "Table" tab.

**Use when:** a component produces one primary plot (or a small plot list) with a summary table — the most common pattern for QC components.

See `inst/quarto/quality_metrics.qmd` (control markers, k-coreness, reads per cell).

### `close_tabset()`

Emits the closing `:::` for a tabset opened by `tabset_figure_table()`. Pair it one-to-one with each `tabset_figure_table()` call.

```r
tabset_figure_table(component$plot, component$table, level = 5)
close_tabset()
```

For multiple independent figure-table pairs in one section, repeat the pair:

```r
#| results: 'asis'
component <- component_denoising(qc_metrics_tables, sample_levels = sample_aliases)

tabset_figure_table(component$plots$removed_umis, component$tables$removed_umis, level = 5, mode = "title")
close_tabset()

tabset_figure_table(component$plots$by_method, component$tables$by_method, level = 5, mode = "title")
close_tabset()

tabset_figure_table(component$plots$isotype_reduction, component$tables$isotype_reduction, level = 5, mode = "title")
close_tabset()
```

Separate markdown `####` headings above each chunk can label the metric groups when using this pattern.

### `title_plotlist()`

Prints a named list of plots under markdown headings **without** creating a tabset. Used internally by `tabset_plotlist()` and as the figure renderer when `tabset_figure_table(..., mode = "title")` receives a plot list.

Rarely called directly from `.qmd` files; prefer `tabset_plotlist()` unless you are building custom tabset markup.

### `section_table()` and `section_intro()`

Lower-level helpers for non-tabbed sections ([`R/render_utils.R`](R/render_utils.R)).

- **`section_intro(title, text, level)`** — prints a heading and body text.
- **`section_table(table, title, level)`** — prints a heading and a `datatables` object.

Often used inside a manually opened tabset:

```r
#| results: 'asis'
cat("::: {.panel-tabset .nav-pills}\n")
section_table(key_tables$pool, "Hash pool metrics", 3)
section_table(key_tables$sample, "Sample metrics", 3)
cat(":::\n")
```

See `inst/quarto/quality_metrics.qmd` (key metrics).

### `style_table()`

Formats a data frame as a non-interactive (or interactive) DT table for report embedding. Components should return tables created with `style_table(caption = "...", interactive = FALSE)` so they render consistently in `tabset_figure_table()` and `section_table()`.

### Adding a new component to the report

1. Implement `component_*()` in [`R/components.R`](R/components.R) returning `ggplot` objects and/or `style_table()` output.
2. In the appropriate `inst/quarto/*.qmd` file, add a chunk with `results: 'asis'`.
3. Choose a layout helper:
   - **One plot + one table** → `tabset_figure_table()` + `close_tabset()`
   - **Several plots, no table** → `tabset_plotlist()`
   - **Several plots + one table** → `tabset_figure_table(..., mode = "tabset_nested")` + `close_tabset()`
   - **Table or intro only** → `section_table()` / `section_intro()`
4. Guard the chunk with `#| eval: !expr ...` when data may be absent (e.g. `!is.null(qc_metrics_tables$denoising)`).

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
