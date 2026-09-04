# Developer Guide

This document contains instructions for developers working on the Proxiome Experiment Summary (`pixelatorES`) package.

## Table of Contents

- [Code Style](#code-style)
- [Linting](#linting)
- [Styler](#styler)
- [Type assertions and error messages](#type-assertions-and-error-messages)
- [Package dev tasks](#package-dev-tasks)
- [The `es_data` ingestion system](#the-es_data-ingestion-system)
- [Quarto report rendering](#quarto-report-rendering)
- [Violin plots with points](#violin-plots-with-points)
- [Updating test data](#updating-test-data)

---

## Code Style

This project uses `lintr` for linting and `styler` for code formatting. The configuration is compatible with `pixelatorR` style conventions.

---

## Linting

To run the linter, you need to install `lintr`. Then you can use one of the following commands:

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

To style the code, you need to install `styler`. Then you can use one of the following commands:

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

## The `es_data` ingestion system

All Experiment Summary preprocessing flows through a single `es_data` object. The report calls `build_es_data(params)` once, and every child `.qmd` and every `component_*()` function reads what it needs from that object. The code lives in [`R/es_data.R`](R/es_data.R) (the object and its extractor implementations), [`R/workflow_registry.R`](R/workflow_registry.R) (the registry), and [`R/workflow_amplicon_demux.R`](R/workflow_amplicon_demux.R) (the built-in workflow's three factories).

### The `es_data` object

`es_data` is a plain `list` with the S3 class `c("es_data", "list")` and a fixed set of slots:

| Slot | Contents |
|---|---|
| `params` | The report parameters (`knit_param_list` or `list`), including `workflow`. |
| `diagnostics` | List of non-fatal problems recorded during the build. |
| `samplesheet` | The parsed experiment samplesheet. |
| `sample_aliases` | Named character vector mapping sample (and pool) ids to aliases; defines sample ordering across the report. |
| `effective_samplesheet` | The samplesheet reduced to samples that actually loaded. |
| `file_paths` | Discovered input file paths, found once and reused. |
| `pxl_data` | Merged raw PXL data (released after processing to save memory). |
| `qc_raw` | Raw QC data used to derive formatted metrics. |
| `qc` | Nested list of formatted QC tables (`qc$read_stats`, `qc$crossing_edges`, ...). |
| `pxl_data_processed` | The processed Seurat object (normalized, clustered, annotated). |
| `proximity` | Proximity scores. |
| `stages` | The workflow's pipeline stage vocabulary (`all`, `pool`, `pxl_preference`), used for input file discovery. |

`new_es_data(params)` builds the empty shell and attaches the workflow's extractors and stage vocabulary; the data slots start `NULL`/empty and are filled during the build.

### Extractors and the build

The object is populated by **extractors** — functions that each take the current `es_data` and return the value for one slot. They are held in a nested named list where top-level names map to slots and names nested under `qc` map to `es_data$qc$...`:

```r
list(
  samplesheet = .extract_samplesheet,
  pxl_data    = .extract_pxl_data,
  qc = list(
    read_stats     = .extract_read_stats,
    crossing_edges = .extract_crossing_edges,
    # ...
  ),
  pxl_data_processed = .extract_pxl_data_processed,
  proximity          = .extract_proximity
)
```

`build_es_data()` orchestrates the build:

1. `new_es_data(params)` creates the shell and looks up the workflow extractors.
2. The samplesheet is read first — this is the **only** hard requirement, and a failure here stops the build.
3. `sample_aliases` is derived from the samplesheet with `.sample_aliases_from_samplesheet()`, and input files are discovered once with `get_file_paths()`, using the object's stage vocabulary.
4. `run_es_data_extractors()` walks the extractor tree depth-first and fills the remaining slots. Once `pxl_data_processed` is set, the raw `pxl_data` slot is cleared to release memory.

### Soft-fail diagnostics

Ingestion is resilient: a broken sample or missing QC file should not lose the whole report. Two mechanisms support this:

- An extractor can return an `es_data_extractor_result` (via the internal `.new_es_data_extractor_result(value, diagnostics)`) to hand back a **partial** value together with diagnostics — for example, PXL loading returns the samples that loaded plus `pxl_load` diagnostics for those that did not.
- If an extractor throws instead, `run_es_data_extractors()` catches the error, records an `extractor` diagnostic, and leaves that slot `NULL` while continuing with the rest.
- After QC groups load, relative stage completeness is checked within samples and within pools. If peers have a stage that another loaded entity lacks, a `qc_load` diagnostic is recorded for that entity while keeping its partial QC data.

Each diagnostic (`.new_es_data_diagnostic()`) has a `type` (`"pxl_load"`, `"qc_load"`, or `"extractor"`), a `target` (a sample alias, pool id, or slot path such as `"qc$crossing_edges"`), and a human-readable `message`. Inspect them via `es_data$diagnostics`.

On the report, diagnostics surface in two places:

- **Samples**: a red callout lists sample/pool loading issues and analysis-step failures. Sample- or pool-targeted loading issues also add a warning marker in the metadata table.
- **Run info**: the Diagnostics section lists every recorded diagnostic when any exist.

### Registering a workflow

Workflows are stored in a package-local registry ([`R/workflow_registry.R`](R/workflow_registry.R)). `params$workflow` selects one and defaults to `"amplicon_demux"`, which is registered from `.onLoad()` in [`R/zzz.R`](R/zzz.R) with the factories from [`R/workflow_amplicon_demux.R`](R/workflow_amplicon_demux.R). Keep one `R/workflow_<id>.R` file per workflow; R does not allow subdirectories under `R/`. Each workflow registers:

- `extractors`: a zero-argument factory returning the nested extractor list
- `report`: a zero-argument factory returning the Quarto report recipe (`preamble` + `sections`)
- `stages`: a zero-argument factory returning the pipeline stage vocabulary (`all` + `pool` + `pxl_preference`)

Extension packages should call `register_es_data_workflow()` from their `.onLoad()` hook:

```r
.register_child <- function(...) {
  system.file("quarto", ..., package = "myPackage", mustWork = TRUE)
}

register_es_data_workflow(
  name = "my_workflow",
  extractors = function() {
    list(
      samplesheet = my_extract_samplesheet,
      pxl_data = my_extract_pxl_data
      # ...
    )
  },
  report = function() {
    list(
      preamble = .register_child("shared", "preprocessing.qmd"),
      sections = list(
        list(
          id = "samples",
          title = "Samples",
          child = .register_child("shared", "samples.qmd")
        ),
        list(
          id = "quality_metrics",
          title = "Quality metrics",
          child = .register_child("workflows", "my_workflow", "quality_metrics.qmd")
        )
      )
    )
  },
  stages = function() {
    list(
      all = c("amplicon", "demux", "graph", "analysis"),
      pool = c("amplicon", "demux"),
      pxl_preference = c("graph", "analysis")
    )
  }
)
```

Built-in workflows use paths relative to `inst/quarto/`. Extension packages should register absolute paths from `system.file()`. All referenced child paths are checked for existence at registration time.

The stage vocabulary is the source of truth for file discovery: `all` lists every stage the pipeline can emit, `pool` marks the stages whose QC files are pool-level, and `pxl_preference` orders the stages when several produce a PXL file for the same sample. `pool` may be empty, but both it and `pxl_preference` must be subsets of `all`.

`find_stage()`, `get_file_paths()`, and `extract_sample_qc_metrics()` all take a **required** `stages` argument — there is no package-level default, so no caller can silently inherit another workflow's vocabulary. Inside a build, pass `es_data$stages` (`find_stage()` takes just `es_data$stages$all`); elsewhere, ask the registry with `get_es_workflow_stages(workflow)`.

Use `list_es_data_workflows()` to see what is registered, `get_es_workflow_report(name)` to inspect a report recipe, and `get_es_workflow_stages(name)` to inspect a stage vocabulary.

### Consuming `es_data`

Report code should treat `es_data` as the single source of truth: `component_*()` functions, `key_metric_table()`, and the `print_*()` helpers all take `es_data` as their first argument and pull the slots they need internally. When adding a new component, accept `es_data` and read from its slots rather than threading individual objects through the `.qmd` files.

### Test fixtures with `test_es_data()`

For unit tests of components and helpers, build lightweight `es_data` objects with [`test_es_data()`](R/es_data.R) instead of hand-rolling `structure(..., class = c("es_data", "list"))`. The helper wraps `new_es_data()` so the class and slot layout stay aligned with the real constructor, then overwrites the slots you pass:

```r
es <- test_es_data(
  samplesheet = sample_sheet,
  qc = qc_metrics_tables,
  pxl_data_processed = pg_data
)
```

Use this for partial or synthetic fixtures. Prefer `build_es_data(params)` when the test needs the full ingestion pipeline.

---

## Quarto report rendering

The ES report is built with Quarto (see `inst/quarto/`). `pixelatorES.qmd` is a thin dispatcher: it loads the registered report recipe for `params$workflow`, knits the preamble children, then emits the panel tabset from the recipe sections.

### Quarto layout

```text
inst/quarto/
  pixelatorES.qmd                 # dispatcher shell
  shared/                         # reused across workflows
    preprocessing.qmd
    samples.qmd
    run_info.qmd
  workflows/
    amplicon_demux/               # workflow-owned sections
      quality_metrics.qmd
      cell_annotation.qmd
      abundance.qmd
      spatial.qmd
```

Components return `ggplot` objects and/or `DT::datatables` objects from [`style_table()`](R/tables.R). To place them in the HTML report, `.qmd` chunks call helper functions that emit Quarto markdown via `cat()`.

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

Report setup calls `register_tabset_chunk_hooks()`. If a knitr chunk errors after `open_tabset()` or a `tabset_*` helper has written an opening fence, the hook appends the missing `:::` so later sections still parse. In `.qmd` chunks, call `open_tabset()` instead of writing the opening fence with `cat()`. Markdown tabsets written outside R chunks are unchanged.

### `tabset_plotlist()`

Renders a **named list of plots** as a tabset — one tab per plot.

```r
#| results: 'asis'
plots <- component_sequencing_saturation_curve(es_data)
tabset_plotlist(plots, level = 5)
```

- **Input:** named `list` of `ggplot` or `datatables` objects.
- **`close`:** defaults to `TRUE`; the function opens and closes the tabset div.
- **Use when:** a component returns multiple plots with no paired summary table, or tables are rendered separately.

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

### `tabset_figure_table()`

Renders a **figure + table** pair as a two-tab set ("Figure" and "Table").

```r
#| results: 'asis'
component <- component_control_markers(es_data)

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
open_tabset()
section_table(key_tables$pool, "Hash pool metrics", 3)
section_table(key_tables$sample, "Sample metrics", 3)
close_tabset()
```

### `style_table()`

Formats a data frame as a non-interactive (or interactive) DT table for report embedding. Components should return tables created with `style_table(caption = "...", interactive = FALSE)` so they render consistently in `tabset_figure_table()` and `section_table()`.

### Adding a new component to the report

1. Implement `component_*()` in [`R/components.R`](R/components.R) returning `ggplot` objects and/or `style_table()` output.
2. In the appropriate section `.qmd` under `inst/quarto/shared/` or `inst/quarto/workflows/<id>/`, add a chunk with `results: 'asis'`.
3. Choose a layout helper:
   - **One plot + one table** → `tabset_figure_table()` + `close_tabset()`
   - **Several plots, no table** → `tabset_plotlist()`
   - **Several plots + one table** → `tabset_figure_table(..., mode = "tabset_nested")` + `close_tabset()`
   - **Table or intro only** → `section_table()` / `section_intro()`
4. Guard the chunk with `#| eval: !expr ...` when data may be absent (e.g. `!is.null(es_data$qc$denoising)`).
5. If the section is new for a workflow, add it to that workflow's registered report recipe.

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
