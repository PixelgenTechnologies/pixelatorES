#' Create a HTML styled table
#'
#' @param df A data frame
#' @param pageLength An integer specifying the number of rows to display per page.
#' If NULL, it defaults to the number of rows in the data frame.
#' @param caption A string to be used as the table caption. If NULL, no caption is displayed.
#' @param interactive A boolean indicating whether the table should contains interactive
#' elements (default is TRUE).
#' @param search A boolean indicating whether to include a search box (default is TRUE).
#' @param lengthChange A boolean indicating whether to allow users to change the number
#' of rows displayed (default is FALSE).
#' @param buttons A boolean indicating whether to include download buttons for CSV and Excel
#' files (default is TRUE).
#' @param tooltips A boolean indicating whether to enable tooltips for the table cells
#' (default is FALSE). If enabled, it will also ensure that MathJax typesets the content
#' of the cells after the tooltip is shown.
#' @param escape A boolean indicating whether to escape HTML characters in the table cells
#' (default is FALSE).
#' @param ... Additional arguments passed to the `DT::datatable` function.
#'
#' @return A HTML table using the DT package.
#'
#' @export
#'
style_table <- function(
  df,
  pageLength = NULL,
  caption = NULL,
  interactive = TRUE,
  search = TRUE,
  lengthChange = FALSE,
  buttons = TRUE,
  tooltips = FALSE,
  escape = FALSE,
  ...) {
  pixelatorR:::assert_class(df, "data.frame")
  pixelatorR:::assert_single_value(pageLength, "integer", allow_null = TRUE)
  pixelatorR:::assert_single_value(caption, "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(interactive, "bool")
  pixelatorR:::assert_single_value(search, "bool")
  pixelatorR:::assert_single_value(lengthChange, "bool")
  pixelatorR:::assert_single_value(buttons, "bool")

  # Options
  opts <- list(pageLength = ifelse(is.null(pageLength), nrow(df), pageLength))

  dom <- ""

  if (lengthChange) {
    dom <- paste0(dom, "l")
  }
  if (buttons) {
    dom <- paste0(dom, "B")
    # Customize the button to make it smaller and apply a specific class
    opts$buttons <- list(
      list(
        extend = "collection",
        buttons = c("csv", "excel"),
        text = "Download"
      )
    )
  }
  if (search) {
    dom <- paste0(dom, "f") # Add search box
  }

  dom <- paste0(dom, "rt")

  if (opts$pageLength < nrow(df)) {
    dom <- paste0(dom, "ip") # Add pagination if pageLength is less than number of rows
  }

  if (!interactive) {
    if (buttons) {
      dom <- "Bt" # Only show the table and buttons, nothing else
    } else {
      dom <- "t" # Only show the table, nothing else
    }
  }

  if(tooltips) {
    opts$initComplete <-
      htmlwidgets::JS("
      $(function () {
        $('[data-bs-toggle=\"tooltip\"]').tooltip();
        document.addEventListener('shown.bs.tooltip', function (e) {
          if (window.MathJax) {
            MathJax.typesetPromise([e.target.nextSibling]);
          }
        });
      });
    ")

    opts$drawCallback <-
      htmlwidgets::JS("
    function(settings) {
      // Reapply tooltips on each redraw
      $('[data-bs-toggle=\"tooltip\"]').tooltip();
    }
  ")

  }

  opts$dom <- dom

  table_widget <- DT::datatable(
    df,
    caption = caption,
    rownames = FALSE,
    escape = escape,
    options = opts,
    extensions = if (buttons) "Buttons" else list(),
    ...
  )

  return(table_widget)
}
