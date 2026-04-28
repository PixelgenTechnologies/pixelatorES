#' Convert PNG images to WebP format
#'
#' This function is used as a hook in knitr to convert PNG images generated during the rendering of a Quarto document
#' to WebP format. It reads the PNG file, converts it to WebP with specified quality, and then removes the original
#' PNG file.
#'
#' @param filepath The path to the PNG file that knitr has just written.
#' @param options The options passed to the knitr hook.
#'
#' @return The path to the WebP file, which knitr will use for embedding in the document.
#'
#' @export
#'
convert_png_to_webp <-
  function(filepath, options) {
    if (str_detect(filepath, "\\.png$")) {
      # Make a WebP file:
      webp_path <- str_replace(filepath, "\\.png$", ".webp")
      img <- image_read(filepath)
      image_write(img, webp_path, format = "webp", compression = "WebP")

      # Remove the original PNG file and replace the filepath with the WebP path
      file.remove(filepath)
      filepath <- webp_path
    }
    # Return the path to the WebP file for knitr to embed
    hook_plot_md(filepath, options)
  }

#' Render a section with a table
#'
#' This function is a helper for rendering a section with a
#' table in the ES report. It takes a datatables object, a title
#' for the section, and an optional level for the heading (default is 3).
#' It then prints the section heading and the table in the appropriate
#' format for the ES report.
#'
#' @param table A datatables object to be rendered in the section.
#' @param title A string representing the title of the section.
#' @param level An integer representing the heading level for the section
#' title (default is 3).
#'
#' @export
#'
section_table <- function(table, title, level = 3) {
  pixelatorR:::assert_class(table, "datatables")
  pixelatorR:::assert_single_value(title, "string")
  pixelatorR:::assert_single_value(level, "integer")
  cat(paste0("\n\n", strrep("#", level), " ", title, "\n\n"))
  print(htmltools::tagList(table))
  cat("\n\n")
}

#' Render a section with introductory text
#'
#' This function is a helper for rendering a section with introductory text
#' in the ES report. It takes a title for the section, an optional level for
#' the heading (default is 3), and the text to be displayed. It then
#' prints the section heading and the introductory text in the appropriate
#' format for the ES report.
#'
#' @param title A string representing the title of the section.
#' @param text A string representing the introductory text to be displayed
#' @param level An integer representing the heading level for the section
#' title (default is 3).
#'
#' @export
#'
section_intro <- function(title, text, level = 3) {
  pixelatorR:::assert_single_value(title, "string")
  pixelatorR:::assert_single_value(text, "string")
  pixelatorR:::assert_single_value(level, "integer")
  cat(paste0("\n\n", strrep("#", level), " ", title, "\n\n"))
  cat(text)
  cat("\n\n")
}
