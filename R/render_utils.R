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
      image_write(img, webp_path, format = "webp", quality = 100)

      # Remove the original PNG file and replace the filepath with the WebP path
      file.remove(filepath)
      filepath <- webp_path
    }
    # Return the path to the WebP file for knitr to embed
    hook_plot_md(filepath, options)
  }
