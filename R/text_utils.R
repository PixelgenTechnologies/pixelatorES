#' Order markers in a consistent way
#'
#' This function orders markers in a consistent way, with markers named "CD.*" first, non-CD markers second
#' in alphabetic order, and control markers last.
#'
#' @param markers A character vector of marker names.
#' @param control_markers A character vector of control marker names.
#'
#' @return A character vector of markers ordered as described above.
#'
#' @export
order_cd_markers <-
  function(markers, control_markers) {
    pixelatorR:::assert_vector(markers, "character", n = 1)
    pixelatorR:::assert_vector(control_markers, "character", n = 1)

    cd_markers <- str_detect(markers, "^CD\\d")
    control_markers <- markers %in% control_markers


    cd_order <-
      tibble(marker = markers[cd_markers]) %>%
      mutate(marker_i = str_remove(marker, "^CD") %>%
        str_extract("^\\d*") %>%
        as.numeric()) %>%
      arrange(marker_i, marker) %>%
      pull(marker)

    non_cd_order <-
      tibble(marker = markers[!cd_markers & !control_markers]) %>%
      arrange(marker) %>%
      pull(marker)

    control_order <-
      tibble(marker = markers[control_markers]) %>%
      arrange(marker) %>%
      pull(marker)

    return(c(cd_order, non_cd_order, control_order))
  }


#' Order sample alias factors in a data frame
#'
#' This function orders the `sample_alias` factor in a data frame according to the specified levels.
#'
#' @param object A data frame containing a `sample_alias` column.
#' @param levels A character vector specifying the levels to order the `sample_alias` factor.
#'
#' @return A data frame with the `sample_alias` factor ordered according to the specified levels.
#'
#' @export
#'
order_sample_alias_factors <-
  function(object, levels) {
    if (inherits(object, "data.frame")) {
      ordered_object <-
        object %>%
        mutate(sample_alias = factor(sample_alias, levels))
    } else if (inherits(object, "list")) {
      ordered_object <-
        object %>%
        lapply(order_sample_alias_factors, levels = levels)
    }

    return(ordered_object)
  }


#' Beautify axis labels
#'
#' This function beautifies axis labels by removing underscores and replacing them with spaces.
#'
#' @param labels A character vector of axis labels.
#'
#' @return A character vector of beautified axis labels.
#'
#' @export
#'
beaut_label <-
  function(labels) {
    pixelatorR:::assert_vector(labels, "character", n = 1)

    labels %>%
      str_replace_all("_", " ") %>%
      toupper()
  }

#' Order reductions in a preferred way
#'
#' This function orders dimensionality reductions in a preferred way, with
#' "harmony_umap" first, followed by "umap", and then "pca".
#'
#' @param reductions A character vector of reduction names.
#' @param keep_undefined A logical value indicating whether to keep reductions
#' that are not in the preferred order.
#'
#' @return A character vector of reductions ordered according to the preferred order.
#'
#' @export
#'
preferred_dimred_order <-
  function(
    reductions,
    keep_undefined = FALSE
  ) {
    pixelatorR:::assert_vector(reductions, "character", n = 1)
    pixelatorR:::assert_class(keep_undefined, "logical")

    preferred_order <-
      c("harmony_umap", "umap", "pca")

    reductions_ordered <-
      preferred_order[preferred_order %in% reductions]

    if (keep_undefined) {
      reductions_ordered <-
        c(
          reductions_ordered,
          reductions[!reductions %in% preferred_order]
        )
    }

    return(reductions_ordered)
  }

#' Compact number formatting
#'
#' This function formats a numeric value into a compact representation, either in thousands (k)
#' or millions (M), with a specified number of decimal places.
#'
#' @param n A numeric value to be formatted.
#' @param accuracy A numeric value specifying the accuracy of the formatting.
#'
#' @return A character string representing the formatted number with the specified size unit.
#'
#' @export
#'
compact_num <-
  function(
    n,
    accuracy = 0.1
    ) {
    pixelatorR:::assert_vector(n, "numeric", n = 1)
    pixelatorR:::assert_single_value(accuracy, "numeric")

    format_func <-
      scales::label_number(
        accuracy = accuracy,
        scale_cut = c(0, "k" = 1000, "M" = 1e+06, "B" = 1e+09, "T" = 1e+12)
      )

    return(format_func(n))
  }

