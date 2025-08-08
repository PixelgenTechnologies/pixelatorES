#' Create sample palette
#'
#' Create a sample paette based on the number of unique conditions in the experiment.
#'
#' @param conditions Character vector of conditions.
#'
#' @return Character vectors of color hex codes.
#'
#' @export
#'
create_sample_palette <-
  function(conditions) {
    pixelatorR:::assert_vector(conditions, "character", n = 1)

    if (n_distinct(conditions) <= 2) {
      sample_palette <-
        c(Pixelgen_accent_colors$blues[5], Pixelgen_accent_colors$greys[4]) %>%
        rep(10)
    } else if (n_distinct(conditions) <= 6) {
      sample_palette <-
        c(
          PixelgenAccentColors(hue = c("greys"), level = c(4)),
          PixelgenAccentColors(hue = c("blues"), level = c(4, 6)),
          PixelgenAccentColors(hue = c("oranges"), level = c(4, 6)),
          PixelgenAccentColors(hue = c("reds"), level = c(6))
        ) %>%
        unname() %>%
        rep(10)
    } else if (n_distinct(conditions) > 6) {
      sample_palette <-
        c(
          PixelgenAccentColors(hue = c("greys", "blues", "cyans"), level = c(4, 6)),
          PixelgenAccentColors(hue = c("pinks"), level = c(5)),
          PixelgenAccentColors(hue = c("reds"), level = c(6)),
          PixelgenAccentColors(hue = c("oranges"), level = c(4, 6)),
          PixelgenAccentColors(hue = c("purples"), level = c(4, 6))
        ) %>%
        unname() %>%
        rep(10)
    }

    return(sample_palette)
  }

#' Cell palette
#'
#' @export
#'
cell_palette <-
  c(
    "CD4 T" = "#6D92D1",
    "CD4 Naive" = "#B9CDED",
    "Naive CD4 T" = "#B9CDED",
    "CD4 TSCM" = "#92B0E0",
    "CD4 TCM" = "#6D92D1",
    "CD4 TEM" = "#4A73C0",
    "CD4 TEFF" = "#224792",
    "CD4 TEMRA" = "#1C3A76",
    "Treg" = "#2955AE",
    "CD4 CTL" = "#2955AE",
    "CD8 T" = "#75C5B5",
    "CD8 Naive" = "#D0EDE6",
    "Naive CD8 T" = "#D0EDE6",
    "CD8 TSCM" = "#A2DACE",
    "CD8 TCM" = "#75C5B5",
    "CD8 TEM" = "#4AAF9D",
    "CD8 TEFF" = "#209785",
    "CD8 TEMRA" = "#1B7E6F",
    "other T" = "#1B7E9F",
    "Other T" = "#1B7E9F",
    "MAIT" = "#156559",
    "DPT" = "#1B7E9F",
    "DNT" = "#7B787F",
    "dnT" = "#7B787F",
    "NK" = "#A28EDB",
    "CD56dim NK" = "#A28EDB",
    "CD56bright NK" = "#866CCD",
    "NK_CD56bright" = "#866CCD",
    "NKT" = "#6C4ABD",
    "B" = "#917557",
    "B naive" = "#F0DAC4",
    "Naive B" = "#F0DAC4",
    "B intermediate" = "#F0DAC4",
    "Memory B" = "#917557",
    "B memory" = "#917557",
    "Plasma cells" = "#DE9982",
    "Plasmablast" = "#DE9982",
    "DC" = "#DA94C1",
    "cDC1" = "#DA94C1",
    "cDC2" = "#BB5391",
    "pDC" = "#9C4579",
    "Mono" = "#F0C966",
    "Classical Mono" = "#F0C966",
    "CD14 Mono" = "#F0C966",
    "Intermediate Mono" = "#DAAC22",
    "CD16 Mono" = "#836714",
    "Non-classical Mono" = "#836714",
    "other" = "#797979",
    "Neutrophils" = "#CDCDCD",
    "Basophils" = "#797979",
    "Platelet" = "#7C2628",
    "Platelets" = "#7C2628",
    "gdT" = "#5A3E9E"
  )

#' Yes/No gradient
#'
#' @export
#'
yes_no_palette <-
  c("#DAD6D7", "#C86584")

#' Cluster palette
#'
#' A palette for clusters, combining various hues and levels from the PixelgenAccentColors.
#'
#' @export
#'
cluster_palette <-
  c(PixelgenAccentColors(
    level = 7,
    hue = c("blues", "reds",
            "cyans", "oranges",
            "pinks", "greens",
            "purples", "yellows",
            "greys", "beiges")
  ),
  PixelgenAccentColors(
    level = 9,
    hue = c("blues", "pinks",
            "reds", "yellows",
            "greys", "beiges")
  ),
  PixelgenAccentColors(
    level = 5,
    hue = c("purples", "blues",
            "cyans", "greens",
            "reds", "oranges")
  ),
  PixelgenAccentColors(
    level = 3,
    hue = c("blues", "reds")
  ),
  PixelgenAccentColors(
    level = 12,
    hue = c("blues", "reds")
  )) %>%
  unname() %>%
  rep(10)

#' Theme for violin plots
#'
#' @param base_size Base font size for the theme.
#'
#' @return A ggplot2 theme object.
#'
#' @export
#'
theme_violin <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}
