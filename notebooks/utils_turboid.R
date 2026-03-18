library(ComplexHeatmap)
library(circlize)
library(readxl)
library(tidyverse)
library(Cairo)
library(ggrepel)
library(extrafont)
library(GetoptLong)
library(tidytext)
library(lemon)
library(scales)
library(httr)
library(ggforce)

# constants for drawing figures
LINE_WIDTH = 0.25 / 2
FONT_FAMILY = "Arial"
FONT_SIZE = 6
FONT_SIZE_MM = 2 + (7 / 60)
LFC_CUTOFF = log2(1.5)
TEXT_ELEMENT = element_text(size = FONT_SIZE,
                            color = "black",
                            family = FONT_FAMILY)
AXIS_LINE = element_line(colour = "black", linewidth = LINE_WIDTH)
REGULATION_COLORS <-
  c(
    "Significant Up" = "indianred2",
    "Significant Down" = "steelblue2",
    "Not Significant" = "lightgrey",
    "Not Significant Up" = "indianred2",
    "Not Significant Down" = "steelblue2",
    "Significant but <1.5 FC" = "lightgrey"
  )

POINT_FILL = "#9999cb"

ALPHA_PALETTE <-
  c(
    "Significant Up" = 1,
    "Significant Down" = 1,
    "Not Significant" = .3,
    "Not Significant Up" = .3,
    "Not Significant Down" = .3,
    "Significant but <1.5 FC" = .3
  )

# ggplot theme used throughout
turboid_theme <- function() {
  theme_classic() +
    theme(
      axis.text.x = TEXT_ELEMENT,
      axis.text.y = TEXT_ELEMENT,
      axis.title.x = TEXT_ELEMENT,
      axis.title.y = TEXT_ELEMENT,
      legend.title = TEXT_ELEMENT,
      legend.text = TEXT_ELEMENT,
      axis.line = element_line(colour = "black", linewidth = 0),
      axis.ticks = element_line(colour = "black", linewidth = LINE_WIDTH),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = LINE_WIDTH
      ),
      strip.text.x = TEXT_ELEMENT,
      plot.title = element_text(
        size = 6,
        color = "black",
        family = FONT_FAMILY,
        hjust = 0.5
      ),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.key.spacing.x = unit(0, "mm"),
      legend.key.spacing = unit(0, "mm"),
      legend.margin = margin(-10, 0, 0, 0),
      legend.spacing = unit(0, "pt"),
      legend.justification = "top",
      strip.clip = "off"
      
    )
}

legend_theme <- function() {
  turboid_theme() +
    theme(
      aspect.ratio = 1,
      legend.position = "top",
      legend.margin = margin(l = -2, b = -5.5),
      legend.box.margin = margin(l = 0, b = -5.7),
      legend.key.spacing.x = unit(-0.1, "cm"),
      legend.key = element_blank()
    )
}

# constants / functions for data management

PIPELINE_RESULTS_PATH = "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/TurboID manuscript/Mass-spectrometry-datasets/03_results/04_turboid_downstream_results/05_results/"
PLOT_RESULTS_PATH = '/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/TurboID manuscript/Data visualization/03_plots/'

read_experiment <- function(experiment, sample_type) {
  suffix <- ""
  if (sample_type == "tissue") {
    sample_type = "01_tissue"
  } else {
    sample_type = "02_serum"
    suffix <- "_enrichment_filtered"
  }
  df <- paste0(
    PIPELINE_RESULTS_PATH,
    sample_type,
    "/",
    experiment,
    "/",
    "final_protein_table_",
    experiment,
    suffix,
    ".csv"
  ) %>% read_csv(show_col_types = FALSE)
  return(df)
}

save_plot <- function(fn, wd, ht) {
  ggsave(
    paste0(PLOT_RESULTS_PATH, fn, ".pdf"),
    width = wd,
    height = ht,
    dpi = 300,
    device = cairo_pdf
  )
}

add_alias <- function(df) {
  aliases <- "reference_data/protein_aliases.csv" %>%
    read_csv(show_col_types = FALSE)
  df <- df %>%
    mutate(protein = gsub(pattern = "_MOUSE", "", `Entry Name`)) %>%
    left_join(aliases, by = "protein") %>%
    mutate(alias = coalesce(alias, protein))
  return(df)
}