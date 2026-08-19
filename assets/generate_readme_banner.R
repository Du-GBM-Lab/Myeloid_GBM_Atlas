#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
  library(svglite)
  library(ragg)
})

font_family <- "Arial"

palette <- list(
  ink = "#263238",
  teal_dark = "#0F5C6E",
  teal = "#2A9D8F",
  blue = "#4C78A8",
  blue_light = "#EAF2F8",
  purple = "#7E57A5",
  purple_light = "#F2ECF8",
  gold = "#E9C46A",
  gold_light = "#FFF7DE",
  coral = "#E76F51",
  coral_light = "#FDEDEA",
  gray = "#7A858B",
  gray_light = "#F5F7F8",
  border = "#D5E0E5",
  white = "#FFFFFF"
)

txt <- function(label, x, y, size = 12, color = palette$ink,
                face = "plain", just = "center", rot = 0) {
  grid.text(
    label,
    x = unit(x, "npc"), y = unit(y, "npc"),
    just = just, rot = rot,
    gp = gpar(fontfamily = font_family, fontsize = size,
              fontface = face, col = color, lineheight = 1.05)
  )
}

round_box <- function(x, y, w, h, fill = palette$white,
                      color = palette$border, lwd = 1.5, radius = 0.018) {
  grid.roundrect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = unit(radius, "npc"),
    gp = gpar(fill = fill, col = color, lwd = lwd)
  )
}

arrow_line <- function(x1, y1, x2, y2, color = palette$ink,
                       lwd = 3.2, lty = 1, head = 0.11) {
  grid.lines(
    x = unit(c(x1, x2), "npc"), y = unit(c(y1, y2), "npc"),
    arrow = arrow(type = "closed", length = unit(head, "inches")),
    gp = gpar(col = color, lwd = lwd, lty = lty,
              lineend = "round", linejoin = "round")
  )
}

pill <- function(label, x, y, w, h, fill, color, size = 9.5,
                 text_color = palette$ink, face = "bold") {
  round_box(x, y, w, h, fill = fill, color = color, lwd = 1.4, radius = 0.012)
  txt(label, x, y, size = size, color = text_color, face = face)
}

draw_myelin <- function(x, y) {
  offsets <- matrix(c(-0.026, 0.022, 12, -0.006,
                       0.000, 0.000, -8, 0.000,
                       0.026, -0.020, 10, 0.006), ncol = 4, byrow = TRUE)
  for (i in seq_len(nrow(offsets))) {
    pushViewport(viewport(x = x + offsets[i, 1], y = y + offsets[i, 2],
                          width = 0.048, height = 0.022,
                          angle = offsets[i, 3]))
    grid.roundrect(
      r = unit(0.35, "snpc"),
      gp = gpar(fill = palette$gold_light, col = palette$gold, lwd = 2)
    )
    grid.lines(x = unit(c(0.35, 0.35), "npc"), y = unit(c(0.1, 0.9), "npc"),
               gp = gpar(col = palette$gold, lwd = 1))
    grid.lines(x = unit(c(0.65, 0.65), "npc"), y = unit(c(0.1, 0.9), "npc"),
               gp = gpar(col = palette$gold, lwd = 1))
    popViewport()
  }
}

draw_tam <- function(x, y, radius = 0.058) {
  for (ang in seq(0, 315, by = 45)) {
    dx <- cos(ang * pi / 180) * radius * 0.78
    dy <- sin(ang * pi / 180) * radius * 0.78
    grid.circle(x = unit(x + dx, "npc"), y = unit(y + dy, "npc"),
                r = unit(radius * 0.32, "npc"),
                gp = gpar(fill = palette$blue_light, col = palette$blue, lwd = 1.6))
  }
  grid.circle(x = unit(x, "npc"), y = unit(y, "npc"),
              r = unit(radius, "npc"),
              gp = gpar(fill = palette$blue_light, col = palette$blue, lwd = 2.4))
  grid.circle(x = unit(x, "npc"), y = unit(y - 0.004, "npc"),
              r = unit(radius * 0.36, "npc"),
              gp = gpar(fill = palette$blue, col = palette$blue, alpha = 0.75))
  grid.circle(x = unit(x - 0.026, "npc"), y = unit(y + 0.021, "npc"),
              r = unit(0.009, "npc"),
              gp = gpar(fill = palette$gold_light, col = palette$gold, lwd = 1.2))
  txt("SPP1 high", x, y + 0.09, size = 10.5, color = palette$coral, face = "bold")
  txt("Phagocytic-suppressive TAM", x, y - 0.096, size = 9.3, face = "bold")
}

draw_gbm <- function(x, y, radius = 0.06) {
  for (ang in seq(22.5, 337.5, by = 45)) {
    dx <- cos(ang * pi / 180) * radius * 0.82
    dy <- sin(ang * pi / 180) * radius * 0.82
    grid.circle(x = unit(x + dx, "npc"), y = unit(y + dy, "npc"),
                r = unit(radius * 0.29, "npc"),
                gp = gpar(fill = palette$purple_light, col = palette$purple, lwd = 1.6))
  }
  grid.circle(x = unit(x, "npc"), y = unit(y, "npc"),
              r = unit(radius, "npc"),
              gp = gpar(fill = palette$purple_light, col = palette$purple, lwd = 2.4))
  grid.circle(x = unit(x + 0.004, "npc"), y = unit(y - 0.004, "npc"),
              r = unit(radius * 0.38, "npc"),
              gp = gpar(fill = palette$purple, col = palette$purple, alpha = 0.72))
  grid.lines(x = unit(c(x - 0.066, x - 0.052, x - 0.043), "npc"),
             y = unit(c(y + 0.006, y + 0.025, y + 0.006), "npc"),
             gp = gpar(col = palette$purple, lwd = 3, lineend = "round"))
  txt("CD44", x - 0.066, y + 0.064, size = 10, color = palette$purple, face = "bold")
  txt("MES-like GBM cell", x, y - 0.096, size = 9.5, face = "bold")
}

draw_banner <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = palette$white, col = NA))

  txt("A myelin-associated SPP1-CD44 axis shapes an immunosuppressive GBM niche",
      0.5, 0.945, size = 23, color = palette$ink, face = "bold")
  txt("Integrated single-cell, spatial, temporal, and proteomic evidence",
      0.5, 0.895, size = 11.5, color = palette$gray)

  # Left: evidence base
  round_box(0.125, 0.49, 0.215, 0.70, fill = palette$gray_light,
            color = palette$teal_dark, lwd = 1.8, radius = 0.02)
  txt("Integrated atlas", 0.125, 0.79, size = 14, color = palette$teal_dark, face = "bold")
  pill("scRNA-seq\n231 samples | 1,135,677 cells", 0.125, 0.675, 0.175, 0.105,
       palette$white, palette$blue, size = 8.8)
  pill("Spatial transcriptomics\n25 tissue sections", 0.125, 0.545, 0.175, 0.098,
       palette$white, palette$teal, size = 8.8)
  pill("Temporal models", 0.125, 0.425, 0.175, 0.075,
       palette$white, palette$purple, size = 9.3)
  pill("Proteomics", 0.125, 0.325, 0.175, 0.075,
       palette$white, palette$coral, size = 9.3)
  txt("502,672 myeloid cells", 0.125, 0.225, size = 9.8,
      color = palette$teal_dark, face = "bold")

  arrow_line(0.235, 0.49, 0.275, 0.49, color = palette$teal_dark, lwd = 3.4)
  txt("Recurrent\nmyeloid states", 0.255, 0.56, size = 8.5,
      color = palette$teal_dark, face = "bold")

  # Center: canonical mechanism
  round_box(0.505, 0.49, 0.45, 0.70, fill = palette$white,
            color = palette$border, lwd = 1.6, radius = 0.02)
  txt("Myelin-associated niche mechanism", 0.505, 0.79, size = 14,
      color = palette$ink, face = "bold")

  draw_myelin(0.335, 0.51)
  txt("Myelin debris", 0.335, 0.405, size = 9.5, face = "bold")
  arrow_line(0.375, 0.51, 0.42, 0.51, color = palette$ink, lwd = 3.8)
  txt("Myelin phagocytosis", 0.397, 0.585, size = 8.4, color = palette$gray)

  draw_tam(0.475, 0.51)

  # SPP1 ligand field and evidence-bounded receptor-axis arrow
  ligand_x <- c(0.545, 0.563, 0.581, 0.555, 0.576, 0.594)
  ligand_y <- c(0.535, 0.558, 0.538, 0.500, 0.508, 0.532)
  for (i in seq_along(ligand_x)) {
    grid.circle(x = unit(ligand_x[i], "npc"), y = unit(ligand_y[i], "npc"),
                r = unit(0.0053, "npc"),
                gp = gpar(fill = palette$coral, col = palette$coral))
  }
  txt("SPP1", 0.57, 0.595, size = 10, color = palette$coral, face = "bold")
  arrow_line(0.535, 0.51, 0.61, 0.51, color = palette$coral,
             lwd = 2.5, lty = 3, head = 0.09)

  draw_gbm(0.665, 0.51)
  txt("Myeloid-malignant crosstalk", 0.57, 0.295, size = 9.4,
      color = palette$teal_dark, face = "bold")
  pill("Immunosuppressive niche", 0.57, 0.225, 0.19, 0.060,
       palette$gray_light, palette$border, size = 9.2)

  # Right: two distinct intervention branches
  round_box(0.865, 0.49, 0.23, 0.70, fill = palette$gray_light,
            color = palette$teal_dark, lwd = 1.8, radius = 0.02)
  txt("Therapeutic implications", 0.865, 0.79, size = 13.5,
      color = palette$teal_dark, face = "bold")

  pill("Tumor-cell SPP1 knockdown", 0.865, 0.68, 0.195, 0.068,
       palette$coral_light, palette$coral, size = 9.5,
       text_color = palette$coral)
  arrow_line(0.865, 0.642, 0.865, 0.60, color = palette$coral, lwd = 2.2, head = 0.075)
  pill("Reduced tumor growth", 0.815, 0.555, 0.092, 0.068,
       palette$white, palette$coral, size = 8.2)
  pill("Prolonged survival", 0.915, 0.555, 0.092, 0.068,
       palette$white, palette$coral, size = 8.2)

  grid.lines(x = unit(c(0.775, 0.955), "npc"), y = unit(c(0.475, 0.475), "npc"),
             gp = gpar(col = palette$border, lwd = 1.4, lty = 2))

  pill("PLX5622 + anti-PD-1", 0.865, 0.39, 0.195, 0.068,
       palette$purple_light, palette$purple, size = 9.5,
       text_color = palette$purple)
  arrow_line(0.865, 0.352, 0.865, 0.31, color = palette$purple, lwd = 2.2, head = 0.075)
  pill("Improved tumor control", 0.815, 0.265, 0.092, 0.068,
       palette$white, palette$purple, size = 8.2)
  pill("Prolonged survival", 0.915, 0.265, 0.092, 0.068,
       palette$white, palette$purple, size = 8.2)

  # Evidence legend
  grid.lines(x = unit(c(0.31, 0.345), "npc"), y = unit(c(0.075, 0.075), "npc"),
             gp = gpar(col = palette$ink, lwd = 2.5))
  txt("experimentally supported transition or perturbation", 0.353, 0.075,
      size = 7.8, color = palette$gray, just = "left")
  grid.lines(x = unit(c(0.58, 0.615), "npc"), y = unit(c(0.075, 0.075), "npc"),
             gp = gpar(col = palette$coral, lwd = 2.2, lty = 3))
  txt("inferred communication or receptor-axis model", 0.623, 0.075,
      size = 7.8, color = palette$gray, just = "left")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "."
asset_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
qa_dir <- file.path(asset_dir, "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
svg_path <- file.path(asset_dir, "myeloid_gbm_axis.svg")
png_path <- file.path(asset_dir, "myeloid_gbm_axis.png")
pdf_path <- file.path(qa_dir, "myeloid_gbm_axis.pdf")
tiff_path <- file.path(qa_dir, "myeloid_gbm_axis.tiff")
width_mm <- 406.4
height_mm <- 162.56
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(svg_path, width = width_in, height = height_in, bg = palette$white,
                 system_fonts = list(Arial = "Arial"))
draw_banner()
invisible(dev.off())

ragg::agg_png(png_path, width = width_in, height = height_in, units = "in",
              res = 300, background = palette$white, scaling = 1)
draw_banner()
invisible(dev.off())

grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in,
                     family = font_family, bg = palette$white)
draw_banner()
invisible(dev.off())

ragg::agg_tiff(tiff_path, width = width_in, height = height_in, units = "in",
               res = 600, background = palette$white, compression = "lzw")
draw_banner()
invisible(dev.off())

message("Wrote: ", svg_path)
message("Wrote: ", png_path)
message("QA export: ", pdf_path)
message("QA export: ", tiff_path)
