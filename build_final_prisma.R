#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(grid)
})

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(root, "3 outputs", "figures")
manuscript_dir <- file.path(root, "Manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#33506e"

main_boxes <- tribble(
  ~x,   ~y,    ~w,   ~h,   ~fill,      ~label,
  0.30, 0.90,  0.56, 0.10, "#eef3f8",  "Records identified (n = 2,356)\nPubMed 507, Scopus 744\nCENTRAL 848, PEDro 257",
  0.30, 0.74,  0.56, 0.10, "#eef3f8",  "After within-source cleaning (n = 2,314)\nCross-database de-duplication removed 981",
  0.30, 0.575, 0.56, 0.11, "#eef3f8",  "Title / abstract records screened (n = 1,033)\nExcluded: 805, retained for full-text review: 228",
  0.30, 0.40,  0.56, 0.10, "#eef3f8",  "Substantive full-text reports assessed (n = 157)\nAgreement: 151/157 (96.2%); Cohen's kappa = 0.92",
  0.30, 0.18,  0.56, 0.14, "#e8f3ea",  "Studies included in review (n = 80)\nReports of included studies (n = 82)\nPrimary PREV meta-analysis: 2 studies\nNo IPC pool or component network"
)

side_boxes <- tribble(
  ~x,   ~y,    ~w,   ~h,   ~fill,      ~label,
  0.82, 0.74,  0.34, 0.115, "#fff2cc", "Registry-only records (n = 300)\ntracked separately as ongoing /\nawaiting results",
  0.82, 0.575, 0.34, 0.12,  "#eef3f8", "Independent calibration sample (n = 200)\nAgreement: 90.5%; Cohen's kappa = 0.73",
  0.82, 0.40,  0.34, 0.12,  "#fff2cc", "Reports not retrieved as substantive\nfull text (n = 71)\nAbstract-only: 8; no matched full text: 63",
  0.82, 0.22,  0.34, 0.22,  "#fce4d6", "Reports excluded after full-text assessment (n = 75)\nWrong intervention/comparator: 30\nNo eligible volume outcome: 18\nNot randomized/ineligible design: 11\nProtocol/no primary results: 5\nDuplicate/secondary without unique data: 4\nWrong population: 3; not upper-limb BCRL: 2\nPrevention: 2"
)

all_boxes <- bind_rows(main_boxes, side_boxes) %>%
  mutate(text_size = ifelse(grepl("^Reports excluded", label), 2.3, 2.75))

figure <- ggplot() +
  geom_tile(
    data = all_boxes,
    aes(x, y, width = w, height = h, fill = fill),
    colour = ink,
    linewidth = 0.5
  ) +
  geom_text(
    data = all_boxes,
    aes(x, y, label = label, size = text_size),
    lineheight = 0.95,
    colour = "#1d2d40"
  ) +
  scale_fill_identity() +
  scale_size_identity() +
  annotate(
    "segment",
    x = 0.30,
    xend = 0.30,
    y = c(0.85, 0.69, 0.52, 0.33),
    yend = c(0.79, 0.625, 0.45, 0.25),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = ink,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 0.58,
    xend = 0.65,
    y = c(0.74, 0.575, 0.40, 0.34),
    yend = c(0.74, 0.575, 0.40, 0.31),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    colour = ink,
    linewidth = 0.5
  ) +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(0.10, 0.99), expand = FALSE) +
  labs(title = "Figure 1. PRISMA 2020 study-selection flow") +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0,
      margin = margin(b = 6)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

output_path <- file.path(out_dir, "Fig1_PRISMA_flow.pdf")
ggsave(output_path, figure, width = 8.6, height = 8.4, device = cairo_pdf)
file.copy(output_path, file.path(manuscript_dir, "Fig1_PRISMA_flow.pdf"), overwrite = TRUE)

cat("Wrote final PRISMA figure to outputs and manuscript folders\n")
