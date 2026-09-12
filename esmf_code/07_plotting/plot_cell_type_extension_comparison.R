# ============================================================================
# Cell-type Extension Comparison Plot
#
# Purpose : Compare the 6-type and 7-type extension experiments against their
#           corresponding original (baseline) benchmarks, to demonstrate that
#           ESMF does not collapse as the number of cell types increases.
#
# Design  : 2 columns x 2 rows.
#             columns = scenario
#                 col 1 = Generalization / Cross-platform consistency  (WARM)
#                 col 2 = Multi-source reference integration             (COOL)
#             rows    = metric (Pearson on top, RMSE below)
#           Each panel carries ONE scenario's baseline -> extended pair:
#             Generalization : 5 types (base) -> 6 types (extended)
#             Multi-source   : 4 types (base) -> 7 types (extended)
#           Three methods per configuration (ESMF / MuSiC / CIBERSORTx).
#           Baseline bars are solid; extended bars are lighter with a dashed
#           outline. Legend shows scenario-grouped method names.
#
# Style   : Follows LSC_plots_20260531.R — Dunhuang-inspired palette,
#           theme_classic + Arial, 600 dpi output.
#
# Output  : plot/CellTypeExtension_Comparison.png / .pdf
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================
RESULTS_ROOT <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/results"
OUTPUT_DIR   <- "/media/desk16/tjn050/LSC/LSC_deconvolution/LSC/plot"

# Four configurations: directory, tissue, scenario, n_types, status, x position
# NOTE: both scenarios use the same tissue (mouse hypothalamus); they differ in
# experimental design (single-reference vs multi-source), not in tissue.
# x positions are packed tightly so the bars claim more of the axis; the wider
# gap between the two pairs is what separates the scenarios.
CONFIG <- tibble::tribble(
  ~dir,                                                        ~tissue,               ~scenario,        ~n_types, ~status,     ~x,
  "GeneralizationFromSingleReference_Mm_Hypothalamus_threshold200", "Mouse hypothalamus", "Generalization", 5,     "baseline",  1.00,
  "GeneralizationFromSingleReference_Mm_Hypothalamus_6CT",          "Mouse hypothalamus", "Generalization", 6,     "extended",  1.98,
  "MultiSourceReferenceIntegration_Mm_Hypothalamus_threshold200",   "Mouse hypothalamus", "Multi-source",   4,     "baseline",  3.03,
  "MultiSourceReferenceIntegration_Mm_Hypothalamus_7CT",            "Mouse hypothalamus", "Multi-source",   7,     "extended",  4.01
)

# In the 2x2 layout each panel carries ONE scenario, so x positions are local to
# the panel (baseline at 1, extended at 2) instead of global. The baseline-vs-
# extended status is carried by a lighter "powder" tone of the same pigment
# plus a dashed border (no separate legend).
CONFIG <- CONFIG %>%
  mutate(
    x_label       = sprintf("%s\n%s, %d Cell Types", tools::toTitleCase(tissue), tools::toTitleCase(scenario), n_types),
    x_local       = ifelse(status == "baseline", 1.0, 2.0),
    x_label_short = sprintf("%d Cell Types\n(%s)", n_types, tools::toTitleCase(status))
  )

METHODS_OF_INTEREST <- c("ESMF", "MuSiC", "CIBERSORTx")

# Dunhuang palette. Every scenario x method pair carries TWO tones:
#   baseline : deep, muted Dunhuang pigment — 沉稳
#   extended : a lightened, powder version of the same pigment — 粉
# Keys are "<scenario>.<method>.<status>", matching fill_grp2 built below.
METHOD_COLORS <- c(
  # --- Generalization: warm Dunhuang pigments ---
  "Generalization.ESMF.baseline"       = "#8C2B22",   # 深朱砂 cinnabar
  "Generalization.MuSiC.baseline"      = "#A5662A",   # 深赭石 ochre
  "Generalization.CIBERSORTx.baseline" = "#B8891F",   # 深藤黄 gamboge
  "Generalization.ESMF.extended"       = "#D9908A",   # 粉朱砂
  "Generalization.MuSiC.extended"      = "#E0AE84",   # 粉赭
  "Generalization.CIBERSORTx.extended" = "#EBD3A0",   # 粉藤黄
  # --- Multi-source: cool Dunhuang pigments ---
  "Multi-source.ESMF.baseline"         = "#1E3A5F",   # 深藏青 navy
  "Multi-source.MuSiC.baseline"        = "#276B45",   # 深石绿 malachite
  "Multi-source.CIBERSORTx.baseline"   = "#3A7070",   # 深石青 azurite
  "Multi-source.ESMF.extended"         = "#9DB4CE",   # 粉藏青
  "Multi-source.MuSiC.extended"        = "#96C4A8",   # 粉石绿
  "Multi-source.CIBERSORTx.extended"   = "#A6C8C8"    # 粉石青
)

# Legend: grouped by scenario. Each scenario is ONE group header (shown once),
# with its three methods listed beneath it. Headers use a blank key (NA colour).
GRP_WARM <- "grp_warm"
GRP_COOL <- "grp_cool"
LEGEND_SEP <- "legend_sep"          # blank separator row between the two groups
LEGEND_BREAKS <- c(
  GRP_WARM,
  "Generalization.ESMF.baseline",
  "Generalization.MuSiC.baseline",
  "Generalization.CIBERSORTx.baseline",
  LEGEND_SEP,
  GRP_COOL,
  "Multi-source.ESMF.baseline",
  "Multi-source.MuSiC.baseline",
  "Multi-source.CIBERSORTx.baseline"
)
# Labels in the SAME order as LEGEND_BREAKS (positional matching — more reliable
# than named vectors for synthetic header keys that don't exist in the data).
# Group headers are wrapped in <b> so they render bold via element_markdown.
# The separator row carries a blank label so it reads as empty space.
LEGEND_LABELS <- c(
  "<b>5 Cell Types To 6 Cell Types</b>",   # GRP_WARM
  "ESMF",                                    # Generalization.ESMF.baseline
  "MuSiC",                                   # Generalization.MuSiC.baseline
  "CIBERSORTx",                              # Generalization.CIBERSORTx.baseline
  "",                                        # LEGEND_SEP (blank gap)
  "<b>4 Cell Types To 7 Cell Types</b>",     # GRP_COOL
  "ESMF",                                    # Multi-source.ESMF.baseline
  "MuSiC",                                   # Multi-source.MuSiC.baseline
  "CIBERSORTx"                               # Multi-source.CIBERSORTx.baseline
)
# Real scenario.method.status colours for the method rows; white (invisible) for
# the two group-header rows and the blank separator so they render as
# text-only / empty legend entries.
LEGEND_VALUES <- c(setNames(c("white", "white", "white"),
                            c(GRP_WARM, GRP_COOL, LEGEND_SEP)),
                   METHOD_COLORS)

# Column headers (above the top panel of each column), wrapped to two lines,
# title case with a colon before Mm_Brain.
COL_TITLE <- c(
  "Generalization" = "Generalization from Single\nReferences: Mm_Brain",
  "Multi-source"   = "Multi-Source Reference\nIntegration: Mm_Brain"
)

# Dodging geometry — three methods per configuration.
# IMPORTANT: position_dodge() defaults to preserve = "total", which divides the
# supplied bar width by the number of groups. BAR_WIDTH must therefore be the
# FULL dodge width — setting it to DODGE_WIDTH / N_METHODS silently produced
# bars one third of the intended width (measured: 44 px bars on 144 px centres).
N_METHODS   <- length(METHODS_OF_INTEREST)
DODGE_WIDTH <- 0.94
BAR_WIDTH   <- DODGE_WIDTH

# ============================================================================
# DATA LOADING
# ============================================================================
# Prefer the CIBERSORTx-corrected "*_updated.txt" when present.
pick_file <- function(dir, pattern) {
  fs  <- list.files(file.path(RESULTS_ROOT, dir), pattern = pattern, full.names = TRUE)
  if (length(fs) == 0) return(NA_character_)
  upd <- fs[grepl("_updated\\.txt$", fs)]
  if (length(upd) > 0) return(upd[1])
  fs[1]
}

read_metric <- function(file) {
  if (is.na(file)) return(setNames(numeric(0), character(0)))
  df <- read.table(file, header = TRUE, check.names = FALSE, row.names = 1)
  stats::setNames(as.numeric(df[[1]]), rownames(df))
}

load_config <- function(dir) {
  p <- read_metric(pick_file(dir, "^H_pearson_results_.*\\.txt$"))
  r <- read_metric(pick_file(dir, "^H_rmse_results_.*\\.txt$"))
  tibble::tibble(
    method  = METHODS_OF_INTEREST,
    pearson = as.numeric(p[METHODS_OF_INTEREST]),
    rmse    = as.numeric(r[METHODS_OF_INTEREST])
  )
}

dat_long <- purrr::map_dfr(CONFIG$dir, function(d) {
  load_config(d) %>% mutate(dir = d)
}) %>%
  left_join(CONFIG %>% select(dir, scenario, n_types, status,
                              x_local, x_label_short), by = "dir") %>%
  pivot_longer(cols = c(pearson, rmse), names_to = "metric", values_to = "value") %>%
  mutate(
    metric  = factor(metric, levels = c("pearson", "rmse"),
                     labels = c("Pearson", "RMSE")),
    method  = factor(method, levels = METHODS_OF_INTEREST),
    status  = factor(status, levels = c("baseline", "extended"))
  ) %>%
  # Fill key = "<scenario>.<method>.<status>", so each scenario gets its own
  # colour family.  The two group-header levels (GRP_WARM / GRP_COOL) are
  # included so that the legend can display them as section titles even though
  # no bar actually carries those values.
  mutate(fill_grp2 = factor(paste(scenario, method, status, sep = "."),
                            levels = c(GRP_WARM, GRP_COOL, names(METHOD_COLORS))))

# Verify that nothing failed to load
stopifnot(!any(is.na(dat_long$value)))
cat("Loaded", nrow(dat_long), "observations across", n_distinct(dat_long$dir), "configurations\n")

# ============================================================================
# DELTA SUMMARY (console only — deliberately not annotated on the plot)
# ============================================================================
esmf <- dat_long %>% filter(method == "ESMF")

delta_tbl <- esmf %>%
  select(scenario, status, metric, value) %>%
  pivot_wider(names_from = status, values_from = value) %>%
  mutate(delta = extended - baseline)

# ============================================================================
# SCI THEME
# ============================================================================
sci_theme <- theme_classic(base_size = 10) +
  theme(
    text               = element_text(family = "Arial", colour = "black"),
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 10.5,
                                      lineheight = 1.05, margin = margin(b = 6)),
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(size = 10, face = "bold", margin = margin(r = 8)),
    # 45-degree x labels, matching LSC_plots_20260531.R
    axis.text.x        = element_text(size = 7, colour = "black",
                                      angle = 45, hjust = 1, vjust = 1.05,
                                      lineheight = 0.9),
    axis.text.y        = element_text(size = 9, colour = "black"),
    axis.line          = element_line(colour = "black", linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    strip.text         = element_text(size = 10, face = "bold", hjust = 0.5),
    strip.background   = element_blank(),
    # Legend on the RIGHT — scenario-grouped method names (2 lines per entry).
    legend.position      = "right",
    legend.justification = "center",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 6.5),
    legend.key.size      = unit(0.35, "cm"),
    legend.key.width     = unit(0.35, "cm"),
    legend.key.height    = unit(0.28, "cm"),
    legend.spacing.x     = unit(0.15, "cm"),
    legend.spacing.y     = unit(0.12, "cm"),
    legend.box.spacing   = unit(0.30, "cm"),
    legend.margin        = margin(t = 0, r = 0, b = 0, l = 0),
    plot.margin          = margin(t = 4, r = 6, b = 4, l = 6)
  )

# ============================================================================
# PANEL BUILDER
# ============================================================================
make_panel <- function(metric_name, scenario_name, is_top = FALSE,
                       y_breaks, y_limits, value_fmt = "%.3f") {

  d     <- dat_long %>% filter(metric == metric_name, scenario == scenario_name)
  d_cfg <- CONFIG %>% filter(scenario == scenario_name)

  ggplot() +
    # bars
    geom_col(data = d,
             aes(x = x_local, y = value, fill = fill_grp2,
                 linetype = status, group = interaction(x_local, method)),
             position = position_dodge(width = DODGE_WIDTH),
             width = BAR_WIDTH, colour = "black", linewidth = 0.22) +

    # numeric labels on top of bars
    geom_text(data = d,
              aes(x = x_local, y = value, label = sprintf(value_fmt, value),
                  group = interaction(x_local, method)),
              position = position_dodge(width = DODGE_WIDTH),
              vjust = -0.35, size = 1.75, family = "Arial") +

    scale_fill_manual(name = NULL, values = LEGEND_VALUES,
                      breaks = LEGEND_BREAKS, labels = LEGEND_LABELS) +
    scale_linetype_manual(name = "Configuration",
                          values = c(baseline = "solid", extended = "dashed")) +
    scale_x_continuous(breaks = d_cfg$x_local, labels = d_cfg$x_label_short,
                       expand = expansion(add = 0.30)) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits,
                       expand = expansion(mult = c(0, 0.02))) +

    # Legend collected by patchwork on the right; baseline vs extended is
    # shown by solid vs dashed outline only (no separate legend entry).
    guides(
      fill = guide_legend(ncol = 1, byrow = TRUE,
                          override.aes = list(linetype = "solid")),
      linetype = "none"
    ) +

    # The scenario names the COLUMN; its title sits above the top panel only.
    # The metric names the ROW and is carried by the y-axis title.
    labs(title = if (is_top) COL_TITLE[[scenario_name]] else NULL,
         y = metric_name, x = NULL) +
    sci_theme
}

# ============================================================================
# COMPOSE  (2 columns = scenario, 2 rows = metric)
# ============================================================================
# Open a throwaway ragg device so every ggplot_build() text measurement (legend +
# panels) resolves the "Arial" font instead of warning under the null device. The
# real figure is written by ggsave() below; this device is closed afterwards.
.tmp_dev <- tempfile(fileext = ".png")
ragg::agg_png(.tmp_dev, width = 10, height = 10)
# Column 1 = Generalization / Cross-platform (WARM),  rows = Pearson / RMSE
# Column 2 = Multi-source integration          (COOL), rows = Pearson / RMSE
# Each column carries its scenario title on the top panel; the right-hand guide
# area shows the grouped Method legend.
p_gen_pearson <- make_panel("Pearson", "Generalization", is_top = TRUE,
                            y_breaks = seq(0, 1.0, by = 0.2),
                            y_limits = c(0, 1.06))
p_gen_rmse    <- make_panel("RMSE", "Generalization", is_top = FALSE,
                            y_breaks = seq(0, 0.16, by = 0.04),
                            y_limits = c(0, 0.175))
p_mul_pearson <- make_panel("Pearson", "Multi-source", is_top = TRUE,
                            y_breaks = seq(0, 1.0, by = 0.2),
                            y_limits = c(0, 1.06))
p_mul_rmse    <- make_panel("RMSE", "Multi-source", is_top = FALSE,
                            y_breaks = seq(0, 0.16, by = 0.04),
                            y_limits = c(0, 0.175))

# Build a DEDICATED legend MANUALLY (not via guide_legend). The guide always
# indents every label to the right of its swatch, which pushed the group headers
# too far in. Drawing it ourselves lets the bold headers start at the swatch's
# left edge (aligned with the colours) while the method names sit to the right
# of their swatches.
#
# Rows use explicit CM coordinates (y grows UPWARD, so the warm group is on top
# and the cool group below). The three methods of each group are packed tightly
# (0.8 cm apart); a larger gap (1.8 cm) separates the two groups.
legend_df <- data.frame(
  y = c(6.6, 5.8, 5.0, 4.2, 2.4, 1.6, 0.8, 0.0),
  kind = c("header", "method", "method", "method",
           "header", "method", "method", "method"),
  label = c("5 Cell Types To 6 Cell Types", "ESMF", "MuSiC", "CIBERSORTx",
            "4 Cell Types To 7 Cell Types", "ESMF", "MuSiC", "CIBERSORTx"),
  color = c(NA, "#8C2B22", "#A5662A", "#B8891F",
            NA, "#1E3A5F", "#276B45", "#3A7070")
)
shared_legend <- ggplot() +
  # coloured swatches for the three methods of each scenario
  geom_rect(data = legend_df %>% filter(kind == "method"),
            aes(xmin = 0, xmax = 0.4, ymin = y - 0.3,
                ymax = y + 0.3, fill = color)) +
  scale_fill_identity() +
  # method names — to the right of their swatch
  ggtext::geom_richtext(data = legend_df %>% filter(kind == "method"),
            aes(x = 0.65, y = y, label = label),
            hjust = 0, vjust = 0.5, size = 2.1,
            colour = "black", fill = NA, label.colour = NA) +
  # bold group headers — left-aligned with the swatch column
  ggtext::geom_richtext(data = legend_df %>% filter(kind == "header"),
            aes(x = 0, y = y, label = paste0("<b>", label, "</b>")),
            hjust = 0, vjust = 0.5, size = 2.1,
            colour = "black", fill = NA, label.colour = NA) +
  scale_x_continuous(limits = c(-0.05, 5.4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.4, 7.0), expand = c(0, 0)) +
  theme_void(base_size = 10) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Suppress per-panel legends so the shared legend is not duplicated.
p_gen_pearson <- p_gen_pearson + theme(legend.position = "none")
p_gen_rmse    <- p_gen_rmse    + theme(legend.position = "none")
p_mul_pearson <- p_mul_pearson + theme(legend.position = "none")
p_mul_rmse    <- p_mul_rmse    + theme(legend.position = "none")

col_gen <- p_gen_pearson / p_gen_rmse          # warm column
col_mul <- p_mul_pearson / p_mul_rmse          # cool column

# Cap the legend height with a spacer so patchwork doesn't stretch it to the
# full panel height (which would spread the within-group items too far apart).
combined <- (col_gen | col_mul | (shared_legend / plot_spacer())) +
  plot_layout(widths = c(2, 2, 1.6))

# ============================================================================
# SAVE
# ============================================================================
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# 2x2 grid: each panel is roughly 10 x 8 cm, so bars stay thick even though
# every panel now carries only two configurations (baseline + extended).
suppressWarnings(ggsave(file.path(OUTPUT_DIR, "CellTypeExtension_Comparison.png"),
       plot = combined, width = 20, height = 16, units = "cm",
       dpi = 600, device = ragg::agg_png, background = "white"))

suppressWarnings(ggsave(file.path(OUTPUT_DIR, "CellTypeExtension_Comparison.pdf"),
       plot = combined, width = 20, height = 16, units = "cm",
       device = cairo_pdf, bg = "white"))

cat("Saved to:", file.path(OUTPUT_DIR, "CellTypeExtension_Comparison.png"), "\n")
cat("          ", file.path(OUTPUT_DIR, "CellTypeExtension_Comparison.pdf"), "\n")

dev.off()   # close the throwaway ragg device opened before COMPOSE

# ============================================================================
# CONSOLE SUMMARY (for verification)
# ============================================================================
cat("\n--- Summary (Pearson / RMSE) ---\n")
print(dat_long %>%
        select(scenario, n_types, status, method, metric, value) %>%
        pivot_wider(names_from = metric, values_from = value) %>%
        arrange(scenario, n_types, method))

cat("\n--- ESMF delta (extended - baseline) ---\n")
print(delta_tbl %>% mutate(across(c(baseline, extended, delta), ~round(., 4))))
