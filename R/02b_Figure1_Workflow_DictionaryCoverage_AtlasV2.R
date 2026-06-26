# =========================================================
# 02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v6_concept_nodes.R
#
# Main Figure 1, clinically focused version with explicit infection-context layer and concept-node labelling.
#
# This revision keeps the corpus/year-count graphing logic and dictionary
# coverage calculations unchanged, but simplifies the dictionary-construction
# panel and explicitly separates drug and infection terms as treatment-context and infection-context layers.
#
# Recommended run order:
#   source("00_setup_AtlasV2.R")
#   source("02_BuildHitsMatrix_AtlasV2.R")
#   source("02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v6_concept_nodes.R")
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)

quiet_install(c(
  "dplyr", "readr", "stringr", "tibble", "tidyr", "ggplot2",
  "patchwork", "scales", "systemfonts", "showtext"
))

cands  <- c("Hiragino Sans", "Hiragino Kaku Gothic ProN", "Yu Gothic", "IPAexGothic",
            "Noto Sans CJK JP", "Arial", "Helvetica", "DejaVu Sans")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
base_font <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

# -------------------------
# Inputs
# -------------------------
DIC <- load_atlas_dictionary(require_dictionary = TRUE)
DIC_PRIMARY <- select_primary_dictionary_terms(DIC)

f_run_sum <- file.path(DIR_TABLE, sprintf("dictionary_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_term_cov <- file.path(DIR_TABLE, sprintf("dictionary_term_coverage_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_year_counts <- file.path(DIR_TABLE, sprintf("year_counts_stacked_by_pubclass_%s__%s.csv", CORPUS_TAG, DIC_TAG))
if (!file.exists(f_year_counts)) {
  f_year_counts <- find_latest_file(DIR_PROC, "^year_counts_stacked_by_pubclass_[0-9]{8}\\.csv$")
}

RUN_SUM <- if (file.exists(f_run_sum)) readr::read_csv(f_run_sum, show_col_types = FALSE) else tibble::tibble()
TERM_COV <- if (file.exists(f_term_cov)) readr::read_csv(f_term_cov, show_col_types = FALSE) else tibble::tibble()
YEAR_COUNTS <- if (!is.na(f_year_counts) && file.exists(f_year_counts)) readr::read_csv(f_year_counts, show_col_types = FALSE) else tibble::tibble()

n_articles <- if (nrow(RUN_SUM) && "n_articles" %in% names(RUN_SUM)) RUN_SUM$n_articles[1] else NA_integer_
n_terms <- if (nrow(RUN_SUM) && "n_dictionary_terms" %in% names(RUN_SUM)) RUN_SUM$n_dictionary_terms[1] else nrow(DIC)
n_primary <- if (nrow(RUN_SUM) && "n_primary_terms" %in% names(RUN_SUM)) RUN_SUM$n_primary_terms[1] else nrow(DIC_PRIMARY)
n_treatment <- if (nrow(RUN_SUM) && "n_treatment_context_terms" %in% names(RUN_SUM)) RUN_SUM$n_treatment_context_terms[1] else sum(DIC$atlas_layer == "treatment_context", na.rm = TRUE)
n_infection <- if (nrow(RUN_SUM) && "n_infection_context_terms" %in% names(RUN_SUM)) RUN_SUM$n_infection_context_terms[1] else if ("atlas_layer" %in% names(DIC)) sum(DIC$atlas_layer == "infection_context", na.rm = TRUE) else 0L
n_hit_terms <- if (nrow(RUN_SUM) && "n_terms_with_hit" %in% names(RUN_SUM)) RUN_SUM$n_terms_with_hit[1] else NA_integer_
coverage_pct <- if (nrow(RUN_SUM) && "abstract_coverage_pct" %in% names(RUN_SUM)) RUN_SUM$abstract_coverage_pct[1] else NA_real_
articles_with_hit <- if (nrow(RUN_SUM) && "articles_with_any_dictionary_hit" %in% names(RUN_SUM)) RUN_SUM$articles_with_any_dictionary_hit[1] else NA_integer_

fmt_n <- function(x) {
  if (is.na(x)) "not rerun" else scales::comma(x)
}
fmt_pct <- function(x) {
  if (is.na(x)) "not rerun" else paste0(round(x, 1), "%")
}

# -------------------------
# Drawing helpers
# -------------------------
box <- function(xmin, xmax, ymin, ymax, fill = "#F8F8F8", colour = "#444444", linewidth = 0.55) {
  ggplot2::annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = fill, colour = colour, linewidth = linewidth)
}
label <- function(x, y, text, size = 3.4, fontface = "plain", hjust = 0.5, vjust = 0.5,
                  colour = "#111111", lineheight = 1.03) {
  ggplot2::annotate("text", x = x, y = y, label = text, size = size, fontface = fontface,
                    hjust = hjust, vjust = vjust, colour = colour, lineheight = lineheight,
                    family = base_font)
}
arrow_seg <- function(x, xend, y, yend, colour = "#555555") {
  ggplot2::annotate("segment", x = x, xend = xend, y = yend,
                    colour = colour, linewidth = 0.65,
                    arrow = grid::arrow(type = "closed", length = grid::unit(0.13, "inches")))
}
panel_base <- function() {
  ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_font),
      plot.margin = ggplot2::margin(3, 5, 3, 5)
    )
}

# -------------------------
# Panel A: dictionary construction / safety / coverage
# -------------------------
p_dict <- ggplot2::ggplot() +
  ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
  panel_base() +
  label(2, 96, "Dictionary construction, safety audit and corpus coverage",
        hjust = 0, fontface = "bold", size = 4.2) +
  box(2, 18, 55, 84, fill = "#F2F5F9") +
  label(10, 76, "Specialist\nliterature review", fontface = "bold", size = 3.55) +
  label(10, 64, "RA-ILD clinicians\nidentified candidate\nclinical concepts", size = 2.85) +
  arrow_seg(19, 25, 69.5, 69.5) +
  box(26, 42, 55, 84, fill = "#F8F8F8") +
  label(34, 76, "Initial concept\ndictionary", fontface = "bold", size = 3.55) +
  label(34, 64, "host / imaging /\nphysiology /\nbiomarker terms", size = 2.85) +
  arrow_seg(43, 49, 69.5, 69.5) +
  box(50, 66, 55, 84, fill = "#F2F8F2") +
  label(58, 76, "Overlap and\nacronym correction", fontface = "bold", size = 3.45) +
  label(58, 64, "regex safety\ncollision checks\nAPRIL / ROS / RF", size = 2.80) +
  arrow_seg(67, 73, 69.5, 69.5) +
  box(74, 88, 55, 84, fill = "#F9F5EE") +
  label(81, 76, "Abstract-informed\nexpansion", fontface = "bold", size = 3.45) +
  label(81, 64, "outcome wording\nfunctional decline\npathway terms", size = 2.80) +
  arrow_seg(89, 92, 69.5, 69.5) +
  box(93, 99, 55, 84, fill = "#EFF3FF") +
  label(96, 76, "Final\natlas", fontface = "bold", size = 3.35) +
  label(96, 63, paste0(fmt_n(n_terms), "\nconcept nodes"), size = 2.85) +
  box(5, 95, 12, 42, fill = "#FFFFFF", linewidth = 0.55) +
  label(8, 35, "Run-level dictionary grounding", fontface = "bold", hjust = 0, size = 3.75) +
  label(8, 27, paste0(
    "Original abstracts: ", fmt_n(n_articles),
    "    Coverage: ", fmt_pct(coverage_pct),
    "    Abstracts with ≥1 hit: ", fmt_n(articles_with_hit),
    "    Hit concept nodes: ", fmt_n(n_hit_terms)
  ), hjust = 0, size = 3.20) +
  label(8, 19, paste0("Primary atlas concept nodes: ", fmt_n(n_primary),
                      "    Treatment-context concept nodes: ", fmt_n(n_treatment),
                      "    Infection-context concept nodes: ", fmt_n(n_infection)),
        hjust = 0, size = 3.05) +
  label(8, 13, "Primary eligibility fixed before ranking; treatment and infection context retained separately.",
        hjust = 0, size = 2.95)

# -------------------------
# Panel B: clinically framed atlas analysis workflow
# -------------------------
p_workflow <- ggplot2::ggplot() +
  ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
  panel_base() +
  label(3, 96, "Outcome-preserving literature-derived analysis workflow",
        hjust = 0, fontface = "bold", size = 4.2) +
  box(4, 25, 30, 76, fill = "#F2F5F9") +
  label(14.5, 67, "PubMed\nRA-ILD literature", fontface = "bold", size = 3.8) +
  label(14.5, 50, "original research\narticles\nabstract-level analysis", size = 3.15) +
  arrow_seg(26, 34, 53, 53) +
  box(35, 67, 20, 86, fill = "#F9F9F9") +
  label(51, 79, "ABC analytic framework", fontface = "bold", size = 3.75) +
  label(39, 68, "A", fontface = "bold", colour = "#2F5597", size = 3.8) +
  label(45, 68, "upstream host / exposure /\ngenetic / serology terms", hjust = 0, size = 3.0) +
  label(39, 55, "B", fontface = "bold", colour = "#2F5597", size = 3.8) +
  label(45, 55, "imaging / physiology /\nbiomarker / pathway bridges", hjust = 0, size = 3.0) +
  label(39, 42, "C", fontface = "bold", colour = "#2F5597", size = 3.8) +
  label(45, 42, "AE-ILD / progression /\nmortality outcomes", hjust = 0, size = 3.0) +
  label(51, 28, "Treatment and infection terms: separate context layers, not primary disease-state bridge claims",
        size = 2.65) +
  arrow_seg(68, 74, 53, 53) +
  box(75, 97, 30, 76, fill = "#EFF3FF") +
  label(86, 67, "RA-ILD evidence\natlas", fontface = "bold", size = 3.8) +
  label(86, 51, "integrated across\noutcomes", size = 3.15) +
  label(86, 38, "shared core +\noutcome-centred signals", size = 3.05)

# -------------------------
# Panel C: corpus composition over time
# Graphing logic unchanged.
# -------------------------
if (nrow(YEAR_COUNTS)) {
  YEAR_COUNTS <- YEAR_COUNTS |>
    dplyr::mutate(
      corpus_role = as.character(corpus_role),
      corpus_role = dplyr::case_when(
        corpus_role == "main_original" ~ "Original",
        corpus_role == "case_report" ~ "Case report",
        corpus_role == "review_guideline" ~ "Review/Guideline",
        TRUE ~ corpus_role
      ),
      corpus_role = factor(corpus_role, levels = c("Original", "Case report", "Review/Guideline"))
    )

  p_year <- ggplot2::ggplot(YEAR_COUNTS, ggplot2::aes(x = year, y = n, fill = corpus_role)) +
    ggplot2::geom_col(width = 0.85, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_fill_manual(
      values = c("Original" = "#4E79A7", "Case report" = "#59A14F", "Review/Guideline" = "#F28E2B"),
      drop = FALSE
    ) +
    ggplot2::labs(title = "Corpus composition over time", x = "Year", y = "Records", fill = "Corpus type") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_font),
      plot.title = ggplot2::element_text(face = "plain", size = 14, margin = ggplot2::margin(b = 8)),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = c(0.23, 0.79),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9.5),
      legend.background = ggplot2::element_rect(fill = "white", colour = "#555555", linewidth = 0.3),
      plot.margin = ggplot2::margin(4, 10, 4, 4)
    )
} else {
  p_year <- ggplot2::ggplot() +
    ggplot2::theme_void(base_size = 12) +
    label(0, 0, "Run 01_FetchPubMed_AtlasV2.R to generate annual corpus counts.", hjust = 0, size = 3.2) +
    ggplot2::labs(title = "Corpus composition over time")
}

# -------------------------
# Panel D: predefined atlas layers
# Uses atlas_layer when available; otherwise falls back to analysis_tier.
# -------------------------
tier_source <- if (nrow(TERM_COV) && "atlas_layer" %in% names(TERM_COV)) TERM_COV else DIC
if (!"atlas_layer" %in% names(tier_source)) {
  tier_source <- tier_source |>
    dplyr::mutate(
      atlas_layer = dplyr::case_when(
        analysis_tier == "main" & include_primary == TRUE ~ "primary_atlas",
        analysis_tier == "treatment_context" ~ "treatment_context",
        analysis_tier == "exploratory" ~ "exploratory",
        analysis_tier == "context" ~ "context_only",
        TRUE ~ "sensitivity"
      )
    )
}

# Normalise dictionary-layer names before plotting.
# In v2.2, infection_context is a legitimate analysis layer.
# It must be listed explicitly; otherwise factor() converts it to NA and
# creates an erroneous rightmost NA bar.
allowed_layers <- c(
  "primary_atlas",
  "treatment_context",
  "infection_context",
  "sensitivity",
  "exploratory",
  "context_only"
)
allowed_node_levels <- c("specific", "aggregate", "modifier", "context")

tier_source <- tier_source |>
  dplyr::mutate(
    atlas_layer = stringr::str_trim(as.character(atlas_layer)),
    atlas_layer = dplyr::na_if(atlas_layer, ""),
    node_level = stringr::str_trim(as.character(node_level)),
    node_level = dplyr::na_if(node_level, "")
  )

unclassified_terms <- tier_source |>
  dplyr::filter(is.na(atlas_layer) | !(atlas_layer %in% allowed_layers) |
                  is.na(node_level) | !(node_level %in% allowed_node_levels)) |>
  dplyr::select(dplyr::any_of(c(
    "term", "preferred_label", "original_term", "abc_role", "abc_subrole",
    "domain", "concept_family", "analysis_tier", "atlas_layer",
    "node_level", "include_primary"
  )))

if (nrow(unclassified_terms) > 0) {
  qc_file <- file.path(
    DIR_TABLE,
    sprintf("Figure1_unclassified_dictionary_terms_v5_%s__%s.csv", CORPUS_TAG, DIC_TAG)
  )
  readr::write_csv(unclassified_terms, qc_file)
  stop(
    "Figure 1 Panel D found unclassified atlas_layer/node_level values. ",
    "A QC file was written: ", qc_file,
    ". Fix the dictionary metadata before drawing the main figure."
  )
}

tier_df <- tier_source |>
  dplyr::mutate(
    atlas_layer = factor(atlas_layer, levels = allowed_layers),
    node_level = factor(node_level, levels = allowed_node_levels)
  ) |>
  dplyr::count(atlas_layer, node_level, name = "concept_nodes") |>
  dplyr::mutate(
    layer_label = dplyr::recode(
      as.character(atlas_layer),
      "primary_atlas" = "Primary\natlas",
      "treatment_context" = "Treatment\ncontext",
      "infection_context" = "Infection\ncontext",
      "sensitivity" = "Sensitivity",
      "exploratory" = "Exploratory",
      "context_only" = "Context\nonly"
    ),
    layer_label = factor(layer_label, levels = c(
      "Primary\natlas", "Treatment\ncontext", "Infection\ncontext",
      "Sensitivity", "Exploratory", "Context\nonly"
    )),
    node_label = dplyr::recode(
      as.character(node_level),
      "specific" = "Specific",
      "aggregate" = "Aggregate",
      "modifier" = "Modifier",
      "context" = "Context"
    ),
    node_label = factor(node_label, levels = c("Specific", "Aggregate", "Modifier", "Context"))
  )

p_tier <- ggplot2::ggplot(tier_df, ggplot2::aes(x = layer_label, y = concept_nodes, fill = node_label)) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::labs(
    title = "Predefined dictionary layers for analysis",
    x = "Dictionary layer", y = "Dictionary concept nodes", fill = "Node level"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    text = ggplot2::element_text(family = base_font),
    plot.title = ggplot2::element_text(face = "plain", size = 14, margin = ggplot2::margin(b = 8)),
    axis.title = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(size = 9.7, angle = 0, hjust = 0.5),
    axis.text.y = ggplot2::element_text(size = 10),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9.5),
    plot.margin = ggplot2::margin(4, 4, 4, 10)
  )

# -------------------------
# Assemble and save
# -------------------------
p_all <- p_dict / p_workflow / (p_year | p_tier) +
  patchwork::plot_layout(heights = c(1.20, 1.04, 1.05)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(family = base_font, face = "bold", size = 16),
    plot.margin = ggplot2::margin(5, 8, 5, 8)
  )

out_stub <- file.path(
  DIR_FIG2,
  sprintf("Figure1_workflow_dictionary_coverage_v6_concept_nodes_%s__%s", CORPUS_TAG, DIC_TAG)
)

ggplot2::ggsave(paste0(out_stub, ".pdf"), p_all, device = cairo_pdf_device, width = 12.8, height = 10.6)
ggplot2::ggsave(paste0(out_stub, ".png"), p_all, width = 12.8, height = 10.6, dpi = 450)

figure1_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  dic_file = basename(DIC_FILE),
  figure_version = "v6_clinical_treatment_and_infection_context_concept_nodes",
  data_change_from_previous_02b = "display-only update: Panel D y-axis and figure text use dictionary concept nodes instead of terms; no scientific graph data changed",
  n_articles = n_articles,
  n_dictionary_terms = n_terms,
  n_primary_terms = n_primary,
  n_treatment_context_terms = n_treatment,
  n_infection_context_terms = n_infection,
  n_terms_with_hit = n_hit_terms,
  articles_with_any_dictionary_hit = articles_with_hit,
  abstract_coverage_pct = coverage_pct,
  year_counts_file = ifelse(is.na(f_year_counts), NA_character_, basename(f_year_counts)),
  dictionary_run_summary_file = ifelse(file.exists(f_run_sum), basename(f_run_sum), NA_character_)
)
readr::write_csv(
  figure1_summary,
  file.path(DIR_TABLE, sprintf("Figure1_input_summary_v6_concept_nodes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
)

log_msg("WROTE:", paste0(out_stub, ".pdf"))
log_msg("WROTE:", paste0(out_stub, ".png"))
log_msg("=== DONE 02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v6_concept_nodes ===")
