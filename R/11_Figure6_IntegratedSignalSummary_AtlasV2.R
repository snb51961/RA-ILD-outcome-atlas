# =========================================================
# 11_Figure6_IntegratedSignalSummary_AtlasV2.R
# RA-ILD atlas v2.2.x manuscript-facing integrated signal summary
#
# Purpose
#   Build a data-driven integrated shared/outcome-facing signal map directly
#   from the final dictionary, PMID-level hit matrix and optional signed-effect
#   table; the script does not read, modify or patch a previous Figure 6 PDF.
#   The compact box layout follows the manuscript-facing integrated-summary style.
#
# Scientific design
#   - Uses current dictionary metadata and the v2.2 hit matrix.
#   - Excludes treatment, infection, comparator/context and measurement layers.
#   - Includes clinical/bridge concepts and biomarker/molecular concepts.
#   - Assigns terms algorithmically into shared/recurrent and outcome-facing
#     boxes from term-outcome co-mention support.
#   - Does not hard-code named results into the figure.
#   - Displays term names only for readability; full counts and provenance are
#     written to CSV tables.
#
# Required prior outputs
#   00_setup_AtlasV2.R
#   02_BuildHitsMatrix_AtlasV2.R output: hits_matrix_*.csv
#   Optional: 06_SignedEffects_AtlasV2.R output for direction-bearing support.
#
# Outputs
#   Figure6_integrated_signal_summary_FINAL_*.pdf/png
#   Figure6_integrated_signal_assignments_ALL_*.csv
#   Figure6_integrated_signal_display_TERMS_*.csv
#   Figure6_integrated_selection_rules_*.csv
#   Figure6_integrated_run_summary_*.csv
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)
quiet_install(c("dplyr", "readr", "stringr", "tibble", "tidyr", "purrr", "ggplot2", "grid"))
if (!exists("write_csv2")) write_csv2 <- function(x, path) readr::write_csv(x, path)

RUN_NAME <- "11_Figure6_IntegratedSignalSummary_AtlasV2"
log_msg("=== START ", RUN_NAME, " ===")

# -------------------------
# 0) Parameters
# -------------------------
OUTCOMES <- PRIMARY_OUTCOME_TERMS
OUTCOME_LABELS <- c("AE-ILD" = "AE-ILD", "progression" = "Progression", "mortality" = "Mortality")
OUTCOME_ORDER <- c("AE-ILD", "progression", "mortality")
OUTCOME_DISPLAY_ORDER <- c("AE-ILD", "Progression", "Mortality")

# Shared/recurrent is intentionally strict for the manuscript summary:
# only terms supported in all 3 primary outcomes are placed in the shared box.
# Terms supported in exactly 2 outcomes are retained in source tables
# (dagger_2of3) but are not marked in the compact manuscript figure.
SHARED_MIN_OUTCOMES <- as.integer(Sys.getenv("RAILD_FIG6_SHARED_MIN_OUTCOMES", unset = "3"))
SUPPORT_MIN_FOR_OUTCOME <- as.integer(Sys.getenv("RAILD_FIG6_SUPPORT_MIN_FOR_OUTCOME", unset = "1"))
MIN_DISPLAY_COHIT <- as.integer(Sys.getenv("RAILD_FIG6_MIN_DISPLAY_COHIT", unset = "2"))
MIN_AE_DISPLAY_COHIT <- as.integer(Sys.getenv("RAILD_FIG6_MIN_AE_DISPLAY_COHIT", unset = "1"))

# Display limits. These are intentionally larger than the first v2.2 Figure 6
# draft, to approximate the information density of the previous manuscript Figure 4.
N_SHARED_CLINICAL <- as.integer(Sys.getenv("RAILD_FIG6_N_SHARED_CLINICAL", unset = "10"))
N_OUTCOME_CLINICAL <- as.integer(Sys.getenv("RAILD_FIG6_N_OUTCOME_CLINICAL", unset = "8"))
N_SHARED_BIOMARKER <- as.integer(Sys.getenv("RAILD_FIG6_N_SHARED_BIOMARKER", unset = "8"))
N_OUTCOME_BIOMARKER <- as.integer(Sys.getenv("RAILD_FIG6_N_OUTCOME_BIOMARKER", unset = "7"))
SHOW_COUNTS_IN_LABEL <- tolower(Sys.getenv("RAILD_FIG6_SHOW_COUNTS", unset = "false")) %in% c("true", "1", "yes")

# For the manuscript figure, do not show dagger symbols by default.
# Two-outcome support is retained in the source tables (dagger_2of3) but
# displaying † in compact top-N boxes can be visually inconsistent when a
# two-outcome term is selected in one column but not in the other after
# per-column ranking.
DISPLAY_DAGGER_IN_FIGURE <- tolower(Sys.getenv("RAILD_FIG6_SHOW_DAGGER", unset = "false")) %in% c("true", "1", "yes")

# Outcome-adjacent terms in the AE-ILD column are usually components or synonyms
# of AE-ILD rather than independent bridge signals.  Keep them in provenance
# tables, but suppress them from the manuscript-facing AE-ILD outcome-facing box
# unless explicitly requested.
DISPLAY_AE_OUTCOME_ADJACENT <- tolower(Sys.getenv("RAILD_FIG6_SHOW_AE_OUTCOME_ADJACENT", unset = "false")) %in% c("true", "1", "yes")

FIG_WIDTH <- as.numeric(Sys.getenv("RAILD_FIG6_WIDTH_IN", unset = "13.2"))
FIG_HEIGHT <- as.numeric(Sys.getenv("RAILD_FIG6_HEIGHT_IN", unset = "7.0"))
FIG_DPI <- as.integer(Sys.getenv("RAILD_FIG6_DPI", unset = "360"))

MEASUREMENT_SUBROLES <- c("B_measurement_modality", "B_measurement_derived_imaging_feature")
MEASUREMENT_DOMAINS <- c("imaging_modality", "measurement_derived_imaging_feature")
MEASUREMENT_FAMILIES <- c("measurement_modality", "measurement_derived_imaging_feature")
EXCLUDED_ATLAS_LAYERS <- c("treatment_context", "infection_context", "context_only")
EXCLUDED_NODE_LEVELS <- c("modifier")

CARD_FILL <- tibble::tibble(
  outcome_display = c("AE-ILD", "Progression", "Mortality"),
  fill = c("#C74B50", "#E08E2B", "#4E79A7")
)

# -------------------------
# 1) Helpers
# -------------------------
safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

safe_int <- function(x, default = 0L) {
  y <- suppressWarnings(as.integer(x))
  ifelse(is.na(y), default, y)
}

safe_num <- function(x, default = NA_real_) {
  y <- suppressWarnings(as.numeric(x))
  ifelse(is.na(y), default, y)
}

first_nonempty <- function(a, b) {
  a <- safe_chr(a); b <- safe_chr(b)
  ifelse(nzchar(a), a, b)
}

hit_col <- function(term) paste0("hit__", term)

get_vec <- function(df, term) {
  nm <- hit_col(term)
  if (!nm %in% names(df)) return(integer(nrow(df)))
  as.integer(ifelse(is.na(df[[nm]]), 0L, df[[nm]] > 0L))
}

cohit_count <- function(df, term1, term2) {
  x <- get_vec(df, term1)
  y <- get_vec(df, term2)
  as.integer(sum(x == 1L & y == 1L, na.rm = TRUE))
}

find_latest_file_local <- function(dir, regex) {
  if (!dir.exists(dir)) return(NA_character_)
  xs <- list.files(dir, pattern = regex, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

read_optional_table <- function(prefix) {
  exact <- file.path(DIR_TABLE, sprintf("%s_%s__%s.csv", prefix, CORPUS_TAG, DIC_TAG))
  if (file.exists(exact)) {
    log_msg("Loaded table:", exact)
    return(readr::read_csv(exact, show_col_types = FALSE))
  }
  regex <- sprintf("^%s_.*__%s\\.csv$", prefix, DIC_TAG)
  alt <- find_latest_file_local(DIR_TABLE, regex)
  if (!is.na(alt) && file.exists(alt)) {
    log_msg("Loaded latest matching table:", alt)
    return(readr::read_csv(alt, show_col_types = FALSE))
  }
  tibble::tibble()
}

short_label <- function(x, max_chars = 34L) {
  x <- safe_chr(x)
  x <- gsub("anti-citrullinated protein antibody/rheumatoid factor", "ACPA/RF", x, ignore.case = TRUE)
  x <- gsub("anti-citrullinated protein antibody", "ACPA", x, ignore.case = TRUE)
  x <- gsub("rheumatoid factor", "RF", x, ignore.case = TRUE)
  x <- gsub("Surfactant protein", "SP", x, ignore.case = TRUE)
  x <- gsub("Neutrophil-to-lymphocyte ratio", "NLR", x, ignore.case = TRUE)
  x <- gsub("Krebs von den Lungen-6", "KL-6", x, ignore.case = TRUE)
  x <- gsub("Cytokine inflammation", "Cytokine/inflammation", x, ignore.case = TRUE)
  x <- gsub("extracellular matrix", "ECM", x, ignore.case = TRUE)
  x <- gsub("transforming growth factor", "TGF", x, ignore.case = TRUE)
  x <- gsub("matrix metalloproteinase", "MMP", x, ignore.case = TRUE)
  x <- gsub("SPP1 osteopontin", "SPP1/osteopontin", x, ignore.case = TRUE)
  x <- gsub("CRP ESR", "CRP/ESR", x, ignore.case = TRUE)
  x <- gsub("\\bMMP7\\b", "MMP-7", x, ignore.case = TRUE)
  x <- gsub("ILD onset incidence", "ILD onset/incidence", x, ignore.case = TRUE)
  x <- gsub("traction bronchiectasis", "traction bronchiectasis", x, ignore.case = TRUE)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 1), "…")
  x
}

is_measurement_node_vec <- function(abc_subrole, domain, concept_family) {
  safe_chr(abc_subrole) %in% MEASUREMENT_SUBROLES |
    safe_chr(domain) %in% MEASUREMENT_DOMAINS |
    safe_chr(concept_family) %in% MEASUREMENT_FAMILIES
}

# Outcome-adjacent terms are retained in the integrated summary, but they are
# categorised as clinical/event descriptors, not biomarkers or molecular signals.
# This preserves data-driven inclusion while preventing outcome synonyms/components
# (e.g. acute respiratory deterioration) from being displayed in biomarker rows.
is_outcome_adjacent_vec <- function(abc_role, abc_subrole, domain, concept_family) {
  role <- safe_chr(abc_role)
  sub <- safe_chr(abc_subrole)
  dom <- safe_chr(domain)
  fam <- safe_chr(concept_family)
  role == "C" |
    stringr::str_detect(sub, "^C_|outcome|AE_ILD_synonym|severe_outcome|progression_component") |
    dom %in% c("clinical_outcome", "severe_outcome", "outcome_component") |
    stringr::str_detect(fam, "AE_ILD|progression|mortality|severe_outcome|respiratory_failure|transplant|oxygen|ventilation")
}

is_biomarker_molecular_vec <- function(term, preferred_label, class, abc_role, abc_subrole, domain, concept_family, atlas_layer, node_level) {
  # Clinical outcomes and outcome-adjacent descriptors must not be classified as
  # biomarker/molecular signals even if their labels contain broad biological words.
  ifelse(
    is_outcome_adjacent_vec(abc_role, abc_subrole, domain, concept_family),
    FALSE,
    {
      blob <- stringr::str_to_lower(paste(term, preferred_label, class, abc_role, abc_subrole, domain, concept_family, atlas_layer, node_level, sep = " | "))
      dom <- safe_chr(domain)
      sub <- safe_chr(abc_subrole)
      dom %in% c("biomarker", "serological_host_factor", "molecular_pathway", "cellular_process") |
    sub %in% c("A_serological_phenotype", "B_circulating_or_tissue_biomarker", "B_cellular_process", "B_molecular_pathway") |
    stringr::str_detect(blob, paste(c(
      "biomarker", "serolog", "autoantibody", "acpa", "rheumatoid factor", "anti_ro", "anti-ro", "ssa", "ana", "seropositive",
      "chemokine", "cytokine", "interleukin", "cxcl", "ccl", "il-", "il[0-9]", "gm.?csf", "baff", "pd1", "pdl1",
      "kl-?6", "surfactant", "mmp", "y.?kl", "he4", "cc16", "epithelial", "cellular_process", "macrophage", "neutrophil", "fibroblast", "myofibroblast", "th17", "b cell", "b_cell",
      "tgf", "wnt", "ecm", "matrix", "remodel", "ctgf", "ccn2", "oxidative", "mda", "mpo", "adma", "endothelin", "spp1", "osteopontin", "periostin", "checkpoint"
    ), collapse = "|"))
    }
  )
}

assign_integrated_group <- function(term, preferred_label, class, abc_role, abc_subrole, domain, concept_family, atlas_layer, node_level, is_biomarker_molecular, is_outcome_adjacent) {
  blob <- stringr::str_to_lower(paste(term, preferred_label, class, abc_role, abc_subrole, domain, concept_family, atlas_layer, node_level, sep = " | "))
  dplyr::case_when(
    is_outcome_adjacent ~ "Outcome-adjacent / severe clinical descriptor",
    is_biomarker_molecular ~ "Biomarker / serology / molecular",
    stringr::str_detect(blob, "host|exposure|environment|smok|sex|age|duration|bmi|hypoalbumin|genetic|muc5b|tert|telomere|comorbidity|copd|eora") ~ "Host / exposure / clinical context",
    stringr::str_detect(blob, "imaging|radiology|pattern|uip|nsip|honeycomb|reticul|bronchiect|ground|dad|dip|organizing|organising|cpfe|fibrosis_extent|extent|subpleural") ~ "Imaging / pattern / extent",
    stringr::str_detect(blob, "physiology|pulmonary_function|fvc|dlco|pft|severity|index|gap|cpi|pulmonary hypertension|oxygen|respiratory failure|mechanical ventilation|icu|transplant") ~ "Physiology / severity / burden",
    TRUE ~ "Other primary / exploratory signal"
  )
}

signed_lookup <- function(signed06, terms, outcomes) {
  if (nrow(signed06) && all(c("A", "C", "articles", "balance") %in% names(signed06))) {
    signed06 |>
      dplyr::filter(A %in% terms, C %in% outcomes) |>
      dplyr::select(term = A, outcome = C, signed_articles = articles, signed_balance = balance) |>
      dplyr::mutate(signed_articles = safe_int(signed_articles, 0L), signed_balance = safe_num(signed_balance))
  } else {
    tibble::tibble(term = character(), outcome = character(), signed_articles = integer(), signed_balance = numeric())
  }
}

make_display_entry <- function(label, n, show_n = SHOW_COUNTS_IN_LABEL) {
  if (show_n) paste0(label, " (n=", n, ")") else label
}

collapse_terms <- function(x, per_line = 3, empty = "—") {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[nzchar(x)]
  if (!length(x)) return(empty)
  idx <- ceiling(seq_along(x) / per_line)
  lines <- split(x, idx)
  paste(vapply(lines, function(z) paste(z, collapse = " / "), character(1)), collapse = "\n")
}

annot_box <- function(p, x, y, label, fill = "white", colour = "black",
                      size = 3.05, fontface = "plain", box.size = 0.33,
                      padding = 0.34, hjust = 0.5, lineheight = 1.0) {
  p + ggplot2::annotate(
    "label", x = x, y = y, label = label,
    fill = fill, colour = colour, size = size, fontface = fontface,
    label.size = box.size, label.padding = grid::unit(padding, "lines"),
    hjust = hjust, lineheight = lineheight
  )
}

# -------------------------
# 2) Load data
# -------------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)
for (col in c("preferred_label", "class", "abc_role", "abc_subrole", "domain", "concept_family", "atlas_layer", "node_level", "analysis_tier", "include_primary", "primary_atlas_eligible")) {
  if (!col %in% names(dic)) dic[[col]] <- ""
}
DF <- load_hit_matrix()
if (!"pmid" %in% names(DF)) stop("Hit matrix missing pmid column")
DF$pmid <- as.character(DF$pmid)
signed06 <- read_optional_table("signed_effects_summary_primary_atlas")

# -------------------------
# 3) Current-policy integrated scope
# -------------------------
dic_scope <- dic |>
  dplyr::mutate(
    preferred_label_clean = first_nonempty(preferred_label, term),
    is_measurement = is_measurement_node_vec(abc_subrole, domain, concept_family),
    is_outcome_adjacent = is_outcome_adjacent_vec(abc_role, abc_subrole, domain, concept_family),
    is_biomarker_molecular = is_biomarker_molecular_vec(term, preferred_label, class, abc_role, abc_subrole, domain, concept_family, atlas_layer, node_level),
    integrated_group = assign_integrated_group(term, preferred_label, class, abc_role, abc_subrole, domain, concept_family, atlas_layer, node_level, is_biomarker_molecular, is_outcome_adjacent),
    display_label = short_label(preferred_label_clean, 34),
    include_primary_bool = tolower(safe_chr(include_primary)) %in% c("true", "1", "yes"),
    primary_eligible_bool = tolower(safe_chr(primary_atlas_eligible)) %in% c("true", "1", "yes")
  ) |>
  dplyr::mutate(
    in_integrated_scope = (
      (atlas_layer %in% c("primary_atlas", "exploratory", "sensitivity")) &
        !(atlas_layer %in% EXCLUDED_ATLAS_LAYERS) &
        !(term %in% OUTCOMES) &
        !(is_measurement) &
        !(safe_chr(abc_subrole) %in% c("A_treatment_context", "context_population_or_comparator", "B_measurement_modality", "B_measurement_derived_imaging_feature", "C_treatment_safety_context", "infection_safety_context", "infection_respiratory_context", "AE_ILD_infection_context", "modifier")) &
        !(safe_chr(domain) %in% c("treatment_context", "infection_safety_context", "context_population", "imaging_modality", "measurement_derived_imaging_feature", "modifier")) &
        !(safe_chr(node_level) %in% EXCLUDED_NODE_LEVELS)
    ),
    broad_panel = dplyr::if_else(is_biomarker_molecular & !is_outcome_adjacent, "Biomarker / molecular signals", "Clinical / bridge signals")
  )

scope_qc <- dic_scope |>
  dplyr::mutate(scope_flag = dplyr::case_when(
    in_integrated_scope ~ "included: integrated summary scope",
    atlas_layer == "treatment_context" ~ "excluded: treatment context",
    atlas_layer == "infection_context" ~ "excluded: infection context",
    atlas_layer == "context_only" ~ "excluded: comparator/context",
    term %in% OUTCOMES ~ "excluded: primary outcome node",
    is_measurement ~ "excluded: measurement/modality",
    TRUE ~ "excluded: outside integrated summary scope"
  )) |>
  dplyr::select(term, preferred_label, abc_role, abc_subrole, domain, concept_family, node_level, analysis_tier, atlas_layer, integrated_group, broad_panel, is_biomarker_molecular, is_outcome_adjacent, in_integrated_scope, scope_flag)
write_csv2(scope_qc, file.path(DIR_TABLE, sprintf("Figure6_integrated_scope_QC_all_terms_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

integrated_terms <- dic_scope |>
  dplyr::filter(in_integrated_scope) |>
  dplyr::filter(hit_col(term) %in% names(DF))
if (!nrow(integrated_terms)) stop("No integrated summary terms selected. Check dictionary metadata.")
write_csv2(integrated_terms, file.path(DIR_TABLE, sprintf("Figure6_integrated_scope_INCLUDED_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

# -------------------------
# 4) Term-outcome evidence
# -------------------------
term_outcome <- tidyr::expand_grid(term = integrated_terms$term, outcome = OUTCOMES) |>
  dplyr::mutate(cohit_n = purrr::map2_int(term, outcome, ~ cohit_count(DF, .x, .y))) |>
  dplyr::left_join(signed_lookup(signed06, integrated_terms$term, OUTCOMES), by = c("term", "outcome")) |>
  dplyr::mutate(
    signed_articles = safe_int(signed_articles, 0L),
    signed_balance = safe_num(signed_balance),
    direction_available = signed_articles > 0 & !is.na(signed_balance),
    outcome_display = OUTCOME_LABELS[outcome],
    outcome_display = factor(outcome_display, levels = OUTCOME_DISPLAY_ORDER)
  )
write_csv2(term_outcome, file.path(DIR_TABLE, sprintf("Figure6_integrated_term_outcome_evidence_ALL_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

dominant <- term_outcome |>
  dplyr::group_by(term) |>
  dplyr::arrange(dplyr::desc(cohit_n), dplyr::desc(signed_articles), factor(outcome, levels = OUTCOMES), .by_group = TRUE) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::select(term, dominant_outcome = outcome, dominant_outcome_display = outcome_display, max_outcome_cohit = cohit_n, dominant_signed_articles = signed_articles)

term_summary <- integrated_terms |>
  dplyr::select(term, preferred_label, display_label, broad_panel, integrated_group, abc_role, abc_subrole, domain, concept_family, node_level, analysis_tier, atlas_layer, is_biomarker_molecular, is_outcome_adjacent) |>
  dplyr::left_join(
    term_outcome |>
      dplyr::group_by(term) |>
      dplyr::summarise(
        total_cohit = sum(cohit_n, na.rm = TRUE),
        outcomes_supported = sum(cohit_n >= SUPPORT_MIN_FOR_OUTCOME, na.rm = TRUE),
        signed_articles_total = sum(signed_articles, na.rm = TRUE),
        mean_abs_balance = mean(abs(signed_balance[direction_available]), na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(mean_abs_balance = ifelse(is.nan(mean_abs_balance), 0, mean_abs_balance)),
    by = "term"
  ) |>
  dplyr::left_join(dominant, by = "term") |>
  dplyr::mutate(
    total_cohit = safe_int(total_cohit, 0L),
    max_outcome_cohit = safe_int(max_outcome_cohit, 0L),
    outcomes_supported = safe_int(outcomes_supported, 0L),
    signed_articles_total = safe_int(signed_articles_total, 0L),
    mean_abs_balance = safe_num(mean_abs_balance, 0),
    # Data-derived display score; no named anchor terms.
    display_score = log1p(total_cohit) + 0.75 * outcomes_supported + 0.20 * log1p(signed_articles_total) + 0.08 * mean_abs_balance
  ) |>
  dplyr::arrange(broad_panel, dplyr::desc(display_score), dplyr::desc(max_outcome_cohit), display_label)


# -------------------------
# 4b) Final manuscript assignment: all-three backbone, progression-mortality axis,
#     and progression-/mortality-facing signals
# -------------------------
# Design rationale
#   The final Figure 6 is intended as a compact clinical summary rather than a
#   three-column repetition of every non-shared term.  Because the AE-ILD row has
#   no robust non-shared displayed signal after suppressing outcome-adjacent AE
#   descriptors, the display is reorganised into:
#     1) shared across AE-ILD, progression and mortality;
#     2) shared/balanced progression-mortality axis without AE-ILD support;
#     3) progression-facing and mortality-facing enriched signals.
#   All assignments are computed from term-outcome support in the hit matrix.
#   The script does not read or edit any previous Figure 6 PDF.

PM_BALANCE_LOW <- as.numeric(Sys.getenv("RAILD_FIG6_PM_BALANCE_LOW", unset = "0.5"))
PM_BALANCE_HIGH <- as.numeric(Sys.getenv("RAILD_FIG6_PM_BALANCE_HIGH", unset = "2.0"))
MIN_PM_SHARED_COHIT <- as.integer(Sys.getenv("RAILD_FIG6_MIN_PM_SHARED_COHIT", unset = "2"))
MIN_FACING_COHIT <- as.integer(Sys.getenv("RAILD_FIG6_MIN_FACING_COHIT", unset = "2"))
DISPLAY_AE_NONSHARED <- tolower(Sys.getenv("RAILD_FIG6_DISPLAY_AE_NONSHARED", unset = "false")) %in% c("true", "1", "yes")

N_PM_SHARED_CLINICAL <- as.integer(Sys.getenv("RAILD_FIG6_N_PM_SHARED_CLINICAL", unset = "8"))
N_PROGRESSION_FACING_CLINICAL <- as.integer(Sys.getenv("RAILD_FIG6_N_PROGRESSION_FACING_CLINICAL", unset = "6"))
N_MORTALITY_FACING_CLINICAL <- as.integer(Sys.getenv("RAILD_FIG6_N_MORTALITY_FACING_CLINICAL", unset = "6"))
N_PM_SHARED_BIOMARKER <- as.integer(Sys.getenv("RAILD_FIG6_N_PM_SHARED_BIOMARKER", unset = "8"))
N_PROGRESSION_FACING_BIOMARKER <- as.integer(Sys.getenv("RAILD_FIG6_N_PROGRESSION_FACING_BIOMARKER", unset = "6"))
N_MORTALITY_FACING_BIOMARKER <- as.integer(Sys.getenv("RAILD_FIG6_N_MORTALITY_FACING_BIOMARKER", unset = "6"))

support_pattern <- term_outcome |>
  dplyr::select(term, outcome, cohit_n) |>
  tidyr::pivot_wider(names_from = outcome, values_from = cohit_n, values_fill = 0) |>
  dplyr::rename(
    ae_cohit = `AE-ILD`,
    progression_cohit = progression,
    mortality_cohit = mortality
  ) |>
  dplyr::mutate(
    has_AE_ILD = ae_cohit >= SUPPORT_MIN_FOR_OUTCOME,
    has_progression = progression_cohit >= SUPPORT_MIN_FOR_OUTCOME,
    has_mortality = mortality_cohit >= SUPPORT_MIN_FOR_OUTCOME,
    pm_ratio = (progression_cohit + 0.5) / (mortality_cohit + 0.5),
    pm_ratio_inverse = (mortality_cohit + 0.5) / (progression_cohit + 0.5)
  )

term_summary2 <- term_summary |>
  dplyr::left_join(support_pattern, by = "term") |>
  dplyr::mutate(
    ae_cohit = safe_int(ae_cohit, 0L),
    progression_cohit = safe_int(progression_cohit, 0L),
    mortality_cohit = safe_int(mortality_cohit, 0L),
    has_AE_ILD = ifelse(is.na(has_AE_ILD), FALSE, has_AE_ILD),
    has_progression = ifelse(is.na(has_progression), FALSE, has_progression),
    has_mortality = ifelse(is.na(has_mortality), FALSE, has_mortality),
    pm_ratio = safe_num(pm_ratio, 1),
    pm_ratio_inverse = safe_num(pm_ratio_inverse, 1),
    support_pattern = dplyr::case_when(
      has_AE_ILD & has_progression & has_mortality ~ "AE-ILD+progression+mortality",
      !has_AE_ILD & has_progression & has_mortality ~ "progression+mortality",
      !has_AE_ILD & has_progression & !has_mortality ~ "progression only",
      !has_AE_ILD & !has_progression & has_mortality ~ "mortality only",
      has_AE_ILD & has_progression & !has_mortality ~ "AE-ILD+progression",
      has_AE_ILD & !has_progression & has_mortality ~ "AE-ILD+mortality",
      has_AE_ILD & !has_progression & !has_mortality ~ "AE-ILD only",
      TRUE ~ "no primary-outcome support"
    ),
    final_display_layer = dplyr::case_when(
      has_AE_ILD & has_progression & has_mortality & total_cohit >= MIN_DISPLAY_COHIT ~ "Shared across all three outcomes",
      !has_AE_ILD & has_progression & has_mortality &
        pm_ratio >= PM_BALANCE_LOW & pm_ratio <= PM_BALANCE_HIGH &
        (progression_cohit + mortality_cohit) >= MIN_PM_SHARED_COHIT ~ "Progression-mortality shared axis",
      !has_AE_ILD & has_progression &
        (!has_mortality | pm_ratio > PM_BALANCE_HIGH) &
        progression_cohit >= MIN_FACING_COHIT ~ "Progression-facing signals",
      !has_AE_ILD & has_mortality &
        (!has_progression | pm_ratio < PM_BALANCE_LOW) &
        mortality_cohit >= MIN_FACING_COHIT ~ "Mortality-facing signals",
      DISPLAY_AE_NONSHARED & has_AE_ILD & !has_progression & !has_mortality &
        ae_cohit >= MIN_AE_DISPLAY_COHIT & !is_outcome_adjacent ~ "AE-ILD-facing signals",
      TRUE ~ "Retained in source tables only"
    ),
    assigned_cohit = dplyr::case_when(
      final_display_layer == "Shared across all three outcomes" ~ max_outcome_cohit,
      final_display_layer == "Progression-mortality shared axis" ~ pmin(progression_cohit, mortality_cohit),
      final_display_layer == "Progression-facing signals" ~ progression_cohit,
      final_display_layer == "Mortality-facing signals" ~ mortality_cohit,
      final_display_layer == "AE-ILD-facing signals" ~ ae_cohit,
      TRUE ~ max_outcome_cohit
    ),
    facing_enrichment_ratio = dplyr::case_when(
      final_display_layer == "Progression-facing signals" ~ pm_ratio,
      final_display_layer == "Mortality-facing signals" ~ pm_ratio_inverse,
      TRUE ~ NA_real_
    ),
    is_progression_mortality_balanced = final_display_layer == "Progression-mortality shared axis",
    is_progression_enriched = final_display_layer == "Progression-facing signals",
    is_mortality_enriched = final_display_layer == "Mortality-facing signals",
    can_display = final_display_layer != "Retained in source tables only",
    display_marker = dplyr::if_else(is_outcome_adjacent, "‡", ""),
    display_entry = dplyr::case_when(
      SHOW_COUNTS_IN_LABEL ~ paste0(display_label, display_marker, " (n=", assigned_cohit, ")"),
      TRUE ~ paste0(display_label, display_marker)
    ),
    final_display_layer = factor(
      final_display_layer,
      levels = c(
        "Shared across all three outcomes",
        "Progression-mortality shared axis",
        "Progression-facing signals",
        "Mortality-facing signals",
        "AE-ILD-facing signals",
        "Retained in source tables only"
      )
    )
  ) |>
  dplyr::arrange(broad_panel, final_display_layer, dplyr::desc(display_score), dplyr::desc(assigned_cohit), display_label)

term_assignments <- term_summary2 |>
  dplyr::mutate(
    assignment_rule = dplyr::case_when(
      final_display_layer == "Shared across all three outcomes" ~ "support in AE-ILD, progression and mortality; placed in all-three shared backbone",
      final_display_layer == "Progression-mortality shared axis" ~ paste0("support in progression and mortality without AE-ILD support; progression/mortality support ratio within [", PM_BALANCE_LOW, ", ", PM_BALANCE_HIGH, "]"),
      final_display_layer == "Progression-facing signals" ~ paste0("progression-only or progression-enriched support; progression/mortality support ratio > ", PM_BALANCE_HIGH, " or mortality support absent"),
      final_display_layer == "Mortality-facing signals" ~ paste0("mortality-only or mortality-enriched support; progression/mortality support ratio < ", PM_BALANCE_LOW, " or progression support absent"),
      final_display_layer == "AE-ILD-facing signals" ~ "AE-ILD-only non-outcome-adjacent support; optional display disabled by default",
      TRUE ~ "retained in source/provenance tables but not displayed in compact manuscript figure"
    )
  )

write_csv2(term_summary2, file.path(DIR_TABLE, sprintf("Figure6_integrated_signal_term_summary_ALL_%s__%s.csv", CORPUS_TAG, DIC_TAG)))
write_csv2(term_assignments, file.path(DIR_TABLE, sprintf("Figure6_integrated_signal_assignments_ALL_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

# -------------------------
# 5) Build display text maps
# -------------------------
select_terms <- function(panel, layer, n, per_line = 3) {
  term_assignments |>
    dplyr::filter(broad_panel == panel, final_display_layer == layer, can_display) |>
    dplyr::arrange(dplyr::desc(display_score), dplyr::desc(assigned_cohit), dplyr::desc(total_cohit), display_label) |>
    dplyr::slice_head(n = n)
}

clinical_all3 <- select_terms("Clinical / bridge signals", "Shared across all three outcomes", N_SHARED_CLINICAL)
clinical_pm_shared <- select_terms("Clinical / bridge signals", "Progression-mortality shared axis", N_PM_SHARED_CLINICAL)
clinical_P_face <- select_terms("Clinical / bridge signals", "Progression-facing signals", N_PROGRESSION_FACING_CLINICAL)
clinical_M_face <- select_terms("Clinical / bridge signals", "Mortality-facing signals", N_MORTALITY_FACING_CLINICAL)

bio_all3 <- select_terms("Biomarker / molecular signals", "Shared across all three outcomes", N_SHARED_BIOMARKER)
bio_pm_shared <- select_terms("Biomarker / molecular signals", "Progression-mortality shared axis", N_PM_SHARED_BIOMARKER)
bio_P_face <- select_terms("Biomarker / molecular signals", "Progression-facing signals", N_PROGRESSION_FACING_BIOMARKER)
bio_M_face <- select_terms("Biomarker / molecular signals", "Mortality-facing signals", N_MORTALITY_FACING_BIOMARKER)

display_terms <- dplyr::bind_rows(
  clinical_all3 |> dplyr::mutate(display_box = "Clinical shared across all three"),
  clinical_pm_shared |> dplyr::mutate(display_box = "Clinical progression-mortality shared"),
  clinical_P_face |> dplyr::mutate(display_box = "Clinical progression-facing"),
  clinical_M_face |> dplyr::mutate(display_box = "Clinical mortality-facing"),
  bio_all3 |> dplyr::mutate(display_box = "Biomarker shared across all three"),
  bio_pm_shared |> dplyr::mutate(display_box = "Biomarker progression-mortality shared"),
  bio_P_face |> dplyr::mutate(display_box = "Biomarker progression-facing"),
  bio_M_face |> dplyr::mutate(display_box = "Biomarker mortality-facing")
)
write_csv2(display_terms, file.path(DIR_TABLE, sprintf("Figure6_integrated_signal_display_TERMS_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

clinical_all3_text <- collapse_terms(clinical_all3$display_entry, per_line = 3)
clinical_pm_shared_text <- collapse_terms(clinical_pm_shared$display_entry, per_line = 2)
clinical_P_face_text <- collapse_terms(clinical_P_face$display_entry, per_line = 2)
clinical_M_face_text <- collapse_terms(clinical_M_face$display_entry, per_line = 2)

bio_all3_text <- collapse_terms(bio_all3$display_entry, per_line = 3)
bio_pm_shared_text <- collapse_terms(bio_pm_shared$display_entry, per_line = 2)
bio_P_face_text <- collapse_terms(bio_P_face$display_entry, per_line = 2)
bio_M_face_text <- collapse_terms(bio_M_face$display_entry, per_line = 2)

# -------------------------
# 6) Figure canvas: all-three backbone spanning AE-ILD/progression/mortality,
#    progression-mortality axis and facing signals
# -------------------------
p <- ggplot2::ggplot() +
  ggplot2::xlim(0, 100) +
  ggplot2::ylim(0, 100) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(8, 10, 12, 10))

# Layout anchors. The all-three shared box is centred across the AE-ILD,
# progression and mortality columns.  The progression-mortality shared axis is
# centred only across the progression and mortality columns.  AE-ILD-facing boxes
# are shown with the same visual grammar as the progression- and mortality-facing
# boxes, but display a dash when no non-shared AE-ILD signal meets the rule.
x_ae <- 27
x_prog <- 55
x_mort <- 83
x_all3 <- mean(c(x_ae, x_mort))
x_pm <- mean(c(x_prog, x_mort))

# Light span indicators clarify what each shared box covers without changing the
# data-driven selection.  These are graphical guides only.
add_span_indicator <- function(p, x1, x2, y, tick = 0.9, colour = "grey65", linewidth = 0.35) {
  p +
    ggplot2::annotate("segment", x = x1, xend = x2, y = y, yend = y, colour = colour, linewidth = linewidth) +
    ggplot2::annotate("segment", x = x1, xend = x1, y = y, yend = y - tick, colour = colour, linewidth = linewidth) +
    ggplot2::annotate("segment", x = x2, xend = x2, y = y, yend = y - tick, colour = colour, linewidth = linewidth)
}

# Title and subtitle
p <- p +
  ggplot2::annotate("text", x = 5, y = 97.3, hjust = 0, label = "Figure 6. Integrated shared and outcome-facing RA-ILD worsening signals", size = 4.6, fontface = "bold") +
  ggplot2::annotate("text", x = 5, y = 94.1, hjust = 0, label = "Algorithmic summary from final dictionary and PMID-level hit matrix; treatment, infection and measurement layers excluded", size = 2.75, colour = "grey25")

# Section labels and separators
p <- p +
  ggplot2::annotate("text", x = 5, y = 85.6, label = "Clinical / bridge\nsignals", hjust = 0, size = 4.1, fontface = "bold") +
  ggplot2::annotate("segment", x = 5, xend = 95, y = 47.4, yend = 47.4, colour = "grey76", linewidth = 0.45) +
  ggplot2::annotate("text", x = 5, y = 39.8, label = "Biomarker / molecular\nsignals", hjust = 0, size = 4.1, fontface = "bold")

# Headers
all3_header <- "Shared across AE-ILD, progression, and mortality"
pm_header <- "Shared / balanced progression–mortality axis"

# Clinical / bridge half: all-three shared backbone spanning all three outcomes
p <- annot_box(p, x_all3, 88.7, all3_header, fill = "#55565B", colour = "white", size = 3.35, fontface = "bold", box.size = 0, padding = 0.22)
p <- add_span_indicator(p, x_ae, x_mort, 86.7)
p <- annot_box(p, x_all3, 81.8, clinical_all3_text, fill = "white", colour = "black", size = 2.62, padding = 0.30, lineheight = 1.00)

# Clinical / bridge half: progression-mortality shared axis spanning only the two chronic-outcome columns
p <- annot_box(p, x_pm, 72.1, pm_header, fill = "#77797F", colour = "white", size = 3.18, fontface = "bold", box.size = 0, padding = 0.22)
p <- add_span_indicator(p, x_prog, x_mort, 70.2)
p <- annot_box(p, x_pm, 64.6, clinical_pm_shared_text, fill = "white", colour = "black", size = 2.52, padding = 0.30, lineheight = 1.00)

# Outcome-facing clinical / bridge signals
p <- annot_box(p, x_ae, 55.6, "AE-ILD-facing", fill = "#C74B50", colour = "white", size = 3.35, fontface = "bold", box.size = 0, padding = 0.22)
p <- annot_box(p, x_prog, 55.6, "Progression-facing", fill = "#E08E2B", colour = "white", size = 3.35, fontface = "bold", box.size = 0, padding = 0.22)
p <- annot_box(p, x_mort, 55.6, "Mortality-facing", fill = "#4E79A7", colour = "white", size = 3.35, fontface = "bold", box.size = 0, padding = 0.22)
p <- annot_box(p, x_ae, 49.0, "—", fill = "white", colour = "black", size = 3.0, padding = 0.28, lineheight = 1.00)
p <- annot_box(p, x_prog, 49.0, clinical_P_face_text, fill = "white", colour = "black", size = 2.40, padding = 0.28, lineheight = 1.00)
p <- annot_box(p, x_mort, 49.0, clinical_M_face_text, fill = "white", colour = "black", size = 2.40, padding = 0.28, lineheight = 1.00)

# Biomarker/molecular half: all-three shared backbone spanning all three outcomes
p <- annot_box(p, x_all3, 40.8, all3_header, fill = "#55565B", colour = "white", size = 3.20, fontface = "bold", box.size = 0, padding = 0.20)
p <- add_span_indicator(p, x_ae, x_mort, 38.9)
p <- annot_box(p, x_all3, 34.4, bio_all3_text, fill = "white", colour = "black", size = 2.45, padding = 0.28, lineheight = 0.98)

# Biomarker/molecular half: progression-mortality shared axis spanning only progression/mortality
p <- annot_box(p, x_pm, 26.7, pm_header, fill = "#77797F", colour = "white", size = 3.05, fontface = "bold", box.size = 0, padding = 0.20)
p <- add_span_indicator(p, x_prog, x_mort, 25.0)
p <- annot_box(p, x_pm, 20.2, bio_pm_shared_text, fill = "white", colour = "black", size = 2.38, padding = 0.28, lineheight = 0.98)

# Outcome-facing biomarker / molecular signals
p <- annot_box(p, x_ae, 13.3, "AE-ILD-facing", fill = "#C74B50", colour = "white", size = 3.15, fontface = "bold", box.size = 0, padding = 0.20)
p <- annot_box(p, x_prog, 13.3, "Progression-facing", fill = "#E08E2B", colour = "white", size = 3.15, fontface = "bold", box.size = 0, padding = 0.20)
p <- annot_box(p, x_mort, 13.3, "Mortality-facing", fill = "#4E79A7", colour = "white", size = 3.15, fontface = "bold", box.size = 0, padding = 0.20)
p <- annot_box(p, x_ae, 7.3, "—", fill = "white", colour = "black", size = 2.85, padding = 0.26, lineheight = 0.96)
p <- annot_box(p, x_prog, 7.3, bio_P_face_text, fill = "white", colour = "black", size = 2.24, padding = 0.26, lineheight = 0.96)
p <- annot_box(p, x_mort, 7.3, bio_M_face_text, fill = "white", colour = "black", size = 2.24, padding = 0.26, lineheight = 0.96)

foot <- "Displayed concept nodes are selected algorithmically. All-three shared = support detected in AE-ILD, progression and mortality. Progression–mortality axis = support in both progression and mortality without AE-ILD support and without strong outcome dominance. Progression-/mortality-facing = one-outcome or outcome-enriched support. ‡ = outcome-adjacent descriptor/component, not interpreted as a biomarker or upstream bridge."
p <- p + ggplot2::annotate("text", x = 50, y = 2.15, label = foot, hjust = 0.5, size = 2.05, colour = "grey35")

# -------------------------
# 7) Outputs
# -------------------------
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)

stub <- file.path(DIR_FIG, sprintf("Figure6_integrated_signal_summary_FINAL_%s__%s", CORPUS_TAG, DIC_TAG))
fig_pdf <- paste0(stub, ".pdf")
fig_png <- paste0(stub, ".png")

ggplot2::ggsave(fig_pdf, p, width = FIG_WIDTH, height = max(FIG_HEIGHT, 7.4), device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(fig_png, p, width = FIG_WIDTH, height = max(FIG_HEIGHT, 7.4), dpi = FIG_DPI, bg = "white")

rules <- tibble::tribble(
  ~parameter, ~value,
  "scope", "Final dictionary terms excluding treatment_context, infection_context, context_only, measurement/modality, modifiers and primary outcome nodes",
  "all_three_shared_definition", "Term supported in AE-ILD, progression and mortality",
  "progression_mortality_axis_definition", "Term supported in progression and mortality, not supported in AE-ILD, and progression/mortality support ratio within configured balanced range",
  "progression_facing_definition", "Term supported in progression and not AE-ILD, with mortality support absent or progression/mortality support ratio above configured high threshold",
  "mortality_facing_definition", "Term supported in mortality and not AE-ILD, with progression support absent or progression/mortality support ratio below configured low threshold",
  "AE_nonshared_display", as.character(DISPLAY_AE_NONSHARED),
  "PM_BALANCE_LOW", as.character(PM_BALANCE_LOW),
  "PM_BALANCE_HIGH", as.character(PM_BALANCE_HIGH),
  "MIN_PM_SHARED_COHIT", as.character(MIN_PM_SHARED_COHIT),
  "MIN_FACING_COHIT", as.character(MIN_FACING_COHIT),
  "display_dagger_in_figure", "FALSE; two-outcome support patterns are encoded in source tables rather than by dagger symbols",
  "display_AE_outcome_adjacent", as.character(DISPLAY_AE_OUTCOME_ADJACENT),
  "display_score", "log1p(total_cohit) + 0.75*outcomes_supported + 0.20*log1p(signed_articles_total) + 0.08*mean_abs_balance",
  "manual_named_term_forcing", "none",
  "show_counts_in_labels", as.character(SHOW_COUNTS_IN_LABEL),
  "clinical_all_three_display_limit", as.character(N_SHARED_CLINICAL),
  "clinical_progression_mortality_shared_display_limit", as.character(N_PM_SHARED_CLINICAL),
  "clinical_progression_facing_display_limit", as.character(N_PROGRESSION_FACING_CLINICAL),
  "clinical_mortality_facing_display_limit", as.character(N_MORTALITY_FACING_CLINICAL),
  "biomarker_all_three_display_limit", as.character(N_SHARED_BIOMARKER),
  "biomarker_progression_mortality_shared_display_limit", as.character(N_PM_SHARED_BIOMARKER),
  "biomarker_progression_facing_display_limit", as.character(N_PROGRESSION_FACING_BIOMARKER),
  "biomarker_mortality_facing_display_limit", as.character(N_MORTALITY_FACING_BIOMARKER)
)
write_csv2(rules, file.path(DIR_TABLE, sprintf("Figure6_integrated_selection_rules_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

run_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  n_integrated_terms = nrow(integrated_terms),
  n_terms_with_any_support = sum(term_summary2$total_cohit > 0, na.rm = TRUE),
  n_all_three_shared_terms = sum(term_summary2$final_display_layer == "Shared across all three outcomes", na.rm = TRUE),
  n_progression_mortality_axis_terms = sum(term_summary2$final_display_layer == "Progression-mortality shared axis", na.rm = TRUE),
  n_progression_facing_terms = sum(term_summary2$final_display_layer == "Progression-facing signals", na.rm = TRUE),
  n_mortality_facing_terms = sum(term_summary2$final_display_layer == "Mortality-facing signals", na.rm = TRUE),
  n_source_only_terms = sum(term_summary2$final_display_layer == "Retained in source tables only", na.rm = TRUE),
  n_clinical_terms_displayed = sum(display_terms$broad_panel == "Clinical / bridge signals"),
  n_biomarker_terms_displayed = sum(display_terms$broad_panel == "Biomarker / molecular signals"),
  figure_pdf = fig_pdf,
  figure_png = fig_png
)
write_csv2(run_summary, file.path(DIR_TABLE, sprintf("Figure6_integrated_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

log_msg("WROTE:", fig_pdf)
log_msg("WROTE:", fig_png)
log_msg("=== DONE ", RUN_NAME, " ===")
