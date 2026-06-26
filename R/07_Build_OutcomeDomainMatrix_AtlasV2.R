# =========================================================
# 07_Build_OutcomeDomainMatrix_AtlasV2.R
# RA-ILD atlas v2.2.2
#
# Purpose
#   Build the final Figure 2 outcome-domain atlas directly from the outputs of
#   03-06, without requiring any previous 07 output.
#
# Scientific policy
#   - Primary disease-state atlas is separated from treatment-context and
#     infection-context layers.
#   - Measurement/modality nodes such as HRCT, quantitative_CT and VRS are not
#     displayed as primary disease-state bridge nodes.
#   - Figure labels are shortened for readability only. Full top-node rankings
#     are written to tables.
#   - Bubble size is based on unique article-level domain-outcome co-mentions
#     when the hit matrix is available. Hit-matrix columns are resolved from the
#     02 convention hit__<term>. Summed term-outcome co-mention counts are
#     retained in tables as a secondary support metric.
#
# Required prior scripts/outputs
#   00_setup_AtlasV2.R
#   02_BuildHitsMatrix_AtlasV2.R
#   03_CoocAndCollocation_AtlasV2.R
#   04_ABC_Rankings_AtlasV2.R
#   05_AC_NPMI_AtlasV2.R
#   06_SignedEffects_AtlasV2.R
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)
quiet_install(c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr", "ggplot2"))

RUN_NAME <- "07_Build_OutcomeDomainMatrix_AtlasV2"
log_msg("=== START ", RUN_NAME, " ===")

# -------------------------
# 0) Display constants
# -------------------------
# Figure labels default to 3 nodes per cell for readability. This is not an
# analytic filter; full ranked nodes are saved in outcome_domain_matrix_top_terms.
TOP_N_LABEL <- as.integer(Sys.getenv("RAILD_FIG2_LABEL_N", unset = "3"))
TOP_N_TABLE <- as.integer(Sys.getenv("RAILD_FIG2_TABLE_TOP_N", unset = "5"))
MIN_COHIT_FOR_LABEL <- as.integer(Sys.getenv("RAILD_FIG2_MIN_COHIT_FOR_LABEL", unset = "1"))
FIG_WIDTH <- as.numeric(Sys.getenv("RAILD_FIG2_WIDTH_IN", unset = "11.2"))
FIG_HEIGHT <- as.numeric(Sys.getenv("RAILD_FIG2_HEIGHT_IN", unset = "7.4"))
FIG_DPI <- as.integer(Sys.getenv("RAILD_FIG2_DPI", unset = "360"))

OUTCOME_ORDER <- c("AE-ILD", "progression", "mortality")
OUTCOME_LABELS <- c(
  "AE-ILD" = "AE-ILD",
  "progression" = "Progression",
  "mortality" = "Mortality"
)

DOMAIN_ORDER <- c(
  "Host/exposure/genetic/serology",
  "Imaging pattern/extent",
  "Physiology/severity",
  "Biomarker/pathway",
  "Treatment context",
  "Infection context"
)
DOMAIN_Y <- setNames(rev(seq_along(DOMAIN_ORDER)), DOMAIN_ORDER)
OUTCOME_X <- setNames(seq_along(OUTCOME_ORDER), OUTCOME_ORDER)

# These are measurement/modality nodes. They remain in the dictionary but are not
# displayed as primary disease-state bridges.
MEASUREMENT_GUARD_TERMS <- c("HRCT", "quantitative_CT", "VRS")
MEASUREMENT_GUARD_SUBROLES <- c("B_measurement_modality", "B_measurement_derived_imaging_feature")
MEASUREMENT_GUARD_DOMAINS <- c("imaging_modality", "measurement_derived_imaging_feature")

# -------------------------
# 1) Helpers
# -------------------------
get_col <- function(df, nm, default = NA) {
  n <- nrow(df)
  if (nm %in% names(df)) return(df[[nm]])
  rep(default, n)
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

read_required_table <- function(prefix) {
  exact <- file.path(DIR_TABLE, sprintf("%s_%s__%s.csv", prefix, CORPUS_TAG, DIC_TAG))
  if (file.exists(exact)) {
    log_msg("Loaded table:", exact)
    return(readr::read_csv(exact, show_col_types = FALSE))
  }
  regex <- sprintf("^%s_.*__%s\\.csv$", prefix, DIC_TAG)
  alt <- find_latest_file(DIR_TABLE, regex)
  if (!is.na(alt) && file.exists(alt)) {
    log_msg("Loaded latest matching table:", alt)
    return(readr::read_csv(alt, show_col_types = FALSE))
  }
  stop("Required table not found for prefix: ", prefix, "\nExpected: ", exact,
       "\nRun 03-06 before 07.")
}

read_optional_table <- function(prefix) {
  exact <- file.path(DIR_TABLE, sprintf("%s_%s__%s.csv", prefix, CORPUS_TAG, DIC_TAG))
  if (file.exists(exact)) {
    log_msg("Loaded table:", exact)
    return(readr::read_csv(exact, show_col_types = FALSE))
  }
  regex <- sprintf("^%s_.*__%s\\.csv$", prefix, DIC_TAG)
  alt <- find_latest_file(DIR_TABLE, regex)
  if (!is.na(alt) && file.exists(alt)) {
    log_msg("Loaded latest matching table:", alt)
    return(readr::read_csv(alt, show_col_types = FALSE))
  }
  log_msg("Optional table not found; continuing with empty table for prefix:", prefix)
  tibble::tibble()
}

standardize_cooc <- function(tab, default_scope, source_table) {
  if (!nrow(tab)) {
    return(tibble::tibble(
      source_table = character(), scope = character(), term = character(), outcome = character(),
      cohit_n = integer(), n10 = integer(), n01 = integer(), n00 = integer(),
      p_x = double(), p_y = double(), p_xy = double(), or = double(), p = double(),
      lift = double(), npmi = double(), q = double()
    ))
  }
  term <- if ("A" %in% names(tab)) get_col(tab, "A") else get_col(tab, "term")
  outcome <- if ("C" %in% names(tab)) get_col(tab, "C") else get_col(tab, "outcome")
  scope <- get_col(tab, "scope", default_scope)
  tibble::tibble(
    source_table = source_table,
    scope = as.character(scope),
    term = as.character(term),
    outcome = as.character(outcome),
    cohit_n = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "n11", 0)), 0)),
    n10 = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "n10", 0)), 0)),
    n01 = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "n01", 0)), 0)),
    n00 = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "n00", 0)), 0)),
    p_x = as.numeric(get_col(tab, "p_x", NA_real_)),
    p_y = as.numeric(get_col(tab, "p_y", NA_real_)),
    p_xy = as.numeric(get_col(tab, "p_xy", NA_real_)),
    or = as.numeric(get_col(tab, "or", NA_real_)),
    p = as.numeric(get_col(tab, "p", NA_real_)),
    lift = as.numeric(get_col(tab, "lift", NA_real_)),
    npmi = as.numeric(get_col(tab, "npmi", NA_real_)),
    q = as.numeric(get_col(tab, "q", NA_real_))
  )
}

standardize_signed <- function(tab, default_scope, source_table) {
  if (!nrow(tab)) {
    return(tibble::tibble(
      source_table_signed = character(), scope = character(), term = character(), outcome = character(),
      signed_articles = integer(), pos_articles = integer(), neg_articles = integer(),
      null_or_mix = integer(), net_score = double(), pos_ratio = double(), balance = double(),
      balance_low = double(), balance_high = double(), interpretation_layer = character()
    ))
  }
  scope <- if ("analysis_scope" %in% names(tab)) get_col(tab, "analysis_scope", default_scope) else get_col(tab, "scope", default_scope)
  tibble::tibble(
    source_table_signed = source_table,
    scope = as.character(scope),
    term = as.character(get_col(tab, "A")),
    outcome = as.character(get_col(tab, "C")),
    signed_articles = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "articles", 0)), 0)),
    pos_articles = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "pos_articles", 0)), 0)),
    neg_articles = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "neg_articles", 0)), 0)),
    null_or_mix = as.integer(dplyr::coalesce(as.numeric(get_col(tab, "null_or_mix", 0)), 0)),
    net_score = as.numeric(get_col(tab, "net_score", NA_real_)),
    pos_ratio = as.numeric(get_col(tab, "pos_ratio", NA_real_)),
    balance = as.numeric(get_col(tab, "balance", NA_real_)),
    balance_low = as.numeric(get_col(tab, "balance_low", NA_real_)),
    balance_high = as.numeric(get_col(tab, "balance_high", NA_real_)),
    interpretation_layer = as.character(get_col(tab, "interpretation_layer", ""))
  )
}

make_current_meta <- function(dic) {
  for (nm in c("preferred_label", "class", "abc_role", "abc_subrole", "domain", "concept_family", "node_level", "analysis_tier", "atlas_layer", "include_primary", "primary_atlas_eligible")) {
    if (!nm %in% names(dic)) dic[[nm]] <- NA
  }
  dic |>
    dplyr::select(dplyr::any_of(c(
      "term", "preferred_label", "class", "abc_role", "abc_subrole", "domain",
      "concept_family", "node_level", "analysis_tier", "atlas_layer",
      "include_primary", "primary_atlas_eligible"
    ))) |>
    dplyr::mutate(
      preferred_label = dplyr::if_else(is.na(preferred_label) | preferred_label == "", term, as.character(preferred_label)),
      include_primary = as_bool(include_primary),
      primary_atlas_eligible = as_bool(primary_atlas_eligible)
    )
}

classify_domain <- function(scope, term, class, abc_role, abc_subrole, domain, concept_family, atlas_layer) {
  blob <- stringr::str_to_lower(paste(
    safe_chr(term), safe_chr(class), safe_chr(abc_role), safe_chr(abc_subrole),
    safe_chr(domain), safe_chr(concept_family), safe_chr(atlas_layer), sep = " | "
  ))
  dplyr::case_when(
    safe_chr(scope) == "treatment_context" | safe_chr(atlas_layer) == "treatment_context" ~ "Treatment context",
    safe_chr(scope) == "infection_context" | safe_chr(atlas_layer) == "infection_context" ~ "Infection context",
    safe_chr(term) %in% MEASUREMENT_GUARD_TERMS |
      safe_chr(abc_subrole) %in% MEASUREMENT_GUARD_SUBROLES |
      safe_chr(domain) %in% MEASUREMENT_GUARD_DOMAINS |
      stringr::str_detect(blob, "measurement_modality|measurement_derived|imaging_modality|quantitative_ct") ~ "Measurement/sensitivity",
    safe_chr(abc_role) == "A" ~ "Host/exposure/genetic/serology",
    stringr::str_detect(blob, "radiographic|pathologic|ild_pattern|pattern|honeycomb|reticulation|traction|ground.glass|organising|organizing|fibrotic_extent|radiology_extent|bronchiectasis|fibroblastic|dad|dip|lip|nsip|uip") ~ "Imaging pattern/extent",
    stringr::str_detect(blob, "pulmonary_function|pft|physiolog|severity_index|gap|cpi|fvc|dlco|percent.predicted") ~ "Physiology/severity",
    stringr::str_detect(blob, "biomarker|pathway|cytokine|chemokine|molecular|oxidative|ecm|tgf|fibrosis_pathway|immune|macrophage|th17|kl.6|sp.d|mmp|ccl|cxcl|il") ~ "Biomarker/pathway",
    TRUE ~ "Other primary disease-state"
  )
}

weighted_balance <- function(balance, weights) {
  ok <- !is.na(balance) & !is.na(weights) & weights > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(balance[ok], weights[ok])
}

label_friendly <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "ACPA/RF serology", "ACPA/RF")
  x <- stringr::str_replace_all(x, "macrophage cell", "macrophage")
  x <- stringr::str_replace_all(x, "neutrophil cell", "neutrophil")
  x <- stringr::str_replace_all(x, "Th17 cell", "Th17")
  x <- stringr::str_replace_all(x, "respiratory infection context", "respiratory infection")
  x <- stringr::str_replace_all(x, "infection safety event", "infection safety")
  x <- stringr::str_replace_all(x, "infection-associated AE-ILD context", "infection-associated AE-ILD")
  x <- stringr::str_replace_all(x, "immunosuppressive therapy", "immunosuppressive therapy")
  x <- stringr::str_replace_all(x, "drug-induced ILD", "drug-induced ILD")
  x <- stringr::str_replace_all(x, "DMARD category", "DMARD category")
  x
}

wrap_label <- function(x, width = 32) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

# Convert a hit-matrix column to logical. Works for TRUE/FALSE, 0/1, yes/no.
as_hit_bool <- function(x) {
  if (is.logical(x)) return(dplyr::coalesce(x, FALSE))
  if (is.numeric(x)) return(dplyr::coalesce(x > 0, FALSE))
  z <- tolower(as.character(x))
  z %in% c("true", "t", "1", "yes", "y")
}

# Resolve a dictionary concept node or outcome name to its hit-matrix column.
# 02_BuildHitsMatrix writes columns as hit__<term>, while several analytic
# tables use raw term names. The previous version looked for raw names only,
# which made all unique article-level counts zero.
resolve_hit_col <- function(name, hit_names) {
  name <- as.character(name)
  candidates <- unique(c(
    name,
    paste0("hit__", name),
    make.names(name),
    paste0("hit__", make.names(name))
  ))
  hit <- candidates[candidates %in% hit_names]
  if (length(hit)) return(hit[[1]])
  NA_character_
}

load_hits_optional <- function() {
  f_exact <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  if (file.exists(f_exact)) {
    log_msg("Loaded hit matrix:", f_exact)
    return(readr::read_csv(f_exact, show_col_types = FALSE))
  }
  alt <- find_latest_file(DIR_TABLE, sprintf("^hits_matrix_.*__%s\\.csv$", DIC_TAG))
  if (!is.na(alt) && file.exists(alt)) {
    log_msg("Loaded latest hit matrix:", alt)
    return(readr::read_csv(alt, show_col_types = FALSE))
  }
  log_msg("Hit matrix not found. Figure bubble size will use summed term-outcome co-mention support instead of unique article-level domain co-mentions.")
  tibble::tibble()
}

# -------------------------
# 2) Load 03-06 output tables and current dictionary metadata
# -------------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)
current_meta <- make_current_meta(dic)

cooc_primary <- read_required_table("cooc_npmi_lift_primary_atlas")
cooc_treat   <- read_optional_table("cooc_npmi_lift_treatment_context")
cooc_inf     <- read_optional_table("cooc_npmi_lift_infection_context")

signed_primary <- read_required_table("signed_effects_summary_primary_atlas")
signed_treat   <- read_optional_table("signed_effects_summary_treatment_context")
signed_inf     <- read_optional_table("signed_effects_summary_infection_context")

abc_ae <- read_optional_table("abc_rankings_primary_AEILD_main_filtered")
hits <- load_hits_optional()

# -------------------------
# 3) Standardise and merge co-occurrence and signed-effect evidence
# -------------------------
cooc_all <- dplyr::bind_rows(
  standardize_cooc(cooc_primary, "primary_atlas", "cooc_npmi_lift_primary_atlas"),
  standardize_cooc(cooc_treat,   "treatment_context", "cooc_npmi_lift_treatment_context"),
  standardize_cooc(cooc_inf,     "infection_context", "cooc_npmi_lift_infection_context")
) |>
  dplyr::filter(outcome %in% OUTCOME_ORDER)

signed_all <- dplyr::bind_rows(
  standardize_signed(signed_primary, "primary_atlas", "signed_effects_summary_primary_atlas"),
  standardize_signed(signed_treat,   "treatment_context", "signed_effects_summary_treatment_context"),
  standardize_signed(signed_inf,     "infection_context", "signed_effects_summary_infection_context")
) |>
  dplyr::filter(outcome %in% OUTCOME_ORDER)

matrix_long <- cooc_all |>
  dplyr::left_join(signed_all, by = c("scope", "term", "outcome")) |>
  dplyr::left_join(current_meta, by = "term") |>
  dplyr::mutate(
    signed_articles = dplyr::coalesce(signed_articles, 0L),
    pos_articles = dplyr::coalesce(pos_articles, 0L),
    neg_articles = dplyr::coalesce(neg_articles, 0L),
    null_or_mix = dplyr::coalesce(null_or_mix, 0L),
    preferred_label = dplyr::if_else(is.na(preferred_label) | preferred_label == "", term, preferred_label),
    display_preferred_label = label_friendly(preferred_label),
    outcome_label = dplyr::recode(outcome, !!!OUTCOME_LABELS),
    outcome_label = factor(outcome_label, levels = OUTCOME_LABELS[OUTCOME_ORDER]),
    clinical_domain = classify_domain(scope, term, class, abc_role, abc_subrole, domain, concept_family, atlas_layer),
    clinical_domain = factor(clinical_domain, levels = c(DOMAIN_ORDER, "Measurement/sensitivity", "Other primary disease-state")),
    measurement_guard = term %in% MEASUREMENT_GUARD_TERMS |
      safe_chr(abc_subrole) %in% MEASUREMENT_GUARD_SUBROLES |
      safe_chr(domain) %in% MEASUREMENT_GUARD_DOMAINS |
      clinical_domain == "Measurement/sensitivity",
    primary_scope_but_not_current_primary = scope == "primary_atlas" &
      !(atlas_layer == "primary_atlas" & include_primary == TRUE & primary_atlas_eligible == TRUE),
    display_score = cohit_n + 0.001 * dplyr::coalesce(signed_articles, 0L) + 0.0001 * dplyr::coalesce(npmi, 0)
  ) |>
  dplyr::arrange(scope, outcome_label, clinical_domain, dplyr::desc(cohit_n), dplyr::desc(signed_articles), term)

# -------------------------
# 4) QC exclusions and matrix main set
# -------------------------
qc_primary_non_current <- matrix_long |>
  dplyr::filter(scope == "primary_atlas", primary_scope_but_not_current_primary) |>
  dplyr::select(scope, term, preferred_label, outcome, cohit_n, atlas_layer, analysis_tier,
                include_primary, primary_atlas_eligible, abc_role, abc_subrole, domain, concept_family,
                measurement_guard)

qc_measurement_excluded <- matrix_long |>
  dplyr::filter(scope == "primary_atlas", measurement_guard) |>
  dplyr::select(scope, term, preferred_label, outcome, cohit_n, atlas_layer, analysis_tier,
                include_primary, primary_atlas_eligible, abc_role, abc_subrole, domain, concept_family)

matrix_main <- matrix_long |>
  dplyr::filter(
    outcome %in% OUTCOME_ORDER,
    clinical_domain %in% DOMAIN_ORDER,
    !(scope == "primary_atlas" & primary_scope_but_not_current_primary),
    !(scope == "primary_atlas" & measurement_guard)
  ) |>
  dplyr::mutate(clinical_domain = factor(as.character(clinical_domain), levels = DOMAIN_ORDER))

# -------------------------
# 5) Display-node selection
# -------------------------
ranked_for_display <- matrix_main |>
  dplyr::filter(cohit_n >= MIN_COHIT_FOR_LABEL) |>
  dplyr::group_by(outcome, outcome_label, clinical_domain) |>
  dplyr::arrange(dplyr::desc(cohit_n), dplyr::desc(signed_articles), dplyr::desc(abs(dplyr::coalesce(balance, 0))), dplyr::desc(dplyr::coalesce(npmi, 0)), term, .by_group = TRUE) |>
  dplyr::mutate(display_rank = dplyr::row_number()) |>
  dplyr::ungroup()

matrix_top_terms <- ranked_for_display |>
  dplyr::filter(display_rank <= TOP_N_TABLE) |>
  dplyr::mutate(
    display_label = display_preferred_label,
    display_cell_label = paste0(display_rank, ". ", display_preferred_label, " (n=", cohit_n, ")")
  ) |>
  dplyr::select(
    outcome, outcome_label, clinical_domain, display_rank, term, preferred_label, display_preferred_label,
    scope, abc_role, abc_subrole, domain, concept_family, node_level, atlas_layer,
    cohit_n, signed_articles, pos_articles, neg_articles, null_or_mix, balance,
    balance_low, balance_high, lift, npmi, q, display_score, interpretation_layer, display_cell_label
  ) |>
  dplyr::arrange(outcome_label, clinical_domain, display_rank)

label_terms_by_cell <- ranked_for_display |>
  dplyr::filter(display_rank <= TOP_N_LABEL) |>
  dplyr::group_by(outcome, outcome_label, clinical_domain) |>
  dplyr::summarise(
    top_concept_nodes_label = paste(display_preferred_label, collapse = "; "),
    top_concept_nodes_label_with_counts = paste(paste0(display_preferred_label, " (n=", cohit_n, ")"), collapse = "; "),
    .groups = "drop"
  )

full_top_terms_by_cell <- matrix_top_terms |>
  dplyr::group_by(outcome, outcome_label, clinical_domain) |>
  dplyr::summarise(
    top_concept_nodes_table = paste(display_preferred_label, collapse = "; "),
    top_concept_nodes_table_with_counts = paste(display_cell_label, collapse = "; "),
    .groups = "drop"
  )

# -------------------------
# 6) Domain-level matrix summary
# -------------------------
domain_summary <- matrix_main |>
  dplyr::group_by(outcome, outcome_label, clinical_domain) |>
  dplyr::summarise(
    concept_nodes_tested = dplyr::n_distinct(term),
    concept_nodes_with_cohit = dplyr::n_distinct(term[cohit_n > 0]),
    summed_term_outcome_cohits = sum(cohit_n, na.rm = TRUE),
    max_single_node_cohit = ifelse(any(cohit_n > 0), max(cohit_n, na.rm = TRUE), 0),
    total_signed_articles = sum(signed_articles, na.rm = TRUE),
    positive_direction_articles = sum(pos_articles, na.rm = TRUE),
    decreasing_or_protective_articles = sum(neg_articles, na.rm = TRUE),
    mixed_or_null_articles = sum(null_or_mix, na.rm = TRUE),
    weighted_direction_balance = weighted_balance(balance, signed_articles),
    median_npmi_nonzero = ifelse(any(cohit_n > 0), median(npmi[cohit_n > 0], na.rm = TRUE), NA_real_),
    max_npmi = ifelse(any(cohit_n > 0), max(npmi[cohit_n > 0], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  dplyr::left_join(label_terms_by_cell, by = c("outcome", "outcome_label", "clinical_domain")) |>
  dplyr::left_join(full_top_terms_by_cell, by = c("outcome", "outcome_label", "clinical_domain")) |>
  dplyr::mutate(
    outcome_label = factor(as.character(outcome_label), levels = OUTCOME_LABELS[OUTCOME_ORDER]),
    clinical_domain = factor(as.character(clinical_domain), levels = DOMAIN_ORDER),
    top_concept_nodes_label = dplyr::coalesce(top_concept_nodes_label, ""),
    top_concept_nodes_label_with_counts = dplyr::coalesce(top_concept_nodes_label_with_counts, ""),
    top_concept_nodes_table = dplyr::coalesce(top_concept_nodes_table, ""),
    top_concept_nodes_table_with_counts = dplyr::coalesce(top_concept_nodes_table_with_counts, "")
  ) |>
  dplyr::arrange(outcome_label, clinical_domain)

complete_grid <- tidyr::expand_grid(
  outcome = OUTCOME_ORDER,
  clinical_domain = DOMAIN_ORDER
) |>
  dplyr::mutate(
    outcome_label = factor(OUTCOME_LABELS[outcome], levels = OUTCOME_LABELS[OUTCOME_ORDER]),
    clinical_domain = factor(clinical_domain, levels = DOMAIN_ORDER)
  )

domain_summary_complete <- complete_grid |>
  dplyr::left_join(domain_summary, by = c("outcome", "outcome_label", "clinical_domain")) |>
  dplyr::mutate(
    concept_nodes_tested = dplyr::coalesce(concept_nodes_tested, 0L),
    concept_nodes_with_cohit = dplyr::coalesce(concept_nodes_with_cohit, 0L),
    summed_term_outcome_cohits = dplyr::coalesce(summed_term_outcome_cohits, 0L),
    max_single_node_cohit = dplyr::coalesce(max_single_node_cohit, 0L),
    total_signed_articles = dplyr::coalesce(total_signed_articles, 0L),
    positive_direction_articles = dplyr::coalesce(positive_direction_articles, 0L),
    decreasing_or_protective_articles = dplyr::coalesce(decreasing_or_protective_articles, 0L),
    mixed_or_null_articles = dplyr::coalesce(mixed_or_null_articles, 0L),
    top_concept_nodes_label = dplyr::coalesce(top_concept_nodes_label, ""),
    top_concept_nodes_label_with_counts = dplyr::coalesce(top_concept_nodes_label_with_counts, ""),
    top_concept_nodes_table = dplyr::coalesce(top_concept_nodes_table, ""),
    top_concept_nodes_table_with_counts = dplyr::coalesce(top_concept_nodes_table_with_counts, "")
  ) |>
  dplyr::arrange(outcome_label, clinical_domain)

# -------------------------
# 7) Unique article-level domain-outcome co-mentions from hit matrix
# -------------------------
unique_domain_counts <- tibble::tibble()
hit_matrix_available <- nrow(hits) > 0
unique_count_fallback_to_summed <- FALSE

if (hit_matrix_available) {
  hit_names <- names(hits)

  # Map analytic term names and outcome names to 02 hit-matrix columns.
  term_domain_map <- matrix_main |>
    dplyr::select(term, clinical_domain) |>
    dplyr::distinct() |>
    dplyr::mutate(hit_col = purrr::map_chr(term, resolve_hit_col, hit_names = hit_names))

  outcome_col_map <- tibble::tibble(
    outcome = OUTCOME_ORDER,
    outcome_hit_col = purrr::map_chr(OUTCOME_ORDER, resolve_hit_col, hit_names = hit_names)
  )

  missing_hit_terms <- term_domain_map |>
    dplyr::filter(is.na(hit_col))

  if (nrow(missing_hit_terms)) {
    f_missing <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_QC_terms_missing_from_hit_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
    write_csv2(missing_hit_terms, f_missing)
    log_msg("WARNING: Some matrix terms were absent from hit matrix; see:", f_missing)
  }

  f_outcome_map <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_QC_outcome_hit_column_map_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_csv2(outcome_col_map, f_outcome_map)

  term_domain_map_ok <- term_domain_map |>
    dplyr::filter(!is.na(hit_col))

  unique_domain_counts <- tidyr::expand_grid(
    outcome = OUTCOME_ORDER,
    clinical_domain = DOMAIN_ORDER
  ) |>
    dplyr::left_join(outcome_col_map, by = "outcome") |>
    dplyr::rowwise() |>
    dplyr::mutate(
      unique_domain_outcome_articles = {
        hit_cols <- term_domain_map_ok$hit_col[as.character(term_domain_map_ok$clinical_domain) == as.character(clinical_domain)]
        hit_cols <- unique(hit_cols[!is.na(hit_cols)])
        if (is.na(outcome_hit_col) || length(hit_cols) == 0) {
          0L
        } else {
          out <- as_hit_bool(hits[[outcome_hit_col]])
          mat <- as.data.frame(lapply(hits[hit_cols], as_hit_bool))
          any_term <- if (ncol(mat) == 1) mat[[1]] else rowSums(mat, na.rm = TRUE) > 0
          as.integer(sum(out & any_term, na.rm = TRUE))
        }
      }
    ) |>
    dplyr::ungroup()

  # Safety fallback: if the hit matrix was read but all unique counts are zero
  # despite non-zero term-outcome evidence, do not silently plot a zero-sized map.
  # This should not occur after hit__ column resolution, but the fallback keeps
  # the figure scientifically interpretable if a future hit-matrix schema changes.
  if (sum(unique_domain_counts$unique_domain_outcome_articles, na.rm = TRUE) == 0 &&
      sum(domain_summary_complete$summed_term_outcome_cohits, na.rm = TRUE) > 0) {
    unique_count_fallback_to_summed <- TRUE
    log_msg("WARNING: unique domain-outcome article counts were all zero; falling back to summed term-outcome co-mentions for Figure 2 bubble size.")
  }
}

if (nrow(unique_domain_counts) && !unique_count_fallback_to_summed) {
  domain_summary_complete <- domain_summary_complete |>
    dplyr::left_join(unique_domain_counts, by = c("outcome", "clinical_domain")) |>
    dplyr::mutate(unique_domain_outcome_articles = dplyr::coalesce(unique_domain_outcome_articles, 0L))
} else {
  domain_summary_complete <- domain_summary_complete |>
    dplyr::mutate(unique_domain_outcome_articles = NA_integer_)
}

# Default figure size metric: unique article-level domain-outcome co-mentions if available.
domain_summary_complete <- domain_summary_complete |>
  dplyr::mutate(
    figure_size_metric = dplyr::if_else(!is.na(unique_domain_outcome_articles), as.numeric(unique_domain_outcome_articles), as.numeric(summed_term_outcome_cohits)),
    figure_size_metric_label = dplyr::if_else(!is.na(unique_domain_outcome_articles), "Unique domain-outcome articles", "Summed term-outcome co-mentions"),
    figure_size_for_plot = dplyr::if_else(figure_size_metric > 0, figure_size_metric, NA_real_),
    y_num = as.numeric(DOMAIN_Y[as.character(clinical_domain)]),
    x_num = as.numeric(OUTCOME_X[outcome]),
    # Labels are governed by term-outcome evidence, not by the unique-size metric.
    plot_label = dplyr::if_else(summed_term_outcome_cohits > 0 & top_concept_nodes_label != "", wrap_label(top_concept_nodes_label, width = 30), "")
  )

# -------------------------
# 8) AE-ILD ABC annotation for Figure 3 / optional Figure 2 note
# -------------------------
abc_ae_std <- tibble::tibble()
if (nrow(abc_ae)) {
  abc_ae_std <- abc_ae |>
    dplyr::select(dplyr::any_of(c(
      "A", "B", "C", "AB_n11", "BC_n11", "AC_n11", "AB_lift", "BC_lift", "AB_npmi", "BC_npmi", "score_q",
      "A_preferred_label", "B_preferred_label", "C_preferred_label",
      "A_domain", "B_domain", "A_concept_family", "B_concept_family"
    ))) |>
    dplyr::filter(C == "AE-ILD") |>
    dplyr::arrange(dplyr::desc(score_q))
}

# -------------------------
# 9) Save output tables
# -------------------------
f_matrix_long <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_long_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_matrix_main <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_main_pairs_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_domain_summary <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_domain_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_top_terms <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_top_terms_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_fig2_display <- file.path(DIR_TABLE, sprintf("Figure2_outcome_domain_atlas_DISPLAY_TABLE_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_abc_ae <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_AEILD_ABC_annotation_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_qc_primary <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_QC_primary_noncurrent_terms_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_qc_measure <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_QC_measurement_nodes_excluded_%s__%s.csv", CORPUS_TAG, DIC_TAG))

write_csv2(matrix_long, f_matrix_long)
write_csv2(matrix_main, f_matrix_main)
write_csv2(domain_summary_complete, f_domain_summary)
write_csv2(matrix_top_terms, f_top_terms)
write_csv2(domain_summary_complete, f_fig2_display)
write_csv2(abc_ae_std, f_abc_ae)
write_csv2(qc_primary_non_current, f_qc_primary)
write_csv2(qc_measurement_excluded, f_qc_measure)

# -------------------------
# 10) Final Figure 2
# -------------------------
plot_df <- domain_summary_complete |>
  dplyr::mutate(
    layer_type = dplyr::if_else(clinical_domain %in% c("Treatment context", "Infection context"), "Context layers", "Primary disease-state layers"),
    plot_label = dplyr::if_else(plot_label != "", plot_label, ""),
    signed_balance_for_plot = weighted_direction_balance
  )

# Slight vertical nudges keep text legible and help separate context rows.
plot_df <- plot_df |>
  dplyr::mutate(
    label_y = y_num + dplyr::case_when(
      clinical_domain == "Host/exposure/genetic/serology" ~ 0.29,
      clinical_domain == "Infection context" ~ 0.30,
      TRUE ~ 0.27
    )
  )

size_title <- ifelse(any(!is.na(plot_df$unique_domain_outcome_articles)),
                     "Unique domain-\noutcome articles",
                     "Summed term-outcome\nco-mentions")

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x_num, y = y_num)) +
  # Light context-layer background. This does not encode data; it separates interpretation layers.
  ggplot2::annotate("rect", xmin = 0.45, xmax = 3.55, ymin = 0.45, ymax = 2.50,
                    fill = "grey95", colour = NA) +
  ggplot2::geom_hline(yintercept = 2.50, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  ggplot2::geom_point(
    data = plot_df |> dplyr::filter(!is.na(figure_size_for_plot)),
    ggplot2::aes(size = figure_size_for_plot, colour = signed_balance_for_plot),
    alpha = 0.88
  ) +
  ggplot2::geom_text(
    data = plot_df |> dplyr::filter(plot_label != ""),
    ggplot2::aes(y = label_y, label = plot_label),
    size = 2.65,
    lineheight = 0.90,
    vjust = 0.5,
    colour = "black"
  ) +
  ggplot2::scale_x_continuous(
    breaks = unname(OUTCOME_X),
    labels = OUTCOME_LABELS[OUTCOME_ORDER],
    limits = c(0.45, 3.55),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    breaks = unname(DOMAIN_Y[DOMAIN_ORDER]),
    labels = DOMAIN_ORDER,
    limits = c(0.45, 6.72),
    expand = ggplot2::expansion(mult = c(0.02, 0.03))
  ) +
  ggplot2::scale_size_area(
    name = size_title,
    max_size = 12.5,
    breaks = scales::pretty_breaks(n = 4)
  ) +
  ggplot2::scale_colour_gradient2(
    name = "Signed-effect\nbalance",
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), na.value = "grey70"
  ) +
  ggplot2::labs(
    title = "Outcome-domain atlas of RA-ILD worsening literature",
    subtitle = "Primary disease-state layers are displayed separately from treatment and infection context layers",
    x = "Clinical worsening outcome",
    y = "Dictionary domain layer"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.major.x = ggplot2::element_line(colour = "grey88", linewidth = 0.35),
    axis.text.y = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 8)),
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(size = 10.5, margin = ggplot2::margin(b = 8)),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 9.5),
    legend.text = ggplot2::element_text(size = 8.5),
    plot.margin = ggplot2::margin(t = 10, r = 18, b = 10, l = 10)
  )

f_pdf <- file.path(DIR_FIG2, sprintf("Figure2_outcome_domain_atlas_FINAL_%s__%s.pdf", CORPUS_TAG, DIC_TAG))
f_png <- file.path(DIR_FIG2, sprintf("Figure2_outcome_domain_atlas_FINAL_%s__%s.png", CORPUS_TAG, DIC_TAG))
f_svg <- file.path(DIR_FIG2, sprintf("Figure2_outcome_domain_atlas_FINAL_%s__%s.svg", CORPUS_TAG, DIC_TAG))

ggplot2::ggsave(f_pdf, p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in")
ggplot2::ggsave(f_png, p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", dpi = FIG_DPI)
# SVG may fail if device support is unavailable; do not stop the pipeline for this.
try(ggplot2::ggsave(f_svg, p, width = FIG_WIDTH, height = FIG_HEIGHT, units = "in"), silent = TRUE)

log_msg("WROTE:", f_pdf)
log_msg("WROTE:", f_png)

# -------------------------
# 11) Run summary
# -------------------------
run_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  top_n_label = TOP_N_LABEL,
  top_n_table = TOP_N_TABLE,
  min_cohit_for_label = MIN_COHIT_FOR_LABEL,
  hit_matrix_available = hit_matrix_available,
  unique_count_fallback_to_summed = unique_count_fallback_to_summed,
  size_metric = ifelse(hit_matrix_available && !unique_count_fallback_to_summed, "unique_domain_outcome_articles", "summed_term_outcome_cohits"),
  unique_domain_outcome_articles_total = sum(domain_summary_complete$unique_domain_outcome_articles, na.rm = TRUE),
  figure_size_metric_total = sum(domain_summary_complete$figure_size_metric, na.rm = TRUE),
  input_primary_pairs = nrow(standardize_cooc(cooc_primary, "primary_atlas", "cooc_npmi_lift_primary_atlas")),
  input_treatment_pairs = nrow(standardize_cooc(cooc_treat, "treatment_context", "cooc_npmi_lift_treatment_context")),
  input_infection_pairs = nrow(standardize_cooc(cooc_inf, "infection_context", "cooc_npmi_lift_infection_context")),
  main_matrix_pairs = nrow(matrix_main),
  main_matrix_pairs_with_cohit = sum(matrix_main$cohit_n > 0, na.rm = TRUE),
  display_table_rows = nrow(domain_summary_complete),
  top_term_rows_written = nrow(matrix_top_terms),
  primary_noncurrent_terms_excluded_rows = nrow(qc_primary_non_current),
  measurement_nodes_excluded_rows = nrow(qc_measurement_excluded),
  context_layers_retained = paste(c("treatment_context", "infection_context"), collapse = ";"),
  output_pdf = f_pdf,
  output_png = f_png,
  interpretation = paste(
    "Final Figure 2 displays an outcome-domain atlas separating primary disease-state layers",
    "from treatment and infection context layers. Labels show top concept nodes only for readability;",
    "full ranked cell contents are saved in outcome_domain_matrix_top_terms."
  )
)

f_summary <- file.path(DIR_TABLE, sprintf("outcome_domain_matrix_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(run_summary, f_summary)

log_msg("07 summary | main pairs=", nrow(matrix_main),
        " top-term rows=", nrow(matrix_top_terms),
        " measurement excluded rows=", nrow(qc_measurement_excluded),
        " primary noncurrent excluded rows=", nrow(qc_primary_non_current),
        " size metric=", ifelse(hit_matrix_available && !unique_count_fallback_to_summed, "unique articles", "summed cohits"),
        " unique total=", sum(domain_summary_complete$unique_domain_outcome_articles, na.rm = TRUE))
log_msg("=== DONE ", RUN_NAME, " ===")
