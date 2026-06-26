# 15_Build_SupplementaryTables_S1_to_S7_AtlasV2.R
# Build human-readable Supplementary Tables S1-S7 from final AtlasV2 outputs.
# This script does not add new terms, rerank figures, or alter primary analyses.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(grid)
})

log_msg <- function(...) cat("[15]", ..., "\n")

# -----------------------------
# Configuration and directories
# -----------------------------
# Public-ready path handling:
# Run this script from the repository root, or set RAILD_ROOT / RAILD_CODE_DIR.
ROOT <- normalizePath(Sys.getenv("RAILD_ROOT", unset = getwd()), mustWork = FALSE)
CODE_DIR <- normalizePath(Sys.getenv("RAILD_CODE_DIR", unset = ROOT), mustWork = FALSE)

setup_path <- file.path(CODE_DIR, "00_setup_AtlasV2.R")
if (file.exists(setup_path)) {
  try(source(setup_path), silent = TRUE)
}

CORPUS_TAG <- if (exists("CORPUS_TAG")) CORPUS_TAG else Sys.getenv("RAILD_CORPUS_TAG", unset = "pm_1980_20251231")
ANALYSIS_TAG <- if (exists("ANALYSIS_TAG")) ANALYSIS_TAG else Sys.getenv("RAILD_ANALYSIS_TAG", unset = "analysis_v2_outcome_preserving_atlas")
DIR_OUTPUT <- if (exists("DIR_OUTPUT")) DIR_OUTPUT else file.path(ROOT, "output", CORPUS_TAG, ANALYSIS_TAG)
DIR_TABLE <- if (exists("DIR_TABLE")) DIR_TABLE else file.path(DIR_OUTPUT, "table")
DIR_FIG <- if (exists("DIR_FIG")) DIR_FIG else file.path(DIR_OUTPUT, "fig")
DIR_SUPP <- file.path(DIR_OUTPUT, "supplementary_tables")
DIR_DATA_PROC <- if (exists("DIR_DATA_PROC")) DIR_DATA_PROC else file.path(ROOT, "data_proc")
DIR_DIC <- if (exists("DIR_DIC")) DIR_DIC else file.path(ROOT, "dic")

dir.create(DIR_SUPP, recursive = TRUE, showWarnings = FALSE)

log_msg("=== START 15_Build_SupplementaryTables_S1_to_S7_AtlasV2 ===")
log_msg("ROOT:", ROOT)
log_msg("TABLE:", DIR_TABLE)
log_msg("SUPP:", DIR_SUPP)

OUTCOMES <- c("AE-ILD", "progression", "mortality")

# -----------------------------
# Helper functions
# -----------------------------
read_csv_safe <- function(path, show_col_types = FALSE) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !file.exists(path)) return(tibble())
  suppressMessages(readr::read_csv(path, show_col_types = show_col_types, progress = FALSE))
}

newest_file <- function(pattern, dirs, recursive = FALSE, required = FALSE) {
  files <- character(0)
  for (d in dirs) {
    if (dir.exists(d)) {
      files <- c(files, list.files(d, pattern = pattern, recursive = recursive, full.names = TRUE))
    }
  }
  files <- unique(files)
  if (length(files) == 0) {
    if (required) stop("No file found for pattern: ", pattern)
    return(NA_character_)
  }
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

pick_col <- function(df, candidates) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit) == 0) NA_character_ else hit[1]
}

col_chr <- function(df, candidates, default = "") {
  cc <- pick_col(df, candidates)
  if (is.na(cc)) rep(default, nrow(df)) else as.character(df[[cc]])
}
col_num <- function(df, candidates, default = NA_real_) {
  cc <- pick_col(df, candidates)
  if (is.na(cc)) rep(default, nrow(df)) else suppressWarnings(as.numeric(df[[cc]]))
}
col_int <- function(df, candidates, default = NA_integer_) {
  cc <- pick_col(df, candidates)
  if (is.na(cc)) rep(default, nrow(df)) else suppressWarnings(as.integer(round(as.numeric(df[[cc]]))))
}

write_supp_csv <- function(df, name) {
  path <- file.path(DIR_SUPP, paste0(name, "_", CORPUS_TAG, "__", ANALYSIS_TAG, ".csv"))
  readr::write_csv(df, path, na = "")
  log_msg("WROTE:", path)
  path
}

safe_pct <- function(x, denom, digits = 1) {
  ifelse(is.na(x) | is.na(denom) | denom == 0, NA_real_, round(100 * x / denom, digits))
}

truncate_text <- function(x, n = 70) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
}

# Human-readable label harmonisation for supplementary tables only.
# This does not alter dictionary matching or analytic term identifiers.
pretty_label <- function(x) {
  x <- as.character(x)
  label_map <- c(
    "CRP ESR" = "CRP/ESR",
    "MMP7" = "MMP-7",
    "Surfactant protein A" = "SP-A",
    "Surfactant protein D" = "SP-D",
    "TGFB axis" = "TGF-beta axis",
    "IL17A" = "IL-17A",
    "GMCSF" = "GM-CSF",
    "YKL40" = "YKL-40",
    "SPP1 osteopontin" = "SPP1/osteopontin",
    "IL8 CXCL8" = "IL-8/CXCL8",
    "IL8CXCL8" = "IL-8/CXCL8",
    "ECM remodeling" = "ECM remodelling",
    "Endothelin1" = "Endothelin-1",
    "sPD1" = "sPD-1",
    "sPDL1" = "sPD-L1"
  )
  idx <- match(x, names(label_map))
  out <- x
  hit <- !is.na(idx)
  out[hit] <- unname(label_map[idx[hit]])
  out
}

# -----------------------------
# Load core data
# -----------------------------
dic_path <- newest_file("^ra_ild_dictionary_analysis_v2_outcome_preserving_atlas\\.csv$", DIR_DIC, required = TRUE)
dic <- read_csv_safe(dic_path)

articles_path <- file.path(DIR_DATA_PROC, "articles_main_original_20260323.csv")
if (!file.exists(articles_path)) {
  articles_path <- newest_file("^articles_main_original_.*\\.csv$", DIR_DATA_PROC, required = TRUE)
}
articles <- read_csv_safe(articles_path)

hit_path <- newest_file("^hits_matrix_.*analysis_v2_outcome_preserving_atlas.*\\.csv$", DIR_TABLE, required = FALSE)
hits <- read_csv_safe(hit_path)

# Dictionary metadata safety
if (!"term" %in% names(dic)) stop("Dictionary must contain a 'term' column.")
for (cc in c("preferred_label", "atlas_layer", "abc_role", "abc_subrole", "domain", "concept_family", "node_level", "analysis_tier", "include_primary", "primary_atlas_eligible")) {
  if (!cc %in% names(dic)) dic[[cc]] <- NA
}

dic <- dic |>
  mutate(
    term = as.character(term),
    preferred_label = ifelse(is.na(preferred_label) | preferred_label == "", term, as.character(preferred_label)),
    preferred_label = pretty_label(preferred_label)
  )

n_articles <- nrow(articles)

# hit matrix term hit counts
# IMPORTANT: some upstream dictionaries or coverage tables may already contain
# corpus_hit_n / corpus_hit_pct.  Remove those columns before joining the
# current hit-matrix-derived counts; otherwise dplyr will create .x/.y columns
# and downstream summaries will not find corpus_hit_n.
term_hit_counts <- tibble(term = dic$term) |>
  mutate(
    hit_col = paste0("hit__", term),
    corpus_hit_n = purrr::map_int(hit_col, ~ if (.x %in% names(hits)) sum(as.numeric(hits[[.x]]) > 0, na.rm = TRUE) else 0L),
    corpus_hit_pct = safe_pct(corpus_hit_n, nrow(hits), 2)
  ) |>
  select(term, corpus_hit_n, corpus_hit_pct)

dic_cov <- dic |>
  dplyr::select(-dplyr::any_of(c("hit_col", "corpus_hit_n", "corpus_hit_pct", "corpus_hit_n.x", "corpus_hit_n.y", "corpus_hit_pct.x", "corpus_hit_pct.y"))) |>
  left_join(term_hit_counts, by = "term")

# Safety guard: make sure coverage columns always exist for Supplementary Table S3.
if (!"corpus_hit_n" %in% names(dic_cov)) dic_cov$corpus_hit_n <- 0L
if (!"corpus_hit_pct" %in% names(dic_cov)) dic_cov$corpus_hit_pct <- 0
log_msg("Dictionary coverage columns:", paste(intersect(c("term", "atlas_layer", "corpus_hit_n", "corpus_hit_pct"), names(dic_cov)), collapse = ", "))

# term-outcome cohit helper
compute_outcome_support <- function(terms) {
  tibble(term = terms) |>
    rowwise() |>
    mutate(
      cohit_AE_ILD = {
        hc1 <- paste0("hit__", term); hc2 <- "hit__AE-ILD";
        if (hc1 %in% names(hits) && hc2 %in% names(hits)) sum(as.numeric(hits[[hc1]]) > 0 & as.numeric(hits[[hc2]]) > 0, na.rm = TRUE) else 0L
      },
      cohit_progression = {
        hc1 <- paste0("hit__", term); hc2 <- "hit__progression";
        if (hc1 %in% names(hits) && hc2 %in% names(hits)) sum(as.numeric(hits[[hc1]]) > 0 & as.numeric(hits[[hc2]]) > 0, na.rm = TRUE) else 0L
      },
      cohit_mortality = {
        hc1 <- paste0("hit__", term); hc2 <- "hit__mortality";
        if (hc1 %in% names(hits) && hc2 %in% names(hits)) sum(as.numeric(hits[[hc1]]) > 0 & as.numeric(hits[[hc2]]) > 0, na.rm = TRUE) else 0L
      }
    ) |>
    ungroup() |>
    mutate(
      total_outcome_cohit = cohit_AE_ILD + cohit_progression + cohit_mortality,
      max_outcome_cohit = pmax(cohit_AE_ILD, cohit_progression, cohit_mortality, na.rm = TRUE),
      dominant_outcome = dplyr::case_when(
        max_outcome_cohit == 0 ~ "none",
        cohit_AE_ILD == max_outcome_cohit ~ "AE-ILD",
        cohit_progression == max_outcome_cohit ~ "progression",
        TRUE ~ "mortality"
      ),
      outcomes_supported = purrr::pmap_chr(list(cohit_AE_ILD, cohit_progression, cohit_mortality), function(ae, pr, mo) {
        out <- OUTCOMES[c(ae, pr, mo) > 0]
        if (length(out) == 0) "none" else paste(out, collapse = "; ")
      })
    )
}

# -----------------------------
# Supplementary Table S1
# -----------------------------
pmid_col <- pick_col(articles, c("pmid", "PMID", "PubMedID", "pubmed_id"))
year_col <- pick_col(articles, c("year", "pub_year", "publication_year", "Year"))

unique_pmids <- if (!is.na(pmid_col)) length(unique(articles[[pmid_col]])) else NA_integer_
years <- if (!is.na(year_col)) suppressWarnings(as.integer(articles[[year_col]])) else integer(0)
min_year <- if (length(years) && any(!is.na(years))) min(years, na.rm = TRUE) else 1980L
max_year <- if (length(years) && any(!is.na(years))) max(years, na.rm = TRUE) else 2025L

year_counts_file <- newest_file("^year_counts_stacked_by_pubclass.*\\.csv$", DIR_DATA_PROC, required = FALSE)
year_counts <- read_csv_safe(year_counts_file)
if (nrow(year_counts) == 0 && !is.na(year_col)) {
  year_counts <- articles |>
    mutate(year = as.integer(.data[[year_col]])) |>
    count(year, name = "total_records") |>
    arrange(year)
}

s1a <- tibble(
  item = c(
    "Search platform",
    "Corpus date / frozen run",
    "Publication years analysed",
    "Main analytic corpus",
    "Unique PubMed identifiers in main corpus",
    "Main article file",
    "Year-count file",
    "Frozen PMID list for repository"
  ),
  value = c(
    "MEDLINE via PubMed",
    "23 March 2026; records published 1980-2025",
    paste0(min_year, "-", max_year),
    paste0(n_articles, " abstract-bearing original research articles"),
    as.character(unique_pmids),
    basename(articles_path),
    ifelse(is.na(year_counts_file), "Not found", basename(year_counts_file)),
    "Dataset_S1_frozen_main_original_pmids.csv"
  ),
  note = c(
    "Search and retrieval metadata should be reported in Supplementary Table S1.",
    "The final analytic corpus is defined by the frozen PMID-level corpus file.",
    "Publication years are read from the article metadata used by the final hit matrix.",
    "Reviews/guidelines, case reports and editorial-like records were excluded from primary analyses.",
    "PMID count is used for reproducibility; PubMed abstract text itself is not redistributed.",
    "Input file used by the final pipeline.",
    "Used for corpus composition plotting when available.",
    "Generated by repository export code; contains PMIDs and minimal bibliographic metadata, not abstracts."
  )
)

s1b <- year_counts |>
  mutate(year = suppressWarnings(as.integer(.data[[pick_col(year_counts, c("year", "Year"))]]))) |>
  filter(!is.na(year)) |>
  arrange(year) |>
  tail(25)
# standardise count column for PDF
if (nrow(s1b) > 0) {
  cnt_col <- pick_col(s1b, c("total_records", "Original", "original", "n", "count", "records"))
  if (is.na(cnt_col)) {
    num_cols <- names(s1b)[vapply(s1b, is.numeric, logical(1))]
    cnt_col <- setdiff(num_cols, "year")[1]
  }
  s1b <- s1b |>
    transmute(year = year, total_records = if (!is.na(cnt_col)) as.integer(.data[[cnt_col]]) else NA_integer_) |>
    group_by(year) |>
    summarise(total_records = sum(total_records, na.rm = TRUE), .groups = "drop") |>
    arrange(year)
}

write_supp_csv(s1a, "Supplementary_Table_S1A_Corpus_Definition")
write_supp_csv(s1b, "Supplementary_Table_S1B_Recent_Annual_Corpus_Counts")

# -----------------------------
# Supplementary Table S2
# -----------------------------
s2a <- tibble(
  item = c(
    "Final dictionary file",
    "Total concept nodes",
    "Primary atlas concept nodes",
    "Treatment-context concept nodes",
    "Infection-context concept nodes",
    "Sensitivity / measurement / supporting nodes",
    "Exploratory concept nodes",
    "Context-only concept nodes"
  ),
  value = c(
    basename(dic_path),
    nrow(dic),
    sum(dic$atlas_layer == "primary_atlas", na.rm = TRUE),
    sum(dic$atlas_layer == "treatment_context", na.rm = TRUE),
    sum(dic$atlas_layer == "infection_context", na.rm = TRUE),
    sum(dic$atlas_layer == "sensitivity", na.rm = TRUE),
    sum(dic$atlas_layer == "exploratory", na.rm = TRUE),
    sum(dic$atlas_layer == "context_only", na.rm = TRUE)
  ),
  note = c(
    "Concept nodes are dictionary entries, not raw tokens.",
    "Final dictionary used by all final figures and tables.",
    "Eligible for primary disease-state atlas analyses unless additional role-specific filters apply.",
    "Retained separately; not interpreted as treatment effects or drug-risk estimates.",
    "Retained separately; not interpreted as primary disease-state bridges.",
    "Includes measurement/modality or supporting concepts not eligible for primary bridge claims.",
    "Hypothesis-generating biomarker/molecular/pathway or lower-support concepts.",
    "Comparator/context concepts retained for corpus interpretation, not primary RA-ILD risk claims."
  )
)

s2b <- dic |>
  count(atlas_layer, node_level, name = "concept_nodes") |>
  arrange(atlas_layer, node_level)

# QC summary, built from current dictionary directly. This avoids an empty QC-index table.
regex_col <- pick_col(dic, c("regex", "pattern", "regular_expression"))
regex_compile_errors <- 0L
if (!is.na(regex_col)) {
  regex_compile_errors <- sum(purrr::map_lgl(dic[[regex_col]], function(rx) {
    if (is.na(rx) || rx == "") return(TRUE)
    inherits(try(grepl(rx, "test abstract text", perl = TRUE), silent = TRUE), "try-error")
  }))
}
mandatory_meta <- c("term", "preferred_label", "atlas_layer", "abc_role", "domain", "node_level", "analysis_tier")
missing_meta_cells <- sum(is.na(dic[intersect(mandatory_meta, names(dic))]) | dic[intersect(mandatory_meta, names(dic))] == "")
duplicate_terms <- sum(duplicated(dic$term))

s2c <- tibble(
  qc_domain = c(
    "Regex compilation",
    "Duplicate canonical terms",
    "Mandatory metadata completeness",
    "Ambiguous acronym / substring safety",
    "Treatment context separation",
    "Infection context separation",
    "Measurement / modality separation",
    "Machine-readable QC/change-log files"
  ),
  summary = c(
    paste0(regex_compile_errors, " regex compile errors detected in the final dictionary check"),
    paste0(duplicate_terms, " duplicated canonical term identifiers detected"),
    paste0(missing_meta_cells, " missing mandatory metadata cells across key fields"),
    "Collision-prone terms are controlled at the dictionary layer; full regex and change-log files should be released in the repository.",
    paste0(sum(dic$atlas_layer == "treatment_context", na.rm = TRUE), " treatment-context nodes retained separately from primary bridge claims"),
    paste0(sum(dic$atlas_layer == "infection_context", na.rm = TRUE), " infection-context nodes retained separately from primary bridge claims"),
    paste0(sum(dic$atlas_layer == "sensitivity", na.rm = TRUE), " sensitivity/supporting nodes include measurement or modality concepts not eligible for primary bridge claims"),
    "Dictionary QC, change-log and collision-test CSV files should be included in the public repository release."
  ),
  interpretation = c(
    "Regex patterns should compile before corpus mapping.",
    "Canonical concept identifiers should be unique.",
    "Metadata fields define analysis layers and figure eligibility.",
    "Examples include short acronyms, substring collisions and generic pneumonia terminology.",
    "Context-layer terms are literature context, not treatment-effect estimates.",
    "Infection terms may represent trigger, complication, differential diagnosis or safety context.",
    "Measurement nodes are retained for coverage/context but not as disease-state bridge claims.",
    "Repository files provide machine-readable audit trail beyond this human-readable summary."
  )
)

write_supp_csv(s2a, "Supplementary_Table_S2A_Dictionary_Layer_Summary")
write_supp_csv(s2b, "Supplementary_Table_S2B_Dictionary_Layer_NodeLevel_Counts")
write_supp_csv(s2c, "Supplementary_Table_S2C_Dictionary_QC_Summary")

# -----------------------------
# Supplementary Table S3
# -----------------------------
dict_run_path <- newest_file("^dictionary_run_summary_.*\\.csv$", DIR_TABLE, required = FALSE)
dict_run <- read_csv_safe(dict_run_path)
if (nrow(dict_run) == 0) {
  dict_run <- tibble(
    metric = c("n_articles", "n_dictionary_terms", "abstract_coverage_pct"),
    value = c(
      nrow(hits),
      nrow(dic),
      {
        hit_cols_tmp <- grep("^hit__", names(hits), value = TRUE)
        if (length(hit_cols_tmp) == 0 || nrow(hits) == 0) NA_real_ else
          safe_pct(sum(rowSums(dplyr::select(hits, dplyr::all_of(hit_cols_tmp)), na.rm = TRUE) > 0), nrow(hits), 2)
      }
    )
  )
} else {
  if (ncol(dict_run) > 1 && nrow(dict_run) == 1) {
    dict_run <- tibble(metric = names(dict_run), value = as.character(as.list(dict_run[1, ])))
  }
}

s3a <- dict_run
s3b <- dic_cov |>
  group_by(atlas_layer) |>
  summarise(
    concept_nodes = n(),
    nodes_with_hit = sum(corpus_hit_n > 0, na.rm = TRUE),
    median_hit_articles = median(corpus_hit_n, na.rm = TRUE),
    max_hit_articles = max(corpus_hit_n, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(concept_nodes))

s3c <- dic_cov |>
  arrange(desc(corpus_hit_n), term) |>
  transmute(term, preferred_label, atlas_layer, domain, node_level, corpus_hit_n, corpus_hit_pct) |>
  head(as.integer(Sys.getenv("RAILD_SUPP_TOP_HIT_TERMS", unset = "40")))

write_supp_csv(s3a, "Supplementary_Table_S3A_Dictionary_Run_Summary")
write_supp_csv(s3b, "Supplementary_Table_S3B_Layer_Coverage_Summary")
write_supp_csv(s3c, "Supplementary_Table_S3C_Highest_Frequency_Concept_Nodes")

# -----------------------------
# Supplementary Table S4 Signed-effect direction
# -----------------------------
read_signed_scope <- function(scope, pattern) {
  p <- newest_file(pattern, DIR_TABLE, required = FALSE)
  x <- read_csv_safe(p)
  if (nrow(x) == 0) return(tibble())
  x |> mutate(scope = scope, source_file = basename(p))
}

signed_all <- bind_rows(
  read_signed_scope("primary_atlas", "^signed_effects_summary_primary_atlas_.*\\.csv$"),
  read_signed_scope("treatment_context", "^signed_effects_summary_treatment_context_.*\\.csv$"),
  read_signed_scope("infection_context", "^signed_effects_summary_infection_context_.*\\.csv$")
)

if (nrow(signed_all) > 0) {
  signed_std <- signed_all |>
    transmute(
      scope = scope,
      term = col_chr(cur_data_all(), c("term", "A", "concept", "term_clean")),
      outcome = col_chr(cur_data_all(), c("outcome", "C")),
      signed_balance = col_num(cur_data_all(), c("signed_balance", "directional_balance", "balance")),
      article_support = col_num(cur_data_all(), c("article_support", "n_articles", "article_n", "articles", "n_article", "total_articles", "n"), default = 0),
      positive_wording = col_num(cur_data_all(), c("positive_wording", "positive_articles", "risk_up_articles", "n_pos", "positive"), default = NA_real_),
      decreasing_or_protective_wording = col_num(cur_data_all(), c("decreasing_or_protective_wording", "negative_articles", "risk_down_articles", "n_neg", "decreasing"), default = NA_real_),
      mixed_or_unclear_wording = col_num(cur_data_all(), c("mixed_or_unclear_wording", "mixed_articles", "no_effect_mixed_articles", "n_mixed", "mixed"), default = NA_real_)
    ) |>
    filter(term != "", outcome %in% OUTCOMES)
} else {
  signed_std <- tibble(scope = character(), term = character(), outcome = character(), signed_balance = numeric(), article_support = numeric(), positive_wording = numeric(), decreasing_or_protective_wording = numeric(), mixed_or_unclear_wording = numeric())
}

s4a <- signed_std |>
  group_by(scope, outcome) |>
  summarise(
    pairs = n(),
    pairs_with_directional_articles = sum(article_support > 0, na.rm = TRUE),
    median_article_support = median(article_support, na.rm = TRUE),
    positive_leaning_pairs = sum(signed_balance > 0, na.rm = TRUE),
    decreasing_or_protective_leaning_pairs = sum(signed_balance < 0, na.rm = TRUE),
    mixed_or_zero_pairs = sum(is.na(signed_balance) | signed_balance == 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(scope, outcome)

s4b <- signed_std |>
  filter(scope == "primary_atlas") |>
  left_join(dic |> select(term, preferred_label, atlas_layer, domain, node_level), by = "term") |>
  arrange(outcome, desc(article_support), desc(abs(signed_balance)), term) |>
  group_by(outcome) |>
  slice_head(n = as.integer(Sys.getenv("RAILD_SUPP_SIGNED_TOP_PER_OUTCOME", unset = "25"))) |>
  ungroup() |>
  transmute(
    outcome,
    term,
    preferred_label = coalesce(preferred_label, term),
    domain,
    node_level,
    signed_balance = round(signed_balance, 3),
    article_support = round(article_support, 1),
    positive_wording,
    decreasing_or_protective_wording,
    mixed_or_unclear_wording
  )

write_supp_csv(s4a, "Supplementary_Table_S4A_Signed_Effect_Direction_Summary")
write_supp_csv(s4b, "Supplementary_Table_S4B_Primary_Atlas_Signed_Effect_TopPairs")

# -----------------------------
# Supplementary Table S5 Context and measurement support
# -----------------------------
context_terms <- dic_cov |>
  # S5 is restricted to prespecified non-primary context/supporting layers.
  # Primary-atlas terms such as disease duration are intentionally excluded here,
  # even if their metadata contain words such as "context".
  filter(atlas_layer %in% c("treatment_context", "infection_context", "sensitivity", "context_only")) |>
  distinct(term, .keep_all = TRUE)
context_support <- compute_outcome_support(context_terms$term)

s5 <- context_terms |>
  left_join(context_support, by = "term") |>
  mutate(
    support_text = tolower(paste(domain, concept_family, abc_subrole, analysis_tier, atlas_layer)),
    layer_class = case_when(
      atlas_layer == "treatment_context" ~ "Treatment context",
      atlas_layer == "infection_context" ~ "Infection context",
      atlas_layer == "context_only" ~ "Comparator/context",
      str_detect(support_text, "measurement|modality|modifier|pct|vrs|hrct|quantitative") ~ "Measurement/sensitivity",
      str_detect(support_text, "clinical_outcome|outcome|hospital|transplant|respiratory failure|ventilation|oxygen|icu|onset") ~ "Outcome-adjacent/supporting descriptor",
      str_detect(support_text, "pathway|fibrosis|cytokine|oxidative|immune") ~ "Supporting pathway/sensitivity",
      TRUE ~ "Sensitivity/supporting"
    ),
    context_interpretation = case_when(
      atlas_layer == "treatment_context" ~ "Treatment-context term; retained separately and not interpreted as treatment effect or drug-risk estimate.",
      atlas_layer == "infection_context" ~ "Infection-context term; retained separately and not interpreted as primary disease-state bridge or infectious-causality estimate.",
      atlas_layer == "context_only" ~ "Comparator/context node; retained for corpus interpretation, not primary RA-ILD worsening claim.",
      layer_class == "Measurement/sensitivity" ~ "Measurement/modality, modifier or measurement-derived node; retained for context/sensitivity, not primary disease-state bridge claim.",
      layer_class == "Outcome-adjacent/supporting descriptor" ~ "Outcome-adjacent or severe-outcome descriptor; retained by prespecified rules but not interpreted as upstream bridge or biomarker.",
      layer_class == "Supporting pathway/sensitivity" ~ "Supporting pathway or aggregate signal; retained for transparency and biomarker/molecular summaries, not primary clinical bridge inference.",
      TRUE ~ "Sensitivity/supporting node; retained outside primary disease-state bridge claim."
    )
  ) |>
  arrange(layer_class, desc(total_outcome_cohit), term) |>
  transmute(
    term,
    preferred_label,
    layer_class,
    atlas_layer,
    domain,
    node_level,
    outcomes_supported,
    total_outcome_cohit,
    cohit_AE_ILD,
    cohit_progression,
    cohit_mortality,
    dominant_outcome,
    context_interpretation
  )

write_supp_csv(s5, "Supplementary_Table_S5_Context_Layer_Support")

# -----------------------------
# Supplementary Table S6 Biomarker/molecular support summary
# -----------------------------
bio_scope <- dic_cov |>
  mutate(scope_text = tolower(paste(abc_role, abc_subrole, domain, concept_family, node_level, analysis_tier, atlas_layer))) |>
  filter(
    str_detect(scope_text, "biomarker|serolog|cell|cytokine|chemokine|molecular|pathway|fibrosis|remodel|oxidative|vascular|stress|immune")
  ) |>
  filter(!atlas_layer %in% c("treatment_context", "infection_context", "context_only")) |>
  filter(!str_detect(scope_text, "measurement|modality|radiology_feature|pulmonary_function|physiology bridge")) |>
  distinct(term, .keep_all = TRUE)

bio_support <- compute_outcome_support(bio_scope$term)

s6 <- bio_scope |>
  left_join(bio_support, by = "term") |>
  mutate(
    concept_group = case_when(
      str_detect(tolower(paste(domain, concept_family)), "serolog|autoantibody") ~ "Serological host phenotype",
      str_detect(tolower(paste(domain, concept_family)), "epithelial|surfactant|biomarker") ~ "Epithelial/fibrosis biomarker",
      str_detect(tolower(paste(domain, concept_family)), "cell|cytokine|chemokine|inflamm") ~ "Inflammatory/cellular process",
      str_detect(tolower(paste(domain, concept_family)), "fibrosis|remodel|tgf|ecm|pathway") ~ "Fibrosis/remodelling pathway",
      str_detect(tolower(paste(domain, concept_family)), "oxidative|vascular|stress") ~ "Oxidative/vascular/stress",
      TRUE ~ "Other biomarker/molecular"
    )
  ) |>
  arrange(desc(total_outcome_cohit), desc(max_outcome_cohit), concept_group, term) |>
  transmute(
    term,
    preferred_label,
    concept_group,
    atlas_layer,
    analysis_tier,
    node_level,
    outcomes_supported,
    total_outcome_cohit,
    cohit_AE_ILD,
    cohit_progression,
    cohit_mortality,
    dominant_outcome,
    corpus_hit_n
  )

write_supp_csv(s6, "Supplementary_Table_S6_Biomarker_Molecular_Support_Summary")

# -----------------------------
# Supplementary Table S7 Internal robustness summaries
# -----------------------------
find_any <- function(pattern) newest_file(pattern, c(DIR_SUPP, DIR_TABLE, DIR_OUTPUT), recursive = TRUE, required = FALSE)

s9_term <- read_csv_safe(find_any("^Supplementary_Table_S9_TemporalSupport_Term_SUMMARY_.*\\.csv$|^Supplementary_Table_S7_Temporal_Term_Summary_.*\\.csv$"))
s9_edge <- read_csv_safe(find_any("^Supplementary_Table_S9_TemporalSupport_Edge_Provenance_.*\\.csv$|^Supplementary_Table_S7_Temporal_Edge_Summary_.*\\.csv$"))
s9_triad <- read_csv_safe(find_any("^Supplementary_Table_S9_TemporalSupport_Triad_Provenance_.*\\.csv$|^Supplementary_Table_S7_Temporal_Triad_Summary_.*\\.csv$"))
s10_term <- read_csv_safe(find_any("^Supplementary_Table_S10_ResamplingStability_Term_SUMMARY_.*\\.csv$|^Supplementary_Table_S7_Resampling_Term_Summary_.*\\.csv$"))
s10_edge <- read_csv_safe(find_any("^Supplementary_Table_S10_ResamplingStability_Edge_Stability_.*\\.csv$|^Supplementary_Table_S7_Resampling_Edge_Summary_.*\\.csv$"))
s10_triad <- read_csv_safe(find_any("^Supplementary_Table_S10_ResamplingStability_Triad_Stability_.*\\.csv$|^Supplementary_Table_S7_Resampling_Triad_Summary_.*\\.csv$"))

s7a <- tibble(
  analysis = c("Publication-era temporal support", "Repeated-resampling stability"),
  design = c(
    "Displayed terms, edges and triads reassessed in 1980-2020 and 2021-2025 publication eras without reselecting terms.",
    "Displayed terms, edges and triads reassessed across fixed-seed repeated 80% subsamples without replacement."
  ),
  interpretation = c(
    "Internal publication-era support check; not external or prospective validation.",
    "Internal sampling-stability check; not external or prospective validation."
  )
)


make_temporal_term_summary <- function(x) {
  if (nrow(x) == 0) return(tibble())
  # If already summarised by a previous S7 run, preserve counts.
  if (all(c("atlas_layer", "temporal_category", "concept_nodes") %in% names(x))) {
    return(x |> dplyr::select(atlas_layer, temporal_category, concept_nodes) |> dplyr::arrange(atlas_layer, temporal_category))
  }
  layer <- col_chr(x, c("atlas_layer", "layer"), default = "unspecified")
  cat <- col_chr(x, c("temporal_support_class", "temporal_category", "support_category"), default = "unspecified")
  tibble(atlas_layer = layer, temporal_category = cat) |>
    count(atlas_layer, temporal_category, name = "concept_nodes") |>
    arrange(atlas_layer, temporal_category)
}
make_temporal_edge_summary <- function(x) {
  if (nrow(x) == 0) return(tibble())
  if (all(c("temporal_category", "edges") %in% names(x))) {
    return(x |> dplyr::select(temporal_category, edges) |> dplyr::arrange(temporal_category))
  }
  cat <- col_chr(x, c("temporal_support_class", "temporal_category", "edge_temporal_class"), default = "unspecified")
  tibble(temporal_category = cat) |> count(temporal_category, name = "edges") |> arrange(temporal_category)
}
make_temporal_triad_summary <- function(x) {
  if (nrow(x) == 0) return(tibble())
  if (all(c("temporal_category", "triads") %in% names(x))) {
    return(x |> dplyr::select(temporal_category, triads) |> dplyr::arrange(temporal_category))
  }
  cat <- col_chr(x, c("AB_BC_temporal_class", "temporal_support_class", "temporal_category"), default = "unspecified")
  tibble(temporal_category = cat) |> count(temporal_category, name = "triads") |> arrange(temporal_category)
}
make_resampling_term_summary <- function(x) {
  if (nrow(x) == 0) return(tibble())
  if (all(c("atlas_layer", "stability_category", "concept_nodes") %in% names(x))) {
    return(x |> dplyr::select(atlas_layer, stability_category, concept_nodes) |> dplyr::arrange(atlas_layer, stability_category))
  }
  layer <- col_chr(x, c("atlas_layer", "layer"), default = "unspecified")
  cat <- col_chr(x, c("resampling_stability_category", "stability_category", "support_stability"), default = "")
  # If no precomputed category exists, derive one from term-outcome support frequency when available.
  if (all(cat == "" | is.na(cat))) {
    p_any <- col_num(x, c("p_any_support_gt0", "p_term_outcome_support_gt0", "p_support_gt0", "selection_frequency"), default = NA_real_)
    cat <- dplyr::case_when(
      !is.na(p_any) & p_any >= 0.80 ~ "stable >=80%",
      !is.na(p_any) & p_any >= 0.50 ~ "moderate 50-79%",
      !is.na(p_any) & p_any > 0 ~ "intermittent <50%",
      TRUE ~ "no_support_in_resamples"
    )
  }
  tibble(atlas_layer = layer, stability_category = cat) |>
    count(atlas_layer, stability_category, name = "concept_nodes") |>
    arrange(atlas_layer, stability_category)
}
make_resampling_edge_summary <- function(x) {
  if (nrow(x) == 0) return(tibble())
  if (all(c("stability_category", "edges") %in% names(x))) {
    return(x |> dplyr::select(stability_category, edges) |> dplyr::arrange(stability_category))
  }
  cat <- col_chr(x, c("resampling_stability_category", "stability_category", "support_stability"), default = "")
  if (all(cat == "" | is.na(cat))) {
    p_any <- col_num(x, c("p_edge_support_gt0", "p_support_gt0", "selection_frequency"), default = NA_real_)
    cat <- dplyr::case_when(
      !is.na(p_any) & p_any >= 0.80 ~ "stable >=80%",
      !is.na(p_any) & p_any >= 0.50 ~ "moderate 50-79%",
      !is.na(p_any) & p_any > 0 ~ "intermittent <50%",
      TRUE ~ "no_support_in_resamples"
    )
  }
  tibble(stability_category = cat) |> count(stability_category, name = "edges") |> arrange(stability_category)
}
make_resampling_triad_summary <- function(x) {
  if (nrow(x) == 0) return(tibble())
  if (all(c("stability_category", "triads") %in% names(x))) {
    return(x |> dplyr::select(stability_category, triads) |> dplyr::arrange(stability_category))
  }
  cat <- col_chr(x, c("resampling_stability_category", "stability_category", "support_stability"), default = "")
  if (all(cat == "" | is.na(cat))) {
    p_any <- col_num(x, c("p_AB_and_BC_support_gt0", "p_AB_BC_support_gt0", "p_support_gt0", "selection_frequency"), default = NA_real_)
    cat <- dplyr::case_when(
      !is.na(p_any) & p_any >= 0.80 ~ "stable >=80%",
      !is.na(p_any) & p_any >= 0.50 ~ "moderate 50-79%",
      !is.na(p_any) & p_any > 0 ~ "intermittent <50%",
      TRUE ~ "no_support_in_resamples"
    )
  }
  tibble(stability_category = cat) |> count(stability_category, name = "triads") |> arrange(stability_category)
}

s7b <- make_temporal_term_summary(s9_term)
s7c <- make_temporal_edge_summary(s9_edge)
s7d <- make_temporal_triad_summary(s9_triad)
s7e <- make_resampling_term_summary(s10_term)
s7f <- make_resampling_edge_summary(s10_edge)
s7g <- make_resampling_triad_summary(s10_triad)

write_supp_csv(s7a, "Supplementary_Table_S7A_Robustness_Design_Note")
write_supp_csv(s7b, "Supplementary_Table_S7B_Temporal_Term_Summary")
write_supp_csv(s7c, "Supplementary_Table_S7C_Temporal_Edge_Summary")
write_supp_csv(s7d, "Supplementary_Table_S7D_Temporal_Triad_Summary")
write_supp_csv(s7e, "Supplementary_Table_S7E_Resampling_Term_Summary")
write_supp_csv(s7f, "Supplementary_Table_S7F_Resampling_Edge_Summary")
write_supp_csv(s7g, "Supplementary_Table_S7G_Resampling_Triad_Summary")

# -----------------------------
# Human-readable PDF rendering
# -----------------------------
# The PDF is intentionally a concise, human-readable summary. Full source
# rows remain in the CSV outputs and repository source-data package.

str_shorten2 <- function(x, n = 42) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "..."), x)
}

fmt_num <- function(x, digits = 2) {
  if (is.numeric(x)) {
    out <- ifelse(is.na(x), "", format(round(x, digits), trim = TRUE, scientific = FALSE))
    return(out)
  }
  as.character(x)
}

compact_table <- function(df, cols = NULL, max_rows = 30, max_chars = 40) {
  df <- as_tibble(df)
  if (!is.null(cols)) {
    cols <- cols[cols %in% names(df)]
    df <- df[, cols, drop = FALSE]
  }
  if (nrow(df) > max_rows) df <- head(df, max_rows)
  if (ncol(df) == 0) return(tibble())
  df[] <- lapply(df, function(x) str_shorten2(fmt_num(x), max_chars))
  df
}

# page management using numeric NPC coordinates
page_no <- 0L
new_pdf_page <- function(title = NULL, subtitle = NULL) {
  page_no <<- page_no + 1L
  grid.newpage()
  if (!is.null(title) && nzchar(title)) {
    grid.text(title, x = unit(0.035, "npc"), y = unit(0.955, "npc"), just = c("left", "top"),
              gp = gpar(fontsize = 13.5, fontface = "bold"))
  }
  if (!is.null(subtitle) && nzchar(subtitle)) {
    grid.text(stringr::str_wrap(subtitle, width = 145), x = unit(0.035, "npc"), y = unit(0.91, "npc"),
              just = c("left", "top"), gp = gpar(fontsize = 8.5))
    y <- 0.855
  } else {
    y <- 0.90
  }
  grid.text(paste0("Page ", page_no), x = unit(0.975, "npc"), y = unit(0.025, "npc"),
            just = c("right", "bottom"), gp = gpar(fontsize = 7, col = "grey35"))
  y
}

draw_text_line <- function(txt, y, x = 0.035, fontsize = 8.2, fontface = "plain", width_chars = 145) {
  grid.text(stringr::str_wrap(txt, width = width_chars), x = unit(x, "npc"), y = unit(y, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = fontsize, fontface = fontface))
}

draw_table2 <- function(df, y_top, x = 0.035, width = 0.93,
                        row_h = 0.026, header_h = 0.032,
                        fontsize = 6.2, header_fontsize = 6.4,
                        col_widths = NULL, header_fill = "#EAF1FB") {
  df <- as_tibble(df)
  if (ncol(df) == 0 || nrow(df) == 0) {
    grid.text("No data available in the current output directory.",
              x = unit(x, "npc"), y = unit(y_top - 0.035, "npc"),
              just = c("left", "top"), gp = gpar(fontsize = fontsize + 1))
    return(y_top - 0.075)
  }
  nr <- nrow(df); nc <- ncol(df)
  if (is.null(col_widths)) col_widths <- rep(1 / nc, nc)
  col_widths <- col_widths / sum(col_widths) * width
  # Header background
  y <- y_top
  grid.rect(x = unit(x + width/2, "npc"), y = unit(y - header_h/2, "npc"),
            width = unit(width, "npc"), height = unit(header_h, "npc"),
            gp = gpar(fill = header_fill, col = "grey55", lwd = 0.45))
  xpos <- x
  for (j in seq_len(nc)) {
    grid.text(str_shorten2(names(df)[j], 26), x = unit(xpos + col_widths[j]/2, "npc"),
              y = unit(y - header_h/2, "npc"), just = "center",
              gp = gpar(fontsize = header_fontsize, fontface = "bold"))
    xpos <- xpos + col_widths[j]
  }
  y <- y - header_h
  for (i in seq_len(nr)) {
    fill <- if (i %% 2 == 0) "#F7F8FA" else "white"
    grid.rect(x = unit(x + width/2, "npc"), y = unit(y - row_h/2, "npc"),
              width = unit(width, "npc"), height = unit(row_h, "npc"),
              gp = gpar(fill = fill, col = "grey85", lwd = 0.25))
    xpos <- x
    for (j in seq_len(nc)) {
      grid.text(as.character(df[[j]][i]), x = unit(xpos + 0.004, "npc"),
                y = unit(y - row_h/2, "npc"), just = c("left", "center"),
                gp = gpar(fontsize = fontsize))
      xpos <- xpos + col_widths[j]
    }
    y <- y - row_h
  }
  y - 0.016
}

draw_section <- function(title, note, df, y, cols = NULL, max_rows = 30, max_chars = 40,
                         row_h = 0.026, fontsize = 6.2, col_widths = NULL,
                         title_size = 10.2, note_size = 7.6) {
  df2 <- compact_table(df, cols = cols, max_rows = max_rows, max_chars = max_chars)
  # Start a new page if there is clearly insufficient vertical space.
  needed <- 0.055 + 0.032 + row_h * max(nrow(df2), 1) + 0.03
  if ((y - needed) < 0.065) {
    y <- new_pdf_page(title = "Supplementary Tables S1-S7 (continued)")
  }
  grid.text(title, x = unit(0.035, "npc"), y = unit(y, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = title_size, fontface = "bold"))
  y <- y - 0.034
  if (!is.null(note) && nzchar(note)) {
    draw_text_line(note, y, fontsize = note_size, width_chars = 145)
    y <- y - 0.034
  }
  draw_table2(df2, y_top = y, row_h = row_h, fontsize = fontsize, header_fontsize = fontsize + 0.2,
              col_widths = col_widths)
}

pdf_path <- file.path(DIR_SUPP, paste0("Supplementary_Tables_S1_to_S7_HumanReadable_Summary_", CORPUS_TAG, "__", ANALYSIS_TAG, ".pdf"))
pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE, family = "Helvetica")

# S1
intro_note <- paste(
  "These human-readable supplementary tables summarise reproducibility, dictionary QC, corpus grounding, signed-effect direction,",
  "context-layer handling, biomarker/molecular support and internal robustness checks. Full machine-readable source tables,",
  "including long-format figure provenance, temporal-support and resampling-stability outputs, should be released in the public repository."
)
y <- new_pdf_page("Supplementary Tables S1-S7. Human-readable analysis summaries", intro_note)
y <- draw_section("Supplementary Table S1A. PubMed retrieval and corpus definition",
                  "Frozen original-article corpus and reproducibility inputs.",
                  s1a, y,
                  cols = c("item", "value", "note"), max_rows = 8, max_chars = 62,
                  row_h = 0.035, fontsize = 6.8,
                  col_widths = c(0.22, 0.26, 0.52))
y <- draw_section("Supplementary Table S1B. Recent annual corpus counts",
                  "Recent publication years shown for readability; full annual summary is provided as CSV.",
                  s1b, y,
                  cols = c("year", "total_records"), max_rows = 25, max_chars = 18,
                  row_h = 0.019, fontsize = 6.6,
                  col_widths = c(0.45, 0.55), title_size = 9.8)

# S2
y <- new_pdf_page("Supplementary Table S2. Final dictionary structure and QC",
                  "Concept nodes are dictionary entries, not raw word tokens; layer assignment was prespecified before downstream atlas displays.")
y <- draw_section("Supplementary Table S2A. Final dictionary and analysis-layer summary", "",
                  s2a, y, cols = c("item", "value", "note"), max_rows = 8, max_chars = 45,
                  row_h = 0.034, fontsize = 6.8, col_widths = c(0.25, 0.30, 0.45))
y <- draw_section("Supplementary Table S2B. Dictionary concept-node counts by layer and node level", "",
                  s2b, y, cols = c("atlas_layer", "node_level", "concept_nodes"), max_rows = 12, max_chars = 28,
                  row_h = 0.026, fontsize = 6.8, col_widths = c(0.40, 0.34, 0.26), title_size = 9.8)
y <- draw_section("Supplementary Table S2C. Dictionary quality-control summary", "",
                  s2c, y, cols = c("qc_domain", "summary", "interpretation"), max_rows = 8, max_chars = 55,
                  row_h = 0.033, fontsize = 6.3, col_widths = c(0.27, 0.36, 0.37), title_size = 9.8)

# S3
y <- new_pdf_page("Supplementary Table S3. Dictionary coverage and corpus grounding",
                  "Coverage was recomputed from the final dictionary and final PMID-level hit matrix.")
y <- draw_section("Supplementary Table S3A. Dictionary run and corpus-grounding summary", "",
                  s3a, y, cols = c("metric", "value"), max_rows = 22, max_chars = 70,
                  row_h = 0.022, fontsize = 6.5, col_widths = c(0.46, 0.54))
y <- draw_section("Supplementary Table S3B. Concept-node coverage by dictionary layer", "",
                  s3b, y, cols = c("atlas_layer", "concept_nodes", "nodes_with_hit", "median_hit_articles", "max_hit_articles"),
                  max_rows = 8, max_chars = 25, row_h = 0.024, fontsize = 6.5,
                  col_widths = c(0.28, 0.18, 0.18, 0.18, 0.18), title_size = 9.8)
y <- new_pdf_page("Supplementary Table S3C. Highest-frequency concept nodes",
                  "Top concept nodes by corpus hit count; full term coverage is repository source data.")
draw_section("Supplementary Table S3C. Highest-frequency concept nodes", "",
             s3c, y, cols = c("term", "preferred_label", "atlas_layer", "domain", "node_level", "corpus_hit_n", "corpus_hit_pct"),
             max_rows = 35, max_chars = 34, row_h = 0.020, fontsize = 5.7,
             col_widths = c(0.16, 0.18, 0.15, 0.19, 0.12, 0.10, 0.10), title_size = 0.1, note_size = 0.1)

# S4
y <- new_pdf_page("Supplementary Table S4. Signed-effect direction summary",
                  "Direction-bearing wording is summarised by scope and outcome; full signed-effect source tables are repository data.")
y <- draw_section("Supplementary Table S4A. Signed-effect direction by scope and outcome", "",
                  s4a, y,
                  cols = c("scope", "outcome", "pairs", "pairs_with_directional_articles", "median_article_support", "positive_leaning_pairs", "decreasing_or_protective_leaning_pairs", "mixed_or_zero_pairs"),
                  max_rows = 12, max_chars = 30, row_h = 0.031, fontsize = 5.9,
                  col_widths = c(0.15, 0.13, 0.08, 0.16, 0.14, 0.12, 0.13, 0.09))
y <- new_pdf_page("Supplementary Table S4B. Primary-atlas signed-effect term-outcome pairs",
                  "Top pairs by article support within each outcome; full signed-effect outputs are repository data.")
draw_section("Supplementary Table S4B. Primary-atlas signed-effect term-outcome pairs", "",
             s4b, y, cols = c("outcome", "preferred_label", "domain", "node_level", "signed_balance", "article_support", "positive_wording", "decreasing_or_protective_wording"),
             max_rows = 36, max_chars = 31, row_h = 0.020, fontsize = 5.6,
             col_widths = c(0.11, 0.22, 0.18, 0.11, 0.10, 0.10, 0.09, 0.09), title_size = 0.1, note_size = 0.1)

# S5
y <- new_pdf_page("Supplementary Table S5. Treatment, infection, measurement and supporting/context-layer support",
                  "Non-primary treatment, infection, comparator, measurement and supporting nodes are retained for transparency but not interpreted as primary disease-state bridge claims.")
s5_display <- s5 |>
  arrange(layer_class, desc(total_outcome_cohit), term)
y <- draw_section("Supplementary Table S5. Context/supporting-layer support summary", "",
                  s5_display, y,
                  cols = c("preferred_label", "layer_class", "domain", "outcomes_supported", "total_outcome_cohit", "cohit_AE_ILD", "dominant_outcome"),
             max_rows = 46, max_chars = 36, row_h = 0.0185, fontsize = 5.5,
             col_widths = c(0.24, 0.21, 0.17, 0.17, 0.08, 0.06, 0.07), title_size = 0.1, note_size = 0.1)

# S6
y <- new_pdf_page("Supplementary Table S6. Biomarker/molecular signal support summary",
                  "Biomarker, serological, cellular, molecular and pathway scope summary for Figures 5-6; full source tables are repository data.")
y <- draw_section("Supplementary Table S6. Biomarker/molecular signal support summary", "",
                  s6, y, cols = c("preferred_label", "concept_group", "atlas_layer", "analysis_tier", "outcomes_supported", "total_outcome_cohit", "cohit_AE_ILD", "cohit_progression", "cohit_mortality"),
             max_rows = 42, max_chars = 34, row_h = 0.019, fontsize = 5.4,
             col_widths = c(0.20, 0.20, 0.11, 0.10, 0.17, 0.07, 0.05, 0.05, 0.05), title_size = 0.1, note_size = 0.1)

# S7
y <- new_pdf_page("Supplementary Table S7. Internal temporal and resampling robustness summaries",
                  "Temporal support and resampling stability are internal robustness checks, not external validation.")
y <- draw_section("Supplementary Table S7A. Internal robustness-check design", "",
                  s7a, y, cols = c("analysis", "design", "interpretation"), max_rows = 4, max_chars = 65,
                  row_h = 0.040, fontsize = 6.4, col_widths = c(0.22, 0.45, 0.33))
y <- draw_section("Supplementary Table S7B. Publication-era temporal support by layer", "Earlier era: 1980-2020; recent era: 2021-2025.",
                  s7b, y, cols = c("atlas_layer", "temporal_category", "concept_nodes"), max_rows = 20, max_chars = 28,
                  row_h = 0.021, fontsize = 6.2, col_widths = c(0.40, 0.38, 0.22), title_size = 9.6)
y <- new_pdf_page("Supplementary Table S7 (continued)",
                  "Edge/triad summaries use the displayed figure-level objects eligible for each internal robustness check.")
y <- draw_section("Supplementary Table S7C. Publication-era temporal support for network edges", "",
                  s7c, y, cols = c("temporal_category", "edges"), max_rows = 8, max_chars = 28,
                  row_h = 0.030, fontsize = 6.7, col_widths = c(0.65, 0.35))
y <- draw_section("Supplementary Table S7D. Publication-era temporal support for ABC triads", "",
                  s7d, y, cols = c("temporal_category", "triads"), max_rows = 8, max_chars = 32,
                  row_h = 0.030, fontsize = 6.7, col_widths = c(0.70, 0.30))
y <- draw_section("Supplementary Table S7E. Repeated-resampling stability by layer", "",
                  s7e, y, cols = c("atlas_layer", "stability_category", "concept_nodes"), max_rows = 14, max_chars = 28,
                  row_h = 0.023, fontsize = 6.3, col_widths = c(0.40, 0.38, 0.22))
y <- draw_section("Supplementary Table S7F. Repeated-resampling stability for network edges", "",
                  s7f, y, cols = c("stability_category", "edges"), max_rows = 6, max_chars = 28,
                  row_h = 0.030, fontsize = 6.7, col_widths = c(0.70, 0.30))
y <- draw_section("Supplementary Table S7G. Repeated-resampling stability for ABC triads", "",
                  s7g, y, cols = c("stability_category", "triads"), max_rows = 6, max_chars = 28,
                  row_h = 0.030, fontsize = 6.7, col_widths = c(0.70, 0.30))

dev.off()
log_msg("WROTE:", pdf_path)
log_msg("=== DONE 15_Build_SupplementaryTables_S1_to_S7_AtlasV2 ===")
