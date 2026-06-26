# =========================================================
# 02_BuildHitsMatrix_AtlasV2.R
# RA-ILD atlas v2 pipeline
#
# Purpose
#   - Load the main original article corpus.
#   - Load the v2 outcome-preserving atlas dictionary.
#   - Build a document-level hit matrix for every dictionary term.
#   - Preserve dictionary-level regex case handling by using
#     ignore.case = FALSE externally.
#   - Write dictionary coverage and composition summaries for Figure 1, including treatment- and infection-context layers.
#
# Important
#   The hit matrix is built for the full v2 dictionary.  Primary downstream
#   analyses should use dictionary metadata columns, especially:
#     include_primary, analysis_tier, node_level, abc_role, abc_subrole,
#     domain, concept_family, atlas_layer, and interpretation_note.
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)

quiet_install(c("dplyr", "readr", "stringr", "tibble", "purrr"))

STAMP_DATE <- stamp_date()

# ------------------ 1) Article loader ------------------
manual_article_file <- Sys.getenv("RAILD_ARTICLE_FILE", unset = "")

make_article_candidates <- function(dir_proc) {
  files <- unique(c(
    Sys.glob(file.path(dir_proc, "articles_main_original_*.csv")),
    Sys.glob(file.path(dir_proc, "articles_original_*.csv")),
    Sys.glob(file.path(dir_proc, "articles_*.csv"))
  ))
  files <- files[file.exists(files)]
  if (!length(files)) return(tibble::tibble())

  bn <- basename(files)
  type <- dplyr::case_when(
    grepl("^articles_main_original_[0-9]{8}\\.csv$", bn) ~ "main_original_named",
    grepl("^articles_original_[0-9]{8}\\.csv$", bn)      ~ "original_named",
    grepl("^articles_[0-9]{8}\\.csv$", bn)               ~ "legacy_main",
    TRUE ~ NA_character_
  )

  tibble::tibble(file = files, basename = bn, type = type) |>
    dplyr::filter(!is.na(type)) |>
    dplyr::mutate(
      file_date = stringr::str_extract(basename, "[0-9]{8}"),
      priority = dplyr::case_when(
        type == "main_original_named" ~ 1L,
        type == "original_named"      ~ 2L,
        type == "legacy_main"         ~ 3L,
        TRUE ~ 99L
      )
    ) |>
    dplyr::arrange(dplyr::desc(file_date), priority, basename)
}

if (nzchar(manual_article_file)) {
  manual_article_file <- path.expand(manual_article_file)
  if (!file.exists(manual_article_file)) stop("Manual RAILD_ARTICLE_FILE does not exist: ", manual_article_file)
  f_articles_csv <- manual_article_file
  log_msg("Using manually specified article file from RAILD_ARTICLE_FILE:", f_articles_csv)
} else {
  cand_tbl <- make_article_candidates(DIR_PROC)
  if (nrow(cand_tbl) == 0) {
    stop(
      "No main original article CSV found in ", DIR_PROC, ".\n",
      "Expected one of: articles_main_original_YYYYMMDD.csv, articles_original_YYYYMMDD.csv, or articles_YYYYMMDD.csv.\n",
      "Run 01_FetchPubMed_AtlasV2.R first, or set Sys.setenv(RAILD_ARTICLE_FILE = '...')."
    )
  }
  log_msg("Article candidates considered:", paste(paste0(cand_tbl$basename, "[", cand_tbl$type, "]"), collapse = " / "))
  f_articles_csv <- cand_tbl$file[1]
}

articles <- readr::read_csv(f_articles_csv, show_col_types = FALSE)
log_msg("02 loaded articles:", f_articles_csv, " n=", nrow(articles))

if ("corpus_role" %in% names(articles)) {
  n_before_role_filter <- nrow(articles)
  articles <- articles |> dplyr::filter(corpus_role == "main_original")
  if (nrow(articles) != n_before_role_filter) {
    log_msg("Filtered by corpus_role == main_original:", n_before_role_filter, "->", nrow(articles))
  }
}

if (!"text" %in% names(articles)) {
  if (all(c("title", "abstract") %in% names(articles))) {
    articles <- articles |>
      dplyr::mutate(text = paste(ifelse(is.na(title), "", title), ifelse(is.na(abstract), "", abstract), sep = " "))
    log_msg("Created text column from title + abstract.")
  } else {
    stop("Article file has no text column and lacks title/abstract columns: ", f_articles_csv)
  }
}

required_article_cols <- c("pmid", "year", "text")
missing_article_cols <- setdiff(required_article_cols, names(articles))
if (length(missing_article_cols)) {
  stop("Article file is missing required columns: ", paste(missing_article_cols, collapse = ", "), "\nFile: ", f_articles_csv)
}
if (!"journal" %in% names(articles)) articles$journal <- NA_character_
if (!"pubtype" %in% names(articles)) articles$pubtype <- NA_character_
if (!"title" %in% names(articles)) articles$title <- NA_character_
if (!"abstract" %in% names(articles)) articles$abstract <- NA_character_

if (nrow(articles) == 0) stop("No main_original articles remain after loading/filtering: ", f_articles_csv)

articles <- articles |>
  dplyr::mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    text = ifelse(is.na(text), "", text)
  )

if (anyDuplicated(articles$pmid)) {
  n_before_dedup <- nrow(articles)
  articles <- articles |> dplyr::distinct(pmid, .keep_all = TRUE)
  log_msg("Removed duplicated PMIDs:", n_before_dedup, "->", nrow(articles))
}
log_msg("Selected main article corpus:", basename(f_articles_csv), " n_main=", nrow(articles))

# ------------------ 2) Dictionary loader and validation ------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)

needed_v2_cols <- c("concept_id", "preferred_label", "abc_role", "abc_subrole", "domain", "concept_family", "node_level", "analysis_tier", "atlas_layer", "include_primary", "interpretation_note")
missing_v2_cols <- setdiff(needed_v2_cols, names(dic))
if (length(missing_v2_cols)) {
  log_msg("WARNING: v2 metadata columns missing:", paste(missing_v2_cols, collapse = ", "))
}

bad_regex <- which(is.na(dic$regex) | !nzchar(dic$regex))
if (length(bad_regex)) stop("Empty regex found for terms: ", paste(dic$term[bad_regex], collapse = ", "))

# Compile-test regex before scanning.
regex_ok <- vapply(seq_len(nrow(dic)), function(i) {
  ok <- tryCatch({ grepl(dic$regex[i], "REGEX_COMPILE_TEST", ignore.case = FALSE, perl = TRUE); TRUE }, error = function(e) FALSE)
  ok
}, logical(1))
if (any(!regex_ok)) {
  stop("Invalid regex found for terms: ", paste(dic$term[!regex_ok], collapse = ", "))
}

log_msg("Using dictionary:", basename(DIC_FILE), " n=", nrow(dic))
log_msg("Primary dictionary terms:", nrow(select_primary_dictionary_terms(dic)))
if ("atlas_layer" %in% names(dic)) log_msg("Dictionary atlas layers:", paste(names(table(dic$atlas_layer, useNA = "ifany")), as.integer(table(dic$atlas_layer, useNA = "ifany")), sep = "=", collapse = " / "))
log_msg("DIC_TAG:", DIC_TAG)
log_msg("DIC head:", paste(head(dic$term, 6), collapse = " / "))

# ------------------ 3) Hit matrix ------------------------
match_one <- function(pattern, texts) {
  ok <- !is.na(texts)
  out <- integer(length(texts))
  if (any(ok)) {
    out[ok] <- as.integer(grepl(pattern, texts[ok], ignore.case = REGEX_IGNORE_CASE_EXTERNALLY, perl = TRUE))
  }
  out
}

for (i in seq_len(nrow(dic))) {
  nm <- paste0("hit__", dic$term[i])
  if (!nm %in% names(articles)) {
    articles[[nm]] <- match_one(dic$regex[i], articles$text)
  }
}

# v2 no longer requires _var aggregation, but keep a defensive compatibility block.
var_childs <- unique(dic$term[grepl("_var$", dic$term)])
term_groups <- tibble::tibble(child = var_childs, parent = sub("_var$", "", var_childs)) |>
  dplyr::distinct(child, .keep_all = TRUE)
if (nrow(term_groups)) {
  for (k in seq_len(nrow(term_groups))) {
    ch <- term_groups$child[k]
    pa <- term_groups$parent[k]
    ch_col <- paste0("hit__", ch)
    pa_col <- paste0("hit__", pa)
    if (ch_col %in% names(articles)) {
      if (!pa_col %in% names(articles)) {
        articles[[pa_col]] <- articles[[ch_col]]
      } else {
        articles[[pa_col]] <- as.integer((articles[[pa_col]] == 1L) | (articles[[ch_col]] == 1L))
      }
    }
  }
}

# Mutual exclusivity / specificity guards retained from prior pipeline.
safe_col <- function(nm) {
  if (nm %in% names(articles)) as.integer(articles[[nm]]) else integer(nrow(articles))
}
z <- function(nm) safe_col(nm)

if (all(c("hit__ACPA_neg", "hit__ACPA_pos") %in% names(articles))) {
  idx_neg <- which(z("hit__ACPA_neg") == 1L)
  if (length(idx_neg)) {
    articles$hit__ACPA_pos[idx_neg] <- 0L
    if ("hit__ACPA_RF_high" %in% names(articles)) articles$hit__ACPA_RF_high[idx_neg] <- 0L
  }
}
if (all(c("hit__CRP_high", "hit__CRP_low") %in% names(articles))) {
  both <- which(z("hit__CRP_high") == 1L & z("hit__CRP_low") == 1L)
  if (length(both)) articles$hit__CRP_low[both] <- 0L
}
if (all(c("hit__ESR_high", "hit__ESR_low") %in% names(articles))) {
  both <- which(z("hit__ESR_high") == 1L & z("hit__ESR_low") == 1L)
  if (length(both)) articles$hit__ESR_low[both] <- 0L
}

hit_cols <- paste0("hit__", dic$term)
hit_cols <- intersect(hit_cols, names(articles))
articles_with_any_hit <- if (length(hit_cols)) rowSums(articles[, hit_cols, drop = FALSE], na.rm = TRUE) > 0 else rep(FALSE, nrow(articles))

# ------------------ 4) Coverage and dictionary summaries ------------------
term_coverage <- tibble::tibble(term = dic$term) |>
  dplyr::mutate(
    hit_col = paste0("hit__", term),
    corpus_hit_n_recomputed = vapply(hit_col, function(col) if (col %in% names(articles)) sum(articles[[col]] == 1L, na.rm = TRUE) else 0L, integer(1)),
    corpus_hit_pct_recomputed = round(100 * corpus_hit_n_recomputed / nrow(articles), 2)
  ) |>
  dplyr::left_join(dic, by = "term") |>
  dplyr::select(
    term, dplyr::any_of(c("concept_id", "preferred_label", "class", "abc_role", "abc_subrole", "domain", "concept_family", "node_level", "analysis_tier", "atlas_layer", "include_primary", "primary_atlas_eligible")),
    corpus_hit_n_recomputed, corpus_hit_pct_recomputed,
    dplyr::any_of(c("corpus_hit_n", "corpus_hit_pct", "regex_qc_action", "retain_reason", "interpretation_note", "notes"))
  ) |>
  dplyr::arrange(dplyr::desc(corpus_hit_n_recomputed), term)

class_coverage <- term_coverage |>
  dplyr::group_by(class, abc_role, abc_subrole, domain, analysis_tier, atlas_layer, node_level, include_primary) |>
  dplyr::summarise(
    terms = dplyr::n(),
    hit_terms = sum(corpus_hit_n_recomputed > 0, na.rm = TRUE),
    total_mentions = sum(corpus_hit_n_recomputed, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(analysis_tier, dplyr::desc(total_mentions), class)

safe_count_layer <- function(layer) {
  if (!"atlas_layer" %in% names(dic)) return(0L)
  sum(dic$atlas_layer == layer, na.rm = TRUE)
}
safe_hit_layer <- function(layer) {
  if (!"atlas_layer" %in% names(term_coverage)) return(0L)
  sum(term_coverage$atlas_layer == layer & term_coverage$corpus_hit_n_recomputed > 0, na.rm = TRUE)
}

primary_terms_tbl <- select_primary_dictionary_terms(dic)

dictionary_run_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  dic_file = basename(DIC_FILE),
  article_file = basename(f_articles_csv),
  n_articles = nrow(articles),
  n_dictionary_terms = nrow(dic),
  n_primary_terms = nrow(primary_terms_tbl),
  n_primary_atlas_terms = safe_count_layer("primary_atlas"),
  n_treatment_context_terms = safe_count_layer("treatment_context"),
  n_infection_context_terms = safe_count_layer("infection_context"),
  n_sensitivity_terms = safe_count_layer("sensitivity"),
  n_exploratory_terms = safe_count_layer("exploratory"),
  n_context_only_terms = safe_count_layer("context_only"),
  n_other_unspecified_terms = if ("atlas_layer" %in% names(dic)) sum(is.na(dic$atlas_layer) | !nzchar(as.character(dic$atlas_layer)), na.rm = TRUE) else NA_integer_,
  n_terms_with_hit = sum(term_coverage$corpus_hit_n_recomputed > 0, na.rm = TRUE),
  n_primary_terms_with_hit = sum(
    term_coverage$include_primary == TRUE &
      term_coverage$analysis_tier == "main" &
      term_coverage$node_level == "specific" &
      term_coverage$atlas_layer == "primary_atlas" &
      term_coverage$corpus_hit_n_recomputed > 0,
    na.rm = TRUE
  ),
  n_treatment_context_terms_with_hit = safe_hit_layer("treatment_context"),
  n_infection_context_terms_with_hit = safe_hit_layer("infection_context"),
  articles_with_any_dictionary_hit = sum(articles_with_any_hit),
  abstract_coverage_pct = round(100 * mean(articles_with_any_hit), 2),
  total_dictionary_mentions = sum(term_coverage$corpus_hit_n_recomputed, na.rm = TRUE),
  regex_case_policy = "dictionary-scoped; external ignore.case=FALSE",
  primary_selection_rule = "include_primary==TRUE & atlas_layer==primary_atlas & analysis_tier==main & node_level==specific & abc_role in A/B/C; treatment- and infection-context terms are retained separately and excluded from primary disease-state bridge claims"
)

f_term_cov <- file.path(DIR_TABLE, sprintf("dictionary_term_coverage_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_class_cov <- file.path(DIR_TABLE, sprintf("dictionary_class_coverage_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_run_sum <- file.path(DIR_TABLE, sprintf("dictionary_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(term_coverage, f_term_cov)
readr::write_csv(class_coverage, f_class_cov)
readr::write_csv(dictionary_run_summary, f_run_sum)
log_msg("WROTE:", f_term_cov)
log_msg("WROTE:", f_class_cov)
log_msg("WROTE:", f_run_sum)

# Also save a copy under data_proc for easy inspection.
readr::write_csv(term_coverage, file.path(DIR_PROC, sprintf("dictionary_term_coverage_%s__%s.csv", CORPUS_TAG, DIC_TAG)))
readr::write_csv(dictionary_run_summary, file.path(DIR_PROC, sprintf("dictionary_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

# ------------------ 5) Save hit matrix ------------------
keep_cols <- c("pmid", "year", "journal", "pubtype", "title", "abstract", "text", hit_cols)
keep_cols <- intersect(keep_cols, names(articles))
hits <- articles |> dplyr::select(dplyr::all_of(keep_cols))

f_hits_proc_tag <- file.path(DIR_PROC, sprintf("hits_matrix_%s__%s.csv", STAMP_DATE, DIC_TAG))
f_hits_tag <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
readr::write_csv(hits, f_hits_proc_tag)
readr::write_csv(hits, f_hits_tag)

log_msg("WROTE:", f_hits_proc_tag)
log_msg("WROTE:", f_hits_tag)
log_msg("Hit matrix n_articles=", nrow(hits), " n_terms=", length(hit_cols))
log_msg("Articles with >=1 dictionary hit=", sum(articles_with_any_hit), " coverage_pct=", dictionary_run_summary$abstract_coverage_pct)
log_msg("=== DONE 02_BuildHitsMatrix_AtlasV2 ===")
