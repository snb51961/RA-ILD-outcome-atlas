# 13_TemporalSupport_AtlasV2.R
# RA-ILD Atlas v2.2: figure-level temporal support analysis
#
# Purpose:
#   Internal robustness / publication-era support check for the final figure-level atlas.
#   This is NOT prospective validation and does NOT reselect the main atlas.
#   It recalculates support for already displayed figure terms, edges and ABC triads
#   in earlier and recent publication eras using the final dictionary and unchanged hit matrix.
#
# Default eras:
#   earlier era: 1980-2020
#   recent era : 2021-2025
#
# Required to have already run:
#   00, 02, 03-06, 07-11, 12
#
# Main inputs:
#   dic/ra_ild_dictionary_analysis_v2_outcome_preserving_atlas.csv
#   data_proc/articles_main_original_20260323.csv
#   output/.../table/hits_matrix_*.csv
#   output/.../table/Supplementary_Table_Figure_Term_Provenance_SUMMARY_*.csv
#   output/.../table/Supplementary_Table_Figure_Edge_Provenance_*.csv  (optional)
#   output/.../table/Supplementary_Table_Figure_Triad_Provenance_*.csv (optional)
#
# Outputs:
#   Supplementary_Table_S9_TemporalSupport_TermOutcome_LONG_*.csv
#   Supplementary_Table_S9_TemporalSupport_Term_SUMMARY_*.csv
#   Supplementary_Table_S9_TemporalSupport_Edge_Provenance_*.csv
#   Supplementary_Table_S9_TemporalSupport_Triad_Provenance_*.csv
#   Supplementary_Table_S9_TemporalSupport_SelectionRules_*.csv
#   TemporalSupport_run_summary_*.csv
#   SupplementaryFigure_TemporalSupport_*.pdf/png

options(stringsAsFactors = FALSE)

# ---------------------------
# Setup
# ---------------------------
CODE_DIR <- Sys.getenv("RAILD_CODE_DIR", unset = getwd())
if (dir.exists(CODE_DIR)) setwd(CODE_DIR)
ROOT <- Sys.getenv("RAILD_ROOT", unset = "")
if (!nzchar(ROOT)) {
  candidate_roots <- unique(c(getwd(), dirname(getwd())))
  has_project_dirs <- vapply(candidate_roots, function(p) {
    dir.exists(file.path(p, "dic")) ||
      dir.exists(file.path(p, "data_proc")) ||
      dir.exists(file.path(p, "data_raw"))
  }, logical(1))
  ROOT <- if (any(has_project_dirs)) candidate_roots[which(has_project_dirs)[1]] else getwd()
}
ROOT <- normalizePath(ROOT, mustWork = FALSE)
PM_RANGE <- Sys.getenv("RAILD_PM_RANGE", unset = "pm_1980_20251231")
DIC_TAG <- Sys.getenv("RAILD_DIC_TAG", unset = "analysis_v2_outcome_preserving_atlas")

if (file.exists("00_setup_AtlasV2.R")) {
  suppressWarnings(source("00_setup_AtlasV2.R"))
  if (exists("ROOT", inherits = TRUE)) ROOT <- get("ROOT", inherits = TRUE)
  if (exists("PM_RANGE", inherits = TRUE)) PM_RANGE <- get("PM_RANGE", inherits = TRUE)
  if (exists("DIC_TAG", inherits = TRUE)) DIC_TAG <- get("DIC_TAG", inherits = TRUE)
}

OUT_BASE <- file.path(ROOT, "output", PM_RANGE, DIC_TAG)
DIR_TABLE <- file.path(OUT_BASE, "table")
DIR_FIG <- file.path(OUT_BASE, "fig")
dir.create(DIR_TABLE, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
stop_if_missing <- function(path, label = path) {
  if (length(path) == 0 || is.na(path) || !file.exists(path)) stop("Missing required file: ", label, call. = FALSE)
  invisible(path)
}

read_csv_safe <- function(path) {
  stop_if_missing(path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv_safe <- function(x, filename) {
  utils::write.csv(x, file.path(DIR_TABLE, filename), row.names = FALSE, na = "")
  invisible(file.path(DIR_TABLE, filename))
}

find_latest <- function(pattern, dirs = c(DIR_TABLE, file.path(ROOT, "data_proc"), ROOT), required = TRUE) {
  files <- unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) character(0) else list.files(d, pattern = pattern, full.names = TRUE)
  }))
  files <- files[file.exists(files)]
  if (!length(files)) {
    if (required) stop("Could not find file matching pattern: ", pattern, call. = FALSE)
    return(NA_character_)
  }
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

sanitize_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

# ---------------------------
# Inputs
# ---------------------------
DIC_FILE <- file.path(ROOT, "dic", "ra_ild_dictionary_analysis_v2_outcome_preserving_atlas.csv")
ARTICLE_FILE <- Sys.getenv(
  "RAILD_ARTICLE_FILE",
  unset = file.path(ROOT, "data_proc", "articles_main_original_20260323.csv")
)
HITS_FILE <- find_latest(paste0("^hits_matrix_.*__", DIC_TAG, "\\.csv$"), dirs = c(DIR_TABLE), required = TRUE)
TERM_PROV_FILE <- find_latest(paste0("^Supplementary_Table_Figure_Term_Provenance_SUMMARY_.*__", DIC_TAG, "\\.csv$"), dirs = c(DIR_TABLE), required = TRUE)
EDGE_PROV_FILE <- find_latest(paste0("^Supplementary_Table_Figure_Edge_Provenance_.*__", DIC_TAG, "\\.csv$"), dirs = c(DIR_TABLE), required = FALSE)
TRIAD_PROV_FILE <- find_latest(paste0("^Supplementary_Table_Figure_Triad_Provenance_.*__", DIC_TAG, "\\.csv$"), dirs = c(DIR_TABLE), required = FALSE)

log_msg("Dictionary:", DIC_FILE)
log_msg("Articles:", ARTICLE_FILE)
log_msg("Hits:", HITS_FILE)
log_msg("Term provenance:", TERM_PROV_FILE)
log_msg("Edge provenance:", EDGE_PROV_FILE)
log_msg("Triad provenance:", TRIAD_PROV_FILE)

dic <- read_csv_safe(DIC_FILE)
articles <- read_csv_safe(ARTICLE_FILE)
hits <- read_csv_safe(HITS_FILE)
term_prov <- read_csv_safe(TERM_PROV_FILE)
edge_prov <- if (!is.na(EDGE_PROV_FILE)) read_csv_safe(EDGE_PROV_FILE) else data.frame()
triad_prov <- if (!is.na(TRIAD_PROV_FILE)) read_csv_safe(TRIAD_PROV_FILE) else data.frame()

# ---------------------------
# Standardise identifiers
# ---------------------------
pick_col <- function(df, candidates) {
  nm <- names(df)
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[1] else NA_character_
}

pmid_col_hits <- pick_col(hits, c("pmid", "PMID", "article_id", "id"))
pmid_col_articles <- pick_col(articles, c("pmid", "PMID", "article_id", "id"))
year_col_articles <- pick_col(articles, c("year", "pub_year", "publication_year", "PublicationYear"))
year_col_hits <- pick_col(hits, c("year", "pub_year", "publication_year", "PublicationYear"))

if (is.na(pmid_col_hits)) stop("Could not identify PMID/article id column in hits matrix.", call. = FALSE)
if (is.na(year_col_hits)) {
  if (is.na(pmid_col_articles) || is.na(year_col_articles)) {
    stop("Hits matrix has no year column and article file lacks usable PMID/year columns.", call. = FALSE)
  }
  article_year <- articles[, c(pmid_col_articles, year_col_articles)]
  names(article_year) <- c("pmid_std", "year_std")
  article_year$pmid_std <- as.character(article_year$pmid_std)
  hits$pmid_std <- as.character(hits[[pmid_col_hits]])
  hits <- merge(hits, article_year, by = "pmid_std", all.x = TRUE)
  hits$year_temporal <- suppressWarnings(as.integer(hits$year_std))
} else {
  hits$year_temporal <- suppressWarnings(as.integer(hits[[year_col_hits]]))
}

if (any(is.na(hits$year_temporal))) {
  warning("Some hit matrix rows have missing publication years; these rows will be included in full support but not era-specific support.")
}

EARLY_START <- as.integer(Sys.getenv("RAILD_TEMPORAL_EARLY_START", unset = "1980"))
EARLY_END <- as.integer(Sys.getenv("RAILD_TEMPORAL_EARLY_END", unset = "2020"))
RECENT_START <- as.integer(Sys.getenv("RAILD_TEMPORAL_RECENT_START", unset = "2021"))
RECENT_END <- as.integer(Sys.getenv("RAILD_TEMPORAL_RECENT_END", unset = "2025"))

mask_full <- rep(TRUE, nrow(hits))
mask_early <- !is.na(hits$year_temporal) & hits$year_temporal >= EARLY_START & hits$year_temporal <= EARLY_END
mask_recent <- !is.na(hits$year_temporal) & hits$year_temporal >= RECENT_START & hits$year_temporal <= RECENT_END

# ---------------------------
# Hit column resolution
# ---------------------------
hit_column_candidates <- function(term) {
  c(
    paste0("hit__", term),
    paste0("hit__", make.names(term)),
    term,
    make.names(term)
  )
}

normalise_key <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]+", "", as.character(x)))
}

dic_term_by_label <- setNames(dic$term, normalise_key(dic$preferred_label))
dic_term_by_term <- setNames(dic$term, normalise_key(dic$term))

resolve_term_to_dictionary_term <- function(x) {
  x <- as.character(x)
  if (length(x) != 1 || is.na(x) || x == "") return(NA_character_)
  if (x %in% dic$term) return(x)
  key <- normalise_key(x)
  if (key %in% names(dic_term_by_term)) return(unname(dic_term_by_term[[key]]))
  if (key %in% names(dic_term_by_label)) return(unname(dic_term_by_label[[key]]))
  NA_character_
}

resolve_hit_col <- function(x) {
  term <- resolve_term_to_dictionary_term(x)
  candidates <- unique(c(hit_column_candidates(x), if (!is.na(term)) hit_column_candidates(term) else character(0)))
  candidates <- candidates[candidates %in% names(hits)]
  if (length(candidates)) candidates[1] else NA_character_
}

as_hit <- function(v) {
  if (is.logical(v)) return(v)
  if (is.numeric(v) || is.integer(v)) return(!is.na(v) & v > 0)
  vv <- tolower(trimws(as.character(v)))
  vv %in% c("true", "t", "1", "yes", "y")
}

get_hit_vec <- function(x) {
  col <- resolve_hit_col(x)
  if (is.na(col)) return(rep(FALSE, nrow(hits)))
  as_hit(hits[[col]])
}

count_pair <- function(a, b, mask) {
  va <- get_hit_vec(a)
  vb <- get_hit_vec(b)
  sum(va & vb & mask, na.rm = TRUE)
}

count_single <- function(a, mask) {
  va <- get_hit_vec(a)
  sum(va & mask, na.rm = TRUE)
}

count_triad <- function(a, b, c, mask) {
  va <- get_hit_vec(a)
  vb <- get_hit_vec(b)
  vc <- get_hit_vec(c)
  sum(va & vb & vc & mask, na.rm = TRUE)
}

OUTCOMES <- c("AE-ILD", "progression", "mortality")
OUTCOME_LABELS <- c("AE-ILD", "Progression", "Mortality")

# ---------------------------
# Term outcome temporal support
# ---------------------------
term_col <- pick_col(term_prov, c("term", "concept", "concept_node", "node", "preferred_label"))
if (is.na(term_col)) stop("Could not identify term column in term provenance summary.", call. = FALSE)

term_meta <- term_prov
term_meta$term_raw_from_provenance <- term_meta[[term_col]]
term_meta$term <- vapply(term_meta$term_raw_from_provenance, resolve_term_to_dictionary_term, character(1))
term_meta$term[is.na(term_meta$term)] <- as.character(term_meta$term_raw_from_provenance[is.na(term_meta$term)])
term_meta <- term_meta[!is.na(term_meta$term) & term_meta$term != "", , drop = FALSE]
term_meta <- term_meta[!duplicated(term_meta$term), , drop = FALSE]

# Attach dictionary metadata, preserving provenance if present
meta_keep <- c("term", "concept_id", "preferred_label", "abc_role", "abc_subrole", "domain", "concept_family", "node_level", "analysis_tier", "atlas_layer", "include_primary", "primary_atlas_eligible")
dic_meta <- dic[, intersect(meta_keep, names(dic)), drop = FALSE]
term_meta2 <- merge(term_meta, dic_meta, by = "term", all.x = TRUE, suffixes = c("_provenance", "_dictionary"))

term_outcome_rows <- list()
row_id <- 1L
for (i in seq_len(nrow(term_meta2))) {
  trm <- term_meta2$term[i]
  for (j in seq_along(OUTCOMES)) {
    out <- OUTCOMES[j]
    n_full <- count_pair(trm, out, mask_full)
    n_early <- count_pair(trm, out, mask_early)
    n_recent <- count_pair(trm, out, mask_recent)
    support_class <- if (n_early > 0 && n_recent > 0) {
      "both_eras"
    } else if (n_early > 0 && n_recent == 0) {
      "earlier_only"
    } else if (n_early == 0 && n_recent > 0) {
      "recent_only"
    } else {
      "no_term_outcome_support"
    }
    term_outcome_rows[[row_id]] <- data.frame(
      term = trm,
      outcome = out,
      outcome_label = OUTCOME_LABELS[j],
      full_support_n = n_full,
      earlier_support_n = n_early,
      recent_support_n = n_recent,
      support_class = support_class,
      stringsAsFactors = FALSE
    )
    row_id <- row_id + 1L
  }
}
term_outcome <- do.call(rbind, term_outcome_rows)
term_outcome <- merge(term_outcome, dic_meta, by = "term", all.x = TRUE)

term_summary <- aggregate(
  cbind(full_support_n, earlier_support_n, recent_support_n) ~ term,
  data = term_outcome,
  FUN = sum
)
term_summary$earlier_any <- term_summary$earlier_support_n > 0
term_summary$recent_any <- term_summary$recent_support_n > 0
term_summary$temporal_support_class <- ifelse(
  term_summary$earlier_any & term_summary$recent_any, "both_eras",
  ifelse(term_summary$earlier_any & !term_summary$recent_any, "earlier_only",
         ifelse(!term_summary$earlier_any & term_summary$recent_any, "recent_only", "no_term_outcome_support"))
)

outcomes_supported_full <- aggregate(full_support_n > 0 ~ term, data = term_outcome, FUN = sum)
names(outcomes_supported_full)[2] <- "outcomes_supported_full"
outcomes_supported_early <- aggregate(earlier_support_n > 0 ~ term, data = term_outcome, FUN = sum)
names(outcomes_supported_early)[2] <- "outcomes_supported_earlier"
outcomes_supported_recent <- aggregate(recent_support_n > 0 ~ term, data = term_outcome, FUN = sum)
names(outcomes_supported_recent)[2] <- "outcomes_supported_recent"
term_summary <- Reduce(function(x, y) merge(x, y, by = "term", all.x = TRUE), list(term_summary, outcomes_supported_full, outcomes_supported_early, outcomes_supported_recent))
term_summary <- merge(term_summary, dic_meta, by = "term", all.x = TRUE)

support_vectors <- by(term_outcome, term_outcome$term, function(df) {
  paste0(df$outcome_label, ": full=", df$full_support_n, ", earlier=", df$earlier_support_n, ", recent=", df$recent_support_n, collapse = "; ")
})
term_summary$temporal_support_vector <- as.character(support_vectors[term_summary$term])

# ---------------------------
# Edge temporal support
# ---------------------------
edge_temporal <- data.frame()
if (nrow(edge_prov) > 0) {
  from_col <- pick_col(edge_prov, c("from", "from_term", "source", "source_term", "node_from", "term_from", "A", "edge_from"))
  to_col <- pick_col(edge_prov, c("to", "to_term", "target", "target_term", "node_to", "term_to", "B", "edge_to"))
  if (!is.na(from_col) && !is.na(to_col)) {
    edge_temporal <- edge_prov
    edge_temporal$from_raw <- edge_temporal[[from_col]]
    edge_temporal$to_raw <- edge_temporal[[to_col]]
    edge_temporal$from_term <- vapply(edge_temporal$from_raw, resolve_term_to_dictionary_term, character(1))
    edge_temporal$to_term <- vapply(edge_temporal$to_raw, resolve_term_to_dictionary_term, character(1))
    edge_temporal$from_term[is.na(edge_temporal$from_term)] <- as.character(edge_temporal$from_raw[is.na(edge_temporal$from_term)])
    edge_temporal$to_term[is.na(edge_temporal$to_term)] <- as.character(edge_temporal$to_raw[is.na(edge_temporal$to_term)])
    edge_temporal$edge_support_full <- mapply(count_pair, edge_temporal$from_term, edge_temporal$to_term, MoreArgs = list(mask = mask_full))
    edge_temporal$edge_support_earlier <- mapply(count_pair, edge_temporal$from_term, edge_temporal$to_term, MoreArgs = list(mask = mask_early))
    edge_temporal$edge_support_recent <- mapply(count_pair, edge_temporal$from_term, edge_temporal$to_term, MoreArgs = list(mask = mask_recent))
    edge_temporal$temporal_support_class <- ifelse(
      edge_temporal$edge_support_earlier > 0 & edge_temporal$edge_support_recent > 0, "both_eras",
      ifelse(edge_temporal$edge_support_earlier > 0 & edge_temporal$edge_support_recent == 0, "earlier_only",
             ifelse(edge_temporal$edge_support_earlier == 0 & edge_temporal$edge_support_recent > 0, "recent_only", "no_edge_support"))
    )
  } else {
    warning("Edge provenance table found, but from/to columns could not be identified. Edge temporal support skipped.")
  }
}

# ---------------------------
# Triad temporal support
# ---------------------------
triad_temporal <- data.frame()
if (nrow(triad_prov) > 0) {
  A_col <- pick_col(triad_prov, c("A", "A_term", "term_A", "upstream", "from_A"))
  B_col <- pick_col(triad_prov, c("B", "B_term", "term_B", "bridge", "bridge_term"))
  C_col <- pick_col(triad_prov, c("C", "C_term", "term_C", "outcome", "outcome_term"))
  if (!is.na(A_col) && !is.na(B_col) && !is.na(C_col)) {
    triad_temporal <- triad_prov
    triad_temporal$A_raw <- triad_temporal[[A_col]]
    triad_temporal$B_raw <- triad_temporal[[B_col]]
    triad_temporal$C_raw <- triad_temporal[[C_col]]
    triad_temporal$A_term <- vapply(triad_temporal$A_raw, resolve_term_to_dictionary_term, character(1))
    triad_temporal$B_term <- vapply(triad_temporal$B_raw, resolve_term_to_dictionary_term, character(1))
    triad_temporal$C_term <- vapply(triad_temporal$C_raw, resolve_term_to_dictionary_term, character(1))
    triad_temporal$A_term[is.na(triad_temporal$A_term)] <- as.character(triad_temporal$A_raw[is.na(triad_temporal$A_term)])
    triad_temporal$B_term[is.na(triad_temporal$B_term)] <- as.character(triad_temporal$B_raw[is.na(triad_temporal$B_term)])
    triad_temporal$C_term[is.na(triad_temporal$C_term)] <- as.character(triad_temporal$C_raw[is.na(triad_temporal$C_term)])
    triad_temporal$AB_support_full <- mapply(count_pair, triad_temporal$A_term, triad_temporal$B_term, MoreArgs = list(mask = mask_full))
    triad_temporal$BC_support_full <- mapply(count_pair, triad_temporal$B_term, triad_temporal$C_term, MoreArgs = list(mask = mask_full))
    triad_temporal$AC_support_full <- mapply(count_pair, triad_temporal$A_term, triad_temporal$C_term, MoreArgs = list(mask = mask_full))
    triad_temporal$ABC_support_full <- mapply(count_triad, triad_temporal$A_term, triad_temporal$B_term, triad_temporal$C_term, MoreArgs = list(mask = mask_full))
    triad_temporal$AB_support_earlier <- mapply(count_pair, triad_temporal$A_term, triad_temporal$B_term, MoreArgs = list(mask = mask_early))
    triad_temporal$BC_support_earlier <- mapply(count_pair, triad_temporal$B_term, triad_temporal$C_term, MoreArgs = list(mask = mask_early))
    triad_temporal$AC_support_earlier <- mapply(count_pair, triad_temporal$A_term, triad_temporal$C_term, MoreArgs = list(mask = mask_early))
    triad_temporal$ABC_support_earlier <- mapply(count_triad, triad_temporal$A_term, triad_temporal$B_term, triad_temporal$C_term, MoreArgs = list(mask = mask_early))
    triad_temporal$AB_support_recent <- mapply(count_pair, triad_temporal$A_term, triad_temporal$B_term, MoreArgs = list(mask = mask_recent))
    triad_temporal$BC_support_recent <- mapply(count_pair, triad_temporal$B_term, triad_temporal$C_term, MoreArgs = list(mask = mask_recent))
    triad_temporal$AC_support_recent <- mapply(count_pair, triad_temporal$A_term, triad_temporal$C_term, MoreArgs = list(mask = mask_recent))
    triad_temporal$ABC_support_recent <- mapply(count_triad, triad_temporal$A_term, triad_temporal$B_term, triad_temporal$C_term, MoreArgs = list(mask = mask_recent))
    triad_temporal$AB_BC_temporal_class <- ifelse(
      triad_temporal$AB_support_earlier > 0 & triad_temporal$BC_support_earlier > 0 & triad_temporal$AB_support_recent > 0 & triad_temporal$BC_support_recent > 0,
      "AB_and_BC_in_both_eras",
      ifelse(triad_temporal$AB_support_recent > 0 & triad_temporal$BC_support_recent > 0, "AB_and_BC_recent",
             ifelse(triad_temporal$AB_support_earlier > 0 & triad_temporal$BC_support_earlier > 0, "AB_and_BC_earlier", "incomplete_temporal_support"))
    )
  } else {
    warning("Triad provenance table found, but A/B/C columns could not be identified. Triad temporal support skipped.")
  }
}

# ---------------------------
# Selection rules / interpretation table
# ---------------------------
selection_rules <- data.frame(
  item = c(
    "Purpose",
    "Interpretation",
    "Earlier era",
    "Recent era",
    "Term-outcome support",
    "Edge support",
    "Triad support",
    "No reselection",
    "Recommended reporting language"
  ),
  value = c(
    "Internal publication-era support check for final displayed figure terms, edges and ABC triads.",
    "Temporal support analysis; not prospective validation, external replication, causality or prediction.",
    paste0(EARLY_START, "-", EARLY_END),
    paste0(RECENT_START, "-", RECENT_END),
    "Unique article-level co-mentions of displayed concept node and outcome in the final hit matrix, recalculated by era.",
    "Unique article-level co-mentions of displayed edge endpoints in the final hit matrix, recalculated by era.",
    "A-B, B-C, A-C and A-B-C support recalculated by era for displayed ABC triads.",
    "The main atlas terms and figures are not reselected in this analysis; support is recalculated only for already displayed terms/edges/triads.",
    "Use 'temporal support' or 'publication-era robustness'; avoid calling this external validation or prospective holdout validation."
  ),
  stringsAsFactors = FALSE
)

# ---------------------------
# Write outputs
# ---------------------------
suffix <- paste0(PM_RANGE, "__", DIC_TAG)
write_csv_safe(term_outcome, paste0("Supplementary_Table_S9_TemporalSupport_TermOutcome_LONG_", suffix, ".csv"))
write_csv_safe(term_summary, paste0("Supplementary_Table_S9_TemporalSupport_Term_SUMMARY_", suffix, ".csv"))
if (nrow(edge_temporal) > 0) write_csv_safe(edge_temporal, paste0("Supplementary_Table_S9_TemporalSupport_Edge_Provenance_", suffix, ".csv"))
if (nrow(triad_temporal) > 0) write_csv_safe(triad_temporal, paste0("Supplementary_Table_S9_TemporalSupport_Triad_Provenance_", suffix, ".csv"))
write_csv_safe(selection_rules, paste0("Supplementary_Table_S9_TemporalSupport_SelectionRules_", suffix, ".csv"))

run_summary <- data.frame(
  metric = c(
    "n_articles_full", "n_articles_earlier", "n_articles_recent",
    "earlier_era", "recent_era",
    "n_figure_terms", "n_term_outcome_rows",
    "n_terms_both_eras", "n_terms_earlier_only", "n_terms_recent_only", "n_terms_no_term_outcome_support",
    "n_edges", "n_triads"
  ),
  value = c(
    nrow(hits), sum(mask_early), sum(mask_recent),
    paste0(EARLY_START, "-", EARLY_END), paste0(RECENT_START, "-", RECENT_END),
    length(unique(term_summary$term)), nrow(term_outcome),
    sum(term_summary$temporal_support_class == "both_eras", na.rm = TRUE),
    sum(term_summary$temporal_support_class == "earlier_only", na.rm = TRUE),
    sum(term_summary$temporal_support_class == "recent_only", na.rm = TRUE),
    sum(term_summary$temporal_support_class == "no_term_outcome_support", na.rm = TRUE),
    nrow(edge_temporal), nrow(triad_temporal)
  ),
  stringsAsFactors = FALSE
)
write_csv_safe(run_summary, paste0("TemporalSupport_run_summary_", suffix, ".csv"))

# ---------------------------
# Supplementary figures
# ---------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  plot_term <- term_summary
  plot_term$atlas_layer <- ifelse(is.na(plot_term$atlas_layer) | plot_term$atlas_layer == "", "unclassified", plot_term$atlas_layer)
  plot_term$temporal_support_class <- factor(
    plot_term$temporal_support_class,
    levels = c("both_eras", "earlier_only", "recent_only", "no_term_outcome_support")
  )
  p1 <- ggplot(plot_term, aes(x = atlas_layer, fill = temporal_support_class)) +
    geom_bar(width = 0.75) +
    coord_flip() +
    labs(
      title = "Temporal support categories for displayed concept nodes",
      subtitle = paste0("Earlier era: ", EARLY_START, "-", EARLY_END, "; recent era: ", RECENT_START, "-", RECENT_END),
      x = "Dictionary / atlas layer",
      y = "Displayed concept nodes",
      fill = "Temporal support"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")

  plot_pair <- term_outcome[term_outcome$full_support_n > 0, , drop = FALSE]
  plot_pair$atlas_layer <- ifelse(is.na(plot_pair$atlas_layer) | plot_pair$atlas_layer == "", "unclassified", plot_pair$atlas_layer)
  p2 <- ggplot(plot_pair, aes(x = earlier_support_n + 0.25, y = recent_support_n + 0.25)) +
    geom_point(aes(size = full_support_n, shape = outcome_label), alpha = 0.65) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Era-specific support for displayed term-outcome pairs",
      subtitle = "Axes show article-level co-mentions + 0.25 on log scale; points are term-outcome pairs shown in final atlas provenance.",
      x = paste0("Earlier-era support (", EARLY_START, "-", EARLY_END, ")"),
      y = paste0("Recent-era support (", RECENT_START, "-", RECENT_END, ")"),
      size = "Full support",
      shape = "Outcome"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")

  fig_pdf <- file.path(DIR_FIG, paste0("SupplementaryFigure_TemporalSupport_", suffix, ".pdf"))
  grDevices::pdf(fig_pdf, width = 10.5, height = 7.2, onefile = TRUE)
  print(p1)
  print(p2)
  grDevices::dev.off()

  fig_png1 <- file.path(DIR_FIG, paste0("SupplementaryFigure_TemporalSupport_A_term_categories_", suffix, ".png"))
  grDevices::png(fig_png1, width = 3000, height = 2100, res = 300)
  print(p1)
  grDevices::dev.off()

  fig_png2 <- file.path(DIR_FIG, paste0("SupplementaryFigure_TemporalSupport_B_term_outcome_pairs_", suffix, ".png"))
  grDevices::png(fig_png2, width = 3000, height = 2100, res = 300)
  print(p2)
  grDevices::dev.off()

  log_msg("WROTE:", fig_pdf)
  log_msg("WROTE:", fig_png1)
  log_msg("WROTE:", fig_png2)
} else {
  warning("ggplot2 is not installed; temporal support tables were written but supplementary figures were skipped.")
}

log_msg("WROTE temporal support tables to", DIR_TABLE)
log_msg("DONE 13_TemporalSupport_AtlasV2.R")
