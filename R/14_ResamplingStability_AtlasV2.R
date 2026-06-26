# 14_ResamplingStability_AtlasV2.R
# Purpose: Internal repeated-resampling stability analysis for displayed atlas terms, edges and ABC triads.
# This script does NOT create new main results and does NOT add new terms.
# It uses figure-level provenance outputs from step 12 as the source of displayed terms/edges/triads.
# Interpretation: internal sampling stability / robustness check, not external validation.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

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

# Load project setup when available
if (file.exists(file.path(CODE_DIR, "00_setup_AtlasV2.R"))) {
  source(file.path(CODE_DIR, "00_setup_AtlasV2.R"))
}

# Fallbacks if setup variables are absent
if (!exists("PM_RANGE")) PM_RANGE <- "pm_1980_20251231"
if (!exists("DIC_TAG")) DIC_TAG <- "analysis_v2_outcome_preserving_atlas"
if (!exists("DIR_OUT")) DIR_OUT <- file.path(ROOT, "output", PM_RANGE, DIC_TAG)
if (!exists("DIR_TAB")) DIR_TAB <- file.path(DIR_OUT, "table")
if (!exists("DIR_FIG")) DIR_FIG <- file.path(DIR_OUT, "fig")
if (!exists("DIR_LOG")) DIR_LOG <- file.path(DIR_OUT, "log")
if (!exists("DIR_PROC")) DIR_PROC <- file.path(ROOT, "data_proc")
if (!exists("DIC_FILE")) DIC_FILE <- file.path(ROOT, "dic", "ra_ild_dictionary_analysis_v2_outcome_preserving_atlas.csv")

dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat("[14]", paste(..., collapse = " "), "\n")
}

newest_file <- function(dir, pattern, required = TRUE) {
  x <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(x) == 0 && required) {
    stop("No file found in ", dir, " matching pattern: ", pattern, call. = FALSE)
  }
  if (length(x) == 0) return(NA_character_)
  x[order(file.info(x)$mtime, decreasing = TRUE)][1]
}

read_csv_safe <- function(path, guess_max = 100000) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE, guess_max = guess_max)
}

# ---------------------------
# Parameters
# ---------------------------
N_ITER <- as.integer(Sys.getenv("RAILD_RESAMPLE_N", unset = "200"))
SAMPLE_FRAC <- as.numeric(Sys.getenv("RAILD_RESAMPLE_FRAC", unset = "0.80"))
SEED <- as.integer(Sys.getenv("RAILD_RESAMPLE_SEED", unset = "51961"))
OUTCOMES <- c("AE-ILD", "progression", "mortality")
STABLE_CUT <- as.numeric(Sys.getenv("RAILD_RESAMPLE_STABLE_CUT", unset = "0.80"))
MODERATE_CUT <- as.numeric(Sys.getenv("RAILD_RESAMPLE_MODERATE_CUT", unset = "0.50"))

set.seed(SEED)

log_msg("Root:", ROOT)
log_msg("Table directory:", DIR_TAB)
log_msg("Figure directory:", DIR_FIG)
log_msg("Iterations:", N_ITER, "sample fraction:", SAMPLE_FRAC, "seed:", SEED)

# ---------------------------
# Input files
# ---------------------------
dict_file <- DIC_FILE
hits_file <- newest_file(DIR_TAB, paste0("^hits_matrix_.*__", DIC_TAG, "\\.csv$"))
term_summary_file <- newest_file(DIR_TAB, paste0("^Supplementary_Table_Figure_Term_Provenance_SUMMARY_.*__", DIC_TAG, "\\.csv$"))
term_long_file <- newest_file(DIR_TAB, paste0("^Supplementary_Table_Figure_Term_Provenance_LONG_.*__", DIC_TAG, "\\.csv$"), required = FALSE)
edge_file <- newest_file(DIR_TAB, paste0("^Supplementary_Table_Figure_Edge_Provenance_.*__", DIC_TAG, "\\.csv$"), required = FALSE)
triad_file <- newest_file(DIR_TAB, paste0("^Supplementary_Table_Figure_Triad_Provenance_.*__", DIC_TAG, "\\.csv$"), required = FALSE)

log_msg("Dictionary:", dict_file)
log_msg("Hits:", hits_file)
log_msg("Term summary:", term_summary_file)
if (!is.na(term_long_file)) log_msg("Term long:", term_long_file)
if (!is.na(edge_file)) log_msg("Edge provenance:", edge_file)
if (!is.na(triad_file)) log_msg("Triad provenance:", triad_file)

dic <- read_csv_safe(dict_file)
hits <- read_csv_safe(hits_file)
term_summary <- read_csv_safe(term_summary_file)
term_long <- if (!is.na(term_long_file)) read_csv_safe(term_long_file) else tibble()
edge_prov <- if (!is.na(edge_file)) read_csv_safe(edge_file) else tibble()
triad_prov <- if (!is.na(triad_file)) read_csv_safe(triad_file) else tibble()

# ---------------------------
# Helpers
# ---------------------------
std_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

first_existing_col <- function(df, candidates) {
  z <- intersect(candidates, names(df))
  if (length(z) == 0) return(NA_character_)
  z[1]
}

hit_col_candidates <- function(term) {
  term <- as.character(term)
  c(
    paste0("hit__", term),
    paste0("hit__", make.names(term)),
    paste0("hit__", gsub("[^A-Za-z0-9_]+", "_", term)),
    paste0("hit__", gsub("[^A-Za-z0-9]+", "_", term)),
    term,
    make.names(term),
    gsub("[^A-Za-z0-9_]+", "_", term)
  )
}

resolve_hit_col <- function(term, hit_names) {
  cand <- hit_col_candidates(term)
  z <- cand[cand %in% hit_names]
  if (length(z) > 0) return(z[1])
  NA_character_
}

as_hit_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), FALSE, x > 0))
  xx <- tolower(as.character(x))
  xx %in% c("1", "true", "t", "yes", "y")
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(std_chr(x))
  x <- x[nzchar(x)]
  if (length(x) == 0) "" else paste(x, collapse = sep)
}

# ---------------------------
# Displayed term set from provenance summary
# ---------------------------
term_col <- first_existing_col(term_summary, c("term", "dictionary_term", "concept_node", "node", "from", "to"))
if (is.na(term_col)) stop("Could not find a term column in term provenance summary.", call. = FALSE)

terms0 <- term_summary |>
  mutate(term = .data[[term_col]]) |>
  filter(!is.na(term), nzchar(as.character(term))) |>
  distinct(term, .keep_all = TRUE)

# Join dictionary metadata; keep provenance metadata if already present.
terms_meta <- terms0 |>
  left_join(dic, by = "term", suffix = c(".prov", ".dict"))

# Prefer dictionary metadata when available; otherwise keep provenance values.
prefer_col <- function(df, base) {
  prov <- paste0(base, ".prov")
  dict <- paste0(base, ".dict")
  if (dict %in% names(df) && prov %in% names(df)) {
    dplyr::coalesce(as.character(df[[dict]]), as.character(df[[prov]]))
  } else if (dict %in% names(df)) {
    as.character(df[[dict]])
  } else if (prov %in% names(df)) {
    as.character(df[[prov]])
  } else if (base %in% names(df)) {
    as.character(df[[base]])
  } else {
    rep(NA_character_, nrow(df))
  }
}

terms_meta <- terms_meta |>
  mutate(
    preferred_label = prefer_col(cur_data_all(), "preferred_label"),
    concept_id = prefer_col(cur_data_all(), "concept_id"),
    atlas_layer = prefer_col(cur_data_all(), "atlas_layer"),
    abc_role = prefer_col(cur_data_all(), "abc_role"),
    abc_subrole = prefer_col(cur_data_all(), "abc_subrole"),
    domain = prefer_col(cur_data_all(), "domain"),
    concept_family = prefer_col(cur_data_all(), "concept_family"),
    node_level = prefer_col(cur_data_all(), "node_level"),
    analysis_tier = prefer_col(cur_data_all(), "analysis_tier"),
    include_primary = prefer_col(cur_data_all(), "include_primary")
  ) |>
  select(term, preferred_label, concept_id, atlas_layer, abc_role, abc_subrole,
         domain, concept_family, node_level, analysis_tier, include_primary, everything())

terms_meta$hit_col <- vapply(terms_meta$term, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits))
terms_meta <- terms_meta |> mutate(has_hit_col = !is.na(hit_col))

if (sum(terms_meta$has_hit_col) == 0) stop("No displayed terms could be resolved to hit matrix columns.", call. = FALSE)
terms_for_resampling <- terms_meta |> filter(has_hit_col)

# Outcomes
outcome_cols <- vapply(OUTCOMES, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits))
if (any(is.na(outcome_cols))) {
  stop("Could not resolve outcome hit columns: ", paste(OUTCOMES[is.na(outcome_cols)], collapse = ", "), call. = FALSE)
}

# Use PMID/year when available
pmid_col <- first_existing_col(hits, c("pmid", "PMID", "uid"))
year_col <- first_existing_col(hits, c("year", "pub_year", "publication_year"))
if (is.na(pmid_col)) hits$pmid_tmp__ <- seq_len(nrow(hits)) else hits$pmid_tmp__ <- hits[[pmid_col]]
if (is.na(year_col)) hits$year_tmp__ <- NA_integer_ else hits$year_tmp__ <- suppressWarnings(as.integer(hits[[year_col]]))

# Matrices
term_mat <- as.data.frame(lapply(terms_for_resampling$hit_col, function(cc) as_hit_bool(hits[[cc]])))
names(term_mat) <- terms_for_resampling$term
term_mat <- as.matrix(term_mat)

out_mat <- as.data.frame(lapply(outcome_cols, function(cc) as_hit_bool(hits[[cc]])))
names(out_mat) <- OUTCOMES
out_mat <- as.matrix(out_mat)

n_articles <- nrow(hits)
sample_n <- max(1L, min(n_articles, round(SAMPLE_FRAC * n_articles)))
n_terms <- ncol(term_mat)
n_out <- ncol(out_mat)

log_msg("Displayed terms resolved:", n_terms, "of", nrow(terms_meta))
log_msg("Articles:", n_articles, "sample_n:", sample_n)

# ---------------------------
# Full support
# ---------------------------
full_counts <- t(term_mat) %*% out_mat
full_counts <- matrix(as.integer(full_counts), nrow = n_terms, ncol = n_out,
                      dimnames = list(colnames(term_mat), colnames(out_mat)))

# ---------------------------
# Repeated subsampling
# ---------------------------
counts_arr <- array(0L, dim = c(n_terms, n_out, N_ITER),
                    dimnames = list(colnames(term_mat), colnames(out_mat), paste0("iter", seq_len(N_ITER))))

for (ii in seq_len(N_ITER)) {
  idx <- sample.int(n_articles, size = sample_n, replace = FALSE)
  cc <- t(term_mat[idx, , drop = FALSE]) %*% out_mat[idx, , drop = FALSE]
  counts_arr[, , ii] <- as.integer(cc)
  if (ii %% 25 == 0) log_msg("Resample", ii, "of", N_ITER)
}

# Term-outcome summary
term_outcome_long <- expand_grid(
  term = colnames(term_mat),
  outcome = colnames(out_mat)
) |>
  mutate(
    full_support = as.integer(full_counts[cbind(term, outcome)]),
    mean_support = purrr::map2_dbl(term, outcome, ~ mean(counts_arr[.x, .y, ])),
    median_support = purrr::map2_dbl(term, outcome, ~ median(counts_arr[.x, .y, ])),
    q25_support = purrr::map2_dbl(term, outcome, ~ as.numeric(quantile(counts_arr[.x, .y, ], 0.25))),
    q75_support = purrr::map2_dbl(term, outcome, ~ as.numeric(quantile(counts_arr[.x, .y, ], 0.75))),
    p_support_gt0 = purrr::map2_dbl(term, outcome, ~ mean(counts_arr[.x, .y, ] > 0)),
    p_support_ge2 = purrr::map2_dbl(term, outcome, ~ mean(counts_arr[.x, .y, ] >= 2)),
    p_support_ge3 = purrr::map2_dbl(term, outcome, ~ mean(counts_arr[.x, .y, ] >= 3))
  ) |>
  left_join(terms_for_resampling |> select(term, preferred_label, concept_id, atlas_layer, abc_role, abc_subrole,
                                           domain, concept_family, node_level, analysis_tier, include_primary), by = "term") |>
  relocate(term, preferred_label, outcome)

# Term-level summary
tterm <- lapply(seq_len(n_terms), function(i) {
  mat_i <- counts_arr[i, , , drop = TRUE] # outcomes x iter
  if (is.vector(mat_i)) mat_i <- matrix(mat_i, nrow = n_out)
  support_outcomes_each_iter <- colSums(mat_i > 0)
  tibble(
    term = colnames(term_mat)[i],
    full_total_support = sum(full_counts[i, ]),
    full_outcomes_supported = sum(full_counts[i, ] > 0),
    full_support_vector = paste(paste0(colnames(out_mat), "=", as.integer(full_counts[i, ])), collapse = "; "),
    mean_total_support = mean(colSums(mat_i)),
    median_total_support = median(colSums(mat_i)),
    p_any_outcome_support = mean(support_outcomes_each_iter >= 1),
    p_two_or_more_outcomes = mean(support_outcomes_each_iter >= 2),
    p_all_three_outcomes = mean(support_outcomes_each_iter == 3),
    p_no_outcome_support = mean(support_outcomes_each_iter == 0)
  )
}) |> bind_rows() |>
  left_join(terms_for_resampling |> select(term, preferred_label, concept_id, atlas_layer, abc_role, abc_subrole,
                                           domain, concept_family, node_level, analysis_tier, include_primary), by = "term") |>
  mutate(
    resampling_stability_category = case_when(
      p_any_outcome_support >= STABLE_CUT ~ "stable_any_support",
      p_any_outcome_support >= MODERATE_CUT ~ "moderate_any_support",
      p_any_outcome_support > 0 ~ "intermittent_any_support",
      TRUE ~ "no_support_in_resamples"
    )
  ) |>
  relocate(term, preferred_label)

# ---------------------------
# Edge stability
# ---------------------------
extract_edge_cols <- function(edge_df) {
  if (nrow(edge_df) == 0) return(tibble())
  from_col <- first_existing_col(edge_df, c("from", "source", "node1", "term1", "A", "term_a"))
  to_col <- first_existing_col(edge_df, c("to", "target", "node2", "term2", "B", "term_b", "C"))
  type_col <- first_existing_col(edge_df, c("edge_type", "relationship", "type", "edge_class"))
  fig_col <- first_existing_col(edge_df, c("figure", "fig"))
  panel_col <- first_existing_col(edge_df, c("panel"))
  if (is.na(from_col) || is.na(to_col)) return(tibble())
  edge_df |>
    transmute(
      figure = if (!is.na(fig_col)) as.character(.data[[fig_col]]) else "",
      panel = if (!is.na(panel_col)) as.character(.data[[panel_col]]) else "",
      edge_type = if (!is.na(type_col)) as.character(.data[[type_col]]) else "edge",
      from = as.character(.data[[from_col]]),
      to = as.character(.data[[to_col]])
    ) |>
    filter(!is.na(from), !is.na(to), nzchar(from), nzchar(to)) |>
    distinct()
}

edge_list <- extract_edge_cols(edge_prov)
if (nrow(edge_list) > 0) {
  edge_list <- edge_list |>
    mutate(
      from_col = vapply(from, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits)),
      to_col = vapply(to, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits)),
      has_cols = !is.na(from_col) & !is.na(to_col)
    )
  edge_stab <- lapply(seq_len(nrow(edge_list)), function(i) {
    ee <- edge_list[i, ]
    if (!ee$has_cols) {
      return(tibble(full_edge_support = NA_integer_, mean_edge_support = NA_real_, median_edge_support = NA_real_,
                    q25_edge_support = NA_real_, q75_edge_support = NA_real_, p_edge_support_gt0 = NA_real_,
                    p_edge_support_ge2 = NA_real_, p_edge_support_ge3 = NA_real_))
    }
    v1 <- as_hit_bool(hits[[ee$from_col]])
    v2 <- as_hit_bool(hits[[ee$to_col]])
    full <- sum(v1 & v2)
    vals <- replicate(N_ITER, {
      idx <- sample.int(n_articles, size = sample_n, replace = FALSE)
      sum(v1[idx] & v2[idx])
    })
    tibble(
      full_edge_support = as.integer(full),
      mean_edge_support = mean(vals),
      median_edge_support = median(vals),
      q25_edge_support = as.numeric(quantile(vals, 0.25)),
      q75_edge_support = as.numeric(quantile(vals, 0.75)),
      p_edge_support_gt0 = mean(vals > 0),
      p_edge_support_ge2 = mean(vals >= 2),
      p_edge_support_ge3 = mean(vals >= 3)
    )
  }) |> bind_rows()
  edge_stability <- bind_cols(edge_list, edge_stab)
} else {
  edge_stability <- tibble()
}

# ---------------------------
# ABC triad stability
# ---------------------------
extract_triad_cols <- function(triad_df) {
  if (nrow(triad_df) == 0) return(tibble())
  A_col <- first_existing_col(triad_df, c("A", "term_A", "a_term"))
  B_col <- first_existing_col(triad_df, c("B", "term_B", "b_term"))
  C_col <- first_existing_col(triad_df, c("C", "outcome", "term_C", "c_term"))
  fig_col <- first_existing_col(triad_df, c("figure", "fig"))
  panel_col <- first_existing_col(triad_df, c("panel"))
  if (is.na(A_col) || is.na(B_col) || is.na(C_col)) return(tibble())
  triad_df |>
    transmute(
      figure = if (!is.na(fig_col)) as.character(.data[[fig_col]]) else "",
      panel = if (!is.na(panel_col)) as.character(.data[[panel_col]]) else "",
      A = as.character(.data[[A_col]]),
      B = as.character(.data[[B_col]]),
      C = as.character(.data[[C_col]])
    ) |>
    filter(!is.na(A), !is.na(B), !is.na(C), nzchar(A), nzchar(B), nzchar(C)) |>
    distinct()
}

triad_list <- extract_triad_cols(triad_prov)
if (nrow(triad_list) > 0) {
  triad_list <- triad_list |>
    mutate(
      A_col = vapply(A, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits)),
      B_col = vapply(B, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits)),
      C_col = vapply(C, resolve_hit_col, FUN.VALUE = character(1), hit_names = names(hits)),
      has_cols = !is.na(A_col) & !is.na(B_col) & !is.na(C_col)
    )
  triad_stab <- lapply(seq_len(nrow(triad_list)), function(i) {
    tt <- triad_list[i, ]
    if (!tt$has_cols) {
      return(tibble(full_AB_support = NA_integer_, full_BC_support = NA_integer_,
                    mean_AB_support = NA_real_, mean_BC_support = NA_real_,
                    median_AB_support = NA_real_, median_BC_support = NA_real_,
                    p_AB_support_gt0 = NA_real_, p_BC_support_gt0 = NA_real_,
                    p_AB_and_BC_support_gt0 = NA_real_, p_AB_and_BC_support_ge2 = NA_real_, p_AB_and_BC_support_ge3 = NA_real_))
    }
    va <- as_hit_bool(hits[[tt$A_col]])
    vb <- as_hit_bool(hits[[tt$B_col]])
    vc <- as_hit_bool(hits[[tt$C_col]])
    full_AB <- sum(va & vb)
    full_BC <- sum(vb & vc)
    vals <- replicate(N_ITER, {
      idx <- sample.int(n_articles, size = sample_n, replace = FALSE)
      c(AB = sum(va[idx] & vb[idx]), BC = sum(vb[idx] & vc[idx]))
    })
    AB <- vals["AB", ]; BC <- vals["BC", ]
    tibble(
      full_AB_support = as.integer(full_AB),
      full_BC_support = as.integer(full_BC),
      mean_AB_support = mean(AB),
      mean_BC_support = mean(BC),
      median_AB_support = median(AB),
      median_BC_support = median(BC),
      p_AB_support_gt0 = mean(AB > 0),
      p_BC_support_gt0 = mean(BC > 0),
      p_AB_and_BC_support_gt0 = mean(AB > 0 & BC > 0),
      p_AB_and_BC_support_ge2 = mean(AB >= 2 & BC >= 2),
      p_AB_and_BC_support_ge3 = mean(AB >= 3 & BC >= 3)
    )
  }) |> bind_rows()
  triad_stability <- bind_cols(triad_list, triad_stab)
} else {
  triad_stability <- tibble()
}

# ---------------------------
# Selection rules and summary tables
# ---------------------------
selection_rules <- tibble::tribble(
  ~component, ~rule,
  "Purpose", "Internal repeated-subsampling stability check for displayed atlas concept nodes, edges and ABC triads; not external validation.",
  "Input terms", "Displayed concept nodes were obtained from Supplementary_Table_Figure_Term_Provenance_SUMMARY and mapped back to the final dictionary and hit matrix.",
  "Sampling", paste0(N_ITER, " repeated subsamples without replacement; sample fraction = ", SAMPLE_FRAC, "; fixed seed = ", SEED, "."),
  "Term-outcome support", "For each displayed term and each primary outcome, support was recalculated as article-level co-mention in each subsample.",
  "Term stability", "Term-level stability summarises whether a displayed concept node retained any, two-or-more, or all-three outcome support across subsamples.",
  "Edge stability", "Displayed network edges from figure-level edge provenance were recalculated as article-level co-mentions in each subsample.",
  "ABC triad stability", "Displayed ABC triads were assessed by recalculating A-B and B-C supports in each subsample; structural support required both A-B and B-C support.",
  "Interpretation", "Results are internal sampling-stability metrics and should not be interpreted as external validation, prospective replication, causal inference or prediction-model validation."
)

run_summary <- tibble::tibble(
  metric = c(
    "n_articles", "sample_fraction", "sample_n", "n_iterations", "seed",
    "displayed_terms_in_summary", "displayed_terms_resolved_to_hits",
    "term_outcome_pairs", "edges_assessed", "triads_assessed",
    "stable_terms_any_support_ge_cut", "moderate_terms_any_support_ge_cut",
    "edges_with_p_support_gt0_ge_0_8", "triads_with_p_AB_and_BC_gt0_ge_0_8"
  ),
  value = c(
    n_articles, SAMPLE_FRAC, sample_n, N_ITER, SEED,
    nrow(terms_meta), n_terms,
    nrow(term_outcome_long), nrow(edge_stability), nrow(triad_stability),
    sum(tterm$p_any_outcome_support >= STABLE_CUT, na.rm = TRUE),
    sum(tterm$p_any_outcome_support >= MODERATE_CUT, na.rm = TRUE),
    if (nrow(edge_stability) > 0) sum(edge_stability$p_edge_support_gt0 >= STABLE_CUT, na.rm = TRUE) else 0,
    if (nrow(triad_stability) > 0) sum(triad_stability$p_AB_and_BC_support_gt0 >= STABLE_CUT, na.rm = TRUE) else 0
  )
)

# ---------------------------
# Write tables
# ---------------------------
prefix <- paste0("pm_1980_20251231__", DIC_TAG)
write_csv(term_outcome_long, file.path(DIR_TAB, paste0("Supplementary_Table_S10_ResamplingStability_TermOutcome_LONG_", prefix, ".csv")))
write_csv(tterm, file.path(DIR_TAB, paste0("Supplementary_Table_S10_ResamplingStability_Term_SUMMARY_", prefix, ".csv")))
write_csv(edge_stability, file.path(DIR_TAB, paste0("Supplementary_Table_S10_ResamplingStability_Edge_Stability_", prefix, ".csv")))
write_csv(triad_stability, file.path(DIR_TAB, paste0("Supplementary_Table_S10_ResamplingStability_Triad_Stability_", prefix, ".csv")))
write_csv(selection_rules, file.path(DIR_TAB, paste0("Supplementary_Table_S10_ResamplingStability_SelectionRules_", prefix, ".csv")))
write_csv(run_summary, file.path(DIR_TAB, paste0("ResamplingStability_run_summary_", prefix, ".csv")))

# ---------------------------
# Figures
# ---------------------------
layer_order <- c("primary_atlas", "treatment_context", "infection_context", "sensitivity", "exploratory", "context_only", "")
stability_order <- c("stable_any_support", "moderate_any_support", "intermittent_any_support", "no_support_in_resamples")

plot_term_layer <- tterm |>
  mutate(
    atlas_layer = if_else(is.na(atlas_layer) | atlas_layer == "", "unknown", atlas_layer),
    atlas_layer = factor(atlas_layer, levels = unique(c(layer_order, sort(unique(atlas_layer))))),
    resampling_stability_category = factor(resampling_stability_category, levels = stability_order)
  ) |>
  count(atlas_layer, resampling_stability_category, name = "n")

pA <- ggplot(plot_term_layer, aes(x = n, y = atlas_layer, fill = resampling_stability_category)) +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = c(
      stable_any_support = "#E15759",
      moderate_any_support = "#F28E2B",
      intermittent_any_support = "#59A14F",
      no_support_in_resamples = "#B07BEF"
    ),
    labels = c(
      stable_any_support = paste0("Stable (p≥", STABLE_CUT, ")"),
      moderate_any_support = paste0("Moderate (p≥", MODERATE_CUT, ")"),
      intermittent_any_support = "Intermittent",
      no_support_in_resamples = "No support"
    ),
    drop = FALSE
  ) +
  labs(
    title = "A. Resampling stability of displayed concept nodes",
    subtitle = paste0(N_ITER, " fixed-seed subsamples without replacement; sample fraction = ", SAMPLE_FRAC),
    x = "Displayed concept nodes",
    y = "Dictionary / atlas layer",
    fill = "Sampling stability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "right"
  )

# Edge/triad panel
edge_plot <- edge_stability |>
  transmute(
    object_type = "Network edge",
    label = paste(from, to, sep = " → "),
    full_support = full_edge_support,
    stability = p_edge_support_gt0
  )
triad_plot <- triad_stability |>
  transmute(
    object_type = "ABC triad",
    label = paste(A, B, C, sep = " → "),
    full_support = pmin(full_AB_support, full_BC_support, na.rm = TRUE),
    stability = p_AB_and_BC_support_gt0
  )
combined_obj <- bind_rows(edge_plot, triad_plot) |>
  filter(!is.na(full_support), !is.na(stability)) |>
  mutate(
    object_type = factor(object_type, levels = c("Network edge", "ABC triad")),
    full_support_plot = full_support + 0.25
  )

pB <- ggplot(combined_obj, aes(x = full_support_plot, y = stability, colour = object_type, shape = object_type)) +
  geom_hline(yintercept = STABLE_CUT, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_point(alpha = 0.72, size = 2.35) +
  scale_x_log10(labels = scales::label_number()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.02)) +
  scale_colour_manual(values = c("Network edge" = "#4E79A7", "ABC triad" = "#F28E2B")) +
  labs(
    title = "B. Resampling stability of displayed edges and ABC triads",
    subtitle = "X-axis shows full-corpus support + 0.25 on log scale; y-axis shows support frequency across subsamples.",
    x = "Full-corpus support",
    y = "Support frequency",
    colour = "Object type",
    shape = "Object type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "right"
  )

fig_pdf <- file.path(DIR_FIG, paste0("SupplementaryFigure_ResamplingStability_AB_combined_", prefix, ".pdf"))
fig_png <- file.path(DIR_FIG, paste0("SupplementaryFigure_ResamplingStability_AB_combined_", prefix, ".png"))

# Save combined figure
save_combined <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") grDevices::pdf(path, width = 9.2, height = 11.0, onefile = TRUE)
  if (device == "png") grDevices::png(path, width = 2760, height = 3300, res = 300)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(0.95, 1.05), "null"))))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.draw(ggplot2::ggplotGrob(pA))
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::grid.draw(ggplot2::ggplotGrob(pB))
  grid::popViewport(2)
  grDevices::dev.off()
}

save_combined(fig_pdf, "pdf")
save_combined(fig_png, "png")

# Also save individual panels for convenience
ggplot2::ggsave(file.path(DIR_FIG, paste0("SupplementaryFigure_ResamplingStability_A_term_categories_", prefix, ".pdf")), pA, width = 8.8, height = 4.8)
ggplot2::ggsave(file.path(DIR_FIG, paste0("SupplementaryFigure_ResamplingStability_B_edges_triads_", prefix, ".pdf")), pB, width = 8.8, height = 5.4)

log_msg("WROTE:", fig_pdf)
log_msg("WROTE:", fig_png)
log_msg("DONE")
