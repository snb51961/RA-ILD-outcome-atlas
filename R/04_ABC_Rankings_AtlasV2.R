# =========================================================
# 04_ABC_Rankings_AtlasV2.R
# RA-ILD atlas v2.2 core analysis
#
# Purpose
#   - Rank A -> B -> C bridge structures using the primary disease-state atlas.
#   - Keep treatment-context and infection-context terms outside primary ABC
#     bridge claims.
#   - Write full candidate rankings plus main-display-filtered rankings.
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)
quiet_install(c("dplyr", "readr", "stringr", "tibble", "purrr"))

RUN_NAME <- "04_ABC_Rankings_AtlasV2"
log_msg("=== START ", RUN_NAME, " ===")

# -------------------------
# 0) Parameters
# -------------------------
MAX_PMIDS_PER_PAIR <- 15L
TOP_EVIDENCE_TRIADS <- 50L
TOP_AE_NETWORK_TRIADS <- 30L
NOVELTY_LAMBDA <- 0.75

# These filters are not used to erase evidence.  The full candidate table is written.
# They define the more conservative main-display table.
MIN_AB_N11 <- MIN_AB_N11_MAIN
MIN_BC_N11 <- MIN_BC_N11_MAIN
LIFT_MIN_AB <- LIFT_MIN_MAIN
LIFT_MIN_BC <- LIFT_MIN_MAIN
NPMI_MIN_AB <- NPMI_MIN_MAIN
NPMI_MIN_BC <- NPMI_MIN_MAIN

# -------------------------
# 1) Load inputs
# -------------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)
df <- load_hit_matrix()
terms_present <- hit_terms_present(df, dic)

A_set <- intersect(select_primary_A_terms(dic), terms_present)
B_set <- intersect(select_primary_B_terms(dic), terms_present)
C_set <- intersect(select_primary_C_terms(dic, primary_outcomes_only = TRUE), terms_present)
C_set <- unique(c(intersect(PRIMARY_OUTCOME_TERMS, C_set), setdiff(C_set, PRIMARY_OUTCOME_TERMS)))

if (!length(A_set)) stop("No primary A terms available for ABC analysis.")
if (!length(B_set)) stop("No primary B terms available for ABC analysis.")
if (!length(C_set)) stop("No primary C outcome terms available for ABC analysis.")

# Defensive: ensure context layers cannot enter primary ABC.
A_set <- setdiff(A_set, c(treatment_context_terms(dic), infection_context_terms(dic)))
B_set <- setdiff(B_set, c(treatment_context_terms(dic), infection_context_terms(dic)))

log_msg(sprintf("04 primary ABC sets |A|=%d |B|=%d |C|=%d [%s]", length(A_set), length(B_set), length(C_set), paste(C_set, collapse = ", ")))

# -------------------------
# 2) Pre-compute B IDF and helpers
# -------------------------
N_docs <- nrow(df)
df_B <- vapply(B_set, function(b) sum(get_hit_vec(df, b) == 1L, na.rm = TRUE), integer(1))
idf_B <- log((N_docs + 1) / (df_B + 1)) + 1
idfB_of <- function(b) as.numeric(idf_B[[b]])

safelog <- function(z) {
  x <- suppressWarnings(log(z))
  x[!is.finite(x)] <- 0
  x
}

# -------------------------
# 3) Full ABC candidate ranking
# -------------------------
res <- list()
for (c in C_set) {
  for (a in A_set) {
    for (b in B_set) {
      AB <- pair_stats(df, a, b)
      BC <- pair_stats(df, b, c)
      if (AB$n11 <= 0 || BC$n11 <= 0) next
      AC <- pair_stats(df, a, c)
      res[[length(res) + 1L]] <- tibble::tibble(
        A = a, B = b, C = c,
        AB_or = AB$or, AB_p = AB$p, AB_n11 = AB$n11, AB_lift = AB$lift, AB_npmi = AB$npmi,
        BC_or = BC$or, BC_p = BC$p, BC_n11 = BC$n11, BC_lift = BC$lift, BC_npmi = BC$npmi,
        AC_or = AC$or, AC_p = AC$p, AC_n11 = AC$n11, AC_lift = AC$lift, AC_npmi = AC$npmi
      )
    }
  }
}

tab_all <- if (!length(res)) {
  tibble::tibble(
    A = character(), B = character(), C = character(),
    AB_or = double(), AB_p = double(), AB_n11 = integer(), AB_lift = double(), AB_npmi = double(),
    BC_or = double(), BC_p = double(), BC_n11 = integer(), BC_lift = double(), BC_npmi = double(),
    AC_or = double(), AC_p = double(), AC_n11 = integer(), AC_lift = double(), AC_npmi = double()
  )
} else {
  dplyr::bind_rows(res)
}

if (!nrow(tab_all)) {
  tab_all$AB_q <- double()
  tab_all$BC_q <- double()
  tab_all$AC_q <- double()
  tab_all$AB_log <- double()
  tab_all$BC_log <- double()
  tab_all$AC_log <- double()
  tab_all$AB_qsig <- double()
  tab_all$BC_qsig <- double()
  tab_all$AC_qsig <- double()
  tab_all$idfB <- double()
  tab_all$support_hmean <- double()
  tab_all$AC_penalty <- double()
  tab_all$score_q <- double()
  tab_all$main_display_eligible <- logical()
}

if (nrow(tab_all)) {
  tab_all <- tab_all |>
    dplyr::mutate(
      AB_q = p.adjust(AB_p, method = "BH"),
      BC_q = p.adjust(BC_p, method = "BH"),
      AC_q = p.adjust(AC_p, method = "BH"),
      AB_log = safelog(AB_or),
      BC_log = safelog(BC_or),
      AC_log = safelog(AC_or),
      AB_qsig = -log10(pmax(AB_q, 1e-300)),
      BC_qsig = -log10(pmax(BC_q, 1e-300)),
      AC_qsig = -log10(pmax(AC_q, 1e-300)),
      idfB = vapply(B, idfB_of, numeric(1)),
      support_hmean = 2 * AB_n11 * BC_n11 / pmax(1, AB_n11 + BC_n11),
      AC_penalty = pmax(AC_log, 0) * as.integer(AC_qsig >= 3),
      score_q = ((pmax(AB_log, 0) * AB_qsig * idfB) + (pmax(BC_log, 0) * BC_qsig * idfB)) * sqrt(pmax(1, support_hmean)) - (NOVELTY_LAMBDA * AC_penalty),
      main_display_eligible = AB_n11 >= MIN_AB_N11 & BC_n11 >= MIN_BC_N11 &
        AB_lift >= LIFT_MIN_AB & BC_lift >= LIFT_MIN_BC &
        AB_npmi >= NPMI_MIN_AB & BC_npmi >= NPMI_MIN_BC
    ) |>
    add_term_metadata(dic, term_col = "A", prefix = "A") |>
    add_term_metadata(dic, term_col = "B", prefix = "B") |>
    add_term_metadata(dic, term_col = "C", prefix = "C") |>
    dplyr::arrange(dplyr::desc(main_display_eligible), dplyr::desc(score_q), dplyr::desc(support_hmean), A, B, C)
}

tab_main <- tab_all |> dplyr::filter(main_display_eligible == TRUE) |> dplyr::arrange(dplyr::desc(score_q))
tab_ae   <- tab_main |> dplyr::filter(C == "AE-ILD") |> dplyr::arrange(dplyr::desc(score_q))

# -------------------------
# 4) Treatment/infection context pair diagnostics, not primary ABC
# -------------------------
context_terms <- c(treatment_context_terms(dic), infection_context_terms(dic))
context_terms <- intersect(context_terms, terms_present)
context_pair_stats <- if (length(context_terms)) {
  purrr::map_dfr(context_terms, function(t) {
    purrr::map_dfr(C_set, function(c) {
      ps <- pair_stats(df, t, c)
      tibble::tibble(context_term = t, outcome = c) |> dplyr::bind_cols(ps)
    })
  }) |>
    dplyr::mutate(q = p.adjust(p, method = "BH")) |>
    add_term_metadata(dic, term_col = "context_term", prefix = "context_term") |>
    add_term_metadata(dic, term_col = "outcome", prefix = "outcome") |>
    dplyr::arrange(context_term_atlas_layer, outcome, dplyr::desc(n11), dplyr::desc(npmi))
} else {
  tibble::tibble(context_term = character(), outcome = character())
}

# -------------------------
# 5) Evidence PMIDs for top triads
# -------------------------
pmids_for <- function(term) df$pmid[get_hit_vec(df, term) == 1L]

join_titles <- function(ev) {
  if (!nrow(ev)) return(ev)
  have <- intersect(c("pmid", "year", "journal", "title"), names(df))
  if (!"pmid" %in% have) return(ev)
  dplyr::left_join(ev, df |> dplyr::select(dplyr::any_of(c("pmid", "year", "journal", "title"))), by = "pmid") |>
    dplyr::arrange(A, B, C, dplyr::desc(year))
}

make_evidence <- function(top_triads) {
  if (!nrow(top_triads)) {
    return(list(AB = tibble::tibble(), BC = tibble::tibble()))
  }
  top_triads <- top_triads |> dplyr::slice_head(n = min(TOP_EVIDENCE_TRIADS, nrow(top_triads)))
  ev_AB <- purrr::map_dfr(seq_len(nrow(top_triads)), function(i) {
    a <- top_triads$A[i]; b <- top_triads$B[i]; c <- top_triads$C[i]
    tibble::tibble(A = a, B = b, C = c, evidence_edge = "AB", pmid = head(intersect(pmids_for(a), pmids_for(b)), MAX_PMIDS_PER_PAIR))
  })
  ev_BC <- purrr::map_dfr(seq_len(nrow(top_triads)), function(i) {
    a <- top_triads$A[i]; b <- top_triads$B[i]; c <- top_triads$C[i]
    tibble::tibble(A = a, B = b, C = c, evidence_edge = "BC", pmid = head(intersect(pmids_for(b), pmids_for(c)), MAX_PMIDS_PER_PAIR))
  })
  list(AB = join_titles(ev_AB), BC = join_titles(ev_BC))
}

ev_all <- make_evidence(tab_main)
ev_ae  <- make_evidence(tab_ae)

# -------------------------
# 6) Network edges for AE-ILD-focused downstream figures
# -------------------------
top_ae_net <- tab_ae |> dplyr::slice_head(n = min(TOP_AE_NETWORK_TRIADS, nrow(tab_ae)))
edge_AB <- if (nrow(top_ae_net)) {
  top_ae_net |> dplyr::transmute(from = A, to = B, w = support_hmean, kind = "AB", C = C, score_q = score_q)
} else tibble::tibble(from = character(), to = character(), w = double(), kind = character(), C = character(), score_q = double())
edge_BC <- if (nrow(top_ae_net)) {
  top_ae_net |> dplyr::transmute(from = B, to = C, w = support_hmean, kind = "BC", C = C, score_q = score_q)
} else tibble::tibble(from = character(), to = character(), w = double(), kind = character(), C = character(), score_q = double())
edges_ae <- dplyr::bind_rows(edge_AB, edge_BC) |> dplyr::distinct()

# -------------------------
# 7) Save outputs
# -------------------------
f_all  <- file.path(DIR_TABLE, sprintf("abc_rankings_primary_all_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_main <- file.path(DIR_TABLE, sprintf("abc_rankings_primary_main_filtered_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_ae   <- file.path(DIR_TABLE, sprintf("abc_rankings_primary_AEILD_main_filtered_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_ctx  <- file.path(DIR_TABLE, sprintf("abc_context_pair_stats_treatment_infection_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_evAB <- file.path(DIR_TABLE, sprintf("abc_evidence_AB_primary_top%d_%s__%s.csv", TOP_EVIDENCE_TRIADS, CORPUS_TAG, DIC_TAG))
f_evBC <- file.path(DIR_TABLE, sprintf("abc_evidence_BC_primary_top%d_%s__%s.csv", TOP_EVIDENCE_TRIADS, CORPUS_TAG, DIC_TAG))
f_evAB_ae <- file.path(DIR_TABLE, sprintf("abc_evidence_AB_AEILD_top%d_%s__%s.csv", TOP_EVIDENCE_TRIADS, CORPUS_TAG, DIC_TAG))
f_evBC_ae <- file.path(DIR_TABLE, sprintf("abc_evidence_BC_AEILD_top%d_%s__%s.csv", TOP_EVIDENCE_TRIADS, CORPUS_TAG, DIC_TAG))
f_edges <- file.path(DIR_TABLE, sprintf("abc_edges_AEILD_top%d_%s__%s.csv", TOP_AE_NETWORK_TRIADS, CORPUS_TAG, DIC_TAG))

write_csv2(tab_all, f_all)
write_csv2(tab_main, f_main)
write_csv2(tab_ae, f_ae)
write_csv2(context_pair_stats, f_ctx)
write_csv2(ev_all$AB, f_evAB)
write_csv2(ev_all$BC, f_evBC)
write_csv2(ev_ae$AB, f_evAB_ae)
write_csv2(ev_ae$BC, f_evBC_ae)
write_csv2(edges_ae, f_edges)

# Compatibility output: AE-ILD main filtered ranking, because old downstream scripts often expected abc_rankings_*.csv.
f_compat <- file.path(DIR_TABLE, sprintf("abc_rankings_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(tab_ae, f_compat)

summary_tab <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  A_terms = length(A_set),
  B_terms = length(B_set),
  C_terms = length(C_set),
  C_list = paste(C_set, collapse = ";"),
  all_candidate_triads_with_AB_and_BC_support = nrow(tab_all),
  main_display_triads = nrow(tab_main),
  AE_ILD_main_display_triads = nrow(tab_ae),
  min_AB_n11 = MIN_AB_N11,
  min_BC_n11 = MIN_BC_N11,
  treatment_context_terms_excluded_from_primary_ABC = length(intersect(treatment_context_terms(dic), terms_present)),
  infection_context_terms_excluded_from_primary_ABC = length(intersect(infection_context_terms(dic), terms_present)),
  interpretation = "Primary ABC uses primary_atlas A/B/C terms only. Treatment and infection-context terms are retained in context pair diagnostics but excluded from primary disease-state bridge claims."
)
f_summary <- file.path(DIR_TABLE, sprintf("abc_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(summary_tab, f_summary)

log_msg("04 ABC rows all=", nrow(tab_all), " main=", nrow(tab_main), " AE=", nrow(tab_ae))
log_msg("=== DONE ", RUN_NAME, " ===")
