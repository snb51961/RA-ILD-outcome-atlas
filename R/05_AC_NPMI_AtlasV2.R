# =========================================================
# 05_AC_NPMI_AtlasV2.R
# RA-ILD atlas v2.2 core analysis
#
# Purpose
#   - Compute A/C or term/outcome coherence using lift and NPMI.
#   - Separate primary disease-state atlas from treatment-context and
#     infection-context layers.
#   - Provide clean inputs for outcome-domain matrix and Supplementary coherence figures.
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)
quiet_install(c("dplyr", "readr", "stringr", "tibble", "purrr"))

RUN_NAME <- "05_AC_NPMI_AtlasV2"
log_msg("=== START ", RUN_NAME, " ===")

# -------------------------
# 1) Load inputs
# -------------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)
df <- load_hit_matrix()
terms_present <- hit_terms_present(df, dic)

C_primary <- intersect(select_primary_C_terms(dic, primary_outcomes_only = TRUE), terms_present)
C_primary <- unique(c(intersect(PRIMARY_OUTCOME_TERMS, C_primary), setdiff(C_primary, PRIMARY_OUTCOME_TERMS)))
if (!length(C_primary)) stop("No primary outcome terms found for AC/NPMI analysis.")

primary_signal_terms <- intersect(select_signed_primary_terms(dic), terms_present)
primary_A_terms <- intersect(select_primary_A_terms(dic), terms_present)
primary_B_terms <- intersect(select_primary_B_terms(dic), terms_present)
treatment_terms <- intersect(treatment_context_terms(dic), terms_present)
infection_terms <- intersect(infection_context_terms(dic), terms_present)
exploratory_terms <- intersect(terms_by_atlas_layer(dic, "exploratory"), terms_present)
sensitivity_terms <- intersect(terms_by_atlas_layer(dic, "sensitivity"), terms_present)
context_only_terms <- intersect(terms_by_atlas_layer(dic, "context_only"), terms_present)
all_non_outcome_terms <- setdiff(terms_present, C_primary)

log_msg("05 outcome terms:", paste(C_primary, collapse = ", "))
log_msg("05 signal sets | primary=", length(primary_signal_terms),
        " A=", length(primary_A_terms),
        " B=", length(primary_B_terms),
        " treatment=", length(treatment_terms),
        " infection=", length(infection_terms))

# -------------------------
# 2) Coherence helper
# -------------------------
make_npmi <- function(terms, outcomes, scope_name) {
  if (!length(terms) || !length(outcomes)) {
    return(tibble::tibble(scope = character(), A = character(), C = character()))
  }
  tab <- purrr::map_dfr(terms, function(a) {
    purrr::map_dfr(outcomes, function(c) {
      ps <- pair_stats(df, a, c)
      tibble::tibble(scope = scope_name, A = a, C = c) |>
        dplyr::bind_cols(ps)
    })
  })
  if (!nrow(tab)) return(tab)
  tab |>
    dplyr::mutate(q = p.adjust(p, method = "BH")) |>
    add_term_metadata(dic, term_col = "A", prefix = "A") |>
    add_term_metadata(dic, term_col = "C", prefix = "C") |>
    dplyr::arrange(scope, C, dplyr::desc(n11), dplyr::desc(npmi), A)
}

npmi_primary <- make_npmi(primary_signal_terms, C_primary, "primary_atlas")
npmi_A       <- make_npmi(primary_A_terms, C_primary, "primary_A_only")
npmi_B       <- make_npmi(primary_B_terms, C_primary, "primary_B_only")
npmi_treat   <- make_npmi(treatment_terms, C_primary, "treatment_context")
npmi_inf     <- make_npmi(infection_terms, C_primary, "infection_context")
npmi_expl    <- make_npmi(exploratory_terms, C_primary, "exploratory")
npmi_sens    <- make_npmi(sensitivity_terms, C_primary, "sensitivity")
npmi_ctx     <- make_npmi(context_only_terms, C_primary, "context_only")
npmi_all     <- make_npmi(all_non_outcome_terms, C_primary, "all_layers")

combined_scoped <- dplyr::bind_rows(npmi_primary, npmi_treat, npmi_inf, npmi_expl, npmi_sens, npmi_ctx)

# Domain-level summary for clinical reading.
domain_summary <- combined_scoped |>
  dplyr::group_by(scope, C, A_domain, A_abc_role, A_abc_subrole) |>
  dplyr::summarise(
    terms_tested = dplyr::n_distinct(A),
    terms_with_cohit = dplyr::n_distinct(A[n11 > 0]),
    total_cohits = sum(n11, na.rm = TRUE),
    median_npmi_nonzero = ifelse(any(n11 > 0), median(npmi[n11 > 0], na.rm = TRUE), NA_real_),
    top_term = A[which.max(n11)],
    top_term_cohit = max(n11, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(scope, C, dplyr::desc(total_cohits))

# -------------------------
# 3) Save outputs
# -------------------------
f_primary <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_primary_atlas_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_A       <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_primary_A_only_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_B       <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_primary_B_only_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_treat   <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_treatment_context_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_inf     <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_infection_context_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_all     <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_all_layers_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_scoped  <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_scoped_layers_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_domain  <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_domain_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))

write_csv2(npmi_primary, f_primary)
write_csv2(npmi_A,       f_A)
write_csv2(npmi_B,       f_B)
write_csv2(npmi_treat,   f_treat)
write_csv2(npmi_inf,     f_inf)
write_csv2(npmi_all,     f_all)
write_csv2(combined_scoped, f_scoped)
write_csv2(domain_summary, f_domain)

# Compatibility output: primary atlas coherence.
f_compat <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(npmi_primary, f_compat)

run_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  primary_signal_terms = length(primary_signal_terms),
  treatment_context_terms = length(treatment_terms),
  infection_context_terms = length(infection_terms),
  outcomes = paste(C_primary, collapse = ";"),
  primary_pairs = nrow(npmi_primary),
  primary_pairs_with_cohit = sum(npmi_primary$n11 > 0, na.rm = TRUE),
  treatment_pairs_with_cohit = sum(npmi_treat$n11 > 0, na.rm = TRUE),
  infection_pairs_with_cohit = sum(npmi_inf$n11 > 0, na.rm = TRUE),
  interpretation = "Primary coherence is computed separately from treatment and infection-context layers to avoid conflating disease-state bridges with safety/context signals."
)
f_summary <- file.path(DIR_TABLE, sprintf("cooc_npmi_lift_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(run_summary, f_summary)

log_msg("=== DONE ", RUN_NAME, " ===")
