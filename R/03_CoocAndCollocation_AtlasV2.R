# =========================================================
# 03_CoocAndCollocation_AtlasV2.R
# RA-ILD atlas v2.2 core analysis
#
# Purpose
#   - Compute document-level co-occurrence between dictionary terms and
#     primary worsening outcomes.
#   - Keep primary disease-state atlas terms separate from treatment-context
#     and infection-context terms.
#   - Write layer-aware tables for Figure 2 planning and supplementary checks.
#   - Optionally compute corpus collocations from abstract text.
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)
quiet_install(c("dplyr", "readr", "stringr", "tibble", "purrr"))

RUN_NAME <- "03_CoocAndCollocation_AtlasV2"
log_msg("=== START ", RUN_NAME, " ===")

# -------------------------
# 1) Load inputs
# -------------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)
df  <- load_hit_matrix()
terms_present <- hit_terms_present(df, dic)

C_primary <- intersect(select_primary_C_terms(dic, primary_outcomes_only = TRUE), terms_present)
C_primary <- unique(c(intersect(PRIMARY_OUTCOME_TERMS, C_primary), setdiff(C_primary, PRIMARY_OUTCOME_TERMS)))
if (!length(C_primary)) stop("No primary outcome terms found in hit matrix. Expected: ", paste(PRIMARY_OUTCOME_TERMS, collapse = ", "))

primary_terms <- intersect(select_primary_dictionary_terms(dic) |>
                             dplyr::filter(abc_role %in% c("A", "B")) |>
                             dplyr::pull(term), terms_present)

treatment_terms <- intersect(treatment_context_terms(dic), terms_present)
infection_terms <- intersect(infection_context_terms(dic), terms_present)
exploratory_terms <- intersect(terms_by_atlas_layer(dic, "exploratory"), terms_present)
sensitivity_terms <- intersect(terms_by_atlas_layer(dic, "sensitivity"), terms_present)
context_only_terms <- intersect(terms_by_atlas_layer(dic, "context_only"), terms_present)

# Full set excludes primary outcome nodes themselves so that AE-ILD vs AE-ILD is not shown.
all_non_primary_outcomes <- setdiff(terms_present, C_primary)

log_msg("03 term sets | primary=", length(primary_terms),
        " treatment=", length(treatment_terms),
        " infection=", length(infection_terms),
        " exploratory=", length(exploratory_terms),
        " sensitivity=", length(sensitivity_terms),
        " context_only=", length(context_only_terms),
        " outcomes=", paste(C_primary, collapse = ", "))

# -------------------------
# 2) Pair table helper
# -------------------------
make_term_outcome_pairs <- function(terms, outcomes, scope_name) {
  if (!length(terms) || !length(outcomes)) {
    return(tibble::tibble(scope = character(), term = character(), outcome = character()))
  }
  tab <- purrr::map_dfr(terms, function(t) {
    purrr::map_dfr(outcomes, function(c) {
      ps <- pair_stats(df, t, c)
      tibble::tibble(scope = scope_name, term = t, outcome = c) |>
        dplyr::bind_cols(ps)
    })
  })
  if (!nrow(tab)) return(tab)
  tab <- tab |>
    dplyr::mutate(q = p.adjust(p, method = "BH")) |>
    add_term_metadata(dic, term_col = "term", prefix = "term") |>
    add_term_metadata(dic, term_col = "outcome", prefix = "outcome") |>
    dplyr::arrange(scope, outcome, dplyr::desc(n11), dplyr::desc(npmi), term)
  tab
}

cooc_primary <- make_term_outcome_pairs(primary_terms, C_primary, "primary_atlas")
cooc_treat   <- make_term_outcome_pairs(treatment_terms, C_primary, "treatment_context")
cooc_inf     <- make_term_outcome_pairs(infection_terms, C_primary, "infection_context")
cooc_expl    <- make_term_outcome_pairs(exploratory_terms, C_primary, "exploratory")
cooc_sens    <- make_term_outcome_pairs(sensitivity_terms, C_primary, "sensitivity")
cooc_ctx     <- make_term_outcome_pairs(context_only_terms, C_primary, "context_only")
cooc_all     <- dplyr::bind_rows(cooc_primary, cooc_treat, cooc_inf, cooc_expl, cooc_sens, cooc_ctx)

# -------------------------
# 3) Layer and domain summaries
# -------------------------
layer_summary <- cooc_all |>
  dplyr::group_by(scope, outcome) |>
  dplyr::summarise(
    terms_tested = dplyr::n_distinct(term),
    terms_with_cohit = dplyr::n_distinct(term[n11 > 0]),
    total_cohits = sum(n11, na.rm = TRUE),
    max_cohit = max(n11, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(outcome, scope)

domain_summary <- cooc_all |>
  dplyr::group_by(scope, outcome, term_domain, term_abc_role, term_abc_subrole) |>
  dplyr::summarise(
    terms_tested = dplyr::n_distinct(term),
    terms_with_cohit = dplyr::n_distinct(term[n11 > 0]),
    total_cohits = sum(n11, na.rm = TRUE),
    top_term = term[which.max(n11)],
    top_term_cohit = max(n11, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(scope, outcome, dplyr::desc(total_cohits))

# -------------------------
# 4) Save outputs
# -------------------------
f_primary <- file.path(DIR_TABLE, sprintf("cooc_primary_atlas_to_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_treat   <- file.path(DIR_TABLE, sprintf("cooc_treatment_context_to_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_inf     <- file.path(DIR_TABLE, sprintf("cooc_infection_context_to_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_all     <- file.path(DIR_TABLE, sprintf("cooc_all_layers_to_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_layer   <- file.path(DIR_TABLE, sprintf("cooc_layer_summary_to_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_domain  <- file.path(DIR_TABLE, sprintf("cooc_domain_summary_to_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))

write_csv2(cooc_primary, f_primary)
write_csv2(cooc_treat,   f_treat)
write_csv2(cooc_inf,     f_inf)
write_csv2(cooc_all,     f_all)
write_csv2(layer_summary, f_layer)
write_csv2(domain_summary, f_domain)

# Backward-compatible generic co-occurrence table for later scripts if needed.
f_compat <- file.path(DIR_TABLE, sprintf("cooc_terms_outcomes_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(cooc_all, f_compat)

# -------------------------
# 5) Optional collocations from abstract text
# -------------------------
load_articles_for_text <- function() {
  if ("text" %in% names(df) && mean(nzchar(ifelse(is.na(df$text), "", df$text))) > 0.90) return(df)

  env_file <- Sys.getenv("RAILD_ARTICLE_FILE", unset = "")
  if (nzchar(env_file) && file.exists(path.expand(env_file))) {
    art_file <- path.expand(env_file)
  } else {
    cand <- unique(c(
      Sys.glob(file.path(DIR_PROC, "articles_main_original_*.csv")),
      Sys.glob(file.path(DIR_PROC, "articles_original_*.csv")),
      Sys.glob(file.path(DIR_PROC, "articles_*.csv"))
    ))
    cand <- cand[file.exists(cand)]
    cand <- cand[!grepl("review|guideline|case_report|excluded|articles_all", basename(cand), ignore.case = TRUE)]
    if (!length(cand)) return(tibble::tibble())
    art_file <- cand[which.max(file.info(cand)$mtime)]
  }
  art <- readr::read_csv(art_file, show_col_types = FALSE)
  if ("corpus_role" %in% names(art)) art <- art |> dplyr::filter(corpus_role == "main_original")
  if (!"text" %in% names(art) && all(c("title", "abstract") %in% names(art))) {
    art <- art |> dplyr::mutate(text = paste(ifelse(is.na(title), "", title), ifelse(is.na(abstract), "", abstract), sep = " "))
  }
  if (!"pmid" %in% names(art) || !"text" %in% names(art)) return(tibble::tibble())
  art |> dplyr::mutate(pmid = as.character(pmid)) |> dplyr::semi_join(df |> dplyr::select(pmid), by = "pmid")
}

articles <- load_articles_for_text()
if (nrow(articles) && "text" %in% names(articles) &&
    requireNamespace("quanteda", quietly = TRUE) &&
    requireNamespace("quanteda.textstats", quietly = TRUE)) {

  sw_general <- quanteda::stopwords("en")
  sw_domain <- c(
    "rheumatoid", "arthritis", "interstitial", "lung", "disease", "ild", "ip", "pulmonary", "fibrosis",
    "patients", "patient", "study", "studies", "review", "case", "cases", "report", "reports",
    "introduction", "conclusion", "methods", "background", "objective", "aim", "result", "results", "purpose",
    "mg", "ml", "day", "days", "week", "weeks", "year", "years"
  )

  corp <- quanteda::corpus(articles, text_field = "text")
  toks <- corp |>
    quanteda::tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE) |>
    quanteda::tokens_tolower() |>
    quanteda::tokens_remove(c(sw_general, sw_domain))

  coll2 <- quanteda.textstats::textstat_collocations(toks, size = 2, min_count = 10, smoothing = 0.5)
  coll3 <- quanteda.textstats::textstat_collocations(toks, size = 3, min_count = 8,  smoothing = 0.5)

  clean_noise <- function(D) {
    if (is.null(D) || !nrow(D)) return(D)
    D |>
      dplyr::mutate(phrase = stringr::str_replace_all(collocation, "_", " ")) |>
      dplyr::filter(!grepl("\\d", phrase), nchar(phrase) >= 5, !grepl("(^[a-z]$|^[-_]+$)", phrase))
  }

  f2 <- file.path(DIR_TABLE, sprintf("collocation_bigram_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  f3 <- file.path(DIR_TABLE, sprintf("collocation_trigram_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  write_csv2(clean_noise(coll2) |> dplyr::arrange(dplyr::desc(lambda)), f2)
  write_csv2(clean_noise(coll3) |> dplyr::arrange(dplyr::desc(lambda)), f3)
} else {
  log_msg("03 collocations skipped: article text not available or quanteda/quanteda.textstats not installed.")
}

log_msg("=== DONE ", RUN_NAME, " ===")
