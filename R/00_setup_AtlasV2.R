# =========================================================
# 00_setup_AtlasV2.R
# RA-ILD literature-derived atlas pipeline
# Main dictionary: analysis_v2_outcome_preserving_atlas
# Dictionary policy: v2.2 final, with separated treatment-context and infection-context layers
# =========================================================

ensure_packages <- function(pkgs) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) install.packages(missing, dependencies = TRUE)
  invisible(suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE)))
}

ensure_packages(c("readr", "dplyr", "stringr", "tibble", "purrr", "lubridate"))

# -------------------------
# 0) ROOT auto-detect + optional override
# -------------------------
# To run outside the project root:
#   Sys.setenv(RAILD_ROOT = "/path/to/RAILD_project")
ROOT_OVERRIDE <- Sys.getenv("RAILD_ROOT", unset = NA_character_)

guess_root <- function() {
  candidates <- c(ROOT_OVERRIDE, getwd(), dirname(getwd()))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) {
    if (dir.exists(p) && (dir.exists(file.path(p, "dic")) || dir.exists(file.path(p, "data_raw")) || dir.exists(file.path(p, "data_proc")))) return(p)
  }
  getwd()
}

ROOT <- guess_root()

# -------------------------
# 1) Corpus window and dictionary identity
# -------------------------
YEAR_MIN <- 1980L
DATE_MAX <- "2025/12/31"
YEAR_MAX <- suppressWarnings(as.integer(substr(DATE_MAX, 1, 4)))
if (!is.finite(YEAR_MAX)) stop("DATE_MAX is invalid: ", DATE_MAX)

CORPUS_TAG <- sprintf("pm_%d_%s", YEAR_MIN, gsub("[^0-9]", "", DATE_MAX))

DIC_TAG <- "analysis_v2_outcome_preserving_atlas"
DIC_LABEL <- "final v2.2 outcome-preserving atlas dictionary with treatment- and infection-context layers"
DIC_FILE_NAME <- "ra_ild_dictionary_analysis_v2_outcome_preserving_atlas.csv"

# -------------------------
# 2) Directories
# -------------------------
DIR_RAW   <- file.path(ROOT, "data_raw")
DIR_PROC  <- file.path(ROOT, "data_proc")
DIR_DIC   <- file.path(ROOT, "dic")
DIR_FIG   <- file.path(ROOT, "fig")
DIR_OUT   <- file.path(ROOT, "output", CORPUS_TAG, DIC_TAG)
DIR_LOG   <- file.path(DIR_OUT, "log")
DIR_TABLE <- file.path(DIR_OUT, "table")
DIR_FIG2  <- file.path(DIR_OUT, "fig")

for (d in c(DIR_RAW, DIR_PROC, DIR_DIC, DIR_FIG, DIR_OUT, DIR_LOG, DIR_TABLE, DIR_FIG2)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

DIC_FILE <- file.path(DIR_DIC, DIC_FILE_NAME)

# Optional QC/provenance files produced during dictionary construction.
DIC_QC_FILE          <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v2_2_final_QC.csv")
DIC_CHANGE_LOG_FILE  <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v2_2_change_log.csv")
DIC_LAYER_SUMMARY    <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v2_2_layer_summary.csv")
DIC_INFECTION_QC     <- file.path(DIR_DIC, "ra_ild_dictionary_analysis_v2_2_infection_context_QC.csv")

# -------------------------
# 3) Dictionary policy constants
# -------------------------
PRIMARY_ANALYSIS_TIERS <- c("main")
PRIMARY_NODE_LEVELS    <- c("specific")
PRIMARY_ABC_ROLES      <- c("A", "B", "C")
PRIMARY_ATLAS_LAYERS   <- c("primary_atlas")
TREATMENT_CONTEXT_LAYERS <- c("treatment_context")
INFECTION_CONTEXT_LAYERS <- c("infection_context")
SUPPORTING_ATLAS_LAYERS  <- c("treatment_context", "infection_context", "sensitivity", "exploratory", "context_only")

PRIMARY_OUTCOME_TERMS <- c("AE-ILD", "progression", "mortality")
PROGRESSION_COMPONENT_TERMS <- c("FVC_decline", "DLCO_decline", "PPF", "PF_ILD")
AE_COMPONENT_TERMS <- c("AE_interstitial_pneumonia", "acute_respiratory_deterioration")
SECONDARY_SEVERE_OUTCOME_TERMS <- c("hospitalization", "respiratory_failure", "LTOT", "oxygen_requirement", "mechanical_ventilation", "ICU_admission", "lung_transplantation")

# Never force ignore_case=TRUE outside the regex. The dictionary contains scoped inline flags
# and selected case-sensitive acronym guards.
REGEX_IGNORE_CASE_EXTERNALLY <- FALSE

# Main ABC display filters. Full candidate tables are also written.
MIN_AB_N11_MAIN <- 3L
MIN_BC_N11_MAIN <- 3L
MIN_AC_N11_CONTEXT <- 1L
LIFT_MIN_MAIN <- 1.00
NPMI_MIN_MAIN <- 0.00

# -------------------------
# 4) Logging
# -------------------------
log_file <- file.path(DIR_LOG, "pipeline.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  message(msg)
  readr::write_lines(msg, log_file, append = TRUE)
}

stamp_date <- function() format(Sys.Date(), "%Y%m%d")
stamp_time <- function() format(Sys.time(), "%Y%m%d_%H%M")

write_csv2 <- function(df, file) {
  readr::write_csv(df, file)
  log_msg("WROTE:", file)
  invisible(file)
}

quiet_install <- ensure_packages

# -------------------------
# 5) Dictionary loader and selectors
# -------------------------
assert_atlas_v2_setup <- function(require_dictionary = TRUE) {
  if (!exists("DIC_TAG") || !identical(as.character(DIC_TAG), "analysis_v2_outcome_preserving_atlas")) {
    stop("Atlas v2 setup guard failed. Current DIC_TAG = ", if (exists("DIC_TAG")) as.character(DIC_TAG) else "<missing>")
  }
  if (require_dictionary && (!exists("DIC_FILE") || !is.character(DIC_FILE) || length(DIC_FILE) != 1 || !file.exists(DIC_FILE))) {
    stop("Atlas v2 dictionary not found: ", if (exists("DIC_FILE")) as.character(DIC_FILE) else "<missing>",
         "\nCopy ", DIC_FILE_NAME, " into DIR_DIC: ", DIR_DIC)
  }
  invisible(TRUE)
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

load_atlas_dictionary <- function(require_dictionary = TRUE) {
  assert_atlas_v2_setup(require_dictionary = require_dictionary)
  dic <- readr::read_csv(DIC_FILE, show_col_types = FALSE)
  required_cols <- c("term", "regex", "class")
  miss <- setdiff(required_cols, names(dic))
  if (length(miss)) stop("Dictionary is missing required columns: ", paste(miss, collapse = ", "))
  if (anyDuplicated(dic$term)) {
    dup_terms <- unique(dic$term[duplicated(dic$term)])
    stop("Duplicate dictionary terms found: ", paste(dup_terms, collapse = ", "))
  }

  default_cols <- list(
    include_primary = TRUE,
    analysis_tier = "main",
    node_level = "specific",
    abc_role = NA_character_,
    abc_subrole = NA_character_,
    domain = NA_character_,
    concept_family = NA_character_,
    atlas_layer = NA_character_,
    interpretation_note = NA_character_,
    primary_atlas_eligible = NA
  )
  for (nm in names(default_cols)) {
    if (!nm %in% names(dic)) dic[[nm]] <- default_cols[[nm]]
  }

  dic <- dic |>
    dplyr::mutate(
      include_primary = as_bool(include_primary),
      analysis_tier = as.character(analysis_tier),
      node_level = as.character(node_level),
      abc_role = as.character(abc_role),
      abc_subrole = as.character(abc_subrole),
      domain = as.character(domain),
      concept_family = as.character(concept_family),
      atlas_layer = dplyr::case_when(
        !is.na(atlas_layer) & nzchar(as.character(atlas_layer)) ~ as.character(atlas_layer),
        analysis_tier == "treatment_context" ~ "treatment_context",
        analysis_tier == "infection_context" ~ "infection_context",
        analysis_tier == "exploratory" ~ "exploratory",
        analysis_tier == "context" ~ "context_only",
        analysis_tier == "sensitivity" ~ "sensitivity",
        include_primary == TRUE & analysis_tier == "main" ~ "primary_atlas",
        TRUE ~ "sensitivity"
      ),
      interpretation_note = as.character(interpretation_note),
      primary_atlas_eligible = include_primary == TRUE &
        analysis_tier %in% PRIMARY_ANALYSIS_TIERS &
        node_level %in% PRIMARY_NODE_LEVELS &
        abc_role %in% PRIMARY_ABC_ROLES &
        atlas_layer %in% PRIMARY_ATLAS_LAYERS
    )
  dic
}

select_primary_dictionary_terms <- function(dic) {
  dic |>
    dplyr::filter(
      include_primary == TRUE,
      analysis_tier %in% PRIMARY_ANALYSIS_TIERS,
      node_level %in% PRIMARY_NODE_LEVELS,
      abc_role %in% PRIMARY_ABC_ROLES,
      atlas_layer %in% PRIMARY_ATLAS_LAYERS
    )
}

select_primary_A_terms <- function(dic) {
  select_primary_dictionary_terms(dic) |>
    dplyr::filter(abc_role == "A") |>
    dplyr::pull(term) |>
    unique()
}

select_primary_B_terms <- function(dic) {
  select_primary_dictionary_terms(dic) |>
    dplyr::filter(abc_role == "B") |>
    dplyr::pull(term) |>
    unique()
}

select_primary_C_terms <- function(dic, primary_outcomes_only = TRUE) {
  d <- select_primary_dictionary_terms(dic) |> dplyr::filter(abc_role == "C")
  if (primary_outcomes_only) {
    d <- d |> dplyr::filter(term %in% PRIMARY_OUTCOME_TERMS | abc_subrole == "C_primary_worsening_outcome")
  }
  d |> dplyr::pull(term) |> unique()
}

select_signed_primary_terms <- function(dic) {
  # Signed-effect panels may use both upstream A terms and disease-state B terms against outcomes.
  select_primary_dictionary_terms(dic) |>
    dplyr::filter(abc_role %in% c("A", "B")) |>
    dplyr::pull(term) |>
    unique()
}

terms_by_atlas_layer <- function(dic, layers) {
  dic |>
    dplyr::filter(atlas_layer %in% layers) |>
    dplyr::pull(term) |>
    unique()
}

treatment_context_terms <- function(dic) terms_by_atlas_layer(dic, TREATMENT_CONTEXT_LAYERS)
infection_context_terms <- function(dic) terms_by_atlas_layer(dic, INFECTION_CONTEXT_LAYERS)

terms_by_abc_role <- function(dic, roles, primary_only = TRUE) {
  d <- dic
  if (primary_only) d <- select_primary_dictionary_terms(d)
  d |> dplyr::filter(abc_role %in% roles) |> dplyr::pull(term) |> unique()
}

terms_by_class <- function(dic, classes, primary_only = FALSE) {
  d <- dic
  if (primary_only) d <- select_primary_dictionary_terms(d)
  d |> dplyr::filter(class %in% classes) |> dplyr::pull(term) |> unique()
}

find_latest_file <- function(dir, regex) {
  if (is.na(dir) || !dir.exists(dir)) return(NA_character_)
  xs <- list.files(dir, pattern = regex, full.names = TRUE)
  if (!length(xs)) return(NA_character_)
  xs[which.max(file.info(xs)$mtime)]
}

find_hits_matrix <- function() {
  f <- file.path(DIR_TABLE, sprintf("hits_matrix_%s__%s.csv", CORPUS_TAG, DIC_TAG))
  if (file.exists(f)) return(f)
  f2 <- find_latest_file(DIR_PROC, sprintf("^hits_matrix_.*__%s\\.csv$", DIC_TAG))
  if (!is.na(f2) && file.exists(f2)) return(f2)
  stop("No tagged hit matrix found. Expected: ", f, "\nRun 02_BuildHitsMatrix_AtlasV2.R first.")
}

load_hit_matrix <- function() {
  f <- find_hits_matrix()
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (!"pmid" %in% names(df)) stop("Hit matrix missing pmid column: ", f)
  df$pmid <- as.character(df$pmid)
  log_msg("Loaded hit matrix:", f, " n=", nrow(df), " cols=", ncol(df))
  df
}

hit_terms_present <- function(df, dic = NULL) {
  terms <- sub("^hit__", "", names(df)[startsWith(names(df), "hit__")])
  if (!is.null(dic)) terms <- intersect(terms, dic$term)
  unique(terms)
}

get_hit_vec <- function(df, term) {
  nm <- paste0("hit__", term)
  if (nm %in% names(df)) as.integer(ifelse(is.na(df[[nm]]), 0L, df[[nm]] > 0L)) else integer(nrow(df))
}

binv <- function(v) as.integer(ifelse(is.na(v), 0L, v > 0L))

pair_stats <- function(df, x_term, y_term) {
  x <- get_hit_vec(df, x_term)
  y <- get_hit_vec(df, y_term)
  n11 <- sum(x == 1L & y == 1L, na.rm = TRUE)
  n10 <- sum(x == 1L & y == 0L, na.rm = TRUE)
  n01 <- sum(x == 0L & y == 1L, na.rm = TRUE)
  n00 <- sum(x == 0L & y == 0L, na.rm = TRUE)
  N <- n11 + n10 + n01 + n00
  OR <- ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))
  p <- tryCatch(
    suppressWarnings(stats::fisher.test(matrix(c(n11, n10, n01, n00), 2, byrow = TRUE))$p.value),
    error = function(e) 1
  )
  pX <- (n11 + n10) / max(1, N)
  pY <- (n11 + n01) / max(1, N)
  pXY <- n11 / max(1, N)
  eps <- 1e-12
  pmi <- if (pXY > 0 && pX > 0 && pY > 0) log((pXY + eps) / (pX * pY + eps)) else 0
  npmi <- if (pXY > 0) pmi / (-log(pXY + eps)) else 0
  lift <- if (pX > 0 && pY > 0) (pXY + eps) / (pX * pY + eps) else 0
  tibble::tibble(
    n11 = as.integer(n11), n10 = as.integer(n10), n01 = as.integer(n01), n00 = as.integer(n00),
    p_x = pX, p_y = pY, p_xy = pXY,
    or = as.numeric(OR), p = as.numeric(p), lift = as.numeric(lift), npmi = as.numeric(npmi)
  )
}

add_term_metadata <- function(tab, dic, term_col = "term", prefix = "term") {
  meta_cols <- c("term", "preferred_label", "class", "abc_role", "abc_subrole", "domain", "concept_family", "node_level", "analysis_tier", "atlas_layer", "include_primary", "primary_atlas_eligible")
  meta <- dic |> dplyr::select(dplyr::any_of(meta_cols))
  names(meta) <- ifelse(names(meta) == "term", term_col, paste0(prefix, "_", names(meta)))
  dplyr::left_join(tab, meta, by = term_col)
}

set.seed(42)
log_msg("Atlas v2 setup loaded | ROOT=", ROOT, " | DIC_LABEL=", DIC_LABEL)
