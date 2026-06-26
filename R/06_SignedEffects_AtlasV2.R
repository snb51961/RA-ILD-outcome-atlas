# =========================================================
# 06_SignedEffects_AtlasV2.R
# RA-ILD atlas v2.2 core analysis
#
# Purpose
#   - Extract directional wording for term -> outcome relations from abstracts.
#   - Use primary disease-state atlas terms as the main scope.
#   - Keep treatment-context and infection-context signed signals separate.
#
# Interpretation
#   - Directional labels describe abstract wording, not causality.
#   - Treatment-context and infection-context rows must not be interpreted as
#     treatment effect, drug risk, infection causality, or clinical prediction.
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = TRUE)
quiet_install(c("dplyr", "readr", "stringr", "tibble", "purrr"))

RUN_NAME <- "06_SignedEffects_AtlasV2"
log_msg("=== START ", RUN_NAME, " ===")

# -------------------------
# 0) Parameters
# -------------------------
OUTCOMES_FOR_SIGNED <- PRIMARY_OUTCOME_TERMS
CHAR_WINDOW <- 180L
MAX_CONTEXTS <- 3L
SAME_SENTENCE_ONLY <- FALSE
WINDOW <- 1L
MAX_TERMS_PER_SCOPE <- NULL  # set to e.g. 100 for debugging only; keep NULL for manuscript runs

# -------------------------
# 1) Load inputs
# -------------------------
dic <- load_atlas_dictionary(require_dictionary = TRUE)
DF <- load_hit_matrix()
terms_present <- hit_terms_present(DF, dic)

# Ensure text is available. 02 normally writes text into the hit matrix.
load_articles_for_join <- function() {
  env_file <- Sys.getenv("RAILD_ARTICLE_FILE", unset = "")
  if (nzchar(env_file) && file.exists(path.expand(env_file))) {
    f_art <- path.expand(env_file)
  } else {
    cand <- unique(c(
      Sys.glob(file.path(DIR_PROC, "articles_main_original_*.csv")),
      Sys.glob(file.path(DIR_PROC, "articles_original_*.csv")),
      Sys.glob(file.path(DIR_PROC, "articles_*.csv"))
    ))
    cand <- cand[file.exists(cand)]
    cand <- cand[!grepl("review|guideline|case_report|excluded|articles_all", basename(cand), ignore.case = TRUE)]
    if (!length(cand)) return(tibble::tibble())
    f_art <- cand[which.max(file.info(cand)$mtime)]
  }
  art <- readr::read_csv(f_art, show_col_types = FALSE)
  if ("corpus_role" %in% names(art)) art <- art |> dplyr::filter(corpus_role == "main_original")
  if (!"text" %in% names(art) && all(c("title", "abstract") %in% names(art))) {
    art <- art |> dplyr::mutate(text = paste(ifelse(is.na(title), "", title), ifelse(is.na(abstract), "", abstract), sep = " "))
  }
  if (!"pmid" %in% names(art) || !"text" %in% names(art)) return(tibble::tibble())
  art |>
    dplyr::mutate(pmid = as.character(pmid)) |>
    dplyr::select(dplyr::any_of(c("pmid", "title", "abstract", "text", "year", "journal", "pubtype"))) |>
    dplyr::distinct(pmid, .keep_all = TRUE)
}

if (!"text" %in% names(DF) || mean(nzchar(ifelse(is.na(DF$text), "", DF$text))) < 0.90) {
  ART <- load_articles_for_join()
  if (!nrow(ART)) stop("Signed-effect extraction requires abstract text. Hit matrix lacks text and no article file could be found.")
  DF <- DF |>
    dplyr::left_join(ART, by = "pmid", suffix = c("", ".article"))
  if ("text.article" %in% names(DF)) {
    DF$text <- ifelse(is.na(DF$text) | !nzchar(DF$text), DF$text.article, DF$text)
  }
}

if (!"text" %in% names(DF)) stop("DF has no text column after loading/joining.")
if (!"year" %in% names(DF)) DF$year <- NA_integer_
if (!"pubtype" %in% names(DF)) DF$pubtype <- "Original research"
if (!"journal" %in% names(DF)) DF$journal <- NA_character_
if (!"title" %in% names(DF)) DF$title <- NA_character_

DF <- DF |>
  dplyr::mutate(
    pmid = as.character(pmid),
    year = suppressWarnings(as.integer(year)),
    text = ifelse(is.na(text), "", text),
    pubtype = ifelse(is.na(pubtype), "Original research", pubtype)
  )

text_nonempty_rate <- mean(nzchar(DF$text))
log_msg(sprintf("06 text non-empty rate: %.3f", text_nonempty_rate))
if (is.na(text_nonempty_rate) || text_nonempty_rate < 0.90) {
  stop("Too few non-empty texts for signed-effect extraction: ", round(text_nonempty_rate, 3))
}

# -------------------------
# 2) Term scopes
# -------------------------
C_terms <- intersect(OUTCOMES_FOR_SIGNED, terms_present)
if (!length(C_terms)) stop("No signed-effect outcome terms found in hit matrix: ", paste(OUTCOMES_FOR_SIGNED, collapse = ", "))

primary_terms <- intersect(select_signed_primary_terms(dic), terms_present)
treatment_terms <- intersect(treatment_context_terms(dic), terms_present)
infection_terms <- intersect(infection_context_terms(dic), terms_present)

# Avoid scanning outcome terms as A.
primary_terms <- setdiff(primary_terms, C_terms)
treatment_terms <- setdiff(treatment_terms, C_terms)
infection_terms <- setdiff(infection_terms, C_terms)

cap_terms <- function(terms, scope) {
  if (is.null(MAX_TERMS_PER_SCOPE) || length(terms) <= MAX_TERMS_PER_SCOPE) return(terms)
  freq <- vapply(terms, function(t) sum(get_hit_vec(DF, t) == 1L, na.rm = TRUE), integer(1))
  kept <- names(sort(freq, decreasing = TRUE))[seq_len(MAX_TERMS_PER_SCOPE)]
  log_msg("06 capped ", scope, " terms to ", length(kept), " by frequency.")
  kept
}

primary_terms <- cap_terms(primary_terms, "primary_atlas")
treatment_terms <- cap_terms(treatment_terms, "treatment_context")
infection_terms <- cap_terms(infection_terms, "infection_context")

scopes <- list(
  primary_atlas = primary_terms,
  treatment_context = treatment_terms,
  infection_context = infection_terms
)
log_msg("06 signed scopes | primary=", length(primary_terms),
        " treatment=", length(treatment_terms),
        " infection=", length(infection_terms),
        " outcomes=", paste(C_terms, collapse = ", "))

# -------------------------
# 3) Pair diagnostics
# -------------------------
make_pair_diag <- function(scope_name, terms) {
  if (!length(terms)) return(tibble::tibble())
  tibble::as_tibble(expand.grid(analysis_scope = scope_name, A = unique(terms), C = unique(C_terms), stringsAsFactors = FALSE)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_A_hits = sum(get_hit_vec(DF, A) == 1L, na.rm = TRUE),
      n_C_hits = sum(get_hit_vec(DF, C) == 1L, na.rm = TRUE),
      n_cohit  = sum(get_hit_vec(DF, A) == 1L & get_hit_vec(DF, C) == 1L, na.rm = TRUE)
    ) |>
    dplyr::ungroup()
}

pair_diag <- dplyr::bind_rows(purrr::imap(scopes, function(terms, scope_name) make_pair_diag(scope_name, terms))) |>
  add_term_metadata(dic, term_col = "A", prefix = "A") |>
  add_term_metadata(dic, term_col = "C", prefix = "C") |>
  dplyr::arrange(analysis_scope, dplyr::desc(n_cohit), A, C)

f_pair_diag <- file.path(DIR_TABLE, sprintf("signed_effects_pair_diagnostics_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(pair_diag, f_pair_diag)
log_msg("06 cohit pairs by scope: ", paste(pair_diag |>
                                            dplyr::group_by(analysis_scope) |>
                                            dplyr::summarise(n = sum(n_cohit > 0), .groups = "drop") |>
                                            dplyr::transmute(x = paste0(analysis_scope, "=", n)) |>
                                            dplyr::pull(x), collapse = " / "))

if (!any(pair_diag$analysis_scope == "primary_atlas" & pair_diag$n_cohit > 0)) {
  stop("No primary_atlas A-C co-hit pairs found for signed-effect extraction. See: ", f_pair_diag)
}

# -------------------------
# 4) Text and direction helpers
# -------------------------
AE_EXPAND_RE <- "(?i)((acute\\s+(exacerbation|worsen(ing|ed)?|deterioration|respiratory\\s+(failure|worsen(ing|ed)?)|decompensation)|acute[-\\s]*on[-\\s]*chronic).{0,80}(interstitial|pulmon|lung|fibros|\\bILD\\b|\\bIP\\b)|(interstitial|pulmon|lung|fibros|\\bILD\\b|\\bIP\\b).{0,80}acute\\s+(exacerbation|worsen(ing|ed)?|deterioration|respiratory\\s+(failure|worsen(ing|ed)?)|decompensation))"
NEG_RE <- "(?i)\\b(no|not|without|neither|never|absence\\s+of|lack\\s+of)\\b"
HEDGE_RE <- "(?i)\\b(may|might|tend\\s+to|trend(s)?\\s+toward|borderline|possibly|suggests?)\\b"

split_sentences <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) return(character(0))
  t <- gsub("[\\r\\n]+", " ", txt)
  xs <- unlist(strsplit(t, "(?<=[\\.!?;:])\\s+|\\s+(?=Results?:|Conclusions?:|Background:|Methods?:|Objectives?:)", perl = TRUE))
  xs <- trimws(xs)
  xs[nzchar(xs)]
}

regex_detect <- function(text, pattern) {
  grepl(pattern, text, perl = TRUE, ignore.case = REGEX_IGNORE_CASE_EXTERNALLY)
}

all_spans <- function(pattern, text) {
  if (is.na(text) || !nzchar(text)) return(list(starts = integer(0), ends = integer(0)))
  m <- gregexpr(pattern, text, perl = TRUE, ignore.case = REGEX_IGNORE_CASE_EXTERNALLY)
  if (length(m) == 0 || m[[1]][1] == -1) return(list(starts = integer(0), ends = integer(0)))
  starts <- as.integer(m[[1]])
  lens <- attr(m[[1]], "match.length")
  list(starts = starts, ends = starts + lens - 1L)
}

span_contexts <- function(text, reA, reC, char_window = CHAR_WINDOW, max_ctx = MAX_CONTEXTS) {
  if (is.na(text) || !nzchar(text)) return(character(0))
  sa <- all_spans(reA, text)
  sc <- all_spans(reC, text)
  if (!length(sa$starts) || !length(sc$starts)) return(character(0))
  comb <- expand.grid(i = seq_along(sa$starts), j = seq_along(sc$starts))
  comb$dist <- mapply(function(i, j) {
    a_mid <- (sa$starts[i] + sa$ends[i]) / 2
    c_mid <- (sc$starts[j] + sc$ends[j]) / 2
    abs(a_mid - c_mid)
  }, comb$i, comb$j)
  comb <- comb[order(comb$dist), , drop = FALSE]
  ctxs <- character(0)
  used <- matrix(FALSE, nrow = length(sa$starts), ncol = length(sc$starts))
  for (k in seq_len(nrow(comb))) {
    i <- comb$i[k]; j <- comb$j[k]
    if (used[i, j]) next
    used[i, ] <- TRUE; used[, j] <- TRUE
    bgn <- max(1L, min(sa$starts[i], sc$starts[j]) - char_window)
    end <- min(nchar(text), max(sa$ends[i], sc$ends[j]) + char_window)
    ctxs <- c(ctxs, substr(text, bgn, end))
    if (length(ctxs) >= max_ctx) break
  }
  unique(ctxs)
}

rx <- list(
  inc = "(?i)(associated\\s+with\\s+(an?\\s+)?increase(d)?|increase(d|s)?\\s+(in|of)|higher\\s+risk|elevat(ed|es)|worsen(ed|ing)|exacerbat(ed|es)|risk\\s+factor|predict(or|ive)\\s+of|trigger(s|ed)?|independent\\s+predictor)",
  dec = "(?i)(associated\\s+with\\s+(a\\s+)?decrease(d)?|decrease(d|s)?\\s+(in|of)|lower\\s+risk|reduc(ed|es)|protective|improv(ed|es|ement)|ameliorat(ed|es)|attenuat(ed|es)|inhibit(ed|s)|stabili[sz](ed|ation))",
  null = "(?i)(not\\s+associated|no\\s+(significant\\s+)?association|no\\s+significant\\s+(difference|effect)|nonsignificant|non[-\\s]?significant|did\\s+not\\s+(find|show)|failed\\s+to\\s+(show|demonstrate)|similar\\s+(rates?|risk)|comparable)",
  meas = "(?i)\\b((?:a\\s*OR|aOR|OR|HR|aHR|RR|IRR|SHR|adj(?:usted)?\\s*(?:odds|hazard|risk)\\s*ratio))\\b\\s*[:=]?\\s*([0-9]+\\.?[0-9]*)",
  ci = "(?i)95%\\s*CI\\s*[:=]?\\s*\\(?\\s*([0-9]+\\.?[0-9]*)\\s*[-–,]\\s*([0-9]+\\.?[0-9]*)\\s*\\)?",
  p = "(?i)\\bp\\s*[<=>]\\s*0\\.?0*([0-9]\\d?)",
  auc = "(?i)\\bAUC\\b\\s*[:=]?\\s*([0-9]\\.?[0-9]*)",
  cut = "(?i)(cut[-\\s]*off|threshold)\\s*[:=]?\\s*([0-9]+\\.?[0-9]*)\\s*(U\\s*/\\s*mL|U/mL|mg/L|ng/mL|pg/mL)?",
  sens = "(?i)(sensitivity)\\s*[:=]?\\s*([0-9]{1,3}\\.?[0-9]*)\\s*%?",
  spec = "(?i)(specificity)\\s*[:=]?\\s*([0-9]{1,3}\\.?[0-9]*)\\s*%?"
)
SERO_NEG_RE <- "(?i)\\b(seronegative|ACPA\\s*negative|anti-?CCP\\s*negative|RF\\s*negative)\\b"

infer_one <- function(s) {
  null_flag <- grepl(rx$null, s, perl = TRUE)
  inc_flag <- grepl(rx$inc, s, perl = TRUE)
  dec_flag <- grepl(rx$dec, s, perl = TRUE)
  neg_flag <- grepl(NEG_RE, s, perl = TRUE)
  if (neg_flag && (inc_flag || dec_flag)) {
    inc_flag <- FALSE; dec_flag <- FALSE; null_flag <- TRUE
  }
  hedge_flag <- grepl(HEDGE_RE, s, perl = TRUE)

  m <- stringr::str_match(s, rx$meas)
  val <- suppressWarnings(as.numeric(if (is.null(m)) NA else m[, 3]))
  meas_name <- if (!is.null(m) && !is.na(m[, 2])) m[, 2] else NA
  ci <- stringr::str_match(s, rx$ci)
  lo <- suppressWarnings(as.numeric(ci[, 2]))
  hi <- suppressWarnings(as.numeric(ci[, 3]))
  p_sig <- grepl(rx$p, s, perl = TRUE)
  m_auc <- stringr::str_match(s, rx$auc); AUC <- suppressWarnings(as.numeric(m_auc[, 2]))
  m_cut <- stringr::str_match(s, rx$cut); CUT <- suppressWarnings(as.numeric(m_cut[, 3])); UNIT <- ifelse(is.na(m_cut[, 4]), NA_character_, m_cut[, 4])
  m_sen <- stringr::str_match(s, rx$sens); SENS <- suppressWarnings(as.numeric(m_sen[, 3]))
  m_spe <- stringr::str_match(s, rx$spec); SPEC <- suppressWarnings(as.numeric(m_spe[, 3]))

  sig_by_ci <- if (!is.na(lo) && !is.na(hi)) if (lo > 1 || hi < 1) TRUE else FALSE else NA
  sig <- if (!is.na(sig_by_ci)) sig_by_ci else if (p_sig) TRUE else NA

  sgn <- NA_real_
  weight <- 1
  if (!is.na(val)) sgn <- if (val > 1) +1 else if (val < 1) -1 else 0
  if (is.na(sgn)) {
    sgn <- if (null_flag) 0 else if (inc_flag && !dec_flag) +1 else if (dec_flag && !inc_flag) -1 else NA
  }
  if (!is.na(sig) && sig) weight <- weight + 1
  if (grepl("(?i)(adjusted|multivaria(te|ted))", s, perl = TRUE)) weight <- weight + 0.5
  if (hedge_flag) weight <- weight * 0.6
  if (is.na(sgn) && grepl(SERO_NEG_RE, s, perl = TRUE)) { sgn <- -1; weight <- weight * 0.8 }

  list(
    sign = sgn, weight = weight,
    inc = as.integer(inc_flag), dec = as.integer(dec_flag), null = as.integer(null_flag),
    measure = if (!is.na(meas_name)) meas_name else NA_character_, value = val,
    ci_low = lo, ci_high = hi, p_sig = as.integer(p_sig), auc = AUC,
    cutoff = CUT, unit = UNIT, sens = SENS, spec = SPEC
  )
}

term_regex <- function(term) {
  if (identical(term, "AE-ILD")) return(AE_EXPAND_RE)
  r <- dic$regex[match(term, dic$term)]
  if (is.na(r) || !nzchar(r)) return(NA_character_)
  r
}

# Design weight is deliberately modest and heuristic; it supports direction extraction only.
design_weight <- function(pubtype, year) {
  year <- suppressWarnings(as.integer(year))
  w <- 1
  if (!is.na(pubtype) && grepl("Randomized|Controlled Clinical Trial", pubtype, ignore.case = TRUE)) w <- w + 0.8
  if (!is.na(pubtype) && grepl("Cohort|Case-Control", pubtype, ignore.case = TRUE)) w <- w + 0.3
  if (!is.na(pubtype) && grepl("Case Reports?", pubtype, ignore.case = TRUE)) w <- w - 0.3
  if (!is.na(year)) w <- w + 0.1 * pmax(0, year - 2015) / 10
  pmax(w, 0.2)
}

analyze_pair <- function(scope_name, a_term, c_term, window = WINDOW, same_sentence_only = SAME_SENTENCE_ONLY) {
  mask <- which(get_hit_vec(DF, a_term) == 1L & get_hit_vec(DF, c_term) == 1L)

  if (identical(c_term, "AE-ILD")) {
    c_extra <- grepl(AE_EXPAND_RE, DF$text, perl = TRUE, ignore.case = FALSE)
    mask <- which(get_hit_vec(DF, a_term) == 1L & (get_hit_vec(DF, c_term) == 1L | c_extra))
  }
  if (!length(mask)) return(list(sent = NULL, art = NULL))

  a_re <- term_regex(a_term)
  c_re <- term_regex(c_term)
  if (is.na(a_re) || is.na(c_re)) return(list(sent = NULL, art = NULL))

  out_sent <- list()
  out_art <- list()

  for (idx in mask) {
    pmid <- DF$pmid[idx]
    yr <- DF$year[idx]
    abs_txt <- DF$text[idx]
    if (is.na(abs_txt) || !nzchar(abs_txt)) next

    ss <- split_sentences(abs_txt)
    contexts <- character(0)
    src_weight <- 1

    Sa <- which(vapply(ss, regex_detect, logical(1), pattern = a_re))
    Sc <- which(vapply(ss, regex_detect, logical(1), pattern = c_re))

    if (length(Sa) && length(Sc)) {
      if (same_sentence_only) {
        S_both <- intersect(Sa, Sc)
        if (length(S_both)) contexts <- c(contexts, ss[S_both])
      } else {
        idxs <- sort(unique(c(Sa, Sc)))
        for (k in idxs) {
          i1 <- max(1, k - window); i2 <- min(length(ss), k + window)
          seg <- paste(ss[i1:i2], collapse = " ")
          if (regex_detect(seg, a_re) && regex_detect(seg, c_re)) contexts <- c(contexts, seg)
        }
      }
    }

    if (!length(contexts)) {
      ctx_span <- span_contexts(abs_txt, a_re, c_re, CHAR_WINDOW, MAX_CONTEXTS)
      if (length(ctx_span)) { contexts <- c(contexts, ctx_span); src_weight <- 0.8 }
    }

    if (!length(contexts) && regex_detect(abs_txt, a_re) && regex_detect(abs_txt, c_re)) {
      contexts <- c(contexts, abs_txt)
      src_weight <- 0.4
    }
    if (!length(contexts)) next

    art_w <- design_weight(DF$pubtype[idx], DF$year[idx])
    S <- dplyr::bind_rows(lapply(unique(contexts), function(s) {
      r <- infer_one(s)
      r$weight <- r$weight * src_weight * art_w
      tibble::tibble(
        analysis_scope = scope_name,
        pmid = pmid, year = yr, A = a_term, C = c_term, sentence = s,
        sign = r$sign, weight = r$weight, inc = r$inc, dec = r$dec, null = r$null,
        measure = r$measure, value = r$value, ci_low = r$ci_low, ci_high = r$ci_high,
        p_sig = r$p_sig, auc = r$auc, cutoff = r$cutoff, unit = r$unit, sens = r$sens, spec = r$spec
      )
    }))
    if (!nrow(S)) next
    out_sent[[length(out_sent) + 1L]] <- S

    sc <- sum(ifelse(is.na(S$sign), 0, S$sign * S$weight), na.rm = TRUE)
    legacy_label <- if (sc > 0) "risk_up" else if (sc < 0) "risk_down" else "no_effect_or_mixed"
    wording_label <- if (sc > 0) "positive_outcome_associated_wording" else if (sc < 0) "decreasing_or_protective_wording" else "mixed_or_no_clear_direction"
    out_art[[length(out_art) + 1L]] <- tibble::tibble(
      analysis_scope = scope_name,
      pmid = pmid, year = yr, A = a_term, C = c_term,
      n_sent = nrow(S),
      pos = sum(S$sign == 1, na.rm = TRUE),
      neg = sum(S$sign == -1, na.rm = TRUE),
      null = sum(S$sign == 0, na.rm = TRUE),
      score = sc,
      label = legacy_label,
      wording_label = wording_label
    )
  }

  list(
    sent = if (length(out_sent)) dplyr::bind_rows(out_sent) else NULL,
    art = if (length(out_art)) dplyr::bind_rows(out_art) else NULL
  )
}

# -------------------------
# 5) Run extraction
# -------------------------
sent_res <- list()
art_res <- list()
for (scope_name in names(scopes)) {
  terms <- unique(scopes[[scope_name]])
  if (!length(terms)) next
  terms_with_cohit <- pair_diag |>
    dplyr::filter(analysis_scope == scope_name, n_cohit > 0) |>
    dplyr::pull(A) |>
    unique()
  terms <- intersect(terms, terms_with_cohit)
  log_msg("06 scanning scope=", scope_name, " terms_with_cohit=", length(terms))
  for (a in terms) {
    for (c in C_terms) {
      if (!any(pair_diag$analysis_scope == scope_name & pair_diag$A == a & pair_diag$C == c & pair_diag$n_cohit > 0)) next
      z <- analyze_pair(scope_name, a, c)
      if (!is.null(z$sent)) sent_res[[length(sent_res) + 1L]] <- z$sent
      if (!is.null(z$art)) art_res[[length(art_res) + 1L]] <- z$art
    }
  }
}

signed_sent <- if (length(sent_res)) dplyr::bind_rows(sent_res) else {
  tibble::tibble(
    analysis_scope = character(), pmid = character(), year = integer(), A = character(), C = character(), sentence = character(),
    sign = double(), weight = double(), inc = integer(), dec = integer(), null = integer(),
    measure = character(), value = double(), ci_low = double(), ci_high = double(), p_sig = integer(),
    auc = double(), cutoff = double(), unit = character(), sens = double(), spec = double()
  )
}

signed_art <- if (length(art_res)) dplyr::bind_rows(art_res) else {
  tibble::tibble(
    analysis_scope = character(), pmid = character(), year = integer(), A = character(), C = character(),
    n_sent = integer(), pos = integer(), neg = integer(), null = integer(), score = double(), label = character(), wording_label = character()
  )
}

if (!any(signed_art$analysis_scope == "primary_atlas")) {
  log_msg("06 WARNING: no primary_atlas signed-effect article rows were extracted. Pair diagnostics exist at: ", f_pair_diag)
}

# -------------------------
# 6) Summary with Wilson CI
# -------------------------
wilson_balance_ci <- function(pos, neg) {
  n <- pos + neg
  if (is.na(n) || n <= 0) return(c(low = NA_real_, high = NA_real_))
  p <- pos / n
  z <- 1.96
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  lo <- center - half
  hi <- center + half
  c(low = 2 * lo - 1, high = 2 * hi - 1)
}

sum_tab <- if (nrow(signed_art)) {
  signed_art |>
    dplyr::group_by(analysis_scope, A, C) |>
    dplyr::summarise(
      articles = dplyr::n(),
      pos_articles = sum(label == "risk_up", na.rm = TRUE),
      neg_articles = sum(label == "risk_down", na.rm = TRUE),
      null_or_mix = sum(label == "no_effect_or_mixed", na.rm = TRUE),
      net_score = sum(score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pos_ratio = pos_articles / pmax(1, pos_articles + neg_articles),
      balance = (pos_articles - neg_articles) / pmax(1, pos_articles + neg_articles)
    )
} else {
  tibble::tibble(
    analysis_scope = character(), A = character(), C = character(), articles = integer(),
    pos_articles = integer(), neg_articles = integer(), null_or_mix = integer(),
    net_score = double(), pos_ratio = double(), balance = double()
  )
}

if (nrow(sum_tab)) {
  ci_mat <- t(mapply(wilson_balance_ci, sum_tab$pos_articles, sum_tab$neg_articles))
  sum_tab$balance_low <- ci_mat[, "low"]
  sum_tab$balance_high <- ci_mat[, "high"]
} else {
  sum_tab$balance_low <- double()
  sum_tab$balance_high <- double()
}

sum_tab <- sum_tab |>
  add_term_metadata(dic, term_col = "A", prefix = "A") |>
  add_term_metadata(dic, term_col = "C", prefix = "C") |>
  dplyr::mutate(
    interpretation_layer = dplyr::case_when(
      analysis_scope == "primary_atlas" ~ "primary disease-state atlas wording signal",
      analysis_scope == "treatment_context" ~ "treatment-context wording signal; not treatment effect or drug risk",
      analysis_scope == "infection_context" ~ "infection-context wording signal; not causal infection attribution",
      TRUE ~ analysis_scope
    )
  ) |>
  dplyr::arrange(analysis_scope, C, dplyr::desc(articles), dplyr::desc(abs(balance)))

# Add metadata to sentence and article rows too.
signed_sent <- signed_sent |>
  add_term_metadata(dic, term_col = "A", prefix = "A") |>
  add_term_metadata(dic, term_col = "C", prefix = "C")
signed_art <- signed_art |>
  add_term_metadata(dic, term_col = "A", prefix = "A") |>
  add_term_metadata(dic, term_col = "C", prefix = "C")

# -------------------------
# 7) Save outputs
# -------------------------
f_sent <- file.path(DIR_TABLE, sprintf("signed_effects_sentence_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_art  <- file.path(DIR_TABLE, sprintf("signed_effects_article_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_sum  <- file.path(DIR_TABLE, sprintf("signed_effects_summary_withCI_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_sum_primary <- file.path(DIR_TABLE, sprintf("signed_effects_summary_primary_atlas_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_sum_treat <- file.path(DIR_TABLE, sprintf("signed_effects_summary_treatment_context_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_sum_inf <- file.path(DIR_TABLE, sprintf("signed_effects_summary_infection_context_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_num <- file.path(DIR_TABLE, sprintf("signed_effects_numeric_%s__%s.csv", CORPUS_TAG, DIC_TAG))
f_biom <- file.path(DIR_TABLE, sprintf("biomarker_metrics_%s__%s.csv", CORPUS_TAG, DIC_TAG))

write_csv2(signed_sent, f_sent)
write_csv2(signed_art, f_art)
write_csv2(sum_tab, f_sum)
write_csv2(sum_tab |> dplyr::filter(analysis_scope == "primary_atlas"), f_sum_primary)
write_csv2(sum_tab |> dplyr::filter(analysis_scope == "treatment_context"), f_sum_treat)
write_csv2(sum_tab |> dplyr::filter(analysis_scope == "infection_context"), f_sum_inf)

num_cols <- c("analysis_scope", "pmid", "year", "A", "C", "measure", "value", "ci_low", "ci_high", "p_sig", "auc", "cutoff", "unit", "sens", "spec", "sentence")
write_csv2(signed_sent |> dplyr::select(dplyr::any_of(num_cols)), f_num)

biomarker_terms <- dic |>
  dplyr::filter(domain == "biomarker" | class == "biomarker" | grepl("biomarker", abc_subrole, ignore.case = TRUE)) |>
  dplyr::pull(term) |>
  unique()
BM <- signed_sent |>
  dplyr::filter(A %in% biomarker_terms) |>
  dplyr::mutate(marker = A) |>
  dplyr::select(dplyr::any_of(c("marker", "analysis_scope", "pmid", "year", "C", "auc", "cutoff", "unit", "sens", "spec", "measure", "value", "ci_low", "ci_high", "sentence"))) |>
  dplyr::arrange(dplyr::desc(auc), dplyr::desc(sens), dplyr::desc(spec))
write_csv2(BM, f_biom)

run_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  dic_tag = DIC_TAG,
  outcomes = paste(C_terms, collapse = ";"),
  primary_terms_scanned = length(primary_terms),
  treatment_terms_scanned = length(treatment_terms),
  infection_terms_scanned = length(infection_terms),
  sentence_rows = nrow(signed_sent),
  article_rows = nrow(signed_art),
  summary_pairs = nrow(sum_tab),
  primary_summary_pairs = sum(sum_tab$analysis_scope == "primary_atlas", na.rm = TRUE),
  treatment_summary_pairs = sum(sum_tab$analysis_scope == "treatment_context", na.rm = TRUE),
  infection_summary_pairs = sum(sum_tab$analysis_scope == "infection_context", na.rm = TRUE),
  regex_case_policy = "dictionary-scoped; external ignore.case=FALSE for dictionary terms",
  interpretation = "Signed-effect labels describe abstract wording only. Treatment and infection-context signals are separated from primary disease-state claims."
)
f_run <- file.path(DIR_TABLE, sprintf("signed_effects_run_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG))
write_csv2(run_summary, f_run)

log_msg("06 signed rows sentence=", nrow(signed_sent), " article=", nrow(signed_art), " summary=", nrow(sum_tab))
log_msg("=== DONE ", RUN_NAME, " ===")
