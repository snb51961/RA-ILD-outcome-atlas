# =========================================================
# 01_FetchPubMed_AtlasV2.R
# RA-ILD atlas v2 pipeline
#
# Purpose
#   - Retrieve the fixed PubMed RA-ILD corpus through 2025-12-31.
#   - Classify records after retrieval into main original articles,
#     reviews/guidelines, case reports and excluded editorial-like records.
#   - Save the main original corpus for downstream dictionary-based mapping.
#   - Save annual corpus-composition counts used by Figure 1.
#
# Dictionary status
#   - This script is dictionary-independent. It sources 00_setup_AtlasV2.R
#     only for shared paths, corpus tags and logging.
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = FALSE)

quiet_install(c(
  "rentrez", "xml2", "tibble", "dplyr", "readr", "ggplot2",
  "changepoint", "systemfonts", "showtext", "stringr", "tidyr", "scales"
))

cands  <- c("Hiragino Sans", "Hiragino Kaku Gothic ProN", "Yu Gothic", "IPAexGothic",
            "Noto Sans CJK JP", "Arial", "Helvetica", "DejaVu Sans")
avail  <- intersect(cands, unique(systemfonts::system_fonts()$family))
base_font <- if (length(avail)) avail[1] else "sans"
showtext::showtext_auto()
cairo_pdf_device <- grDevices::cairo_pdf

STAMP_DATE <- stamp_date()

# ------------------ 1) PubMed search -----------------------
QUERY <- paste0(
  '(',
    '"Arthritis, Rheumatoid"[MH] OR "rheumatoid arthritis"[Title/Abstract]',
  ') AND (',
    '"interstitial lung disease"[Title/Abstract] OR ILD[Title/Abstract] OR ',
    '"interstitial pneumonia"[Title/Abstract] OR "pulmonary fibrosis"[Title/Abstract] OR ',
    '"diffuse parenchymal lung disease"[Title/Abstract]',
  ') AND hasabstract[text]',
  ' AND ("', YEAR_MIN, '/01/01"[pdat] : "', DATE_MAX, '"[pdat])'
)

log_msg("Searching PubMed for fixed RA-ILD corpus ...")
s0 <- rentrez::entrez_search(db = "pubmed", term = QUERY, use_history = TRUE, retmax = 0)
TOTAL <- as.integer(s0$count)
log_msg("PubMed total abstract-bearing records:", TOTAL)

pmids <- character(0)
step <- 5000L
if (TOTAL > 0) {
  for (start in seq(0L, TOTAL - 1L, by = step)) {
    ids <- rentrez::entrez_search(
      db = "pubmed", term = QUERY, retstart = start, retmax = step
    )$ids
    pmids <- c(pmids, ids)
    log_msg(sprintf("PMID retrieval: %d / %d", length(pmids), TOTAL))
    Sys.sleep(0.34)
  }
}
pmids <- unique(pmids)

readr::write_csv(
  tibble::tibble(pmid = pmids),
  file.path(DIR_RAW, sprintf("pmids_%s.csv", CORPUS_TAG))
)

# ------------------ 2) EFetch / XML parsing ----------------
extract_one <- function(n){
  z   <- function(xpath) xml2::xml_text(xml2::xml_find_first(n, xpath))
  zx  <- function(xpath) paste(xml2::xml_text(xml2::xml_find_all(n, xpath)), collapse = " ")
  zxs <- function(xpath) paste(xml2::xml_text(xml2::xml_find_all(n, xpath)), collapse = ";")

  yr <- z(".//PubDate/Year")
  if (is.na(yr) || yr == "") yr <- z(".//ArticleDate/Year")

  tibble::tibble(
    pmid             = z(".//PMID"),
    title            = z(".//ArticleTitle"),
    abstract         = zx(".//Abstract/AbstractText"),
    year             = suppressWarnings(as.integer(yr)),
    pubtype          = zxs(".//PublicationType"),
    mesh             = zxs(".//MeshHeading/DescriptorName"),
    journal          = z(".//Journal/Title"),
    affiliation_last = z(".//AuthorList/Author[last()]/AffiliationInfo/Affiliation")
  )
}

BATCH <- 200L
rows <- list()
if (length(pmids) > 0) {
  for (i in seq(1, length(pmids), by = BATCH)) {
    j <- min(i + BATCH - 1L, length(pmids))
    xmltxt <- rentrez::entrez_fetch(db = "pubmed", id = pmids[i:j], rettype = "xml", parsed = FALSE)
    doc    <- xml2::read_xml(xmltxt)
    nodes  <- xml2::xml_find_all(doc, ".//PubmedArticle")
    rows[[length(rows) + 1L]] <- dplyr::bind_rows(lapply(nodes, extract_one))
    log_msg(sprintf("EFetch: %d / %d", j, length(pmids)))
    Sys.sleep(0.34)
  }
}

articles_raw <- dplyr::bind_rows(rows)

# ------------------ 3) Corpus classification ----------------
collapse_na <- function(x) ifelse(is.na(x), "", x)

is_valid_abstract <- function(x, min_nchar = 80L){
  txt <- trimws(collapse_na(x))
  txt != "" && nchar(txt) >= min_nchar
}

classify_record <- function(pubtype_string, title_string, abstract_string){
  pt <- tolower(collapse_na(pubtype_string))
  ti <- tolower(collapse_na(title_string))
  ab <- tolower(collapse_na(abstract_string))
  txt <- paste(ti, ab)

  if (stringr::str_detect(pt, "editorial|comment|news|interview|historical article|biography|published erratum")) {
    return("exclude_editorial_like")
  }
  if (stringr::str_detect(ti, "^comment on\\b|^reply to\\b|^response to\\b|^editorial\\b")) {
    return("exclude_editorial_like")
  }

  if (stringr::str_detect(pt, "guideline|practice guideline|consensus development conference|consensus")) {
    return("review_guideline")
  }
  if (stringr::str_detect(pt, "systematic review|meta-analysis|review")) {
    return("review_guideline")
  }
  review_title_re <- paste(
    c("\\bsystematic review\\b", "\\bmeta-analysis\\b", "\\bnarrative review\\b",
      "^review\\b", "\\ba review\\b", "\\bguideline\\b", "\\bconsensus\\b", "\\bdelphi\\b",
      "appraisal of .* guideline", "recommendation", "position statement"),
    collapse = "|"
  )
  if (stringr::str_detect(ti, review_title_re)) {
    return("review_guideline")
  }

  if (stringr::str_detect(pt, "case reports?")) {
    return("case_report")
  }

  if (stringr::str_detect(
      pt,
      paste(
        c("clinical trial", "randomized controlled trial", "controlled clinical trial",
          "observational study", "comparative study", "multicenter study",
          "validation study", "evaluation study"),
        collapse = "|"
      )
    )) {
    return("original_or_other")
  }

  if (stringr::str_detect(
      txt,
      paste(
        c("we aimed", "we investigated", "we evaluated", "we analyzed", "we analysed", "we assessed",
          "retrospective", "prospective", "cohort", "case-control", "cross-sectional",
          "multicentre", "multicenter", "patients were", "subjects were",
          "odds ratio", "hazard ratio", "confidence interval", "sensitivity", "specificity",
          "receiver operating characteristic", "auc"),
        collapse = "|"
      )
    )) {
    if (!stringr::str_detect(ti, review_title_re)) return("original_or_other")
  }

  if (stringr::str_detect(pt, "journal article") &&
      !stringr::str_detect(pt, "review|guideline|consensus|editorial|comment|news|interview|historical article|biography|published erratum|case reports?")) {
    return("original_or_other")
  }

  if (stringr::str_detect(pt, "letter")) {
    return("exclude_editorial_like")
  }

  return("original_or_other")
}

articles <- articles_raw |>
  dplyr::filter(!is.na(year), year >= YEAR_MIN, year <= YEAR_MAX) |>
  dplyr::mutate(
    pmid = as.character(pmid),
    title = trimws(collapse_na(title)),
    abstract = trimws(collapse_na(abstract)),
    text = paste(title, abstract, sep = " "),
    pubclass = mapply(classify_record, pubtype, title, abstract, USE.NAMES = FALSE),
    has_valid_abstract = vapply(abstract, is_valid_abstract, logical(1)),
    corpus_role = dplyr::case_when(
      pubclass == "original_or_other" & has_valid_abstract ~ "main_original",
      pubclass == "review_guideline"  & has_valid_abstract ~ "review_guideline",
      pubclass == "case_report"       & has_valid_abstract ~ "case_report",
      TRUE ~ "excluded"
    )
  )

if (nrow(articles) && any(articles$year > YEAR_MAX, na.rm = TRUE)) {
  stop("ERROR: year > ", YEAR_MAX, " remains after filtering.")
}

articles_main <- articles |> dplyr::filter(corpus_role == "main_original")
articles_review <- articles |> dplyr::filter(corpus_role == "review_guideline")
articles_case <- articles |> dplyr::filter(corpus_role == "case_report")
articles_excluded <- articles |> dplyr::filter(corpus_role == "excluded")

# ------------------ 4) Save outputs ------------------------
saveRDS(articles_main, file.path(DIR_RAW,  sprintf("articles_%s.rds", STAMP_DATE)))
readr::write_csv(articles_main, file.path(DIR_PROC, sprintf("articles_%s.csv", STAMP_DATE)))

readr::write_csv(articles,          file.path(DIR_PROC, sprintf("articles_all_%s.csv", STAMP_DATE)))
readr::write_csv(articles_main,     file.path(DIR_PROC, sprintf("articles_main_original_%s.csv", STAMP_DATE)))
readr::write_csv(articles_review,   file.path(DIR_PROC, sprintf("articles_review_guideline_%s.csv", STAMP_DATE)))
readr::write_csv(articles_case,     file.path(DIR_PROC, sprintf("articles_case_report_%s.csv", STAMP_DATE)))
readr::write_csv(articles_excluded, file.path(DIR_PROC, sprintf("articles_excluded_%s.csv", STAMP_DATE)))

pubclass_counts <- articles |>
  dplyr::count(pubclass, corpus_role, name = "n") |>
  dplyr::arrange(dplyr::desc(n))
readr::write_csv(pubclass_counts, file.path(DIR_PROC, sprintf("pubclass_counts_%s.csv", STAMP_DATE)))

corpus_summary <- tibble::tibble(
  corpus_tag = CORPUS_TAG,
  date_max = DATE_MAX,
  pubmed_records = TOTAL,
  all_records_after_year_filter = nrow(articles),
  main_original = nrow(articles_main),
  review_guideline = nrow(articles_review),
  case_report = nrow(articles_case),
  excluded = nrow(articles_excluded)
)
readr::write_csv(corpus_summary, file.path(DIR_PROC, sprintf("corpus_summary_%s.csv", CORPUS_TAG)))
readr::write_csv(corpus_summary, file.path(DIR_TABLE, sprintf("corpus_summary_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

log_msg("Saved all records n=", nrow(articles))
log_msg("Saved main_original n=", nrow(articles_main))
log_msg("Saved review_guideline n=", nrow(articles_review))
log_msg("Saved case_report n=", nrow(articles_case))
log_msg("Saved excluded n=", nrow(articles_excluded))

# ------------------ 5) Annual counts and Figure 1 inputs ----------
year_counts <- articles_main |>
  dplyr::count(year, name = "n") |>
  dplyr::arrange(year)
readr::write_csv(year_counts, file.path(DIR_PROC, sprintf("year_counts_%s.csv", STAMP_DATE)))

p_year <- ggplot2::ggplot(year_counts, ggplot2::aes(year, n)) +
  ggplot2::geom_col() +
  ggplot2::labs(title = "Annual publication counts: main original RA-ILD corpus",
                x = "Year", y = "Count") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(text = ggplot2::element_text(family = base_font))

ggplot2::ggsave(file.path(DIR_FIG, "year_counts.pdf"), p_year,
                device = cairo_pdf_device, width = 8, height = 5)

cpt <- tryCatch(
  changepoint::cpt.meanvar(year_counts$n, method = "PELT", penalty = "BIC"),
  error = function(e) NULL
)
if (!is.null(cpt) && nrow(year_counts) > 0) {
  idx <- changepoint::cpts(cpt)
  change_years <- year_counts$year[idx]
  p_year_cp <- p_year +
    ggplot2::geom_vline(xintercept = change_years, linetype = "dashed", color = "red")
  ggplot2::ggsave(file.path(DIR_FIG, "year_counts_changepoint.pdf"), p_year_cp,
                  device = cairo_pdf_device, width = 8, height = 5)
  log_msg("Changepoints:", paste(change_years, collapse = ", "))
}

plot_df <- articles |>
  dplyr::filter(corpus_role %in% c("main_original", "review_guideline", "case_report")) |>
  dplyr::filter(!is.na(year)) |>
  dplyr::mutate(
    corpus_role = factor(corpus_role, levels = c("main_original", "case_report", "review_guideline"))
  )

year_type_counts <- plot_df |>
  dplyr::count(year, corpus_role, name = "n") |>
  tidyr::complete(year, corpus_role, fill = list(n = 0)) |>
  dplyr::arrange(year, corpus_role)

f_year_stacked <- file.path(DIR_PROC, sprintf("year_counts_stacked_by_pubclass_%s.csv", STAMP_DATE))
readr::write_csv(year_type_counts, f_year_stacked)
readr::write_csv(year_type_counts, file.path(DIR_TABLE, sprintf("year_counts_stacked_by_pubclass_%s__%s.csv", CORPUS_TAG, DIC_TAG)))

p_stacked <- ggplot2::ggplot(
  year_type_counts,
  ggplot2::aes(x = year, y = n, fill = corpus_role)
) +
  ggplot2::geom_col(width = 0.85, position = ggplot2::position_stack(reverse = TRUE)) +
  ggplot2::labs(
    title = "Annual publication counts by corpus type",
    subtitle = "Fixed RA-ILD PubMed corpus",
    x = "Year", y = "Count", fill = "Corpus type"
  ) +
  ggplot2::scale_fill_manual(
    values = c("main_original" = "#4E79A7", "case_report" = "#59A14F", "review_guideline" = "#F28E2B"),
    breaks = c("main_original", "case_report", "review_guideline"),
    labels = c("main_original" = "Original", "case_report" = "Case report", "review_guideline" = "Review/Guideline")
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(reverse = FALSE)) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(text = ggplot2::element_text(family = base_font), panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(file.path(DIR_FIG, "year_counts_stacked_by_pubclass.pdf"), p_stacked,
                device = cairo_pdf_device, width = 9, height = 5.5)

year_type_pct <- year_type_counts |>
  dplyr::group_by(year) |>
  dplyr::mutate(pct = ifelse(sum(n) > 0, n / sum(n), 0)) |>
  dplyr::ungroup()
readr::write_csv(year_type_pct, file.path(DIR_PROC, sprintf("year_counts_stacked_by_pubclass_percent_%s.csv", STAMP_DATE)))

log_msg("WROTE Figure 1 corpus-composition input:", f_year_stacked)
log_msg("=== DONE 01_FetchPubMed_AtlasV2 ===")
