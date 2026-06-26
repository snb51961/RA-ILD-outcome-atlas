# =========================================================
# 01b_FetchPubMed_ByFrozenPMIDs_AtlasV2.R
# Public-release helper for RA-ILD Atlas V2.
#
# Purpose:
#   Fetch PubMed title/abstract metadata for the frozen PMID list included in
#   data/corpus/RAILD_AtlasV2_frozen_pmids.csv, without redistributing PubMed
#   text in the repository. The output is written to data_proc/ in the format
#   expected by downstream scripts.
#
# Usage:
#   setwd("/path/to/RAILD_AtlasV2_repository_code_data/R")
#   Sys.setenv(RAILD_ROOT = normalizePath(".."))
#   source("01b_FetchPubMed_ByFrozenPMIDs_AtlasV2.R")
# =========================================================

source(file.path(getwd(), "00_setup_AtlasV2.R"))
assert_atlas_v2_setup(require_dictionary = FALSE)

quiet_install(c("rentrez", "xml2", "tibble", "dplyr", "readr", "stringr"))

STAMP_DATE <- stamp_date()
PMID_FILE <- Sys.getenv(
  "RAILD_FROZEN_PMID_FILE",
  unset = file.path(ROOT, "data", "corpus", "RAILD_AtlasV2_frozen_pmids.csv")
)
if (!file.exists(PMID_FILE)) stop("Frozen PMID file not found: ", PMID_FILE)

pmid_tbl <- readr::read_csv(PMID_FILE, show_col_types = FALSE)
if (!"pmid" %in% names(pmid_tbl)) stop("Frozen PMID file must contain a pmid column: ", PMID_FILE)
pmids <- unique(as.character(pmid_tbl$pmid))
pmids <- pmids[!is.na(pmids) & nzchar(pmids)]
if (!length(pmids)) stop("No PMIDs found in frozen PMID file: ", PMID_FILE)

log_msg("Fetching PubMed records for frozen PMID list:", length(pmids), "PMIDs")

extract_one <- function(n) {
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

BATCH <- as.integer(Sys.getenv("RAILD_ENTREZ_BATCH", unset = "200"))
rows <- list()
for (i in seq(1, length(pmids), by = BATCH)) {
  j <- min(i + BATCH - 1L, length(pmids))
  xmltxt <- rentrez::entrez_fetch(db = "pubmed", id = pmids[i:j], rettype = "xml", parsed = FALSE)
  doc <- xml2::read_xml(xmltxt)
  nodes <- xml2::xml_find_all(doc, ".//PubmedArticle")
  rows[[length(rows) + 1L]] <- dplyr::bind_rows(lapply(nodes, extract_one))
  log_msg(sprintf("Frozen PMID EFetch: %d / %d", j, length(pmids)))
  Sys.sleep(0.34)
}

articles <- dplyr::bind_rows(rows) |>
  dplyr::mutate(
    pmid = as.character(pmid),
    title = trimws(ifelse(is.na(title), "", title)),
    abstract = trimws(ifelse(is.na(abstract), "", abstract)),
    text = paste(title, abstract, sep = " "),
    corpus_role = "main_original"
  ) |>
  dplyr::right_join(tibble::tibble(pmid = pmids, frozen_order = seq_along(pmids)), by = "pmid") |>
  dplyr::arrange(frozen_order)

if ("publication_year" %in% names(pmid_tbl)) {
  articles <- articles |>
    dplyr::left_join(
      pmid_tbl |> dplyr::mutate(pmid = as.character(pmid)) |> dplyr::select(pmid, publication_year),
      by = "pmid"
    ) |>
    dplyr::mutate(year = dplyr::coalesce(year, suppressWarnings(as.integer(publication_year)))) |>
    dplyr::select(-publication_year)
}

missing_text <- articles |> dplyr::filter(is.na(abstract) | !nzchar(abstract))
if (nrow(missing_text)) {
  log_msg("WARNING: frozen PMIDs with missing/empty abstracts after EFetch:", nrow(missing_text))
}

out_file <- file.path(DIR_PROC, sprintf("articles_main_original_%s.csv", STAMP_DATE))
readr::write_csv(articles, out_file)
Sys.setenv(RAILD_ARTICLE_FILE = normalizePath(out_file, mustWork = FALSE))
log_msg("WROTE frozen PMID article text file:", out_file)
log_msg("Set RAILD_ARTICLE_FILE for current R session:", Sys.getenv("RAILD_ARTICLE_FILE"))
