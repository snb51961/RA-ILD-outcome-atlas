# R code notes

This directory contains the R analysis scripts associated with the RA-ILD Atlas V2 manuscript.

## Important limitations

PubMed title/abstract text is not redistributed in this repository. The public data package includes a frozen PMID list, dictionary, and selected compact supporting summaries, but not the full article text corpus, all intermediate hit matrices, generated figures, or human-readable journal Supplementary Tables S1-S7. Therefore, a complete rerun requires either:

1. fetching PubMed records for the included frozen PMID list with `01b_FetchPubMed_ByFrozenPMIDs_AtlasV2.R` **recommended for public reproducibility**,
2. re-running the original query-based fetch script `01_FetchPubMed_AtlasV2.R`, or
3. providing a local article text CSV via `RAILD_ARTICLE_FILE`.

## Recommended public run using the frozen PMID list

From R, run from this `R/` directory:

```r
setwd("/path/to/RAILD_AtlasV2_repository_code_data/R")
Sys.setenv(RAILD_ROOT = normalizePath(".."))

# Fetch PubMed title/abstract metadata for the frozen PMID list. Requires internet access.
source("01b_FetchPubMed_ByFrozenPMIDs_AtlasV2.R")

# Run the analysis scripts excluding PubMed fetch.
source("run_all_no01_public_AtlasV2.R")
```

If using a local article text file instead:

```r
setwd("/path/to/RAILD_AtlasV2_repository_code_data/R")
Sys.setenv(RAILD_ROOT = normalizePath(".."))
Sys.setenv(RAILD_ARTICLE_FILE = "/path/to/articles_main_original_YYYYMMDD.csv")
source("run_all_no01_public_AtlasV2.R")
```

`run_all_no01_AtlasV2.R` is retained as the internal analysis runner and includes a manuscript-era hard-coded pre-flight check for `articles_main_original_20260323.csv`. The public helper `run_all_no01_public_AtlasV2.R` avoids that hard-coded date and searches for a local article file.

## Package dependencies

See `PACKAGE_DEPENDENCIES.txt`. The major non-base packages used include readr, dplyr, stringr, tibble, purrr, lubridate, rentrez, xml2, ggplot2, quanteda, quanteda.textstats, igraph, patchwork, ggrepel, changepoint, systemfonts, showtext, tidyr, and scales.
