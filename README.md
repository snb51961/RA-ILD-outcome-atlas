# RA-ILD Outcome Atlas repository: public data and R code

This repository contains the public processed datasets and R analysis scripts associated with the manuscript:

**Outcome-preserving literature-derived atlas of rheumatoid arthritis-associated interstitial lung disease (RA-ILD).**

The repository is intended to support reproducibility of the published analyses while complying with PubMed redistribution policies.

---

## Scope of this repository

This repository accompanies, but does not replace, the journal Supplementary Materials.

Human-readable Supplementary Methods, Supplementary Figures, and Supplementary Tables submitted with the manuscript are not duplicated here.

Instead, this repository provides the frozen manuscript corpus identifiers, the final curated dictionary, machine-readable provenance and robustness summaries, and the R scripts required to reproduce the published analyses after retrieving PubMed metadata for the frozen PMID list.

---

## Repository contents

### `data/`

Public processed datasets used in the manuscript.

- `corpus/`
  - Frozen PMID list defining the manuscript corpus.
- `dictionary/`
  - Final curated dictionary.
  - Dictionary layer summary.
  - Dictionary quality-control summary.
- `supporting_summaries/`
  - Figure-level provenance summaries.
  - Figure selection rules.
  - Context-layer quality-control summaries.
  - Internal temporal and resampling robustness summaries.

### `dic/`

Runtime-compatible copy of the final dictionary using the filename expected by the R analysis scripts.

### `R/`

Complete R analysis pipeline, including

- PubMed retrieval scripts
- Dictionary mapping
- Atlas construction
- Figure generation
- Provenance analysis
- Robustness analyses
- Public rerun scripts

### `data_proc/`

Documentation describing intermediate processing files used during manuscript development.

### `data_raw/`

Documentation explaining why raw PubMed text is not redistributed.

### `dataset_manifest.csv`

Manifest describing all publicly released files included in this repository.

---

## Not included

The following materials are intentionally **not** redistributed.

- PubMed title and abstract text.
- Intermediate local processing files.
- Local hit matrices.
- Generated figure image files.
- Journal Supplementary Methods, Supplementary Figures, and Supplementary Tables submitted alongside the manuscript.

---

## Reproducibility

The frozen PMID list together with the final curated dictionary defines the manuscript corpus and concept mapping framework.

Because PubMed title/abstract text is not redistributed, a public rerun first retrieves PubMed records for the included frozen PMID list (or alternatively uses a locally supplied article file).

Recommended public rerun:

```r
setwd("/path/to/RA-ILD-outcome-atlas/R")
Sys.setenv(RAILD_ROOT = normalizePath(".."))
source("01b_FetchPubMed_ByFrozenPMIDs_AtlasV2.R")
source("run_all_no01_public_AtlasV2.R")
```

Alternatively,

```r
setwd("/path/to/RA-ILD-outcome-atlas/R")
Sys.setenv(RAILD_ROOT = normalizePath(".."))
Sys.setenv(RAILD_ARTICLE_FILE="/path/to/articles_main_original_YYYYMMDD.csv")
source("run_all_no01_public_AtlasV2.R")
```

See `R/README_code.md` for package requirements and execution details.

---

## License

Please select an appropriate repository license during GitHub/Zenodo deposition according to institutional and journal requirements.

---

## Citation

If you use this repository, please cite both the associated manuscript and the archived Zenodo record.

**Zenodo DOI:** *(to be added after deposition)*