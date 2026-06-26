# =========================================================
# run_all_no01_AtlasV2.R
# RA-ILD Atlas v2.2.x end-to-end run script, excluding PubMed fetch (01).
#
# Purpose:
#   Run 00, 02, 02b, 03-06, 07, 08, 09, 10, 11 in order.
#   This assumes the frozen corpus and final dictionary already exist.
#
# Usage in R:
#   source("/path/to/RAILD_project/R/run_all_no01_AtlasV2.R")
# =========================================================

# ---- 0) User paths ----
# ---- 0) User paths ----
CODE_DIR <- Sys.getenv("RAILD_CODE_DIR", unset = "")
if (!nzchar(CODE_DIR)) CODE_DIR <- getwd()
CODE_DIR <- normalizePath(CODE_DIR, mustWork = FALSE)

ROOT_DIR <- Sys.getenv("RAILD_ROOT", unset = "")
if (!nzchar(ROOT_DIR)) {
  candidate_roots <- unique(c(CODE_DIR, dirname(CODE_DIR), getwd(), dirname(getwd())))
  has_project_dirs <- vapply(candidate_roots, function(p) {
    dir.exists(file.path(p, "dic")) ||
      dir.exists(file.path(p, "data_proc")) ||
      dir.exists(file.path(p, "data_raw"))
  }, logical(1))
  ROOT_DIR <- if (any(has_project_dirs)) candidate_roots[which(has_project_dirs)[1]] else dirname(CODE_DIR)
}
ROOT_DIR <- normalizePath(ROOT_DIR, mustWork = FALSE)

setwd(CODE_DIR)
Sys.setenv(RAILD_CODE_DIR = CODE_DIR)
Sys.setenv(RAILD_ROOT = ROOT_DIR)

message("=========================================================")
message("RA-ILD Atlas v2.2.x run: 00, 02, 02b, 03-06, 07-11")
message("CODE_DIR: ", CODE_DIR)
message("RAILD_ROOT: ", ROOT_DIR)
message("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("=========================================================")

# ---- 1) Pre-flight checks ----
must_exist <- function(path, label = path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", label, "\nExpected at: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

must_exist(file.path(CODE_DIR, "00_setup_AtlasV2.R"), "00_setup_AtlasV2.R")
must_exist(file.path(CODE_DIR, "02_BuildHitsMatrix_AtlasV2.R"), "02_BuildHitsMatrix_AtlasV2.R")

must_exist(
  file.path(ROOT_DIR, "dic", "ra_ild_dictionary_analysis_v2_outcome_preserving_atlas.csv"),
  "final v2 atlas dictionary"
)

must_exist(
  file.path(ROOT_DIR, "data_proc", "articles_main_original_20260323.csv"),
  "frozen main original article corpus"
)

yc_files <- list.files(
  file.path(ROOT_DIR, "data_proc"),
  pattern = "^year_counts_stacked_by_pubclass_.*\\.csv$",
  full.names = TRUE
)
if (length(yc_files) == 0) {
  warning("No year_counts_stacked_by_pubclass_*.csv found in data_proc. Figure 1 panel C may be incomplete.")
}

# ---- 2) Helper for robust step execution ----
source_step <- function(script, label = script, required = TRUE) {
  path <- file.path(CODE_DIR, script)
  if (!file.exists(path)) {
    msg <- paste0("Script not found: ", script)
    if (required) stop(msg, call. = FALSE)
    warning(msg, immediate. = TRUE)
    return(invisible(FALSE))
  }

  message("\n---------------------------------------------------------")
  message("START: ", label)
  message("SCRIPT: ", script)
  message("TIME:  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  message("---------------------------------------------------------")

  t0 <- Sys.time()
  tryCatch(
    {
      source(path, local = globalenv())
      dt <- difftime(Sys.time(), t0, units = "mins")
      message("DONE:  ", label, "  [", round(as.numeric(dt), 2), " min]")
      gc()
      invisible(TRUE)
    },
    error = function(e) {
      message("\n*** ERROR in step: ", label, " ***")
      message(conditionMessage(e))
      message("Stopped at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      stop(e)
    }
  )
}

# ---- 3) Choose the latest Figure 1 builder available ----
fig1_candidates <- c(
  "02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v6_concept_nodes.R",
  "02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v5_noNA.R",
  "02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v4_infection_context.R",
  "02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v3_clinical.R",
  "02b_Figure1_Workflow_DictionaryCoverage_AtlasV2_v2_readable.R",
  "02b_Figure1_Workflow_DictionaryCoverage_AtlasV2.R"
)
FIG1_SCRIPT <- fig1_candidates[file.exists(file.path(CODE_DIR, fig1_candidates))][1]
if (is.na(FIG1_SCRIPT) || !nzchar(FIG1_SCRIPT)) {
  stop("No 02b Figure 1 builder script found in CODE_DIR.", call. = FALSE)
}
message("Using Figure 1 builder: ", FIG1_SCRIPT)

# ---- 4) Run pipeline excluding 01_FetchPubMed ----

# 00 + 02 + 02b
source_step("00_setup_AtlasV2.R", "00 setup")
source_step("02_BuildHitsMatrix_AtlasV2.R", "02 build hit matrix and dictionary coverage")
source_step(FIG1_SCRIPT, "02b Figure 1 dictionary/workflow/corpus figure")

# 03-06 core analytic tables
source_step("03_CoocAndCollocation_AtlasV2.R", "03 co-occurrence and collocation")
source_step("04_ABC_Rankings_AtlasV2.R", "04 ABC rankings")
source_step("05_AC_NPMI_AtlasV2.R", "05 outcome coherence / NPMI")
source_step("06_SignedEffects_AtlasV2.R", "06 signed-effect extraction")

# 07-11 manuscript figures
source_step("07_Build_OutcomeDomainMatrix_AtlasV2.R", "07 Figure 2 outcome-domain atlas")
source_step("08_Figure3_AEILD_BridgeAtlas_AtlasV2.R", "08 Figure 3 AE-ILD bridge atlas")
source_step("09_Figure4_ProgressionMortality_BridgeAtlas_AtlasV2.R", "09 Figure 4 progression/mortality bridge atlas")
source_step("10_Figure5_BiomarkerMolecularFreeAtlas_AtlasV2.R", "10 Figure 5 biomarker/molecular free-layout atlas")
source_step("11_Figure6_IntegratedSignalSummary_AtlasV2.R", "11 Figure 6 integrated signal summary")

message("\n=========================================================")
message("ALL REQUESTED STEPS COMPLETED")
message("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Output root should be under:")
message(file.path(ROOT_DIR, "output", "pm_1980_20251231", "analysis_v2_outcome_preserving_atlas"))
message("=========================================================")
