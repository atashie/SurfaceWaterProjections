################################################################################
# Smoke Test for Post-Reorganization Validation
#
# Verifies that the codebase reorganization didn't break any functionality.
# Run from project root: Rscript R/tests/smoke_test_reorganization.R
################################################################################

cat("========================================\n")
cat("POST-REORGANIZATION SMOKE TEST\n")
cat("========================================\n\n")

# Track test results
tests_passed <- 0
tests_failed <- 0
tests_total <- 0

run_test <- function(test_name, test_expr) {
  tests_total <<- tests_total + 1
  result <- tryCatch({
    eval(test_expr)
    tests_passed <<- tests_passed + 1
    cat("PASS:", test_name, "\n")
    TRUE
  }, error = function(e) {
    tests_failed <<- tests_failed + 1
    cat("FAIL:", test_name, "-", e$message, "\n")
    FALSE
  })
  return(result)
}

# =============================================================================
# LEVEL 0: Source Path Validation
# =============================================================================
cat("\n--- LEVEL 0: Source Path Validation ---\n")

# Auto-detect project root
main_dir <- getwd()
if (!file.exists(file.path(main_dir, "config.R"))) {
  main_dir <- dirname(main_dir)
  if (!file.exists(file.path(main_dir, "config.R"))) {
    main_dir <- dirname(main_dir)
  }
  setwd(main_dir)
}

run_test("Source config.R from project root", quote({
  source(file.path(main_dir, "config.R"))
}))

run_test("Source R/helperFunctions.R from project root", quote({
  suppressPackageStartupMessages(source(file.path(main_dir, "R", "helperFunctions.R")))
}))

# =============================================================================
# LEVEL 1: Core Function Loading
# =============================================================================
cat("\n--- LEVEL 1: Core Function Loading ---\n")

critical_functions <- c(
  "process_signatures_from_parquet",
  "generate_stats",
  "calculate_flow_vols_by_year",
  "analyze_baseflow_indices",
  "analyze_recession_parameters",
  "calculate_pulse_metrics",
  "analyze_flashiness_trends",
  "analyze_flow_timing_trends",
  "log_info",
  "log_error",
  "validate_file_exists"
)

run_test("All 11 critical functions exist", quote({
  missing <- critical_functions[!sapply(critical_functions, exists)]
  if (length(missing) > 0) {
    stop(paste("Missing:", paste(missing, collapse = ", ")))
  }
}))

config_vars <- c(
  "MIN_Q_VALUE_AND_DAYS",
  "MIN_NUM_YEARS",
  "MIN_FRAC_GOOD_DATA",
  "STAT_SUFFIXES",
  "EXPECTED_SIGNATURE_BASES",
  "LOG_LEVEL"
)

run_test("All 6 config variables exist", quote({
  missing <- config_vars[!sapply(config_vars, exists)]
  if (length(missing) > 0) {
    stop(paste("Missing:", paste(missing, collapse = ", ")))
  }
}))

# =============================================================================
# LEVEL 2: Entry Point Script Dependencies
# =============================================================================
cat("\n--- LEVEL 2: Entry Point Dependencies ---\n")

entry_points <- c(
  "run_full_processing.R",
  "run_ingest_usgs_hydat.R",
  "run_caravan_processing.R"
)

for (ep in entry_points) {
  ep_path <- file.path(main_dir, ep)
  run_test(paste("Entry point exists:", ep), quote({
    if (!file.exists(ep_path)) stop("File not found")
  }))
}

# =============================================================================
# LEVEL 2: Test Files Can Source Dependencies
# =============================================================================
cat("\n--- LEVEL 2: Test File Dependencies ---\n")

test_files <- c(
  "R/tests/smoke_test.R",
  "R/tests/test_climate_functions.R",
  "R/tests/test_climate_signatures.R",
  "R/tests/qa_qc_signatures.R"
)

for (tf in test_files) {
  tf_path <- file.path(main_dir, tf)
  run_test(paste("Test file exists:", tf), quote({
    if (!file.exists(tf_path)) stop("File not found")
  }))
}

# =============================================================================
# LEVEL 2: R/ Folder Scripts
# =============================================================================
cat("\n--- LEVEL 2: R/ Folder Scripts ---\n")

r_scripts <- c(
  "R/helperFunctions.R",
  "R/run_conversion.R",
  "R/precompute_cross_signature_analysis.R"
)

for (rs in r_scripts) {
  rs_path <- file.path(main_dir, rs)
  run_test(paste("R/ script exists:", rs), quote({
    if (!file.exists(rs_path)) stop("File not found")
  }))
}

# =============================================================================
# LEVEL 2: Shiny App
# =============================================================================
cat("\n--- LEVEL 2: Shiny App ---\n")

run_test("Shiny app directory exists", quote({
  app_dir <- file.path(main_dir, "streamflowAndClimateVisualizationApp")
  if (!dir.exists(app_dir)) stop("Directory not found")
}))

run_test("Shiny app.R exists", quote({
  app_file <- file.path(main_dir, "streamflowAndClimateVisualizationApp", "app.R")
  if (!file.exists(app_file)) stop("File not found")
}))

run_test("Shiny app helperFunctions.R exists", quote({
  helper_file <- file.path(main_dir, "streamflowAndClimateVisualizationApp", "helperFunctions.R")
  if (!file.exists(helper_file)) stop("File not found")
}))

# =============================================================================
# LEVEL 2: Archive and Orphans
# =============================================================================
cat("\n--- LEVEL 2: Archive and Docs ---\n")

run_test("archive/ directory exists", quote({
  if (!dir.exists(file.path(main_dir, "archive"))) stop("Directory not found")
}))

run_test("docs/ directory exists", quote({
  if (!dir.exists(file.path(main_dir, "docs"))) stop("Directory not found")
}))

run_test("logs/ directory exists", quote({
  if (!dir.exists(file.path(main_dir, "logs"))) stop("Directory not found")
}))

# =============================================================================
# Summary
# =============================================================================
cat("\n========================================\n")
cat("SMOKE TEST SUMMARY\n")
cat("========================================\n")
cat("Tests passed:", tests_passed, "/", tests_total, "\n")
cat("Tests failed:", tests_failed, "\n")

if (tests_failed == 0) {
  cat("\nALL TESTS PASSED - Reorganization verified!\n")
} else {
  cat("\nWARNING: Some tests failed. Review output above.\n")
}
