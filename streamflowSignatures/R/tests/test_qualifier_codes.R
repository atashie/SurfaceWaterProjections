# test_qualifier_codes.R — unit tests for the USGS daily-value acceptance rule
# (R/usgs_qualifiers.R; decision 1d of the 2026-10-21 qualifier-mask rerun plan).
#
# Usage (from the repository root; needs base R only):
#   Rscript R/tests/test_qualifier_codes.R
# Exit status 1 when any check fails.
#
# The toy frame covers every qualification code of the census (§1 of the plan) with and
# without a value, the NWIS no-data sentinel, non-numeric values, infinities, codes without
# a status letter and the Canadian-style bare symbols that must never reach this rule with
# a kept status (they are masked because their status letter is not A or P, or kept when it
# happens to be "A" — the rule is for US rows only, as in generate_streamflow_dt).
# ==============================================================================

src <- if (file.exists("R/usgs_qualifiers.R")) "R/usgs_qualifiers.R" else "../usgs_qualifiers.R"
source(src)

cat("=== USGS qualifier-code acceptance rule tests ===\n\n")
n_pass <- 0L
n_fail <- 0L
check <- function(label, cond) {
  if (isTRUE(cond)) {
    n_pass <<- n_pass + 1L
    cat("  PASS:", label, "\n")
  } else {
    n_fail <<- n_fail + 1L
    cat("  FAIL:", label, "\n")
  }
}

# --- 1. Every code of census §1 with a finite value: approved and provisional codes are
#        kept whatever the remark; the status letter decides --------------------------
approved_remarks <- c("A", "A e", "A R", "A <", "A >", "A [4]", "A R e", "A e R")
provisional_remarks <- c("P", "P e", "P Ice", "P Dis", "P Rat", "P ***", "P Bkw", "P Eqp",
                         "P Fld", "P [4]", "P <", "P >", "P ZFL", "P Ssn", "P Mnt")
codes <- c(approved_remarks, provisional_remarks)
vals <- seq_along(codes) * 1.5
cat("1. Status letter decides, remarks never mask\n")
check("every approved-status code with a finite value is kept",
      all(usgs_qualifier_keep(approved_remarks, vals[seq_along(approved_remarks)])))
check("every provisional-status code with a finite value is kept",
      all(usgs_qualifier_keep(provisional_remarks, vals[length(approved_remarks) + seq_along(provisional_remarks)])))
check("kept values are returned unchanged by usgs_qualifier_mask",
      identical(usgs_qualifier_mask(codes, vals), vals))
check("the four approved codes the old rule masked (A R, A <, A >, A [4]) are now kept",
      all(usgs_qualifier_keep(c("A R", "A <", "A >", "A [4]"), c(10, 0.5, 100, 3))))
check("zero flow under P ZFL is a value, kept",
      usgs_qualifier_keep("P ZFL", 0))
check("the old rule's four exact strings are still kept",
      all(usgs_qualifier_keep(c("A", "A e", "P", "P e"), c(1, 2, 3, 4))))

# --- 2. No usable value: NA, the sentinel, non-numeric, infinite ----------------------
cat("2. A code never manufactures a value\n")
check("NA value under an approved code is not kept",
      !usgs_qualifier_keep("A", NA_real_))
check("the NWIS no-data sentinel -999999 is not kept under A",
      !usgs_qualifier_keep("A", -999999))
check("the sentinel is not kept under A R either",
      !usgs_qualifier_keep("A R", -999999))
check("a different sentinel can be configured",
      !usgs_qualifier_keep("A", -1, no_data_value = -1) && usgs_qualifier_keep("A", -999999, no_data_value = -1))
check("a non-numeric value string (\"Ice\") is not kept",
      !usgs_qualifier_keep("P Ice", "Ice"))
check("Inf and -Inf are not kept",
      !any(usgs_qualifier_keep(c("A", "A"), c(Inf, -Inf))))
check("NaN is not kept",
      !usgs_qualifier_keep("A", NaN))
check("masked entries are NA_real_ in the output",
      identical(usgs_qualifier_mask(c("A", "A", "A"), c(NA, -999999, 7)), c(NA_real_, NA_real_, 7)))
check("a numeric string value is parsed (\"12.5\" under A kept as 12.5)",
      identical(usgs_qualifier_mask("A", "12.5"), 12.5))

# --- 3. Status letter missing or not A / P ---------------------------------------------
cat("3. The status letter must be A or P\n")
check("NA code with a value is not kept (decision 4: no-code days are left alone)",
      !usgs_qualifier_keep(NA_character_, 5))
check("empty code with a value is not kept",
      !usgs_qualifier_keep("", 5))
check("a bare remark without a status letter (\"e\") is not kept",
      !usgs_qualifier_keep("e", 5))
check("lower-case status letters are not accepted (NWIS serves upper case)",
      !any(usgs_qualifier_keep(c("a", "p e"), c(1, 2))))
check("HYDAT-style symbols B and E are not kept by this rule (US rows only)",
      !any(usgs_qualifier_keep(c("B", "E"), c(1, 2))))
check("a status letter followed by ':'-separated remarks is recognised (\"A:e\", \"P:Ice\")",
      all(usgs_qualifier_keep(c("A:e", "P:Ice"), c(1, 2))))
check("a ','-separated code is recognised (\"A,e\")", usgs_qualifier_keep("A,e", 1))
check("the sentinel served as a string (\"-999999\") is not kept", !usgs_qualifier_keep("A", "-999999"))
check("leading / trailing whitespace around the code is ignored",
      all(usgs_qualifier_keep(c(" A ", "P e  "), c(1, 2))))
check("usgs_qualifier_status returns the first token",
      identical(usgs_qualifier_status(c("A R", "P", NA, "", " A e")), c("A", "P", "", "", "A")))

# --- 4. Vector contract -----------------------------------------------------------------
cat("4. Vector contract\n")
check("length mismatch is an error",
      inherits(tryCatch(usgs_qualifier_keep(c("A", "A"), 1), error = function(e) e), "error"))
check("zero-length input gives zero-length output",
      length(usgs_qualifier_keep(character(0), numeric(0))) == 0)
check("negative values are kept by the rule (rejection of negative Q is the preprocessor's decision)",
      usgs_qualifier_keep("A", -0.5))

# --- 5. The audit tally -------------------------------------------------------------------
cat("5. Per-code tally\n")
tf <- c("A", "A", "A R", "A R", "P Ice", "P Ice", NA, "e")
tv <- c(1, NA, 2, -999999, 3, NA, NA, 9)
tl <- usgs_qualifier_tally(tf, tv)
row <- function(code) tl[tl$flag == code, ]
check("tally has one row per distinct code incl. <NA>",
      nrow(tl) == 5 && all(c("A", "A R", "P Ice", "<NA>", "e") %in% tl$flag))
check("A: 2 days, 1 kept, 0 masked with value, 1 without value",
      with(row("A"), days == 2 && kept == 1 && masked_with_value == 0 && no_value == 1))
check("A R: the sentinel counts as no value",
      with(row("A R"), days == 2 && kept == 1 && masked_with_value == 0 && no_value == 1))
check("e: a finite value masked by the rule is counted as masked_with_value",
      with(row("e"), days == 1 && kept == 0 && masked_with_value == 1 && no_value == 0))
check("<NA> code: no value",
      with(row("<NA>"), days == 1 && kept == 0 && no_value == 1))
check("tally days sum to the input length", sum(tl$days) == length(tf))

# --- 6. Old rule vs new rule on a mixed frame: exactly the approved-remark days change ---
cat("6. Old rule vs new rule\n")
old_rule <- function(flag, value) {
  v <- as.numeric(value)
  v[!(flag %in% c("A", "A e", "P", "P e"))] <- NA_real_
  v
}
mf <- data.frame(
  flag = c("A", "A e", "P", "P e", "A R", "A <", "A >", "A [4]", "P Ice", "P Rat", NA, "A"),
  value = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, NA, NA),
  stringsAsFactors = FALSE
)
q_old <- old_rule(mf$flag, mf$value)
q_new <- usgs_qualifier_mask(mf$flag, mf$value)
changed <- which(!(is.na(q_old) & is.na(q_new)) & !(!is.na(q_old) & !is.na(q_new) & q_old == q_new))
check("the old rule's kept days are unchanged",
      identical(q_new[1:4], q_old[1:4]))
check("the approved-remark days and the provisional-remark days with a value are restored",
      identical(mf$flag[changed], c("A R", "A <", "A >", "A [4]", "P Ice", "P Rat")))
check("days without a value stay NA under both rules",
      all(is.na(q_new[11:12])) && all(is.na(q_old[11:12])))

# --- 7. Integration: generate_streamflow_dt on a stubbed NWIS frame ---------------------------
# Needs lubridate (the function calls lubridate::year/month/yday); skipped where it is absent
# (this Mac) and run on the ingestion machine. The function is taken from R/helperFunctions.R
# without sourcing that file (its package requirements are not installed here).
cat("7. Integration (generate_streamflow_dt)\n")
hf <- if (file.exists("R/helperFunctions.R")) "R/helperFunctions.R" else "../helperFunctions.R"
if (requireNamespace("lubridate", quietly = TRUE) && file.exists(hf)) {
  exprs <- parse(hf)
  is_def <- vapply(exprs, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("generate_streamflow_dt")), logical(1))
  env <- new.env()
  env$usgs_qualifier_mask <- usgs_qualifier_mask
  env$log_warn <- function(...) invisible(NULL)
  env$calculate_water_year_info <- function(d) { m <- as.integer(format(d, "%m")); y <- as.integer(format(d, "%Y"))
    wy <- ifelse(m >= 10, y + 1L, y); list(water_year = wy, dowy = as.integer(d - as.Date(paste0(wy - 1, "-10-01")) + 1L)) }
  env$as.data.table <- function(x) x
  dates <- seq(as.Date("1979-10-01"), by = "day", length.out = 365 * 21)
  codes <- rep(c("A", "A e", "P", "P e", "A R", "A <", "A >", "A [4]", "P Ice", NA), length.out = length(dates))
  vals <- seq_along(dates) * 0.5; vals[codes %in% "P Ice"] <- NA; vals[is.na(codes)] <- NA; vals[5] <- -999999
  env$readNWISdv <- function(siteNumber, parameterCd, startDate, endDate) data.frame(
    agency_cd = "USGS", site_no = siteNumber, Date = dates, X_00060_00003 = vals, X_00060_00003_cd = codes, stringsAsFactors = FALSE)
  eval(exprs[[which(is_def)]], env)
  res <- env$generate_streamflow_dt(data.frame(STAID = "01234567", DRAIN_SQKM = 250, stringsAsFactors = FALSE), "USGS",
                                    min_num_years = 20, start_date = as.Date("1979-10-01"), end_date = as.Date("2025-10-01"))
  conv <- 60 * 60 * 24 / (250 * 3280.84^3) * 1e6
  check("integration: the flag column is the served code string, verbatim, for every row", identical(res$flag, codes))
  check("integration: Q = served value x conversion on every kept day",
        all(abs(res$Q[!is.na(res$Q)] - vals[!is.na(res$Q)] * conv) <= 1e-12 * abs(res$Q[!is.na(res$Q)])))
  check("integration: the approved-remark days (A R, A <, A >, A [4]) carry values", all(!is.na(res$Q[codes %in% c("A R", "A <", "A >", "A [4]")])))
  check("integration: NA codes, NA values and the sentinel give NA Q", all(is.na(res$Q[is.na(codes) | is.na(vals) | vals == -999999])))
  check("integration: the row count and dates are the served ones", nrow(res) == length(dates) && identical(res$Date, dates))
} else {
  cat("  SKIP: lubridate not installed (run on the ingestion machine)\n")
}

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
