#' Load GAGES-II human interference metadata
#'
#' Reads GAGES-II text files (CONUS + AKHIPR) and returns a data.table of
#' dam, land cover, and classification variables keyed by STAID.
#'
#' @param gages_dir Path to directory containing GAGES-II files.
#'   Defaults to config value.
#' @return A data.table with STAID and interference columns, or NULL if
#'   directory not found.
#' @export
#' @importFrom data.table fread rbindlist
load_gages_ii_interference <- function(gages_dir = NULL) {
  if (is.null(gages_dir)) gages_dir <- pkg_env$metadata_gages_ii_dir
  if (!dir.exists(gages_dir)) {
    warning("GAGES-II directory not found: ", gages_dir)
    return(NULL)
  }

  conus_files <- list(
    hydromod_dams  = "conterm_hydromod_dams.txt",
    pop_infrastr   = "conterm_pop_infrastr.txt",
    hydromod_other = "conterm_hydromod_other.txt",
    bas_classif    = "conterm_bas_classif.txt",
    lc06_basin     = "conterm_lc06_basin.txt"
  )

  akhipr_files <- list(
    hydromod_dams  = "AKHIPR_hydromod_dams.txt",
    pop_infrastr   = "AKHIPR_pop_infrastr.txt",
    hydromod_other = "AKHIPR_hydromod_other.txt",
    bas_classif    = "AKHIPR_bas_classif.txt"
  )

  col_map <- list(
    hydromod_dams  = c("STAID", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009"),
    pop_infrastr   = c("STAID", "IMPNLCD06"),
    hydromod_other = c("STAID", "FRESHW_WITHDRAWAL"),
    bas_classif    = c("STAID", "CLASS", "HYDRO_DISTURB_INDX"),
    lc06_basin     = c("STAID", "DEVNLCD06")
  )

  read_gages_file <- function(filepath, cols) {
    if (!file.exists(filepath)) return(NULL)
    tryCatch({
      dt <- data.table::fread(filepath, sep = ",", header = TRUE,
                              colClasses = list(character = "STAID"),
                              encoding = "Latin-1")
      present <- intersect(cols, names(dt))
      if (!"STAID" %in% present) return(NULL)
      dt[, present, with = FALSE]
    }, error = function(e) NULL)
  }

  merge_files <- function(file_list, dir_path) {
    result <- NULL
    for (type in names(file_list)) {
      fp <- file.path(dir_path, file_list[[type]])
      dt <- read_gages_file(fp, col_map[[type]])
      if (!is.null(dt)) {
        if (is.null(result)) {
          result <- dt
        } else {
          result <- merge(result, dt, by = "STAID", all = TRUE)
        }
      }
    }
    result
  }

  conus_data  <- merge_files(conus_files, gages_dir)
  akhipr_data <- merge_files(akhipr_files, gages_dir)

  if (!is.null(conus_data) && !is.null(akhipr_data)) {
    # Do NOT intersect the CONUS and AKHIPR column sets. The AKHIPR files carry a
    # narrower schema — none of them has IMPNLCD06, DEVNLCD06, FRESHW_WITHDRAWAL
    # or HYDRO_DISTURB_INDX, and there is no AKHIPR_lc06_basin.txt at all — so
    # intersecting dropped all four columns for the ~9,000 CONUS gages that DO
    # have them, silently costing the output 4 metadata columns against canonical
    # Julia (caught 2026-08-26 by the strict schema gate). rbindlist(fill = TRUE)
    # already reconciles differing schemas by filling NA, which is the correct
    # result: those attributes genuinely do not exist for AK/HI/PR gages.
    combined <- data.table::rbindlist(list(conus_data, akhipr_data), fill = TRUE)
  } else if (!is.null(conus_data)) {
    combined <- conus_data
  } else if (!is.null(akhipr_data)) {
    combined <- akhipr_data
  } else {
    return(NULL)
  }

  # Replace -999 sentinel with NA
  num_cols <- names(combined)[vapply(combined, is.numeric, logical(1))]
  for (col in num_cols) {
    combined[get(col) == -999, (col) := NA]
  }

  # Classify
  combined[, human_interference_class := "unknown"]
  if ("CLASS" %in% names(combined)) {
    combined[CLASS == "Ref",     human_interference_class := "reference"]
    combined[CLASS == "Non-ref", human_interference_class := "non-reference"]
  }

  combined
}
