#!/usr/bin/env Rscript
#' Export Canadian HYDAT station metadata to CSV for cross-language use.
#'
#' Extracts RHBN designation and REGULATED status from the tidyhydat database
#' and writes a CSV that Julia/Python can read directly.
#'
#' Output: metadata/canadian_hydat_interference.csv
#'   Columns: STATION_NUMBER, RHBN, REGULATED

library(tidyhydat)
library(data.table)

cat("Exporting HYDAT interference metadata...\n")

# Get station information (includes RHBN flag)
stations <- as.data.table(hy_stations())
cat(sprintf("  Loaded %d stations from hy_stations()\n", nrow(stations)))

# Get regulation information
regulation <- tryCatch(
  as.data.table(hy_stn_regulation()),
  error = function(e) {
    cat("  Warning: Failed to load hy_stn_regulation():", e$message, "\n")
    data.table(STATION_NUMBER = character(), REGULATED = logical())
  }
)
cat(sprintf("  Loaded %d regulation records\n", nrow(regulation)))

# Keep most recent regulation status per station (any TRUE = regulated)
if (nrow(regulation) > 0 && "REGULATED" %in% names(regulation)) {
  regulation <- regulation[, .(REGULATED = any(REGULATED, na.rm = TRUE)),
                           by = STATION_NUMBER]
}

# Select STATION_NUMBER and RHBN from stations
station_cols <- intersect(c("STATION_NUMBER", "RHBN"), names(stations))
stations <- stations[, ..station_cols]

# Merge with regulation status
if (nrow(regulation) > 0) {
  stations <- merge(stations, regulation, by = "STATION_NUMBER", all.x = TRUE)
} else {
  stations[, REGULATED := NA]
}

# Ensure RHBN column exists
if (!"RHBN" %in% names(stations)) {
  stations[, RHBN := NA]
}

cat(sprintf("  Final: %d stations, RHBN non-NA: %d, REGULATED non-NA: %d\n",
            nrow(stations),
            sum(!is.na(stations$RHBN)),
            sum(!is.na(stations$REGULATED))))
cat(sprintf("  RHBN TRUE: %d, FALSE: %d\n",
            sum(stations$RHBN == TRUE, na.rm = TRUE),
            sum(stations$RHBN == FALSE, na.rm = TRUE)))
cat(sprintf("  REGULATED TRUE: %d, FALSE: %d\n",
            sum(stations$REGULATED == TRUE, na.rm = TRUE),
            sum(stations$REGULATED == FALSE, na.rm = TRUE)))

# Write CSV — use working directory (expected to be project root)
output_path <- file.path(getwd(), "metadata", "canadian_hydat_interference.csv")

fwrite(stations, output_path)
cat(sprintf("  Written to: %s\n", output_path))
cat("Done.\n")
