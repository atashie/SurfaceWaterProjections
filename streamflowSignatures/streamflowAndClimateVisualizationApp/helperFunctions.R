# Helper Functions for Streamflow and Climate Visualization App
# Extended from streamflowVisualizationApp/helperFunctions.R

# === LEGEND HELPER ===

addLegend_decreasing_f <- function(map, position = c("topright", "bottomright", "bottomleft", "topleft"),
                                   pal, values, na.label = "NA", bins = 7, colors,
                                   opacity = 0.5, labels = NULL, labFormat = labelFormat(),
                                   title = NULL, className = "info legend", layerId = NULL,
                                   group = NULL, data = getMapData(map), decreasing = FALSE) {
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors)) {
      stop("You must provide either 'pal' or 'colors' (not both)")
    }
    if (missing(title) && inherits(values, "formula")) {
      title <- deparse(values[[2]])
    }
    values <- evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] == 0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins)) {
      warning("'bins' is ignored because the palette type is not numeric")
    }
    if (type == "numeric") {
      cuts <- if (length(bins) == 1) pretty(values, bins) else bins
      if (length(bins) > 2) {
        if (!all(abs(diff(bins, differences = 2)) <= sqrt(.Machine$double.eps))) {
          stop("The vector of breaks 'bins' must be equally spaced")
        }
      }
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1]) / (r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE) {
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      } else {
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
    } else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n]) / 2
      if (decreasing == TRUE) {
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      } else {
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
    } else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n]) / 2, na.rm = TRUE)
      if (decreasing == TRUE) {
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      } else {
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    } else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE) {
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      } else {
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    } else {
      stop("Palette function not supported")
    }
    if (!any(is.na(values))) {
      na.color <- NULL
    }
  } else {
    if (length(colors) != length(labels)) {
      stop("'colors' and 'labels' must be of the same length")
    }
  }
  legend <- list(
    colors = I(unname(colors)), labels = I(unname(labels)),
    na_color = na.color, na_label = na.label, opacity = opacity,
    position = position, type = type, title = title, extra = extra,
    layerId = layerId, className = className, group = group
  )
  invokeMethod(map, data, "addLegend", legend)
}

# === AWS CREDENTIALS ===

check_aws_credentials <- function() {
  required_vars <- c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_DEFAULT_REGION")
  missing_vars <- required_vars[!nzchar(Sys.getenv(required_vars))]

  if (length(missing_vars) > 0) {
    stop(
      "Missing required AWS environment variables: ",
      paste(missing_vars, collapse = ", ")
    )
  }
}

# === S3 DATA READING FUNCTIONS ===

# Function to read CSV directly from S3 using data.table::fread
# @param bucket S3 bucket name
# @param object_key Path to the CSV file in the bucket
# @param select Optional vector of column names to read (reduces memory)
read_csv_from_s3_direct <- function(bucket, object_key, select = NULL) {
  check_aws_credentials()
  require(aws.s3)
  require(data.table)

  # Get the object directly as a raw vector
  s3_object <- aws.s3::get_object(
    object = object_key,
    bucket = bucket,
    region = Sys.getenv("AWS_DEFAULT_REGION")
  )

  # Convert raw vector to character
  csv_text <- rawToChar(s3_object)

  # Read CSV using fread with text input and optional column selection
  csv_data <- if (!is.null(select)) {
    fread(text = csv_text, select = select)
  } else {
    fread(text = csv_text)
  }

  return(csv_data)
}

# Function to read .tif files from S3 with region specification
read_tif_from_s3 <- function(bucket, object_key) {
  check_aws_credentials()
  temp_file <- tempfile(fileext = ".tif")

  # Download the S3 object to the temporary file
  success <- aws.s3::save_object(
    object = object_key,
    bucket = bucket,
    file = temp_file,
    region = Sys.getenv("AWS_DEFAULT_REGION"),
    use_https = TRUE
  )

  # Read the raster from the temporary file
  raster_obj <- raster(temp_file)
  return(raster_obj)
}

# Function to read raster brick from S3
read_rasterbrick_from_s3 <- function(bucket, object_key) {
  check_aws_credentials()
  temp_file <- tempfile(fileext = ".tif")

  # Download the S3 object to the temporary file
  success <- aws.s3::save_object(
    object = object_key,
    bucket = bucket,
    file = temp_file,
    region = Sys.getenv("AWS_DEFAULT_REGION"),
    use_https = TRUE
  )

  # Read as raster brick instead of single raster
  brick_obj <- raster::brick(temp_file)
  return(brick_obj)
}

# === GENERIC PARQUET READER ===

#' Read parquet data from S3 for a specific ID value
#' Handles leading zero variations in a single query
#' Includes retry logic for S3 connection errors (SIGPIPE, broken pipe, etc.)
#'
#' @param id_value The ID to filter on
#' @param id_column Column name containing the ID
#' @param bucket S3 bucket name
#' @param parquet_key Path to the parquet file in the bucket
#' @param select_cols Vector of column names to select (reduces data transfer)
#' @param data_name Display name for messages
#' @param dataset Optional pre-opened Arrow dataset (for connection reuse)
#' @return data.table with filtered data
read_parquet_by_id <- function(id_value, id_column, bucket, parquet_key,
                                select_cols = NULL, data_name = "data",
                                dataset = NULL) {
  require(arrow)
  require(dplyr)
  require(data.table)

  # Helper function to perform the actual query
  perform_query <- function(ds) {
    ids_to_try <- unique(c(as.character(id_value), paste0("0", id_value)))
    query <- ds %>% filter(!!sym(id_column) %in% ids_to_try)
    if (!is.null(select_cols)) {
      query <- query %>% select(all_of(select_cols))
    }
    query %>% collect()
  }

  # Build S3 URL for fresh connections
  s3_parquet_url <- paste0("s3://", bucket, "/", parquet_key)

  # Try with provided dataset first, retry with fresh connection on error

  filtered_data <- tryCatch({
    if (is.null(dataset)) {
      message(paste("Reading", data_name, "from S3 for", id_column, ":", id_value))
      dataset <- open_dataset(s3_parquet_url)
    }
    perform_query(dataset)
  }, error = function(e) {
    # Check if this is a connection-related error
    if (grepl("SIGPIPE|connection|reset|broken pipe|ignoring", e$message, ignore.case = TRUE)) {
      message(paste("Connection error for", data_name, "- retrying with fresh dataset..."))
      tryCatch({
        fresh_dataset <- open_dataset(s3_parquet_url)
        perform_query(fresh_dataset)
      }, error = function(e2) {
        warning(paste("Retry failed for", data_name, ":", e2$message))
        return(data.table())  # Return empty data.table on failure
      })
    } else {
      warning(paste("Error reading", data_name, ":", e$message))
      return(data.table())  # Return empty data.table on non-connection errors
    }
  })

  if (nrow(filtered_data) == 0) {
    message(paste("No", data_name, "found for", id_column, ":", id_value))
  }

  return(as.data.table(filtered_data))
}

# Function to read a subset of parquet data from S3 for a specific gage_id
# Handles cases where leading zeros may be missing from the metadata
#' @param gage_id The gage ID to query
#' @param bucket S3 bucket name
#' @param parquet_key Path to the parquet file
#' @param dataset Optional pre-opened Arrow dataset (for connection reuse)
read_streamflow_by_gage_id <- function(gage_id,
                                       bucket = "climate-ai-data-science-shiny-app-data",
                                       parquet_key = "streamflow/combined_streamflow_data_09feb2026.parquet",
                                       dataset = NULL) {
  read_parquet_by_id(
    id_value = gage_id,
    id_column = "gage_id",
    bucket = bucket,
    parquet_key = parquet_key,
    select_cols = c("gage_id", "Date", "Q"),  # Only needed columns
    data_name = "streamflow data",
    dataset = dataset
  )
}

# === DAYMET CLIMATE DATA READER ===

#' Read Daymet climate data for a specific gage from S3 parquet
#' Includes retry logic for S3 connection errors (SIGPIPE, broken pipe, etc.)
#'
#' @param gage_id The gage ID (site_id in the parquet)
#' @param bucket S3 bucket name
#' @param parquet_key Path to the parquet file in the bucket
#' @param dataset Optional pre-opened Arrow dataset (for connection reuse)
#' @return data.table with columns: site_id, Date, prcp, tmin, tmax, swe, vp, srad
read_daymet_by_gage_id <- function(gage_id,
                                   bucket = "climate-ai-data-science-shiny-app-data",
                                   parquet_key = "streamflow/daymet_1980_2023.parquet",
                                   dataset = NULL) {
  require(arrow)
  require(dplyr)
  require(data.table)

  s3_parquet_url <- paste0("s3://", bucket, "/", parquet_key)

  # Helper function to perform query on a dataset
  perform_daymet_query <- function(ds) {
    col_names <- names(ds)
    date_col <- if ("Date" %in% col_names) "Date" else if ("date" %in% col_names) "date" else NULL

    if (is.null(date_col)) {
      message("Warning: No date column found in Daymet data")
      return(data.table())
    }

    needed_cols <- c("site_id", date_col, "prcp", "tmin", "tmax", "swe", "vp", "srad")
    available_cols <- intersect(needed_cols, col_names)
    ids_to_try <- unique(c(as.character(gage_id), paste0("0", gage_id)))

    result <- ds %>%
      filter(site_id %in% ids_to_try) %>%
      select(all_of(available_cols)) %>%
      collect() %>%
      as.data.table()

    if (nrow(result) > 0 && date_col == "date") {
      setnames(result, "date", "Date")
    }
    if (nrow(result) > 0 && "Date" %in% names(result)) {
      result[, Date := as.Date(Date)]
    }

    return(result)
  }

  # Try with provided dataset first, retry with fresh connection on error
  result <- tryCatch({
    if (is.null(dataset)) {
      message(paste("Reading Daymet climate data from S3 for site_id:", gage_id))
      dataset <- open_dataset(s3_parquet_url)
    }
    perform_daymet_query(dataset)
  }, error = function(e) {
    if (grepl("SIGPIPE|connection|reset|broken pipe|ignoring", e$message, ignore.case = TRUE)) {
      message(paste("Connection error for Daymet data - retrying with fresh dataset..."))
      tryCatch({
        fresh_dataset <- open_dataset(s3_parquet_url)
        perform_daymet_query(fresh_dataset)
      }, error = function(e2) {
        warning(paste("Retry failed for Daymet data:", e2$message))
        return(data.table())
      })
    } else {
      warning(paste("Error reading Daymet data:", e$message))
      return(data.table())
    }
  })

  if (nrow(result) == 0) {
    message(paste("No Daymet climate data found for site_id:", gage_id))
  }

  return(result)
}

# === TREND ANALYSIS FUNCTIONS ===

#' Calculate Theil-Sen slope for two vectors
#'
#' @param x Independent variable (numeric vector)
#' @param y Dependent variable (numeric vector)
#' @return List with slope and intercept
calculate_theil_sen <- function(x, y) {
  require(zyp)

  # Remove NA values
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) < 3) {
    return(list(slope = NA, intercept = NA))
  }

  tryCatch({
    fit <- zyp.sen(y ~ x)
    list(
      slope = fit$coefficients["x"],
      intercept = fit$coefficients["Intercept"]
    )
  }, error = function(e) {
    list(slope = NA, intercept = NA)
  })
}

#' Calculate Mann-Kendall test for a vector
#'
#' @param y Numeric vector (time series)
#' @return List with tau and p-value
calculate_mann_kendall <- function(y) {
  require(Kendall)

  # Remove NA values
  y <- y[!is.na(y)]

  if (length(y) < 4) {
    return(list(tau = NA, pval = NA))
  }

  tryCatch({
    mk_result <- MannKendall(y)
    list(
      tau = as.numeric(mk_result$tau),
      pval = as.numeric(mk_result$sl)
    )
  }, error = function(e) {
    list(tau = NA, pval = NA)
  })
}

#' Calculate Spearman correlation between two vectors
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @return List with rho and p-value
calculate_spearman <- function(x, y) {
  # Remove pairs where either is NA
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) < 4) {
    return(list(rho = NA, pval = NA))
  }

  tryCatch({
    result <- cor.test(x, y, method = "spearman", exact = FALSE)
    list(
      rho = as.numeric(result$estimate),
      pval = as.numeric(result$p.value)
    )
  }, error = function(e) {
    list(rho = NA, pval = NA)
  })
}

# === WEATHER VARIABLE METADATA ===

#' Get display information for weather variables
#'
#' @return Named list with variable info: unit, color, label
get_weather_var_info <- function() {
  list(
    prcp = list(
      unit = "mm",
      color = "rgba(139, 90, 43, 0.7)",  # Sienna brown
      label = "Precipitation"
    ),
    tmin = list(
      unit = "\u00B0C",
      color = "rgb(100, 149, 237)",  # Light blue
      label = "Min Temperature"
    ),
    tmax = list(
      unit = "\u00B0C",
      color = "rgb(220, 20, 60)",  # Crimson red
      label = "Max Temperature"
    ),
    temp = list(
      unit = "\u00B0C",
      color = "rgb(220, 20, 60)",  # Used for combined temp (actual traces use tmin/tmax colors)
      label = "Temperature"
    ),
    swe = list(
      unit = "mm",
      color = "rgb(169, 169, 169)",  # Dark gray
      label = "Snow Water Equivalent"
    ),
    vp = list(
      unit = "Pa",
      color = "rgb(60, 179, 113)",  # Medium sea green
      label = "Vapor Pressure"
    ),
    srad = list(
      unit = "W/m\u00B2",
      color = "rgb(255, 165, 0)",  # Orange
      label = "Solar Radiation"
    )
  )
}

# === FAST THEIL-SEN FOR BATCH USE ===

#' Lightweight Theil-Sen slope estimator for batch pairwise computations
#'
#' Avoids zyp.sen formula overhead. For n > 500, samples pairs to keep runtime
#' manageable while preserving slope estimate accuracy.
#'
#' @param x Numeric vector (independent variable)
#' @param y Numeric vector (dependent variable)
#' @return Numeric scalar: median pairwise slope, or NA if insufficient data
fast_theil_sen_slope <- function(x, y) {
  valid <- complete.cases(x, y)
  x <- x[valid]; y <- y[valid]
  n <- length(x)
  if (n < 3) return(NA_real_)
  if (n > 500) {
    set.seed(42)
    idx_i <- sample.int(n, 50000, replace = TRUE)
    idx_j <- sample.int(n, 50000, replace = TRUE)
    keep <- idx_i != idx_j & (x[idx_i] != x[idx_j])
    slopes <- (y[idx_i[keep]] - y[idx_j[keep]]) / (x[idx_i[keep]] - x[idx_j[keep]])
  } else {
    idx <- combn(n, 2)
    dx <- x[idx[2, ]] - x[idx[1, ]]
    nonzero <- dx != 0
    slopes <- (y[idx[2, nonzero]] - y[idx[1, nonzero]]) / dx[nonzero]
  }
  median(slopes, na.rm = TRUE)
}
