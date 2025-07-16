
addLegend_decreasing_f <- function (map, position = c("topright", "bottomright", "bottomleft","topleft"),
                                    pal, values, na.label = "NA", bins = 7, colors, 
                                    opacity = 0.5, labels = NULL, labFormat = labelFormat(), 
                                    title = NULL, className = "info legend", layerId = NULL, 
                                    group = NULL, data = getMapData(map), decreasing = FALSE) {
  
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors)) 
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(values, "formula")) 
      title <- deparse(values[[2]])
    values <- evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] == 
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins)) 
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1) 
        pretty(values, bins)
      else bins   
      if (length(bins) > 2) 
        if (!all(abs(diff(bins, differences = 2)) <= 
                 sqrt(.Machine$double.eps))) 
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n])/2, na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(values))) 
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels)) 
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)), 
                 na_color = na.color, na_label = na.label, opacity = opacity, 
                 position = position, type = type, title = title, extra = extra, 
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}


# Function to read CSV directly from S3 using data.table::fread
read_csv_from_s3_direct <- function(bucket, object_key) {
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
  
  # Read CSV using fread with text input
  csv_data <- fread(text = csv_text)  # Note the 'text=' parameter
  
  return(csv_data)
}


check_aws_credentials <- function() {
  required_vars <- c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_DEFAULT_REGION")
  missing_vars <- required_vars[!nzchar(Sys.getenv(required_vars))]
  
  if (length(missing_vars) > 0) {
    stop("Missing required AWS environment variables: ", 
         paste(missing_vars, collapse = ", "))
  }
}


# Function to read .tif files from S3 with region specification
read_tif_from_s3 <- function(bucket, object_key) {
  # Create a temporary file
  check_aws_credentials()
  temp_file <- tempfile(fileext = ".tif")
  
  # Download the S3 object to the temporary file
  # Specify use_https = TRUE and region
  success <- aws.s3::save_object(
    object = object_key,
    bucket = bucket,
    file = temp_file,
    region = Sys.getenv("AWS_DEFAULT_REGION"),
    use_https = TRUE
  )
  
  #    if (!success) {
  #      stop("Failed to download the file from S3.")
  #    }
  
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



# Function to read a subset of parquet data from S3 for a specific gage_id
# Handles cases where leading zeros may be missing from the metadata
read_streamflow_by_gage_id <- function(gage_id,
                                       bucket = "climate-ai-data-science-shiny-app-data",
                                       parquet_key = "streamflow/combined_streamflow_data.parquet") {
  
  # Load required libraries
  require(arrow)
  require(dplyr)
  
  # Construct S3 URL for the parquet file
  s3_parquet_url <- paste0("s3://", bucket, "/", parquet_key)
  
  message(paste("Reading streamflow data from S3 for gage_id:", gage_id))
  
  # Open the Parquet file as a dataset directly from S3
  dataset <- open_dataset(s3_parquet_url)
  
  # First attempt: try with the gage_id as provided
  filtered_data <- dataset %>%
    filter(gage_id == !!gage_id) %>%
    collect()  # Collect the filtered data into memory as a data frame
  
  # Check if any data was returned
  if (nrow(filtered_data) == 0) {
    # Try adding a leading zero
    gage_id_with_zero <- paste0("0", gage_id)
    message(paste("No data found for gage_id:", gage_id))
    message(paste("Trying with leading zero:", gage_id_with_zero))
    
    # Second attempt: try with leading zero
    filtered_data <- dataset %>%
      filter(gage_id == !!gage_id_with_zero) %>%
      collect()
    
    # If still no data, provide informative error
    if (nrow(filtered_data) == 0) {
      warning(paste("No data found for gage_id:", gage_id, "or", gage_id_with_zero))
    } else {
      message(paste("Success! Found data using gage_id:", gage_id_with_zero))
    }
  }
  
  return(filtered_data)
}
