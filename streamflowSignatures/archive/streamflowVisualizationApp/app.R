# Load required libraries
library(shiny)
library(shinythemes)
library(leaflet)
library(plotly)
library(dplyr)
library(data.table)
library(arrow)
library(lubridate)
library(sf)
library(shinycssloaders) # For loading spinners

# Source helper functions
source("helperFunctions.R")

# Set S3 bucket name
s3_bucket_name <- "climate-ai-data-science-shiny-app-data"

# === AUTHENTICATION AND DATA LOADING ON APP START ===
# Check AWS credentials once
readRenviron(".Renviron")
check_aws_credentials()

# Read and filter metadata once
message("Loading metadata from S3...")
metadata_raw <- read_csv_from_s3_direct(bucket = s3_bucket_name, 
                                        object_key = "streamflow/combined_watershed_metadata.csv")
goodGages <- subset(metadata_raw, processing_status == 'success')
message(paste("Loaded", nrow(goodGages), "gages with successful processing"))

# Load watershed boundaries
message("Loading watershed boundaries...")
watersheds <- st_read("unified_watershedBoundaries_simplified.gpkg", quiet = TRUE)
# Transform to WGS84 for leaflet
watersheds <- st_transform(watersheds, 4326)

# Extract first basin ID from basin_ids column for matching
watersheds$first_basin_id <- sapply(watersheds$basin_ids, function(x) {
  if (is.na(x)) return(NA)
  strsplit(as.character(x), ";")[[1]][1]
})

message(paste("Loaded", nrow(watersheds), "watershed boundaries"))

# Check matching
n_matched <- sum(goodGages$Downstream_HB_ID %in% watersheds$first_basin_id)
message(paste("Matched", n_matched, "out of", nrow(goodGages), "gages to watershed boundaries"))

# Load streamflow signature summary data
message("Loading streamflow signature summary data from S3...")
signature_data <- read_csv_from_s3_direct(
  bucket = s3_bucket_name,
  object_key = "streamflow/streamflowSignature_summaryData_OCT2025.csv"
)
message(paste("Loaded signature data for", nrow(signature_data), "gages"))

# Extract metric names from columns ending in "_mean"
mean_columns <- grep("_mean$", names(signature_data), value = TRUE)
metric_names <- gsub("_mean$", "", mean_columns)
message(paste("Found", length(metric_names), "metrics:", paste(head(metric_names, 3), collapse = ", "), "..."))

# Define UI
ui <- fluidPage(
  theme = shinytheme("superhero"),
  
  titlePanel("Streamflow Gage Dashboard"),
  
  # Selection input
  fluidRow(
    column(12,
           selectizeInput(
             inputId = "gage_selector",
             label = "Select Gage ID:",
             choices = goodGages$gage_id,
             selected = goodGages$gage_id[1],
             options = list(
               placeholder = "Type or select a gage ID",
               searchField = "label"
             )
           )
    )
  ),
  
  # Map without spinner
  fluidRow(
    column(12,
           leafletOutput("gage_map", height = "500px")
    )
  ),
  
  # Time series plot with scale toggle
  fluidRow(
    column(12,
           h3("Streamflow Time Series"),
           radioButtons("scale_type", 
                        label = "Y-Axis Scale:", 
                        choices = list("Linear" = "linear", "Logarithmic" = "log"),
                        selected = "log",
                        inline = TRUE),
           withSpinner(
             plotlyOutput("streamflow_plot", height = "400px"),
             type = 4,
             color = "#3498db"
           )
    )
  ),
  
  # New section for metric visualization
  fluidRow(
    column(6,
           h3("Streamflow Signature Metrics"),
           selectInput(
             inputId = "metric_selector",
             label = "Select Metric to Plot:",
             choices = metric_names,
             selected = metric_names[1]
           )
    ),
    column(6,
           h3("Point Size Control"),
           sliderInput(
             inputId = "point_size",
             label = "Marker Radius:",
             min = 1,
             max = 5,
             value = 2,
             step = 0.5
           )
    )
  ),
  
  # First row of metric maps (Mean and Median) - clockwise from top left
  fluidRow(
    column(6,
           h4("Mean"),
           withSpinner(
             leafletOutput("metric_map_mean", height = "400px"),
             type = 4,
             color = "#3498db"
           )
    ),
    column(6,
           h4("Median"),
           withSpinner(
             leafletOutput("metric_map_median", height = "400px"),
             type = 4,
             color = "#3498db"
           )
    )
  ),
  
  # Second row of metric maps (Slope and P-value)
  fluidRow(
    column(6,
           h4("Slope (Trend)"),
           withSpinner(
             leafletOutput("metric_map_slp", height = "400px"),
             type = 4,
             color = "#3498db"
           )
    ),
    column(6,
           h4("P-value"),
           withSpinner(
             leafletOutput("metric_map_pval", height = "400px"),
             type = 4,
             color = "#3498db"
           )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Create the base leaflet map
  output$gage_map <- renderLeaflet({
    # Create color palette for num_years
    pal <- colorNumeric(
      palette = c("red", "blue"),
      domain = goodGages$num_years
    )
    
    # Create hover text with all metadata
    hover_text <- apply(goodGages, 1, function(row) {
      info_items <- paste0("<b>", names(row), ":</b> ", row)
      paste(info_items, collapse = "<br>")
    })
    
    # Create base map with ONLY the default basemap
    leaflet() %>%
      addTiles(group = "OpenStreetMap") %>%
      setView(lng = -98.5, lat = 45, zoom = 4) %>%  # Center on US
      # Add gage points
      addCircleMarkers(
        data = goodGages,
        lng = ~longitude,
        lat = ~latitude,
        radius = 6,
        color = ~pal(num_years),
        fillColor = ~pal(num_years),
        fillOpacity = 0.8,
        stroke = TRUE,
        weight = 1,
        popup = hover_text,
        layerId = ~gage_id,
        label = ~paste("Gage ID:", gage_id, "| Years:", num_years),
        group = "Gages"
      ) %>%
      addLegend(
        position = "bottomright",
        pal = pal,
        values = goodGages$num_years,
        title = "Years of Record",
        opacity = 1
      )
  })
  
  # Add additional basemaps after initial render
  observe({
    leafletProxy("gage_map") %>%
      addProviderTiles(providers$Esri.WorldTopoMap, group = "Topographic") %>%
      addProviderTiles(providers$OpenTopoMap, group = "Topographic - 2") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
      addLayersControl(
        baseGroups = c("OpenStreetMap", "Topographic", "Topographic - 2", "Satellite"),
        overlayGroups = c("Gages", "Selected Watershed"),
        options = layersControlOptions(collapsed = FALSE),
        position = "topleft"
      )
  })
  
  # Update map to show ONLY the selected watershed
  observe({
    selected_gage <- input$gage_selector
    
    if (!is.null(selected_gage) && nchar(selected_gage) > 0) {
      # Get the Downstream_HB_ID for this gage
      selected_metadata <- goodGages[goodGages$gage_id == selected_gage, ]
      
      # First, clear any existing watershed
      leafletProxy("gage_map") %>%
        clearGroup("Selected Watershed")
      
      if (nrow(selected_metadata) > 0 && !is.na(selected_metadata$Downstream_HB_ID)) {
        downstream_hb_id <- selected_metadata$Downstream_HB_ID[1]
        
        # Find watershed with matching first_basin_id
        selected_watershed <- watersheds[watersheds$first_basin_id == downstream_hb_id, ]
        
        if (nrow(selected_watershed) > 0) {
          # Calculate bounds for zooming
          bounds <- st_bbox(selected_watershed)
          
          leafletProxy("gage_map") %>%
            # Add ONLY the selected watershed with basin ID in label
            addPolygons(
              data = selected_watershed,
              fillColor = "royalblue",
              fillOpacity = 0.3,
              color = "royalblue",
              weight = 4,
              opacity = 1,
              group = "Selected Watershed",
              label = paste("Basin ID:", selected_watershed$first_basin_id[1])
            ) %>%
            # Zoom to the selected watershed with some padding
            fitBounds(
              lng1 = bounds["xmin"] - 0.5,
              lat1 = bounds["ymin"] - 0.5,
              lng2 = bounds["xmax"] + 0.5,
              lat2 = bounds["ymax"] + 0.5
            )
          
          # Debug message in console
          message(paste("Displayed watershed for gage", selected_gage, 
                        "with Downstream_HB_ID", downstream_hb_id))
        }
      }
    }
  })
  
  # Handle map click events
  observeEvent(input$gage_map_marker_click, {
    click <- input$gage_map_marker_click
    if (!is.null(click)) {
      updateSelectizeInput(session, "gage_selector", selected = click$id)
    }
  })
  
  # Read streamflow data for selected gage
  streamflow_data <- reactive({
    req(input$gage_selector)
    
    # Just load the data without showing a notification
    data <- read_streamflow_by_gage_id(gage_id = input$gage_selector)
    
    return(data)
  })
  
  # Create plotly time series plot with gap detection
  output$streamflow_plot <- renderPlotly({
    data <- streamflow_data()
    
    # Get metadata for selected gage from pre-loaded data
    selected_metadata <- goodGages[goodGages$gage_id == input$gage_selector, ]
    
    plot_title <- paste("Streamflow for Gage", input$gage_selector)
    if (nrow(selected_metadata) > 0 && "station_nm" %in% names(selected_metadata)) {
      plot_title <- paste(plot_title, "-", selected_metadata$station_nm[1])
    }
    
    # Ensure Date column is Date type and sort data by date
    data$Date <- as.Date(data$Date)
    data <- data[order(data$Date), ]
    
    # Remove any duplicate dates (keep first occurrence)
    data <- data[!duplicated(data$Date), ]
    
    # Calculate max Q value for positioning annotations
    max_q <- max(data$Q, na.rm = TRUE)
    if (!is.finite(max_q)) {
      max_q <- 1
    }
    
    # Create complete date sequence
    date_range <- seq(min(data$Date, na.rm = TRUE), 
                      max(data$Date, na.rm = TRUE), 
                      by = "day")
    
    # Find missing dates (gaps in the time series)
    missing_dates <- date_range[!date_range %in% data$Date]
    
    # Identify gap periods (consecutive missing dates)
    gap_periods <- data.frame()
    if (length(missing_dates) > 0) {
      date_diff <- c(1, diff(missing_dates))
      gap_groups <- cumsum(date_diff > 1)
      
      gap_periods <- data.frame(missing_dates, gap_groups) %>%
        group_by(gap_groups) %>%
        summarise(
          start_date = min(missing_dates),
          end_date = max(missing_dates),
          n_days = n()
        ) %>%
        filter(n_days > 1)
    }
    
    # Create the main plot
    p <- plot_ly() %>%
      add_trace(
        x = data$Date, 
        y = data$Q, 
        type = 'scatter', 
        mode = 'lines',
        name = 'Streamflow',
        line = list(color = 'blue', width = 1),
        connectgaps = FALSE,
        hovertemplate = 'Date: %{x}<br>Q: %{y}<extra></extra>'
      )
    
    # Add markers for NA values in the existing data
    na_indices <- which(is.na(data$Q))
    if (length(na_indices) > 0) {
      na_data <- data[na_indices, ]
      min_val <- min(data$Q, na.rm = TRUE)
      placeholder_val <- ifelse(is.finite(min_val), min_val * 0.1, 0.1)
      
      p <- p %>%
        add_trace(
          x = na_data$Date,
          y = rep(placeholder_val, nrow(na_data)),
          type = 'scatter',
          mode = 'markers',
          name = 'NA Values',
          marker = list(
            color = 'orange',
            size = 8,
            symbol = 'x'
          ),
          hovertemplate = 'Date: %{x}<br>Status: NA Value<extra></extra>'
        )
    }
    
    # Create shapes for gap periods
    shapes_list <- list()
    if (nrow(gap_periods) > 0) {
      # Add a dummy trace for legend entry
      p <- p %>%
        add_trace(
          x = as.Date(character(0)),
          y = numeric(0),
          type = 'scatter',
          mode = 'markers',
          name = 'Missing days',
          marker = list(
            color = 'rgba(255, 0, 0, 0.2)',
            size = 15,
            symbol = 'square'
          ),
          showlegend = TRUE
        )
      
      # Create shapes for each gap
      for (i in 1:nrow(gap_periods)) {
        gap <- gap_periods[i, ]
        shapes_list[[i]] <- list(
          type = "rect",
          x0 = gap$start_date - 0.5,
          x1 = gap$end_date + 0.5,
          y0 = 0,
          y1 = 1,
          yref = "paper",
          fillcolor = "rgba(255, 0, 0, 0.2)",
          line = list(color = "rgba(255, 0, 0, 0)"),
          layer = "below"
        )
      }
    }
    
    # Configure layout
    p <- p %>%
      layout(
        title = plot_title,
        xaxis = list(
          title = "Date",
          type = "date",
          rangeslider = list(type = "date")
        ),
        yaxis = list(
          title = "Discharge (Q)",
          type = input$scale_type,
          rangemode = "tozero"
        ),
        hovermode = "x unified",
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          yanchor = "bottom",
          y = 1.02,
          xanchor = "right",
          x = 1
        ),
        shapes = shapes_list
      )
    
    # Add text annotations for long gaps
    if (nrow(gap_periods) > 0) {
      long_gaps <- gap_periods[gap_periods$n_days > 30, ]
      if (nrow(long_gaps) > 0) {
        annotations <- lapply(1:nrow(long_gaps), function(i) {
          gap <- long_gaps[i, ]
          list(
            x = gap$start_date + (gap$end_date - gap$start_date) / 2,
            y = 0.95,
            yref = "paper",
            text = paste(gap$n_days, "days missing"),
            showarrow = FALSE,
            bgcolor = "rgba(255, 255, 255, 0.8)",
            bordercolor = "red",
            borderwidth = 1
          )
        })
        p <- p %>% layout(annotations = annotations)
      }
    }
    
    return(p)
  })
  
  # === NEW METRIC MAPS ===
  
  # Helper function to create metric map
  create_metric_map <- function(metric_name, stat_type) {
    req(input$metric_selector, input$point_size)
    
    # Construct column name
    col_name <- paste0(input$metric_selector, "_", stat_type)
    
    # Check if column exists
    if (!col_name %in% names(signature_data)) {
      return(leaflet() %>% 
               addTiles() %>% 
               setView(lng = -98.5, lat = 45, zoom = 3))
    }
    
    # Get p-value column name for this metric
    pval_col_name <- paste0(input$metric_selector, "_pval")
    
    # Merge signature data with coordinates and p-values
    if (pval_col_name %in% names(signature_data)) {
      map_data <- signature_data %>%
        select(gage_id, latitude, longitude, all_of(col_name), all_of(pval_col_name)) %>%
        filter(!is.na(latitude), !is.na(longitude))
    } else {
      map_data <- signature_data %>%
        select(gage_id, latitude, longitude, all_of(col_name)) %>%
        filter(!is.na(latitude), !is.na(longitude))
    }
    
    # Separate data with valid values from NA values
    map_data_valid <- map_data %>%
      filter(!is.na(!!sym(col_name)))
    
    map_data_na <- map_data %>%
      filter(is.na(!!sym(col_name)))
    
    # Create base map
    m <- leaflet() %>%
      addTiles() %>%
      setView(lng = -98.5, lat = 45, zoom = 3)
    
    # If we have valid data, add colored markers with palette
    if (nrow(map_data_valid) > 0) {
      # Extract the column values
      col_values <- map_data_valid[[col_name]]
      
      # Determine color limits and clamp values based on stat type
      if (stat_type == "slp") {
        # For slope: use max of absolute value of 5th and 95th percentiles
        p05 <- quantile(col_values, 0.05, na.rm = TRUE)
        p95 <- quantile(col_values, 0.95, na.rm = TRUE)
        max_abs <- max(abs(p05), abs(p95))
        
        # Set domain limits
        color_min <- -max_abs
        color_max <- max_abs
        
        # Clamp values to these limits
        col_values_clamped <- pmax(pmin(col_values, color_max), color_min)
        
        # Diverging palette for slope
        pal <- colorNumeric(
          palette = "RdBu",
          domain = c(color_min, color_max),
          reverse = TRUE
        )
        
      } else if (stat_type == "pval") {
        # For p-value: always 0 to 1
        color_min <- 0
        color_max <- 1
        
        # Clamp values (though p-values should already be 0-1)
        col_values_clamped <- pmax(pmin(col_values, color_max), color_min)
        
        # Sequential palette for p-value
        pal <- colorNumeric(
          palette = "YlOrRd",
          domain = c(color_min, color_max),
          reverse = FALSE
        )
        
      } else {
        # For mean and median: use 5th and 95th percentiles
        color_min <- quantile(col_values, 0.05, na.rm = TRUE)
        color_max <- quantile(col_values, 0.95, na.rm = TRUE)
        
        # Clamp values to these limits
        col_values_clamped <- pmax(pmin(col_values, color_max), color_min)
        
        # Sequential palette for mean and median
        pal <- colorNumeric(
          palette = "viridis",
          domain = c(color_min, color_max)
        )
      }
      
      # Create hover text showing ORIGINAL (unclamped) values
      hover_text_valid <- paste0(
        "<b>Gage ID:</b> ", map_data_valid$gage_id, "<br>",
        "<b>", gsub("_", " ", col_name), ":</b> ", 
        round(col_values, 4)
      )
      
      # Determine marker styling based on p-value (applies to ALL plots)
      if (pval_col_name %in% names(map_data_valid)) {
        # Get p-values for these gages
        pval_values <- map_data_valid[[pval_col_name]]
        
        # Black outline (weight=2) if p-value < 0.05, otherwise colored outline (weight=1)
        marker_colors <- ifelse(!is.na(pval_values) & pval_values < 0.05, 
                                "black", 
                                pal(col_values_clamped))
        marker_weights <- ifelse(!is.na(pval_values) & pval_values < 0.05, 2, 1)
      } else {
        # If p-value column doesn't exist, use default styling
        marker_colors <- pal(col_values_clamped)
        marker_weights <- rep(1, length(col_values))
      }
      
      # Add markers for valid data using CLAMPED values for colors
      m <- m %>%
        addCircleMarkers(
          data = map_data_valid,
          lng = ~longitude,
          lat = ~latitude,
          radius = input$point_size,
          color = marker_colors,
          fillColor = pal(col_values_clamped),
          fillOpacity = 0.8,
          stroke = TRUE,
          weight = marker_weights,
          popup = hover_text_valid,
          label = lapply(1:nrow(map_data_valid), function(i) {
            HTML(paste0("Gage: ", map_data_valid$gage_id[i], 
                        " | Value: ", round(col_values[i], 4)))
          }),
          group = "Valid Data"
        ) %>%
        addLegend(
          position = "bottomright",
          pal = pal,
          values = c(color_min, color_max),
          title = gsub("_", " ", col_name),
          opacity = 1
        )
    }
    
    # Add light grey markers for NA values
    if (nrow(map_data_na) > 0) {
      hover_text_na <- paste0(
        "<b>Gage ID:</b> ", map_data_na$gage_id, "<br>",
        "<b>", gsub("_", " ", col_name), ":</b> NA"
      )
      
      m <- m %>%
        addCircleMarkers(
          data = map_data_na,
          lng = ~longitude,
          lat = ~latitude,
          radius = input$point_size,
          color = "lightgrey",
          fillColor = "lightgrey",
          fillOpacity = 0.6,
          stroke = TRUE,
          weight = 1,
          popup = hover_text_na,
          label = ~paste0("Gage: ", gage_id, " | Value: NA"),
          group = "NA Values"
        )
    }
    
    return(m)
  }
  
  # Render the four metric maps
  output$metric_map_mean <- renderLeaflet({
    create_metric_map(input$metric_selector, "mean")
  })
  
  output$metric_map_median <- renderLeaflet({
    create_metric_map(input$metric_selector, "median")
  })
  
  output$metric_map_pval <- renderLeaflet({
    create_metric_map(input$metric_selector, "pval")
  })
  
  output$metric_map_slp <- renderLeaflet({
    create_metric_map(input$metric_selector, "slp")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
