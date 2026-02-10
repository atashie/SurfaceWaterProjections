# Streamflow and Climate Visualization Dashboard
# Enhanced version with:
# - Updated data sources (FEB2026 signatures)
# - 2x3 signature map grid
# - Weather data overlays on hydrograph
# - Custom signature scatter plots with trend analysis

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
library(shinycssloaders)
library(zyp)      # Theil-Sen slope
library(Kendall)  # Mann-Kendall test

# Source helper functions
source("helperFunctions.R")

# === CONSTANTS ===
S3_BUCKET_NAME <- "climate-ai-data-science-shiny-app-data"
SPINNER_COLOR <- "#3498db"
SPINNER_TYPE <- 4
MAP_HEIGHT_SMALL <- "350px"
MAP_HEIGHT_LARGE <- "500px"
DEFAULT_MAP_CENTER <- list(lng = -98.5, lat = 45, zoom = 4)
STAT_TYPES <- c("mean", "median", "linear_slp", "senn_slp", "mk_pval", "spearman_pval")
STAT_LABELS <- c("Mean", "Median", "Linear Slope", "Theil-Sen Slope",
                 "Mann-Kendall P-value", "Spearman P-value")

# Pre-computed p-value palette (fixed 0-1 domain)
PVAL_PALETTE <- colorNumeric(palette = "YlOrRd", domain = c(0, 1), reverse = FALSE)

# Legacy reference (for backwards compatibility)
s3_bucket_name <- S3_BUCKET_NAME

# === AUTHENTICATION AND DATA LOADING ON APP START ===
readRenviron(".Renviron")
check_aws_credentials()

# Read and filter metadata
message("Loading metadata from S3...")
metadata_raw <- read_csv_from_s3_direct(
  bucket = s3_bucket_name,
  object_key = "streamflow/combined_watershed_metadata_09feb2026.csv"
)
goodGages <- subset(metadata_raw, processing_status == "success")
message(paste("Loaded", nrow(goodGages), "gages with successful processing"))

# Load watershed boundaries
message("Loading watershed boundaries...")
watersheds <- st_read("unified_watershedBoundaries_simplified.gpkg", quiet = TRUE)
watersheds <- st_transform(watersheds, 4326)
# Vectorized extraction: sub() returns NA for NA input, so this is safe
watersheds$first_basin_id <- sub(";.*", "", watersheds$basin_ids)
message(paste("Loaded", nrow(watersheds), "watershed boundaries"))

# Check matching
n_matched <- sum(goodGages$Downstream_HB_ID %in% watersheds$first_basin_id)
message(paste("Matched", n_matched, "out of", nrow(goodGages), "gages to watershed boundaries"))

# Load streamflow signature summary data (FEB2026)
message("Loading streamflow signature summary data from S3...")
signature_data <- read_csv_from_s3_direct(
  bucket = s3_bucket_name,
  object_key = "streamflow/streamflow_signatures_full_10feb2026.csv"
)
message(paste("Loaded signature data for", nrow(signature_data), "gages"))

# Pre-filter signature data for maps (removes rows with missing coordinates)
# Keep original signature_data for scatter plot (doesn't need coordinates)
signature_data_for_maps <- signature_data[!is.na(latitude) & !is.na(longitude) &
                                           (is.na(area_normalized) | area_normalized == TRUE)]
message(paste("Filtered to", nrow(signature_data_for_maps), "gages with valid coordinates for maps"))

# Pre-index watershed lookup for faster access (handles NA values)
watersheds_with_id <- watersheds[!is.na(watersheds$first_basin_id), ]
watershed_lookup <- split(watersheds_with_id, watersheds_with_id$first_basin_id)
message(paste("Created watershed lookup with", length(watershed_lookup), "basins"))

# Extract metric names from columns ending in "_mean"
mean_columns <- grep("_mean$", names(signature_data), value = TRUE)
metric_names <- gsub("_mean$", "", mean_columns)
message(paste("Found", length(metric_names), "metrics:", paste(head(metric_names, 3), collapse = ", "), "..."))

# Define statistic choices for dropdowns
stat_choices <- list(
  "Mean" = "mean",
  "Median" = "median",
  "Linear Slope" = "linear_slp",
  "Theil-Sen Slope" = "senn_slp",
  "Mann-Kendall P-val" = "mk_pval",
  "Spearman P-val" = "spearman_pval"
)

# Weather variable info
weather_info <- get_weather_var_info()

# Primary data quality flag for scatter plot filtering
# Remove gages not normalized to watershed area (Q not in mm/day)
scatter_filter_flag <- "area_normalized"

# === DEFINE UI ===
ui <- fluidPage(
  theme = shinytheme("superhero"),

  titlePanel("Streamflow and Climate Dashboard"),

  # Gage selection
  fluidRow(
    column(
      12,
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

  # Main map
  fluidRow(
    column(12, leafletOutput("gage_map", height = "500px"))
  ),

  # Hydrograph section with weather controls
  fluidRow(
    column(12, h3("Streamflow Time Series"))
  ),
  fluidRow(
    column(
      4,
      radioButtons(
        "scale_type",
        label = "Y-Axis Scale:",
        choices = list("Linear" = "linear", "Logarithmic" = "log"),
        selected = "linear",
        inline = TRUE
      )
    ),
    column(
      8,
      checkboxGroupInput(
        "weather_layers",
        label = "Weather Data Overlays:",
        choices = list(
          "Precipitation (hyetograph)" = "prcp",
          "Temperature (min/max)" = "temp",
          "Snow Water Equivalent" = "swe",
          "Vapor Pressure" = "vp",
          "Solar Radiation" = "srad"
        ),
        selected = NULL,
        inline = TRUE
      )
    )
  ),
  fluidRow(
    column(
      12,
      withSpinner(
        plotlyOutput("streamflow_plot", height = MAP_HEIGHT_LARGE),
        type = SPINNER_TYPE,
        color = SPINNER_COLOR
      )
    )
  ),

  # Signature metrics section
  fluidRow(
    column(6, h3("Streamflow Signature Metrics")),
    column(6, h3("Point Size Control"))
  ),
  fluidRow(
    column(
      6,
      selectInput(
        inputId = "metric_selector",
        label = "Select Metric to Plot:",
        choices = metric_names,
        selected = metric_names[1]
      )
    ),
    column(
      6,
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

  # 2x3 metric maps grid - programmatically generated
  # Row 1: Mean and Median
  fluidRow(
    column(6, h4(STAT_LABELS[1]),
      withSpinner(leafletOutput(paste0("metric_map_", STAT_TYPES[1]), height = MAP_HEIGHT_SMALL),
                  type = SPINNER_TYPE, color = SPINNER_COLOR)),
    column(6, h4(STAT_LABELS[2]),
      withSpinner(leafletOutput(paste0("metric_map_", STAT_TYPES[2]), height = MAP_HEIGHT_SMALL),
                  type = SPINNER_TYPE, color = SPINNER_COLOR))
  ),
  # Row 2: Linear Slope and Theil-Sen Slope
  fluidRow(
    column(6, h4(STAT_LABELS[3]),
      withSpinner(leafletOutput(paste0("metric_map_", STAT_TYPES[3]), height = MAP_HEIGHT_SMALL),
                  type = SPINNER_TYPE, color = SPINNER_COLOR)),
    column(6, h4(STAT_LABELS[4]),
      withSpinner(leafletOutput(paste0("metric_map_", STAT_TYPES[4]), height = MAP_HEIGHT_SMALL),
                  type = SPINNER_TYPE, color = SPINNER_COLOR))
  ),
  # Row 3: Mann-Kendall P-value and Spearman P-value
  fluidRow(
    column(6, h4(STAT_LABELS[5]),
      withSpinner(leafletOutput(paste0("metric_map_", STAT_TYPES[5]), height = MAP_HEIGHT_SMALL),
                  type = SPINNER_TYPE, color = SPINNER_COLOR)),
    column(6, h4(STAT_LABELS[6]),
      withSpinner(leafletOutput(paste0("metric_map_", STAT_TYPES[6]), height = MAP_HEIGHT_SMALL),
                  type = SPINNER_TYPE, color = SPINNER_COLOR))
  ),

  # Custom scatter plot section
  fluidRow(
    column(12, hr(), h3("Custom Scatterplots"))
  ),
  fluidRow(
    column(
      6,
      h4("X-Axis"),
      fluidRow(
        column(
          6,
          selectInput(
            "scatter_x_signature",
            "Signature:",
            choices = metric_names,
            selected = metric_names[1]
          )
        ),
        column(
          6,
          selectInput(
            "scatter_x_stat",
            "Statistic:",
            choices = stat_choices,
            selected = "mean"
          )
        )
      )
    ),
    column(
      6,
      h4("Y-Axis"),
      fluidRow(
        column(
          6,
          selectInput(
            "scatter_y_signature",
            "Signature:",
            choices = metric_names,
            selected = if (length(metric_names) > 1) metric_names[2] else metric_names[1]
          )
        ),
        column(
          6,
          selectInput(
            "scatter_y_stat",
            "Statistic:",
            choices = stat_choices,
            selected = "mean"
          )
        )
      )
    )
  ),
  fluidRow(
    column(
      12,
      withSpinner(
        plotlyOutput("signature_scatter_plot", height = MAP_HEIGHT_LARGE),
        type = SPINNER_TYPE,
        color = SPINNER_COLOR
      )
    )
  ),
  fluidRow(
    column(
      12,
      tags$p(
        textOutput("scatter_removed_note", inline = TRUE),
        style = "color: #888; font-size: 12px; margin-top: 5px;"
      )
    )
  )
)

# === DEFINE SERVER ===
server <- function(input, output, session) {

  # === BASE MAP ===
  output$gage_map <- renderLeaflet({
    pal <- colorNumeric(
      palette = c("red", "blue"),
      domain = goodGages$num_years
    )

    # Vectorized hover text: avoids apply() data.frame-to-matrix copy
    # Preserves ALL columns for popup content (identical output)
    col_names <- names(goodGages)
    hover_text <- do.call(paste, c(
      lapply(col_names, function(cn) {
        paste0("<b>", cn, ":</b> ", goodGages[[cn]], "<br>")
      }),
      list(sep = "")
    ))

    leaflet() %>%
      addTiles(group = "OpenStreetMap") %>%
      setView(lng = DEFAULT_MAP_CENTER$lng, lat = DEFAULT_MAP_CENTER$lat,
              zoom = DEFAULT_MAP_CENTER$zoom) %>%
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

  # Add additional basemaps - run once on session start
  observeEvent(TRUE, {
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
  }, once = TRUE)

  # Update map for selected watershed (using pre-indexed lookup for speed)
  observe({
    selected_gage <- input$gage_selector

    if (!is.null(selected_gage) && nchar(selected_gage) > 0) {
      selected_metadata <- goodGages[goodGages$gage_id == selected_gage, ]

      leafletProxy("gage_map") %>%
        clearGroup("Selected Watershed")

      if (nrow(selected_metadata) > 0 && !is.na(selected_metadata$Downstream_HB_ID)) {
        downstream_hb_id <- selected_metadata$Downstream_HB_ID[1]
        hb_id_str <- as.character(downstream_hb_id)

        # Use pre-indexed watershed lookup (handles NA edge case)
        # Note: as.character(NA) returns "NA" string, so check for both
        if (!is.na(hb_id_str) && hb_id_str != "NA" && hb_id_str %in% names(watershed_lookup)) {
          selected_watershed <- watershed_lookup[[hb_id_str]]

          if (nrow(selected_watershed) > 0) {
            bounds <- st_bbox(selected_watershed)

            leafletProxy("gage_map") %>%
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
              fitBounds(
                lng1 = bounds["xmin"] - 0.5,
                lat1 = bounds["ymin"] - 0.5,
                lng2 = bounds["xmax"] + 0.5,
                lat2 = bounds["ymax"] + 0.5
              )
          }
        }
      }
    }
  })

  # Handle map clicks
  observeEvent(input$gage_map_marker_click, {
    click <- input$gage_map_marker_click
    if (!is.null(click)) {
      updateSelectizeInput(session, "gage_selector", selected = click$id)
    }
  })

  # === REACTIVE DATA LOADING ===

  # Single-gage streamflow cache (prevents refetch when toggling weather layers)
  cached_streamflow <- reactiveVal(list(gage_id = NULL, data = NULL))

  # Streamflow data with caching
  streamflow_data <- reactive({
    req(input$gage_selector)
    cache <- cached_streamflow()

    # Return cached data if same gage
    if (!is.null(cache$gage_id) && cache$gage_id == input$gage_selector) {
      return(cache$data)
    }

    # Fetch and cache new data
    data <- read_streamflow_by_gage_id(gage_id = input$gage_selector)
    cached_streamflow(list(gage_id = input$gage_selector, data = data))
    return(data)
  })

  # Daymet climate data (loaded when weather layers are selected)
  daymet_data <- reactive({
    req(input$gage_selector)

    # Only load if weather layers are selected
    if (length(input$weather_layers) == 0) {
      return(data.table())
    }

    read_daymet_by_gage_id(gage_id = input$gage_selector)
  })

  # === HYDROGRAPH WITH WEATHER OVERLAYS ===
  output$streamflow_plot <- renderPlotly({
    data <- streamflow_data()

    selected_metadata <- goodGages[goodGages$gage_id == input$gage_selector, ]
    plot_title <- paste("Streamflow for Gage", input$gage_selector)
    if (nrow(selected_metadata) > 0 && "station_nm" %in% names(selected_metadata)) {
      plot_title <- paste(plot_title, "-", selected_metadata$station_nm[1])
    }

    # Prepare streamflow data (in-place operations for efficiency)
    setDT(data)
    data[, Date := as.Date(Date)]
    setorder(data, Date)
    data <- unique(data, by = "Date")

    max_q <- max(data$Q, na.rm = TRUE)
    if (!is.finite(max_q)) max_q <- 1

    # Date range and gaps
    date_range <- seq(min(data$Date, na.rm = TRUE), max(data$Date, na.rm = TRUE), by = "day")
    missing_dates <- date_range[!date_range %in% data$Date]

    gap_periods <- data.frame()
    if (length(missing_dates) > 0) {
      date_diff <- c(1, diff(missing_dates))
      gap_groups <- cumsum(date_diff > 1)
      gap_periods <- data.frame(missing_dates, gap_groups) %>%
        group_by(gap_groups) %>%
        summarise(
          start_date = min(missing_dates),
          end_date = max(missing_dates),
          n_days = n(),
          .groups = "drop"
        ) %>%
        filter(n_days > 1)
    }

    # Create base plot with streamflow
    p <- plot_ly() %>%
      add_trace(
        x = data$Date,
        y = data$Q,
        type = "scatter",
        mode = "lines",
        name = "Streamflow",
        line = list(color = "blue", width = 1),
        connectgaps = FALSE,
        yaxis = "y",
        hovertemplate = "Date: %{x}<br>Q: %{y:.2f} mm/day<extra></extra>"
      )

    # NA value markers
    na_indices <- which(is.na(data$Q))
    if (length(na_indices) > 0) {
      na_data <- data[na_indices, ]
      min_val <- min(data$Q, na.rm = TRUE)
      placeholder_val <- ifelse(is.finite(min_val), min_val * 0.1, 0.1)

      p <- p %>%
        add_trace(
          x = na_data$Date,
          y = rep(placeholder_val, nrow(na_data)),
          type = "scatter",
          mode = "markers",
          name = "NA Values",
          marker = list(color = "orange", size = 8, symbol = "x"),
          hovertemplate = "Date: %{x}<br>Status: NA Value<extra></extra>"
        )
    }

    # Gap period shapes
    shapes_list <- list()
    if (nrow(gap_periods) > 0) {
      p <- p %>%
        add_trace(
          x = as.Date(character(0)),
          y = numeric(0),
          type = "scatter",
          mode = "markers",
          name = "Missing days",
          marker = list(color = "rgba(255, 0, 0, 0.2)", size = 15, symbol = "square"),
          showlegend = TRUE
        )

      for (i in 1:nrow(gap_periods)) {
        gap <- gap_periods[i, ]
        shapes_list[[length(shapes_list) + 1]] <- list(
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

    # === ADD WEATHER OVERLAYS ===
    weather_layers <- input$weather_layers
    yaxis_list <- list(
      yaxis = list(
        title = "Discharge (mm/day)",
        type = input$scale_type,
        rangemode = if (input$scale_type == "log") "normal" else "tozero",
        side = "left"
      )
    )

    if (length(weather_layers) > 0) {
      daymet <- daymet_data()

      if (nrow(daymet) > 0 && "Date" %in% names(daymet)) {
        # Join by date range of streamflow
        daymet <- daymet[Date >= min(data$Date) & Date <= max(data$Date)]

        axis_position <- 1.0
        axis_offset <- 0.08

        for (i in seq_along(weather_layers)) {
          var_name <- weather_layers[i]
          var_info <- weather_info[[var_name]]

          yaxis_num <- i + 1
          yaxis_name <- paste0("y", yaxis_num)
          yaxis_key <- paste0("yaxis", yaxis_num)

          if (var_name == "prcp" && var_name %in% names(daymet)) {
            # Precipitation as inverted hyetograph (bars from top)
            p <- p %>%
              add_trace(
                x = daymet$Date,
                y = daymet[[var_name]],
                type = "bar",
                name = var_info$label,
                marker = list(color = var_info$color),
                yaxis = yaxis_name,
                hovertemplate = paste0("Date: %{x}<br>", var_info$label, ": %{y:.1f} ", var_info$unit, "<extra></extra>")
              )

            yaxis_list[[yaxis_key]] <- list(
              title = paste0(var_info$label, " (", var_info$unit, ")"),
              overlaying = "y",
              side = "right",
              autorange = "reversed",  # Bars descend from top
              showgrid = FALSE,
              position = axis_position
            )

            axis_position <- axis_position - axis_offset

          } else if (var_name == "temp") {
            # Temperature: plot both tmin and tmax on shared y-axis
            if ("tmin" %in% names(daymet) && "tmax" %in% names(daymet)) {
              # Compute shared range for y-axis limits
              temp_min <- min(daymet$tmin, na.rm = TRUE)
              temp_max <- max(daymet$tmax, na.rm = TRUE)

              # Add tmin trace (blue)
              p <- p %>%
                add_trace(
                  x = daymet$Date,
                  y = daymet$tmin,
                  type = "scatter",
                  mode = "lines",
                  name = "Min Temperature",
                  line = list(color = "rgb(100, 149, 237)", width = 1),
                  yaxis = yaxis_name,
                  hovertemplate = "Date: %{x}<br>Min Temp: %{y:.1f} \u00B0C<extra></extra>"
                )

              # Add tmax trace (red)
              p <- p %>%
                add_trace(
                  x = daymet$Date,
                  y = daymet$tmax,
                  type = "scatter",
                  mode = "lines",
                  name = "Max Temperature",
                  line = list(color = "rgb(220, 20, 60)", width = 1),
                  yaxis = yaxis_name,
                  hovertemplate = "Date: %{x}<br>Max Temp: %{y:.1f} \u00B0C<extra></extra>"
                )

              yaxis_list[[yaxis_key]] <- list(
                title = "Temperature (\u00B0C)",
                overlaying = "y",
                side = "right",
                showgrid = FALSE,
                range = c(temp_min, temp_max),
                position = axis_position
              )

              axis_position <- axis_position - axis_offset
            }

          } else if (var_name %in% names(daymet)) {
            # Other variables as lines
            p <- p %>%
              add_trace(
                x = daymet$Date,
                y = daymet[[var_name]],
                type = "scatter",
                mode = "lines",
                name = var_info$label,
                line = list(color = var_info$color, width = 1),
                yaxis = yaxis_name,
                hovertemplate = paste0("Date: %{x}<br>", var_info$label, ": %{y:.1f} ", var_info$unit, "<extra></extra>")
              )

            yaxis_list[[yaxis_key]] <- list(
              title = paste0(var_info$label, " (", var_info$unit, ")"),
              overlaying = "y",
              side = "right",
              showgrid = FALSE,
              position = axis_position
            )

            axis_position <- axis_position - axis_offset
          }
        }
      }
    }

    # Layout configuration
    layout_args <- c(
      list(
        title = plot_title,
        xaxis = list(
          title = "Date",
          type = "date",
          rangeslider = list(type = "date"),
          domain = c(0, 0.85)  # Leave room for multiple y-axes
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
      ),
      yaxis_list
    )

    p <- do.call(layout, c(list(p), layout_args))

    # Gap annotations
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

  # === METRIC MAPS (2x3 GRID) WITH LEAFLETPROXY ===

  # Helper function to create palette for a metric
  # Returns list with: pal (palette function), color_min, color_max, col_values_clamped
  create_metric_palette <- function(stat_type, col_values) {
    if (stat_type %in% c("mk_pval", "spearman_pval")) {
      # Use pre-computed p-value palette (fixed 0-1 domain)
      return(list(
        pal = PVAL_PALETTE,
        color_min = 0,
        color_max = 1,
        col_values_clamped = pmax(pmin(col_values, 1), 0)
      ))
    } else if (stat_type %in% c("linear_slp", "senn_slp")) {
      # Diverging palette for slopes
      p05 <- quantile(col_values, 0.05, na.rm = TRUE)
      p95 <- quantile(col_values, 0.95, na.rm = TRUE)
      max_abs <- max(abs(p05), abs(p95))
      color_min <- -max_abs
      color_max <- max_abs
      return(list(
        pal = colorNumeric(palette = "RdBu", domain = c(color_min, color_max), reverse = TRUE),
        color_min = color_min,
        color_max = color_max,
        col_values_clamped = pmax(pmin(col_values, color_max), color_min)
      ))
    } else {
      # Sequential palette for mean/median
      color_min <- quantile(col_values, 0.05, na.rm = TRUE)
      color_max <- quantile(col_values, 0.95, na.rm = TRUE)
      return(list(
        pal = colorNumeric(palette = "viridis", domain = c(color_min, color_max)),
        color_min = color_min,
        color_max = color_max,
        col_values_clamped = pmax(pmin(col_values, color_max), color_min)
      ))
    }
  }

  # Initialize all 6 base maps with tiles only (no markers yet)
  lapply(STAT_TYPES, function(stat) {
    output[[paste0("metric_map_", stat)]] <- renderLeaflet({
      leaflet() %>%
        addTiles() %>%
        setView(lng = DEFAULT_MAP_CENTER$lng, lat = DEFAULT_MAP_CENTER$lat, zoom = 3)
    })
  })

  # Update all 6 maps when metric or point size changes (via leafletProxy)
  observe({
    req(input$metric_selector, input$point_size)

    base_col <- input$metric_selector
    point_size <- input$point_size

    # Get p-value column for significance highlighting
    mk_pval_col <- paste0(base_col, "_mk_pval")
    spearman_pval_col <- paste0(base_col, "_spearman_pval")
    pval_col <- if (mk_pval_col %in% names(signature_data_for_maps)) mk_pval_col else spearman_pval_col

    lapply(STAT_TYPES, function(stat_type) {
      col_name <- paste0(base_col, "_", stat_type)
      col_name_display <- gsub("_", " ", col_name)  # Reuse gsub result
      map_id <- paste0("metric_map_", stat_type)

      # Handle column-not-found with control message
      if (!col_name %in% names(signature_data_for_maps)) {
        leafletProxy(map_id) %>%
          clearGroup("Valid Data") %>%
          clearGroup("NA Values") %>%
          clearControls() %>%
          addControl(
            html = paste0("<div style='background: white; padding: 5px;'>Column not found: ", col_name, "</div>"),
            position = "topright"
          )
        return(NULL)
      }

      # Build map data - use pre-filtered signature_data_for_maps
      select_cols <- c("gage_id", "latitude", "longitude", col_name)
      if (pval_col %in% names(signature_data_for_maps)) {
        select_cols <- c(select_cols, pval_col)
      }

      map_data <- signature_data_for_maps %>%
        select(all_of(select_cols))

      map_data_valid <- map_data %>%
        filter(!is.na(!!sym(col_name)))

      map_data_na <- map_data %>%
        filter(is.na(!!sym(col_name)))

      # Start with cleared map
      proxy <- leafletProxy(map_id) %>%
        clearGroup("Valid Data") %>%
        clearGroup("NA Values") %>%
        clearControls()

      if (nrow(map_data_valid) > 0) {
        col_values <- map_data_valid[[col_name]]

        # Create palette (must be computed per-metric, except p-values)
        pal_result <- create_metric_palette(stat_type, col_values)
        pal <- pal_result$pal
        color_min <- pal_result$color_min
        color_max <- pal_result$color_max
        col_values_clamped <- pal_result$col_values_clamped

        # Significance highlighting
        if (pval_col %in% names(map_data_valid)) {
          pval_values <- map_data_valid[[pval_col]]
          marker_colors <- ifelse(!is.na(pval_values) & pval_values < 0.05, "black", pal(col_values_clamped))
          marker_weights <- ifelse(!is.na(pval_values) & pval_values < 0.05, 2, 1)
        } else {
          marker_colors <- pal(col_values_clamped)
          marker_weights <- rep(1, nrow(map_data_valid))
        }

        # Pre-compute hover text and labels (vectorized)
        hover_text_valid <- paste0(
          "<b>Gage ID:</b> ", map_data_valid$gage_id, "<br>",
          "<b>", col_name_display, ":</b> ", round(col_values, 4)
        )
        label_text <- paste0("Gage: ", map_data_valid$gage_id, " | Value: ", round(col_values, 4))

        proxy <- proxy %>%
          addCircleMarkers(
            data = map_data_valid,
            lng = ~longitude,
            lat = ~latitude,
            radius = point_size,
            color = marker_colors,
            fillColor = pal(col_values_clamped),
            fillOpacity = 0.8,
            stroke = TRUE,
            weight = marker_weights,
            popup = hover_text_valid,
            label = label_text,  # Vectorized labels (no lapply)
            group = "Valid Data"
          ) %>%
          addLegend(
            position = "bottomright",
            pal = pal,
            values = c(color_min, color_max),
            title = col_name_display,
            opacity = 1
          )
      }

      # NA markers
      if (nrow(map_data_na) > 0) {
        hover_text_na <- paste0(
          "<b>Gage ID:</b> ", map_data_na$gage_id, "<br>",
          "<b>", col_name_display, ":</b> NA"
        )

        proxy <- proxy %>%
          addCircleMarkers(
            data = map_data_na,
            lng = ~longitude,
            lat = ~latitude,
            radius = point_size,
            color = "lightgrey",
            fillColor = "lightgrey",
            fillOpacity = 0.6,
            stroke = TRUE,
            weight = 1,
            popup = hover_text_na,
            label = paste0("Gage: ", map_data_na$gage_id, " | Value: NA"),
            group = "NA Values"
          )
      }
    })
  })

  # === CUSTOM SCATTER PLOT ===


  # Reactive to compute scatter plot data with flagged gages removed
  # Note: goodGages filtering is NOT applied here - that filter is only used
  # during data pre-processing to select which gages to fetch from HYDAT/USGS.
  # All gages in signature_data have already passed pre-processing QC.
  scatter_plot_data <- reactive({
    req(input$scatter_x_signature, input$scatter_x_stat)
    req(input$scatter_y_signature, input$scatter_y_stat)

    # Build column names
    x_col <- paste0(input$scatter_x_signature, "_", input$scatter_x_stat)
    y_col <- paste0(input$scatter_y_signature, "_", input$scatter_y_stat)

    # Check columns exist
    if (!x_col %in% names(signature_data) || !y_col %in% names(signature_data)) {
      return(list(data = NULL, x_col = x_col, y_col = y_col,
                  n_removed = 0, n_total = 0, removal_reasons = list()))
    }

    # Start with all signature data (non-NA on selected columns)
    all_data <- signature_data %>%
      select(gage_id, all_of(c(x_col, y_col)),
             any_of(scatter_filter_flag)) %>%
      filter(!is.na(!!sym(x_col)), !is.na(!!sym(y_col)))

    n_total <- nrow(all_data)
    removal_reasons <- list()
    plot_data <- all_data

    # Filter out gages not normalized to watershed area
    if (scatter_filter_flag %in% names(plot_data)) {
      before <- nrow(plot_data)
      plot_data <- plot_data %>%
        filter(is.na(!!sym(scatter_filter_flag)) | !!sym(scatter_filter_flag) == TRUE)
      n_removed_qann <- before - nrow(plot_data)
      if (n_removed_qann > 0) {
        removal_reasons$not_area_normalized <- n_removed_qann
      }
    }

    n_removed_total <- n_total - nrow(plot_data)

    list(
      data = plot_data,
      x_col = x_col,
      y_col = y_col,
      n_removed = n_removed_total,
      n_total = n_total,
      removal_reasons = removal_reasons
    )
  })

  output$signature_scatter_plot <- renderPlotly({
    result <- scatter_plot_data()

    x_col <- result$x_col
    y_col <- result$y_col
    plot_data <- result$data

    # Check columns exist
    if (is.null(plot_data)) {
      return(plot_ly() %>%
               layout(title = paste("Column not found:", x_col, "or", y_col)))
    }

    if (nrow(plot_data) < 3) {
      return(plot_ly() %>%
               layout(title = "Insufficient data for scatter plot"))
    }

    x_vals <- plot_data[[x_col]]
    y_vals <- plot_data[[y_col]]

    # Calculate Theil-Sen slope
    sen_result <- calculate_theil_sen(x_vals, y_vals)
    slope <- sen_result$slope
    intercept <- sen_result$intercept

    # Calculate Spearman correlation between x and y
    spearman_result <- calculate_spearman(x_vals, y_vals)

    # Create scatter plot
    p <- plot_ly() %>%
      add_trace(
        x = x_vals,
        y = y_vals,
        type = "scatter",
        mode = "markers",
        text = plot_data$gage_id,
        marker = list(
          color = "rgba(52, 152, 219, 0.7)",
          size = 8,
          line = list(color = "white", width = 1)
        ),
        name = "Gages",
        hovertemplate = paste0(
          "<b>Gage:</b> %{text}<br>",
          "<b>", gsub("_", " ", x_col), ":</b> %{x:.4f}<br>",
          "<b>", gsub("_", " ", y_col), ":</b> %{y:.4f}",
          "<extra></extra>"
        )
      )

    # Add Theil-Sen trend line if valid
    if (!is.na(slope) && !is.na(intercept)) {
      x_range <- range(x_vals)
      y_line <- intercept + slope * x_range

      p <- p %>%
        add_trace(
          x = x_range,
          y = y_line,
          type = "scatter",
          mode = "lines",
          line = list(color = "red", dash = "dash", width = 2),
          name = paste0("Theil-Sen (slope: ", round(slope, 4), ")"),
          hoverinfo = "skip"
        )
    }

    # Build title with statistics
    title_text <- paste0(
      gsub("_", " ", y_col), " vs ", gsub("_", " ", x_col), "<br>",
      "<span style='font-size:12px;'>",
      "Theil-Sen slope: ", ifelse(is.na(slope), "N/A", round(slope, 4)),
      " | Spearman rho: ", ifelse(is.na(spearman_result$rho), "N/A", round(spearman_result$rho, 3)),
      " | p-value: ", ifelse(is.na(spearman_result$pval), "N/A", format(spearman_result$pval, digits = 3, scientific = TRUE)),
      "</span>"
    )

    p <- p %>%
      layout(
        title = list(text = title_text, x = 0.5),
        xaxis = list(title = gsub("_", " ", x_col)),
        yaxis = list(title = gsub("_", " ", y_col)),
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          yanchor = "bottom",
          y = 1.02,
          xanchor = "right",
          x = 1
        )
      )

    return(p)
  })

  # Text output showing number of removed flagged gages
  output$scatter_removed_note <- renderText({
    result <- scatter_plot_data()
    n_plotted <- if (!is.null(result$data)) nrow(result$data) else 0
    reasons <- result$removal_reasons

    if (length(reasons) > 0) {
      parts <- c()
      if (!is.null(reasons$not_area_normalized)) {
        parts <- c(parts, paste0(reasons$not_area_normalized, " not area-normalized"))
      }
      paste0("Showing ", n_plotted, " watersheds. Removed: ",
             paste(parts, collapse = ", "), ".")
    } else {
      paste0("Showing ", n_plotted, " watersheds.")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
