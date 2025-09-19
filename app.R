
# app.R

######################################################
#                 GLOBAL & HELPER FUNCTIONS          #
######################################################
# --- Load necessary libraries ---
# Shiny Core & UI Enhancements
library(shiny)
library(shinydashboard) # Although not used directly, good for potential future expansion
library(shinythemes)
library(shinyWidgets) # For potentially nicer inputs if needed
library(shinyjs)      # To control UI elements dynamically (used for conditionalPanel fix)

# Data Manipulation & Reading
library(readr)
library(data.table) # Useful but potentially replaceable by dplyr/tidyr if preferred
library(dplyr)
library(tidyr)      # Specifically for pivot_longer/wider
library(stringr)
library(reshape2)   # Used in network plot, consider replacing with tidyr if possible

# Plotting & Visualization
library(ggplot2)
library(scales)     # For plot scales
library(plotly)     # Interactive plots
library(heatmaply)  # Interactive heatmaps
library(ggridges)   # Ridgeline plots
library(RColorBrewer) # For clustering plot colors
library(viridis)    # For color scales (used in ridgeline)

# Specific Analysis Libraries
library(igraph)     # Network analysis
library(umap)       # UMAP dimensionality reduction
library(jsonlite)   # Potentially for AG Grid interactions if needed

# AG Grid & DT
library(aggrid)     # Interactive tables
library(DT)         # For DataTables (used in Top Peptides)


#########################
#  Global Functions    #
#########################

# --- Z-score Range Helper ---
get_z_range <- function(slider_val) {
  val <- suppressWarnings(as.numeric(slider_val)) # Avoid warning for NA coercion
  if (is.na(val) || val == 0) {
    return(c(-1, 1))
  } else {
    # Ensure positive value for range
    abs_val <- abs(val)
    return(c(-abs_val, abs_val))
  }
}

# --- Clean Uploaded AG Grid Data ---
# NOTE: This function is highly specific to an expected input format.
# It might need adjustment if the source format changes.
clean_uploaded_data <- function(uploaded_data) {
  # Ensure it's a data frame
  if (!is.data.frame(uploaded_data)) {
    uploaded_data <- as.data.frame(uploaded_data)
  }
  
  # Check dimensions before indexing
  if (nrow(uploaded_data) < 11 || ncol(uploaded_data) < 2) {
    stop("Uploaded data for AG Grid has unexpected dimensions. Required at least 11 rows and 2 columns for cleaning.")
  }
  
  # Perform cleaning steps with basic checks
  tryCatch({
    uploaded_data <- uploaded_data[-c(1, 6:11), ] # Remove specific rows
    if (nrow(uploaded_data) < 4 || ncol(uploaded_data) < 2) {
      stop("Data dimensions incorrect after removing rows.")
    }
    uploaded_data[1:4, 1] <- uploaded_data[1:4, 2] # Copy metadata
    uploaded_data <- uploaded_data[, -2]             # Remove original column 2
    uploaded_data <- as.data.frame(t(uploaded_data)) # Transpose
    
    if (nrow(uploaded_data) < 2) {
      stop("Data dimensions incorrect after transpose.")
    }
    
    colnames(uploaded_data) <- as.character(uploaded_data[1, ]) # Set column names from first row
    uploaded_data <- uploaded_data[-1, ]                     # Remove header row
    rownames(uploaded_data) <- NULL                          # Reset row names
    
    # Identify potential numeric columns (starting from col 3)
    start_col_num <- 3
    if (ncol(uploaded_data) >= start_col_num) {
      numeric_cols_indices <- start_col_num:ncol(uploaded_data)
      # Convert numeric columns, handling potential errors
      # Ensure columns exist before trying to modify them
      valid_indices <- numeric_cols_indices[numeric_cols_indices <= ncol(uploaded_data)]
      if(length(valid_indices) > 0){
        uploaded_data[, valid_indices] <- lapply(uploaded_data[, valid_indices, drop = FALSE], function(x) {
          suppressWarnings(as.numeric(as.character(x))) # Coerce safely
        })
      }
      
    } else {
      warning("No columns found to convert to numeric starting from column 3.")
    }
    
    
    # Rename 'Exposure Time' column if it exists, otherwise check index 4
    exposure_col_name <- 'Exposure_time' # Target name
    current_names <- colnames(uploaded_data)
    
    if ('Exposure Time' %in% current_names) {
      colnames(uploaded_data)[current_names == 'Exposure Time'] <- exposure_col_name
    } else if (length(current_names) >= 4 && !(exposure_col_name %in% current_names)) {
      # Fallback: Rename the 4th column ONLY if 'Exposure_time' doesn't already exist
      colnames(uploaded_data)[4] <- exposure_col_name
      warning("Column 'Exposure Time' not found by name, renamed 4th column to 'Exposure_time'. Check if this is correct.")
    } else if (!(exposure_col_name %in% current_names)) {
      warning("Could not find or rename 'Exposure_time' column (expected at index 4 or named 'Exposure Time').")
    }
    
    return(uploaded_data)
  }, error = function(e) {
    # Provide more context in the error message
    stop(paste("Error during AG Grid data cleaning (check input format & column names):", e$message))
  })
}


# --- Calculate Slopes and Update Signals ---
# Extrapolates/Interpolates signals based on the slope between 2nd and 3rd time points
calculate_slope_and_update_signal <- function(data) {
  required_cols <- c("Sample", "Array", "Cycle", "Exposure_time")
  if (!all(required_cols %in% colnames(data))) {
    missing_cols <- required_cols[!required_cols %in% colnames(data)]
    stop("Missing required columns for slope calculation: ", paste(missing_cols, collapse=", "))
  }
  
  # Ensure Exposure_time is numeric
  data <- data %>% mutate(Exposure_time = as.numeric(as.character(Exposure_time)))
  
  # Identify peptide columns (exclude known metadata)
  known_metadata <- c("Sample", "Array", "Cycle", "Exposure_time", "SampleID", "group_size", # Include potential temp cols
                      "Time2", "Time3", "Signal2", "Signal3", "Slope", "UpdatedSignal")
  peptide_columns <- colnames(data)[!colnames(data) %in% known_metadata]
  if (length(peptide_columns) == 0) {
    warning("No peptide columns detected for slope calculation.")
    return(data) # Return original data if no peptides found
  }
  
  original_col_names <- colnames(data) # Store original column names
  
  # Group and arrange, add group size
  data <- data %>%
    group_by(Sample, Array, Cycle) %>%
    arrange(Exposure_time, .by_group = TRUE) %>%
    mutate(group_size = n())
  
  # Separate groups
  data_large_groups <- data %>% filter(group_size >= 3)
  data_small_groups <- data %>% filter(group_size < 3)
  
  # Process large groups if any exist
  if(nrow(data_large_groups) > 0) {
    # Loop through each peptide column
    for (peptide in peptide_columns) {
      # Ensure peptide column is numeric
      # Using .data[[peptide]] is safer within dplyr context
      data_large_groups <- data_large_groups %>%
        mutate(!!peptide := as.numeric(as.character(.data[[peptide]])))
      
      # Calculate slope and update signal within each group
      data_large_groups <- data_large_groups %>%
        mutate(
          # Use lag/lead or indexing safely within mutate
          Time2 = nth(Exposure_time, 2), # Get 2nd Exposure_time in the group
          Time3 = nth(Exposure_time, 3), # Get 3rd Exposure_time in the group
          Signal2 = nth(.data[[peptide]], 2), # Get 2nd Signal in the group
          Signal3 = nth(.data[[peptide]], 3), # Get 3rd Signal in the group
          # Calculate slope only if valid times/signals exist and time difference is not zero
          Slope = ifelse(
            !is.na(Time2) & !is.na(Time3) & !is.na(Signal2) & !is.na(Signal3) & (Time3 - Time2) != 0,
            (Signal3 - Signal2) / (Time3 - Time2),
            0 # Default slope to 0 if calculation is not possible
          ),
          # Apply update: Use original value if calculation components are missing or slope is 0
          UpdatedSignal = ifelse(Slope == 0 | is.na(Signal2) | is.na(Time2),
                                 .data[[peptide]],
                                 Signal2 + Slope * (Exposure_time - Time2))
        ) %>%
        # Replace the original peptide column with the updated values
        mutate(!!peptide := UpdatedSignal)
    }
    # Remove temporary calculation columns after the loop
    data_large_groups <- data_large_groups %>%
      select(-any_of(c("Time2", "Time3", "Signal2", "Signal3", "Slope", "UpdatedSignal")))
    
  } else {
    warning("No groups found with 3 or more time points for slope calculation.")
  }
  
  # Combine back and clean up
  data_combined <- bind_rows(data_large_groups, data_small_groups) %>%
    ungroup() %>%
    select(-any_of("group_size")) # Remove temporary group size column safely
  
  # Select only original columns to ensure consistency
  # Identify any extra columns that might have been created
  extra_cols <- setdiff(colnames(data_combined), original_col_names)
  if(length(extra_cols) > 0){
    warning("Removing unexpected columns created during slope calculation: ", paste(extra_cols, collapse=", "))
    data_combined <- data_combined %>% select(all_of(original_col_names))
  }
  
  # Round peptide columns (ensure they are numeric first)
  numeric_peptide_cols <- intersect(peptide_columns, colnames(data_combined))
  if (length(numeric_peptide_cols) > 0) {
    data_combined <- data_combined %>%
      mutate(across(all_of(numeric_peptide_cols), ~ round(as.numeric(.), digits = 0)))
  }
  
  return(data_combined)
}


# --- Z-score Calculation ---
z_score <- function(x) {
  # Ensure input is numeric, handle non-numeric gracefully
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  mean_x <- mean(x_num, na.rm = TRUE)
  sd_x <- sd(x_num, na.rm = TRUE)
  
  # Avoid division by zero or NA standard deviation
  if (is.na(sd_x) || sd_x == 0) {
    # Return vector of 0s with same length as input
    return(rep(0, length(x_num)))
  } else {
    return((x_num - mean_x) / sd_x)
  }
}

# --- Percentize Data (Scale 0-100 per column) ---
percentize <- function(data) {
  # Operate on columns using lapply
  data[] <- lapply(data, function(x) {
    # Ensure input is numeric
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    rng <- range(x_num, na.rm = TRUE)
    
    # Handle cases where range is invalid (all NA) or zero (constant value)
    if (any(!is.finite(rng)) || diff(rng) == 0) {
      # Return vector of 0s or NAs matching original input structure
      scaled_x <- rep(0, length(x_num))
      scaled_x[is.na(x_num)] <- NA # Preserve original NAs
      return(scaled_x)
    } else {
      # Perform scaling: (value - min) / (max - min) * 100
      scaled_x <- (x_num - rng[1]) / (rng[2] - rng[1]) * 100
      # NAs in input should remain NA in output
      # scaled_x[is.na(x_num)] <- NA # This should already happen if x_num had NAs
      return(scaled_x)
    }
  })
  return(data) # Return the modified data frame
}

# --- Process Main Kinase Data File ---
process_data <- function(main_path) {
  tryCatch({
    data_raw <- read_csv(main_path, col_names = TRUE, show_col_types = FALSE)
    
    # Basic validation
    if (ncol(data_raw) < 5) {
      stop("Main data file must have at least 5 columns (metadata + >=1 peptide).")
    }
    required_metadata <- c("Sample", "Array", "Cycle", "Exposure_time")
    if (!all(required_metadata %in% colnames(data_raw)[1:min(4, ncol(data_raw))])) {
      warning("Expected metadata columns (Sample, Array, Cycle, Exposure_time) not found in the first 4 columns or insufficient columns. Processing continues, but check format.")
    }
    
    # Ensure correct types for metadata (factors/characters initially are fine)
    # Use mutate + across for safety and clarity
    data_raw <- data_raw %>%
      mutate(across(any_of(required_metadata), as.character)) # Ensure character type for consistency
    
    # Identify peptide columns (assuming they start from column 5)
    peptide_ids <- character(0) # Initialize empty vector
    if (ncol(data_raw) >= 5) {
      peptide_ids <- colnames(data_raw)[5:ncol(data_raw)]
    }
    
    if(length(peptide_ids) == 0) {
      stop("No peptide columns found (expected starting from column 5). Check file format.")
    }
    
    # Convert peptide columns to numeric, coercing errors to NA
    data_raw <- data_raw %>%
      mutate(across(all_of(peptide_ids), ~ suppressWarnings(as.numeric(as.character(.)))))
    
    # Data for heatmap/clustering (data2 is the primary format)
    data2 <- data_raw
    
    # Data for some legacy plots/summaries? (data1: peptides as rows, samples as cols)
    # Consider if this format is truly needed or if data2 + pivoting is sufficient
    # Create data1 only if peptides exist
    data1 <- NULL
    if(length(peptide_ids) > 0){
      # Create a unique sample identifier for column names in data1 if needed
      sample_names_for_data1 <- make.unique(paste(data_raw$Sample, data_raw$Array, data_raw$Cycle, data_raw$Exposure_time, sep="_"))
      if(length(sample_names_for_data1) != nrow(data_raw)){
        warning("Could not create unique sample names for data1 transpose. Using generic names.")
        sample_names_for_data1 <- paste("Sample", 1:nrow(data_raw), sep="_")
      }
      
      data1_matrix <- t(data_raw[, peptide_ids, drop = FALSE])
      colnames(data1_matrix) <- sample_names_for_data1
      data1 <- as.data.frame(data1_matrix) # Convert to data frame
    }
    
    
    # Extract choices for UI controls
    Sample_choices <- sort(unique(data2$Sample))
    array_choices <- sort(unique(data2$Array))
    # Ensure numeric choices are sorted, handle potential NAs/non-numerics
    exposure_choices <- sort(unique(suppressWarnings(as.numeric(data2$Exposure_time))), na.last = TRUE)
    exposure_choices <- exposure_choices[!is.na(exposure_choices)] # Remove NAs if any
    cycle_choices <- sort(unique(suppressWarnings(as.numeric(data2$Cycle))), na.last = TRUE)
    cycle_choices <- cycle_choices[!is.na(cycle_choices)] # Remove NAs
    
    peptide_choices <- peptide_ids # Already extracted
    
    # Return processed data and choices
    list(
      raw = data_raw, # Keep original raw read
      data1 = data1, # Peptides as rows (or NULL if no peptides)
      data2 = data2, # Main data format (samples as rows)
      Sample_choices = Sample_choices,
      array_choices = array_choices,
      exposure_choices = exposure_choices,
      cycle_choices = cycle_choices,
      peptide_choices = peptide_choices
    )
  }, error = function(e) {
    showModal(modalDialog(
      title = "Error Processing Main Data",
      paste("Could not process the main data file:", e$message, "\nPlease check the file format, column names, and data types."),
      easyClose = TRUE, footer = modalButton("Dismiss")
    ))
    return(NULL) # Return NULL on error
  })
}

# --- Generate Network Plot ---
# Based on peptide-peptide correlations
generate_network_plot <- function(data, layout_type = "fr", threshold = 0, seed = 123, top_n_peptides = NULL) {
  set.seed(seed)
  
  if (is.null(data) || ncol(data) < 5) {
    stop("Network plot requires data with metadata and peptide columns (at least 5 total).")
  }
  
  # Extract peptide data (columns 5 onward)
  peptide_cols_indices <- 5:ncol(data)
  if (length(peptide_cols_indices) == 0) stop("No peptide columns found in data for network.")
  peptide_data <- data[, peptide_cols_indices, drop = FALSE]
  
  # Ensure all peptide columns are numeric, handle potential NAs
  peptide_data <- peptide_data %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
  
  # --- Data Filtering & Preparation ---
  # 1. Remove peptides with zero variance or all NAs
  col_variances <- apply(peptide_data, 2, var, na.rm = TRUE)
  cols_to_keep <- !is.na(col_variances) & col_variances > 1e-8 # Use small tolerance
  if(sum(!cols_to_keep) > 0) {
    original_peptides <- colnames(peptide_data)
    removed_peptides <- original_peptides[!cols_to_keep]
    warning("Removing peptides with zero/low variance or all NAs: ", paste(removed_peptides, collapse=", "))
    peptide_data <- peptide_data[, cols_to_keep, drop = FALSE]
  }
  
  if (ncol(peptide_data) < 2) {
    stop("Insufficient peptides remaining (< 2) after removing low/zero variance columns.")
  }
  
  # 2. Select top N peptides based on mean signal if specified
  current_peptide_names <- colnames(peptide_data)
  if (!is.null(top_n_peptides) && is.numeric(top_n_peptides) && top_n_peptides > 0 && top_n_peptides < ncol(peptide_data)) {
    peptide_means_all <- colMeans(peptide_data, na.rm = TRUE)
    # Sort peptides by mean, handle NAs in means if any
    sorted_peptides <- names(sort(peptide_means_all, decreasing = TRUE, na.last = TRUE))
    # Select top N valid peptide names
    top_peptides_names <- head(sorted_peptides[!is.na(peptide_means_all[sorted_peptides])], top_n_peptides)
    
    if (length(top_peptides_names) < 2) {
      stop("Insufficient top peptides found (< 2) to build a network.")
    }
    peptide_data <- peptide_data[, top_peptides_names, drop = FALSE]
    current_peptide_names <- top_peptides_names # Update list of peptides being used
  }
  
  # --- Correlation and Graph Creation ---
  # Compute correlation matrix using pairwise complete observations
  peptide_correlation_matrix <- cor(peptide_data, use = "pairwise.complete.obs")
  
  # Replace NA correlations (can happen with low overlap or zero variance pairs) with 0
  peptide_correlation_matrix[is.na(peptide_correlation_matrix)] <- 0
  
  # Create graph from adjacency matrix
  G_peptide_based <- graph_from_adjacency_matrix(
    peptide_correlation_matrix,
    mode = "undirected", # Correlation is symmetric
    weighted = TRUE,    # Use correlation value as weight
    diag = FALSE        # Ignore self-loops
  )
  
  # Use absolute correlation for edge weights and filtering
  E(G_peptide_based)$weight_absolute <- abs(E(G_peptide_based)$weight)
  
  # Remove edges below the absolute threshold
  G_peptide_filtered <- delete_edges(G_peptide_based, E(G_peptide_based)[weight_absolute < threshold])
  
  # Check if any edges remain
  if (ecount(G_peptide_filtered) == 0) {
    warning("No edges remain after applying the correlation threshold (", threshold, "). Network will only show nodes. Consider lowering the threshold.")
    # Proceed with nodes only
  }
  
  # --- Layout Calculation ---
  layout_func <- switch(
    layout_type,
    fr = layout_with_fr,
    kk = layout_with_kk,
    drl = layout_with_drl,
    nicely = layout_nicely, # Added 'nicely' option
    layout_nicely # Default fallback
  )
  
  # Use the filtered graph for layout if it has edges, otherwise use the full graph structure
  # Ensure the graph passed to layout has the correct vertices even if edges were removed
  graph_for_layout <- G_peptide_based # Always use base graph structure for layout nodes
  V(graph_for_layout)$name <- current_peptide_names # Ensure names are correct
  
  if (vcount(graph_for_layout) == 0) stop("No nodes (peptides) available for layout.")
  
  current_layout <- layout_func(graph_for_layout)
  
  # Ensure layout has row names corresponding to nodes (peptides)
  rownames(current_layout) <- V(graph_for_layout)$name
  colnames(current_layout) <- c('X', 'Y')
  
  # Normalize layout coordinates (0.05 to 0.95 range)
  normalize_coord <- function(coords) {
    # Handle cases with NA or single point
    finite_coords <- coords[is.finite(coords)]
    if (length(finite_coords) == 0) return(rep(0.5, length(coords)))
    min_c <- min(finite_coords)
    max_c <- max(finite_coords)
    if (max_c == min_c) return(rep(0.5, length(coords))) # Center if no range
    # Apply normalization only to finite values
    norm <- coords # Initialize with original values (including potential NAs)
    norm[is.finite(coords)] <- (finite_coords - min_c) / (max_c - min_c)
    return(0.9 * norm + 0.05)
  }
  current_layout[, 1] <- normalize_coord(current_layout[, 1])
  current_layout[, 2] <- normalize_coord(current_layout[, 2])
  
  # Extract layout coordinates, ensuring names are preserved
  Xn <- current_layout[, 1]
  Yn <- current_layout[, 2]
  # Ensure names match the peptides included
  names(Xn) <- rownames(current_layout)
  names(Yn) <- rownames(current_layout)
  
  
  # --- Node Properties ---
  # Compute average expression for each peptide *actually in the network*
  peptides_in_network <- names(Xn)
  peptide_means <- colMeans(peptide_data[, peptides_in_network, drop = FALSE], na.rm = TRUE)
  
  # Node size based on average expression
  min_size <- 10
  max_size <- 30
  # Handle case where all means might be NA or constant
  valid_means <- peptide_means[!is.na(peptide_means)]
  if(length(valid_means) == 0 || diff(range(valid_means)) == 0) {
    node_size_values <- rep(mean(c(min_size, max_size)), length(peptide_means))
  } else {
    node_size_values <- scales::rescale(peptide_means, to = c(min_size, max_size), from = range(valid_means))
  }
  node_size <- setNames(node_size_values, names(peptide_means))
  node_size[is.na(node_size)] <- min_size # Default size for NAs or excluded peptides
  
  # Node color based on Z-score of average expression (Regulation)
  peptide_z <- z_score(peptide_means) # Use our safe z_score function
  max_abs_z <- max(abs(peptide_z), na.rm = TRUE)
  # Normalize Z-score to 0-1 range for color mapping (centering at 0.5)
  z_norm <- if(is.na(max_abs_z) || max_abs_z == 0) {
    rep(0.5, length(peptide_z))
  } else {
    (peptide_z + max_abs_z) / (2 * max_abs_z)
  }
  z_norm[is.na(z_norm)] <- 0.5 # Default color (white) for NAs or constant values
  
  # Create a diverging color scale: Blue (low) -> White (mid) -> Red (high)
  color_palette_func <- colorRampPalette(c("deepskyblue", "white", "tomato"))
  node_color_values <- color_palette_func(100)[cut(z_norm, breaks = 100, include.lowest = TRUE, labels = FALSE)]
  node_color_values[is.na(node_color_values)] <- "#FFFFFF" # White for any remaining NAs
  node_color <- setNames(node_color_values, names(peptide_means))
  
  # --- Edge Properties ---
  edge_plotly_shapes <- list() # Use shapes for edges
  if (ecount(G_peptide_filtered) > 0) {
    edges_df <- as_data_frame(G_peptide_filtered, what = "edges")
    edge_weights_abs <- E(G_peptide_filtered)$weight_absolute
    
    # Scale edge width
    min_width <- 0.5
    max_width <- 4
    valid_weights <- edge_weights_abs[!is.na(edge_weights_abs)]
    if(length(valid_weights) == 0 || diff(range(valid_weights)) == 0){
      edge_widths_values <- rep(mean(c(min_width, max_width)), length(edge_weights_abs))
    } else {
      edge_widths_values <- scales::rescale(edge_weights_abs, to = c(min_width, max_width), from = range(valid_weights))
    }
    edge_widths_values[is.na(edge_widths_values)] <- min_width # Default width
    
    # Create list for plotly shapes (lines)
    for (i in 1:nrow(edges_df)) {
      node0 <- edges_df$from[i]
      node1 <- edges_df$to[i]
      
      # Check if nodes exist in layout (they should, but safety first)
      if (node0 %in% names(Xn) && node1 %in% names(Xn)) {
        edge_plotly_shapes[[length(edge_plotly_shapes) + 1]] <- list(
          type = "line",
          line = list(color = "rgba(180,180,180,0.6)", width = edge_widths_values[i]), # Semi-transparent gray
          x0 = Xn[node0],
          y0 = Yn[node0],
          x1 = Xn[node1],
          y1 = Yn[node1],
          layer = "below" # Draw edges below nodes
        )
      }
    }
  }
  
  
  # --- Build Plotly Graph ---
  p <- plot_ly() %>%
    # Add Nodes
    add_trace(
      type = 'scatter',
      mode = 'markers', # Start with only markers for clarity if text overlaps
      x = Xn[peptides_in_network], # Ensure ordering/selection matches
      y = Yn[peptides_in_network],
      marker = list(
        size = node_size[peptides_in_network],
        color = node_color[peptides_in_network],
        line = list(color = 'black', width = 0.5) # Thinner outline
      ),
      hoverinfo = 'text', # Text displayed on HOVER
      # Provide comprehensive hover text
      text = paste(
        "Peptide:", peptides_in_network,
        "<br>Avg Signal:", round(peptide_means[peptides_in_network], 2),
        "<br>Z-score:", round(peptide_z[peptides_in_network], 2)
      ),
      # Add text labels separately if needed, maybe smaller font
      # Or rely on hover + selection/zooming
      showlegend = FALSE,
      key = peptides_in_network # Add key for potential linking/brushing
    ) %>%
    # Optionally add text labels
    add_trace(
      type = 'scatter',
      mode = 'text',
      x = Xn[peptides_in_network],
      y = Yn[peptides_in_network],
      text = peptides_in_network,
      textposition = 'bottom center',
      textfont = list(size = 8, color = '#555555'), # Smaller, grayish text
      hoverinfo = 'none', # Don't duplicate hover
      showlegend = FALSE
    ) %>%
    # Layout configuration
    layout(
      title = paste("Peptide Correlation Network (Threshold:", threshold, ", Layout:", layout_type, ")"),
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 1)), # Fix range 0-1
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 1)), # Fix range 0-1
      shapes = edge_plotly_shapes, # Add edges as shapes
      hovermode = 'closest',
      paper_bgcolor = 'white',
      plot_bgcolor = 'white'
      # Consider adding dragmode = 'select' or 'lasso' for interactivity
    )
  
  # Return layout coordinates, plot object, and the filtered peptide data used
  # Ensure layout is a data frame with peptide names as rownames
  layout_df <- as.data.frame(current_layout)
  rownames(layout_df) <- peptides_in_network
  
  list(layout = layout_df, plot = p, data = peptide_data)
}


###################################################
#                     UI                          #
###################################################
ui <- fluidPage(
  theme = shinytheme("flatly"), # Or try "lumen", "yeti", "spacelab"
  shinyjs::useShinyjs(), # Initialize shinyjs
  tags$head(
    tags$style(HTML("
      /* Notification Styling */
      .shiny-notification { position:fixed; top: calc(10%); left: calc(50%); width: auto; max-width: 400px; transform: translateX(-50%); z-index: 9999; }
      
      /* Tab Content Padding */
      .tab-content { padding-top: 15px; }
      
      /* Ensure sidebar scrolls if content overflows */
      .well { overflow-y: auto; } /* SidebarPanel has class 'well' */

      /* --- Box Styling (Mimic shinydashboard boxes) --- */
      .custom-box {
        position: relative;
        border-radius: 3px;
        background: #ffffff;
        border-top: 3px solid #d2d6de; /* Default border color */
        margin-bottom: 20px;
        width: 100%;
        box-shadow: 0 1px 1px rgba(0, 0, 0, 0.1);
      }
      .custom-box-header {
        color: #444;
        display: block;
        padding: 10px;
        position: relative;
        border-bottom: 1px solid #f4f4f4;
      }
      .custom-box-header .box-title {
        display: inline-block;
        font-size: 18px;
        margin: 0;
        line-height: 1;
        font-weight: bold;
      }
      .custom-box-body {
        border-top-left-radius: 0;
        border-top-right-radius: 0;
        border-bottom-right-radius: 3px;
        border-bottom-left-radius: 3px;
        padding: 10px;
      }
       /* Colored box borders */
      .custom-box.box-primary { border-top-color: #3c8dbc; }
      .custom-box.box-info { border-top-color: #00c0ef; }
      .custom-box.box-success { border-top-color: #00a65a; }
      .custom-box.box-warning { border-top-color: #f39c12; }
      .custom-box.box-danger { border-top-color: #dd4b39; }
      
      .custom-box-body ol li { margin-bottom: 10px; }
      .custom-box-body ul { margin-top: 5px; margin-left: 20px;} /* Indent sub-lists */
      .custom-box-body ul li { margin-bottom: 5px; }
      
      /* Ensure icons are aligned with text */
      .box-title .fa, .box-title .fas, .box-title .glyphicon { margin-right: 5px; } 
      .list-unstyled { padding-left: 0; list-style: none; }
      .action-button { margin-top: 5px; margin-bottom: 5px; } /* Add space around buttons */
    "))
  ),
  # --- Improved Application Title Bar ---
  div(style="padding: 10px 15px; font-size:24px; background-color:#2c3e50; color:white; margin-bottom:10px;",
      HTML(paste(icon("dna"), "<b>KinoViz</b>: A User-Friendly Web Application for High-Throughput Kinome Profiling Analysis and Visualization in Cancer Research"))
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3, # Maintain width or adjust as needed
      shinyjs::hidden(numericInput("data_is_loaded_flag", "Data Loaded?", value = 0)),
      
      h4(tagList(icon("database"), "1. Load Data")), # Icon added
      fileInput('file_main', 'Upload Main CSV File', accept = c('.csv', '.txt', '.tsv')),
      helpText("Select your dataset or use the default ('clean_kinome.csv') by clicking 'Load & Process Data' directly."),
      actionButton("btn_process", "Load & Process Data", icon = icon("play"), class="btn-primary btn-block"), # Block level button
      hr(),
      
      h4(tagList(icon("filter"), "2. Filter Data")), # Icon added
      helpText("Select conditions to analyze. Filters update plots/tables dynamically."),
      # Conditional panels remain the same logic
      conditionalPanel(
        condition = "input.data_is_loaded_flag > 0",
        # 'Sample' label is okay as it matches the column name used
        uiOutput("ui_Sample"),
        uiOutput("ui_array"),
        uiOutput("ui_exposure"),
        uiOutput("ui_cycle"),
        tags$small(icon("info-circle"), "Note: Specific filter requirements apply to certain tabs (check help text in each tab).")
      ),
      conditionalPanel(
        condition = "input.data_is_loaded_flag == 0",
        tags$p(tags$em("Load data using the button above to enable filters."))
      ),
      hr()
    ),
    mainPanel(
      width=9,
      tabsetPanel(
        id = "main_tabs",
        
        # ---- NEW WELCOME TAB ----
        tabPanel(tagList(icon("home"), "Welcome"), value = "tab_welcome",
                 fluidRow(
                   # Use styled divs to mimic boxes
                   tags$div(class = "col-sm-12", # Full width column
                            tags$div(class="custom-box box-primary", # Box with primary color border
                                     tags$div(class="custom-box-header",
                                              tags$h3(class="box-title", tagList(icon("info-circle"), "Welcome to KinoViz!"))
                                     ),
                                     tags$div(class="custom-box-body",
                                              h4("Analyze and Visualize Your Kinome Profiling Data"),
                                              p("KinoViz provides an interactive interface to explore peptide phosphorylation patterns from high-throughput kinome profiling experiments."),
                                              p("Load your dataset (or use the sample), filter by experimental conditions, perform analyses, and generate visualizations to gain insights into kinase activity and signaling pathways.")
                                     ) # end box-body
                            ) # end custom-box
                   ), # end col-sm-12
                   tags$div(class = "col-sm-12",
                            tags$div(class="custom-box box-info", # Box with info color border
                                     tags$div(class="custom-box-header",
                                              tags$h3(class="box-title", tagList(icon("tasks"), "Getting Started: Recommended Workflow"))
                                     ),
                                     tags$div(class="custom-box-body",
                                              tags$ol(
                                                tags$li(
                                                  strong(tagList(icon("upload"), "Load Data:")),
                                                  " Use the sidebar ('1. Load Data').",
                                                  tags$ul(
                                                    tags$li("Click ", strong("Browse..."), " to upload your main kinome CSV."),
                                                    tags$li("Or click ", strong("Load & Process Data"), " to use the default ", code("clean_kinome.csv"), "(if available)."),
                                                    tags$li(strong("Expected Format:"), " Columns: ", code("Sample"), ", ", code("Array"), ", ", code("Cycle"), ", ", code("Exposure_time"), ", followed by peptide signal columns.")
                                                  )
                                                ),
                                                tags$li(
                                                  strong(tagList(icon("filter"), "Apply Filters:")),
                                                  " Use the sidebar ('2. Filter Data') once data is loaded.",
                                                  tags$ul(
                                                    tags$li("Select ", strong("Sample(s)"), ", ", strong("Array(s)"), ", ", strong("Cycle"), ", and ", strong("Exposure Time"), "."),
                                                    tags$li(em("Note:"), " Filters apply dynamically. Check specific tab requirements.")
                                                  )
                                                ),
                                                tags$li(
                                                  strong(tagList(icon("chart-pie"), "Explore Analysis Tabs:")),
                                                  " Navigate through the tabs above to visualize:",
                                                  tags$ul(
                                                    tags$li(code("Explore & Clean"), ": View/clean data, apply slope correction."),
                                                    tags$li(code("Top Peptides"), ": Identify high-signal peptides."),
                                                    tags$li(code("Kinetics"), ": Track single peptide signals over time."),
                                                    tags$li(code("Comparisons"), ": Compare multiple peptides."),
                                                    tags$li(code("Heatmap"), ": Visualize signal patterns."),
                                                    tags$li(code("Dim Reduction"), ": Cluster samples (PCA/UMAP)."),
                                                    tags$li(code("Network"), ": Analyze peptide correlations."),
                                                    tags$li(code("Summary"), ": Dataset overview.")
                                                  )
                                                ),
                                                tags$li(
                                                  strong(tagList(icon("mouse-pointer"), "Interact:")),
                                                  tags$ul(
                                                    tags$li("Hover over plots for details."),
                                                    tags$li("Use plot tools (zoom, pan, save)."),
                                                    tags$li("Sort/filter tables.")
                                                  )
                                                )
                                              ) # End ordered list
                                     ) # end box-body
                            )# end custom-box
                   ) # end col-sm-12
                 ) # End fluidRow
        ), # End tabPanel("Welcome")
        
        # ---- EXISTING TABS (WITH IMPROVEMENTS) ----
        tabPanel(tagList(icon("table"), "Explore & Clean"), value="tab_aggrid",
                 h3("Interactive Data Table & Optional Cleaning"),
                 fluidRow(
                   column(6,
                          fileInput("file_upload_aggrid", "Optional: Upload Raw CSV for Cleaning", accept = c(".csv")),
                          helpText("Upload a differently formatted file here for automated cleaning if needed.")
                   ),
                   column(6, style = "padding-top: 25px;", # Align button
                          actionButton("apply_slope", "Apply Slope Correction", icon = icon("calculator"), class="btn-info btn-sm"),
                          helpText("Corrects signals based on trend. Requires >= 3 time points per Sample/Array/Cycle.")
                   )
                 ),
                 hr(),
                 helpText("Main table displays currently loaded data ('Sample' shown). Use controls above to load/clean data specifically for this view or apply slope correction. Sort by clicking headers, filter using column menus."),
                 aggridOutput("aggrid_table", height = "600px"),
                 br(),
                 downloadButton("download_aggrid_data", "Download Grid Data", icon = icon("download"), class="btn-success btn-sm")
        ),
        
        tabPanel(tagList(icon("award"), "Top Peptides"), value="tab_top_peptides",
                 sidebarLayout(
                   sidebarPanel(width = 3, # Narrower sidebar within the tab
                                h4("Select Peptides"),
                                radioButtons("peptide_selection_mode", "Selection Method:",
                                             choices = c("Top N by Max Signal" = "top_n",
                                                         "Select Specific Peptides" = "specific"),
                                             selected = "top_n"),
                                hr(),
                                
                                conditionalPanel(
                                  condition = "input.peptide_selection_mode == 'top_n'",
                                  selectInput("num_peptides", "Number of Top Peptides:",
                                              choices = c(5, 10, 15, 20, 30, 50), selected = 10),
                                  helpText(icon("info-circle"), "Finds peptides with highest max signal within the main filtered 'Sample(s)'.")
                                ),
                                
                                conditionalPanel(
                                  condition = "input.peptide_selection_mode == 'specific'",
                                  uiOutput("ui_specific_peptides_selector"), # Dynamic selector
                                  helpText(icon("info-circle"), "Choose specific peptides from the list.")
                                ),
                                
                                hr(),
                                tags$small(icon("filter"), "Uses 'Sample(s)' filter from the main sidebar.")
                   ),
                   mainPanel(width=9,
                             h4(textOutput("top_peptides_table_title")),
                             DTOutput("table_important_peptides"),
                             hr(),
                             h4(textOutput("top_peptides_plot_title")),
                             helpText(icon("chart-density"), "Shows signal density across different 'Arrays' for the selected peptides and filtered 'Sample(s)'. Useful for comparing distribution shapes."),
                             plotOutput("plot_ridgeline", height = "600px")
                   )
                 )
        ),
        
        tabPanel(tagList(icon("chart-line"), "Kinetics"), value="tab_kinetics",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                h4("Analysis Options"),
                                radioButtons("phospho_mode", "Comparison Mode:",
                                             choices = c("Single Sample Focus" = "single",
                                                         "Compare Across Samples" = "all"),
                                             selected = "single"),
                                uiOutput("ui_peptide_single"), # Peptide selection always needed
                                hr(),
                                helpText(icon("info-circle"), "Mode 'Single Sample Focus' uses Sample/Array filters. Mode 'Compare Across Samples' requires fixed Cycle/Exposure filters from main sidebar.")
                   ),
                   mainPanel(width=9,
                             h4("Kinetics vs. Cycle"),
                             helpText("Signal over cycles. Color: Exposure (Single) or Sample (All). Line: Array. Req: Fixed Array/Exposure (All Mode)."),
                             plotlyOutput("plot_cycle_kinetics", height = 400),
                             hr(),
                             h4("Kinetics vs. Exposure Time"),
                             helpText("Signal over exposure. Color: Cycle (Single) or Sample (All). Line: Array. Req: Fixed Array/Cycle (All Mode)."),
                             plotlyOutput("plot_exposure_kinetics", height = 400)
                   )
                 )
        ),
        
        tabPanel(tagList(icon("exchange-alt"), "Comparisons"), value="tab_peptide_details",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                h4("Select Peptides"),
                                uiOutput("ui_peptide_multi"),
                                hr(),
                                helpText(icon("info-circle"), "Select multiple peptides to compare. Requires fixed Cycle or Exposure filter from main sidebar.")
                   ),
                   mainPanel(width=9,
                             h4("Selected Peptides vs. Cycle"),
                             helpText("Compares selected peptides over Cycle. Requires a fixed Exposure Time filter."),
                             plotlyOutput("plot_peptides_vs_cycle", height = 400),
                             hr(),
                             h4("Selected Peptides vs. Exposure Time"),
                             helpText("Compares selected peptides over Exposure Time. Requires a fixed Cycle filter."),
                             plotlyOutput("plot_peptides_vs_exposure", height = 400)
                   )
                 )
        ),
        
        tabPanel(tagList(icon("th"), "Heatmap"), value="tab_heatmap",
                 h4("Peptide Signal Intensity Heatmap"),
                 helpText(icon("info-circle"), " Visualizes signals for the filtered data at a ", strong("fixed Cycle and Exposure Time"),
                          ". Heatmaps are shown per selected 'Array'. Samples (rows) are ordered by Sample then Array name. Peptides (columns) are ordered by global hierarchical clustering across all selected arrays. Color scale (Blue=Low, Red=High) is normalized per peptide (0-100) across all samples shown."),
                 plotlyOutput("plot_heatmap", height = "1000px") # Increased height
        ),
        
        tabPanel(tagList(icon("project-diagram"), "Dim Reduction"), value="tab_dim_red",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                h4("Reduction Options"),
                                selectInput("clustering_method", "Method:", choices = c("UMAP", "PCA"), selected="UMAP"),
                                actionButton("btn_clustering", "Run Reduction Analysis", icon=icon("play"), class="btn-primary btn-block"),
                                hr(),
                                helpText(icon("info-circle"), "Performs dimensionality reduction (UMAP or PCA) on samples. ", strong("Requires fixed Cycle and Exposure Time filters"), " set in the main sidebar. Points represent Sample/Array combinations.")
                   ),
                   mainPanel(width=9,
                             h4(textOutput("dim_reduction_plot_title")), # Dynamic title added
                             plotlyOutput("plot_clustering", height = "600px")
                   )
                 )
        ),
        
        tabPanel(tagList(icon("share-alt"), "Network"), value="tab_network",
                 sidebarLayout(
                   sidebarPanel(width=3,
                                h4("Network Parameters"),
                                numericInput("network_top_peptides", "Use Top N Peptides:",
                                             value = 50, min = 10, max = 200, step = 10),
                                helpText("Build network using N peptides with highest average signal across the entire dataset."),
                                hr(),
                                numericInput("threshold", "Corr. Threshold (|ρ|):", value = 0.5, min=0, max=1, step = 0.05),
                                selectInput("layout_type", "Layout Algorithm:",
                                            choices = c("Fruchterman-Reingold" = "fr",
                                                        "Kamada-Kawai" = "kk",
                                                        "DrL (Large Graphs)" = "drl",
                                                        "Automatic" = "nicely"),
                                            selected = "fr"),
                                numericInput("seed", "Layout Random Seed:", value = 123, min = 1),
                                actionButton("btn_layout", "Generate Network", icon = icon("play"), class="btn-primary btn-block"),
                                hr(),
                                downloadButton("download_layout", "Download Layout Coords", icon=icon("download"), class="btn-success btn-sm btn-block")
                   ),
                   mainPanel(width=9,
                             h4("Peptide Correlation Network"),
                             helpText(icon("info-circle"), " Nodes = peptides (size ~ avg signal, color ~ Z-score: blue=low, red=high). Edges = absolute correlation > threshold. Based on entire dataset (or top N)."),
                             plotlyOutput("network_plot", height = "700px")
                   )
                 )
        ),
        
        tabPanel(tagList(icon("file-alt"), "Summary"), value="tab_summary",
                 h4(tagList(icon("list-ul"), "Dataset Overview")),
                 helpText("Basic counts from the loaded main dataset."),
                 verbatimTextOutput("text_data_summary"),
                 hr(),
                 h4(tagList(icon("signal"), "Top Peptides - Avg Intensity by Array")),
                 helpText("Average signal intensity for the top 10 peptides (by overall mean) across different arrays in the dataset."),
                 plotlyOutput("plot_top10_peptides_by_array", height = "500px")
        ),
        
        # ---- CONTACT TAB ----
        tabPanel(tagList(icon("envelope"), "About / Contact"), value = "tab_contact",
                 fluidRow(
                   tags$div(class = "col-sm-12", # Full width column
                            tags$div(class="custom-box box-primary",
                                     tags$div(class="custom-box-header",
                                              tags$h3(class="box-title", tagList(icon("users"), "Development Team & Support"))
                                     ),
                                     tags$div(class="custom-box-body",
                                              p("KinoViz was developed to facilitate kinome data analysis. We welcome feedback and contributions."),
                                              # Updated Developer Info - Structure using paragraphs or list
                                              h5(strong("Team:")),
                                              p(tags$b("Ehsan Saghapour"), tags$sup("1,2"), "(Developer), ", tags$b("Joshua C. Anderson"), tags$sup("3"), "(Data Collection), ", tags$b("Jake Chen"), tags$sup("1,2"), "(PI), ", tags$b("Christopher D. Willey"), tags$sup("3*"), "(PI)" ),
                                              
                                              h5(strong("Affiliations:")),
                                              tags$ol(
                                                tags$li("Department of Biomedical Informatics and Data Science, The University of Alabama at Birmingham, AL, US"),
                                                tags$li("Systems Pharmacology AI Research Center (SPARC), University of Alabama at Birmingham, AL, US"),
                                                tags$li("Department of Radiation Oncology, The University of Alabama at Birmingham, AL, US")
                                              ),
                                              hr(),
                                              h5(strong("Correspondence (*):")),
                                              p(strong("Christopher D. Willey, MD, PhD")),
                                              tags$address( # Use address tag for contact info
                                                "Department of Radiation Oncology, UAB", br(),
                                                "HSROC 2220B", br(),
                                                "619 19th St. South", br(),
                                                "Birmingham, AL 35249", br(),
                                                "Email:", tags$a(href = "mailto:cwilley@uabmc.edu?subject=KinoViz Inquiry", "cwilley@uabmc.edu")
                                              )
                                     ) # end box-body
                            ) # end custom-box
                   ), # end col
                   tags$div(class = "col-sm-12",
                            tags$div(class="custom-box box-success",
                                     tags$div(class="custom-box-header",
                                              tags$h3(class="box-title", tagList(icon("quote-right"), " Citation"))
                                     ),
                                     tags$div(class="custom-box-body",
                                              h5(strong("How to Cite KinoViz:")),
                                              p("If you use this application in your research, please cite:"),
                                              p(em("Saghapour, E., Anderson, J. C., Chen, J., & Willey, C. D. (2025). KinoViz: A User-Friendly Web Application for High-Throughput Kinome Profiling Analysis and Visualization in Cancer Research. [Provide URL or Journal Reference when available]."))
                                     ) # end box-body
                            ) # end custom-box
                   ) # end col
                 ) # End fluidRow
        ) # End tabPanel("Contact Us")
        
      ) # End tabsetPanel
    ) # End mainPanel
  ) # End sidebarLayout
) # End fluidPage


######################################################
#                    SERVER LOGIC                    #
######################################################
server <- function(input, output, session) {
  
  # --- Reactive Values ---
  kinaseData <- reactiveVal(NULL)
  agGridData <- reactiveVal(NULL)
  
  # --- Observers for Conditional UI ---
  observe({
    if (!is.null(kinaseData())) {
      updateNumericInput(session, "data_is_loaded_flag", value = 1)
    } else {
      updateNumericInput(session, "data_is_loaded_flag", value = 0)
    }
  }) # This observer depends on kinaseData()
  
  ### PART 1: Main Data Loading and Processing ###
  observeEvent(input$btn_process, {
    data_source <- NULL
    
    if (is.null(input$file_main) || input$file_main$datapath == "") {
      default_path <- "clean_kinome.csv"
      if (!file.exists(default_path)) {
        showModal(modalDialog(
          title = "Error",
          paste("Default dataset '", default_path, "' not found. Please upload a file or place it in the app directory."),
          easyClose = TRUE, footer = modalButton("Dismiss")
        ))
        kinaseData(NULL)
        agGridData(NULL)
        return()
      }
      data_source <- default_path
      showNotification("Loading default dataset...", type = "message", duration = 2)
    } else {
      data_source <- input$file_main$datapath
      showNotification(paste("Loading uploaded file:", input$file_main$name), type = "message", duration = 2)
    }
    
    processed <- process_data(data_source)
    
    if (!is.null(processed)) {
      kinaseData(processed)
      agGridData(processed$data2) # Initialize AG Grid with processed data
      showNotification("Data loaded and processed successfully!", type = "message", duration = 4)
    } else {
      kinaseData(NULL)
      agGridData(NULL)
      showNotification("Data loading/processing failed. Check console/logs and file format.", type = "error", duration = 5)
    }
  })
  
  ### PART 2: AG Grid Table Logic ###
  
  # Observer for AG Grid specific file upload
  observeEvent(input$file_upload_aggrid, {
    req(input$file_upload_aggrid)
    showNotification(paste("Processing AG Grid file:", input$file_upload_aggrid$name), type="message", duration=3)
    tryCatch({
      df_raw <- read.csv(input$file_upload_aggrid$datapath, header = FALSE, stringsAsFactors = FALSE)
      df_cleaned <- clean_uploaded_data(df_raw)
      agGridData(df_cleaned) 
      showNotification("AG Grid file cleaned and loaded.", type = "message", duration=3)
    }, error = function(e) {
      showNotification(paste("Error processing AG Grid file:", e$message), type = "error", duration = 10)
    })
  })
  
  # Observer for applying slope calculation to AG Grid data
  observeEvent(input$apply_slope, {
    current_data_in_grid <- agGridData() 
    req(current_data_in_grid) 
    
    showNotification("Applying slope correction...", type = "message", duration=2)
    
    tryCatch({
      updated_data <- calculate_slope_and_update_signal(current_data_in_grid)
      agGridData(updated_data) 
      showNotification("Slope correction applied successfully!", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Slope Correction Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # Render the AG Grid table
  output$aggrid_table <- renderAggrid({
    grid_data <- agGridData()
    req(grid_data)
    
    gridOptions <- list(
      defaultColDef = list(
        sortable = TRUE,
        filter = TRUE,
        resizable = TRUE,
        editable = FALSE 
      ),
      rowSelection = "multiple", 
      pagination = TRUE,         
      paginationPageSize = 50   
    )
    
    aggrid(
      data = grid_data,
      gridOptions = gridOptions,
      theme = "ag-theme-balham", 
      height = "600px",
      width = "100%"
    )
  })
  
  # Download Handler for AG Grid Data
  output$download_aggrid_data <- downloadHandler(
    filename = function() {
      paste0("aggrid_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(agGridData())
      write.csv(agGridData(), file, row.names = FALSE, quote = TRUE)
    }
  )
  
  
  ### PART 3: Dynamic UI Generation (based on loaded data) ###
  
  observe({
    kdata <- kinaseData() 
    req(kdata) 
    
    output$ui_Sample <- renderUI({
      req(kdata$Sample_choices)
      selectInput("Sample", "Select Sample(s)", choices = kdata$Sample_choices, multiple = TRUE, selected = kdata$Sample_choices[1])
    })
    
    output$ui_array <- renderUI({
      req(kdata$array_choices)
      selectInput("array", "Select Array(s)", choices = kdata$array_choices, multiple = TRUE, selected = kdata$array_choices[1])
    })
    
    output$ui_exposure <- renderUI({
      req(kdata$exposure_choices)
      selectInput("exposure_time", "Select Exposure Time",
                  choices = kdata$exposure_choices,
                  selected = kdata$exposure_choices[1], 
                  multiple = FALSE) 
    })
    
    output$ui_cycle <- renderUI({
      req(kdata$cycle_choices)
      selectInput("cycle", "Select Cycle",
                  choices = kdata$cycle_choices,
                  selected = kdata$cycle_choices[1], 
                  multiple = FALSE) 
    })
    
    # Peptide selection for Single Peptide plots (Kinetics)
    output$ui_peptide_single <- renderUI({
      req(kdata$peptide_choices)
      selectizeInput("peptide_single", "Select Peptide", choices = kdata$peptide_choices, selected=kdata$peptide_choices[1], options = list(placeholder = 'Type to search...'))
    })
    
    # Peptide selection for Multi-Peptide plots (Comparisons)
    output$ui_peptide_multi <- renderUI({
      req(kdata$peptide_choices)
      selectizeInput("peptide_multi", "Select Peptides",
                     choices = kdata$peptide_choices,
                     selected = kdata$peptide_choices[1:min(5, length(kdata$peptide_choices))], 
                     multiple = TRUE, options = list(placeholder = 'Type to search...'))
    })
    
    # ---- NEW: Peptide selection for Specific Peptides in Top Peptides Tab ----
    output$ui_specific_peptides_selector <- renderUI({
      req(kdata$peptide_choices) # Ensure peptide choices are loaded
      selectizeInput("specific_peptides_selected", "Select Specific Peptides:",
                     choices = kdata$peptide_choices,
                     selected = NULL, # Start with none selected
                     multiple = TRUE,
                     options = list(
                       placeholder = 'Type or click to select peptides...',
                       plugins = list('remove_button'), # Add clear button
                       maxItems = 50 # Limit for performance/visual clarity
                     ) 
      )
    })
    
  }) # This observer depends on kinaseData()
  
  
  ### PART 4: Server Logic for Analysis Tabs ###
  
  # --- Reactive Data Filtering ---
  filteredKinaseDataFixed <- reactive({
    kdata_list <- kinaseData()
    req(kdata_list, input$Sample, input$array, input$cycle, input$exposure_time) 
    
    kdata <- kdata_list$data2 
    
    kdata %>%
      filter(
        Sample %in% input$Sample,
        Array %in% input$array,
        as.numeric(Cycle) == as.numeric(input$cycle),
        as.numeric(Exposure_time) == as.numeric(input$exposure_time)
      )
  })
  
  filteredKinaseDataFlexible <- reactive({
    kdata_list <- kinaseData()
    req(kdata_list, input$Sample, input$array) 
    
    kdata <- kdata_list$data2
    
    kdata %>%
      filter(
        Sample %in% input$Sample,
        Array %in% input$array
      )
  })
  
  
  # --- Tab: Top Signal Peptides ---
  
  # Reactive to determine the list of peptides to show based on mode
  peptidesToShow <- reactive({
    kdata_list <- kinaseData()
    # Require data, Sample selection, and the mode selector
    req(kdata_list, input$Sample, input$peptide_selection_mode) 
    
    selected_Samples <- input$Sample
    peptide_cols_available <- kdata_list$peptide_choices
    kdata <- kdata_list$data2
    
    # Filter data once by selected Sample(s)
    filtered_data <- kdata %>% filter(Sample %in% selected_Samples)
    
    if (nrow(filtered_data) == 0) {
      showNotification("No data found for the selected Sample(s).", type="warning")
      return(list(peptides = character(0), data = filtered_data)) # Return empty list and filtered data
    }
    
    peptides_to_display <- character(0) # Initialize
    
    # --- Determine peptide list based on mode ---
    if (input$peptide_selection_mode == "top_n") {
      req(input$num_peptides) # Require num_peptides only in this mode
      num_top <- as.numeric(input$num_peptides)
      
      # Calculate max signal per peptide ACROSS selected Samples' data
      peptide_cols_in_data <- intersect(peptide_cols_available, colnames(filtered_data))
      if (length(peptide_cols_in_data) == 0) {
        showNotification("No peptide columns found in the filtered data.", type="warning")
        return(list(peptides = character(0), data = filtered_data))
      }
      
      max_values_per_peptide <- sapply(filtered_data[, peptide_cols_in_data, drop = FALSE], function(col) {
        numeric_col <- suppressWarnings(as.numeric(col))
        if(all(is.na(numeric_col))) return(NA_real_) else return(max(numeric_col, na.rm = TRUE))
      })
      
      max_values_per_peptide <- max_values_per_peptide[!is.na(max_values_per_peptide)]
      if(length(max_values_per_peptide) == 0) {
        showNotification("No valid signal values found for peptide ranking.", type="warning")
        return(list(peptides = character(0), data = filtered_data))
      }
      
      # Sort and get top N names
      sorted_peptides <- names(sort(max_values_per_peptide, decreasing = TRUE))
      peptides_to_display <- head(sorted_peptides, num_top)
      
    } else { # mode == "specific"
      # Use req() for specific_peptides_selected here for cleaner flow later
      req(input$specific_peptides_selected) 
      
      selected_specific <- input$specific_peptides_selected
      
      # Check if selection is NULL or empty (can happen if user clears the input)
      if (is.null(selected_specific) || length(selected_specific) == 0) {
        # Return empty list, let table/plot handle the message
        return(list(peptides = character(0), data = filtered_data)) 
      }
      
      # Validate selected peptides against available columns
      peptides_to_display <- intersect(selected_specific, peptide_cols_available)
      missing_peptides <- setdiff(selected_specific, peptides_to_display)
      if (length(missing_peptides) > 0) {
        showNotification(paste("Ignoring specifically selected peptides not found in data:", paste(missing_peptides, collapse=", ")), type="warning", duration=4)
      }
      
      if (length(peptides_to_display) == 0) {
        # Return empty list, let table/plot handle message
        return(list(peptides = character(0), data = filtered_data)) 
      }
    }
    
    return(list(peptides = peptides_to_display, data = filtered_data))
  })
  
  # Dynamic title for the table
  output$top_peptides_table_title <- renderText({
    req(input$peptide_selection_mode)
    ptp_info <- peptidesToShow() # Get results from the reactive
    # Need to handle the case where specific peptides are required but not yet selected
    # Need input$specific_peptides_selected to exist, even if NULL
    if(input$peptide_selection_mode == 'specific' && !("specific_peptides_selected" %in% names(input))) {
      return("Loading peptide selector...")
    }
    req(ptp_info) 
    
    num_peptides_found <- length(ptp_info$peptides)
    
    if (input$peptide_selection_mode == "top_n") {
      paste("Top", num_peptides_found, "Peptides by Max Signal")
    } else { # specific mode
      if (num_peptides_found > 0) {
        paste("Table for", num_peptides_found, "Selected Specific Peptide(s)")
      } else {
        # Check if the input exists and is empty/null vs not existing yet
        if (is.null(input$specific_peptides_selected) || length(input$specific_peptides_selected) == 0) {
          "Select Specific Peptides to Display Table"
        } else { # Input exists, had peptides, but none were valid/found
          "No Data for Specifically Selected Peptides"
        }
      }
    }
  })
  
  # Render the table
  output$table_important_peptides <- renderDT({
    ptp_info <- peptidesToShow() # Get results from the reactive
    # Need input$specific_peptides_selected to exist, even if NULL, in specific mode
    if(input$peptide_selection_mode == 'specific' && !("specific_peptides_selected" %in% names(input))) {
      return(datatable(data.frame(Message = "Loading peptide selector..."), options = list(searching = FALSE, paging = FALSE, info = FALSE)))
    }
    req(ptp_info) 
    
    peptides_to_display <- ptp_info$peptides
    filtered_data <- ptp_info$data # Use the pre-filtered data
    
    # Handle cases where no peptides should be displayed
    if (length(peptides_to_display) == 0) {
      msg <- if (input$peptide_selection_mode == "specific" && (is.null(input$specific_peptides_selected) || length(input$specific_peptides_selected) == 0)) {
        "Please select specific peptides using the control above."
      } else if (input$peptide_selection_mode == "specific") {
        "None of the selected specific peptides were found or had data for the selected Sample(s)."
      } else { # Top N mode
        "No peptides found based on the criteria (e.g., no data for Sample(s) or no valid signals)."
      }
      return(datatable(data.frame(Message = msg), options = list(searching = FALSE, paging = FALSE, info = FALSE)))
    }
    
    # Calculate max signal for the peptides_to_display within the filtered data
    peptide_cols_in_data <- intersect(peptides_to_display, colnames(filtered_data))
    if (length(peptide_cols_in_data) == 0) {
      return(datatable(data.frame(Message = "Selected peptides not found as columns in filtered data."), options = list(searching = FALSE, paging = FALSE, info = FALSE)))
    }
    
    max_values_per_peptide <- sapply(filtered_data[, peptide_cols_in_data, drop = FALSE], function(col) {
      numeric_col <- suppressWarnings(as.numeric(col))
      if(all(is.na(numeric_col))) return(NA_real_) else return(max(numeric_col, na.rm = TRUE))
    })
    
    # Create data frame for the table
    display_df <- data.frame(
      Peptide_ID = names(max_values_per_peptide),
      Max_Signal_in_Filter = max_values_per_peptide # Rename slightly for clarity
    ) 
    
    # Sort based on the mode
    if (input$peptide_selection_mode == "top_n") {
      display_df <- display_df %>% arrange(desc(Max_Signal_in_Filter))
    } else {
      # For specific peptides, order them as they appear in the input selection
      display_df <- display_df %>% 
        mutate(Peptide_ID = factor(Peptide_ID, levels = peptides_to_display)) %>% 
        arrange(Peptide_ID) %>%
        mutate(Peptide_ID = as.character(Peptide_ID)) # Convert back to character for DT
    }
    
    rownames(display_df) <- NULL
    display_df <- display_df %>% mutate(Max_Signal_in_Filter = round(Max_Signal_in_Filter, 1)) # Round for display
    
    datatable(display_df, options = list(pageLength = 10, searching=FALSE), rownames = FALSE)
  })
  
  # Dynamic title for the plot
  output$top_peptides_plot_title <- renderText({
    req(input$peptide_selection_mode)
    ptp_info <- peptidesToShow() # Get results from the reactive
    # Need input$specific_peptides_selected to exist, even if NULL
    if(input$peptide_selection_mode == 'specific' && !("specific_peptides_selected" %in% names(input))) {
      return("Loading peptide selector...")
    }
    req(ptp_info)
    
    num_peptides_found <- length(ptp_info$peptides)
    
    if (input$peptide_selection_mode == "top_n") {
      paste("Signal Distribution for Top", num_peptides_found, "Peptides")
    } else {
      if (num_peptides_found > 0) {
        paste("Signal Distribution for", num_peptides_found, "Selected Specific Peptide(s)")
      } else {
        if (is.null(input$specific_peptides_selected) || length(input$specific_peptides_selected) == 0) {
          "Select Specific Peptides to Display Plot"
        } else {
          "No Data for Specifically Selected Peptides"
        }
      }
    }
  })
  
  # Render the plot
  output$plot_ridgeline <- renderPlot({
    ptp_info <- peptidesToShow() # Get results from the reactive
    # Need input$specific_peptides_selected to exist, even if NULL
    if(input$peptide_selection_mode == 'specific' && !("specific_peptides_selected" %in% names(input))) {
      return(ggplot() + labs(title="Loading peptide selector...") + theme_void())
    }
    req(ptp_info) 
    
    top_peptides_names <- ptp_info$peptides
    filtered_data <- ptp_info$data # Use the pre-filtered data from the reactive
    selected_Samples <- isolate(input$Sample) # Get samples used for filtering 
    
    # Handle cases where no peptides should be displayed
    if (length(top_peptides_names) == 0) {
      title_msg <- if (input$peptide_selection_mode == "specific" && (is.null(input$specific_peptides_selected) || length(input$specific_peptides_selected) == 0)) {
        "Please select specific peptides using the control above."
      } else if (input$peptide_selection_mode == "specific") {
        "None of the selected specific peptides were found or had data for the selected Sample(s)."
      } else { # Top N mode
        "No peptides found based on the criteria (e.g., no data for Sample(s) or no valid signals)."
      }
      return(ggplot() + labs(title=title_msg) + theme_void())
    }
    
    # Ensure 'Array' exists
    if (!"Array" %in% colnames(filtered_data)) {
      return(ggplot() + labs(title="Error: 'Array' column missing in data.") + theme_void())
    }
    
    # Proceed with plotting logic
    peptide_cols_available <- intersect(top_peptides_names, colnames(filtered_data))
    if (length(peptide_cols_available) == 0) {
      return(ggplot() + labs(title="Selected peptides not found as columns in filtered data") + theme_void())
    }
    
    # Convert the filtered data from wide to long format
    df_long <- filtered_data %>%
      select(Sample, Array, all_of(peptide_cols_available)) %>%
      pivot_longer(cols = all_of(peptide_cols_available),
                   names_to = "Peptide",
                   values_to = "Signal",
                   values_transform = list(Signal = ~ suppressWarnings(as.numeric(.)))) %>%
      filter(!is.na(Signal)) # Remove NA signals before plotting
    
    # Check if data remains for plotting
    if (nrow(df_long) == 0) {
      return(ggplot() + labs(title="No signal data available for the selected peptides and samples") + theme_void())
    }
    
    # Ensure Peptide factor levels match the desired order for plotting
    # For specific peptides, use the order from top_peptides_names (which matches input order)
    # For top N, top_peptides_names is already sorted by max signal
    df_long <- df_long %>%
      mutate(Peptide = factor(Peptide, levels = top_peptides_names))
    
    # Determine plot title based on mode
    plot_main_title <- if (input$peptide_selection_mode == "top_n") {
      paste("Signal Distribution for Top", length(top_peptides_names), "Peptides")
    } else {
      paste("Signal Distribution for", length(top_peptides_names), "Selected Specific Peptides")
    }
    plot_subtitle <- paste("Filtered by Sample(s):", paste(selected_Samples, collapse=", "))
    
    # Plot density curves using standard geom_density and facet_wrap
    ggplot(df_long, aes(x = Signal, color = Array, fill = Array)) + 
      geom_density(alpha = 0.4, na.rm = TRUE) + 
      facet_wrap(~ Peptide, ncol = 2, scales = "free") + 
      labs(
        title = plot_main_title,
        subtitle = plot_subtitle,
        x = "Signal Intensity",
        y = "Density"
      ) +
      scale_fill_viridis_d(name = "Array") + # Add legend title
      scale_color_viridis_d(name = "Array") + # Add legend title
      theme_minimal(base_size = 11) + 
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"), 
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1) 
      )
    
  }, res = 100) # Control resolution
  
  
  # --- Tab: Kinetic Phosphorylation ---
  # Reactive expression to prepare data for kinetics plots
  kineticsPlotData <- reactive({
    kdata_list <- kinaseData()
    # Require necessary inputs depending on the mode
    req(kdata_list, input$peptide_single, input$phospho_mode)
    if(input$phospho_mode == "single") req(input$Sample, input$array)
    if(input$phospho_mode == "all") req(input$array)
    
    peptide <- input$peptide_single
    mode <- input$phospho_mode
    kdata <- kdata_list$data2
    
    # Ensure the selected peptide exists
    if (!peptide %in% colnames(kdata)) {
      showNotification(paste("Selected peptide", peptide, "not found in data."), type="error")
      return(NULL)
    }
    
    # Base filtering and selection
    # Handle potential errors during numeric conversion
    base_data <- kdata %>%
      select(Sample, Array, Cycle, Exposure_time, Signal = all_of(peptide)) %>%
      mutate(
        Cycle = suppressWarnings(as.numeric(as.character(Cycle))),
        Exposure_time = suppressWarnings(as.numeric(as.character(Exposure_time))),
        Signal = suppressWarnings(as.numeric(as.character(Signal)))
      ) %>%
      # Remove rows where conversion failed if necessary (optional)
      filter(!is.na(Cycle), !is.na(Exposure_time)) # Keep only rows with valid time/cycle
    
    # Mode-specific filtering
    if (mode == "single") {
      # Filter by selected Sample(s) and array(s)
      filtered <- base_data %>%
        filter(Sample %in% input$Sample, Array %in% input$array)
      
      return(list(data = filtered, mode = "single", peptide = peptide))
      
    } else { # mode == "all"
      # Filter only by selected array(s)
      # Cycle/Exposure filters from the main sidebar will be applied in the plot rendering
      filtered <- base_data %>% filter(Array %in% input$array)
      
      return(list(data = filtered, mode = "all", peptide = peptide))
    }
  })
  
  # Plot: Kinetics vs Cycle
  output$plot_cycle_kinetics <- renderPlotly({
    plot_info <- kineticsPlotData()
    req(plot_info)
    
    df <- plot_info$data
    mode <- plot_info$mode
    peptide <- plot_info$peptide
    
    if(nrow(df) == 0) return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = "No data for selected filters"))
    
    plot_title <- paste("Kinetic vs Cycle:", peptide)
    
    if (mode == "single") {
      # Single Sample Mode: Plot Cycle vs Signal, color by Exposure_time
      p <- plot_ly(df, x = ~Cycle, y = ~Signal, color = ~factor(Exposure_time), # Color by exposure
                   type = 'scatter', mode = 'lines+markers',
                   # Create a combined group identifier for lines if multiple Samples/arrays selected
                   linetype = ~paste(Sample, Array, sep="-"),
                   hoverinfo = 'text',
                   text = ~paste("Sample:", Sample, "<br>Array:", Array, "<br>Cycle:", Cycle,
                                 "<br>Exposure:", Exposure_time, "<br>Signal:", round(Signal, 2))) %>%
        layout(
          title = plot_title,
          xaxis = list(title = "Cycle"),
          yaxis = list(title = "Signal"),
          legend = list(title = list(text = 'Exposure Time'))
        )
      return(p)
      
    } else { # mode == "all"
      req(input$exposure_time) # Requires fixed exposure time for comparison across Samples
      
      fixed_exposure <- as.numeric(input$exposure_time)
      df_filtered <- df %>% filter(Exposure_time == fixed_exposure)
      
      if(nrow(df_filtered) == 0) return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = paste("No data for Exposure Time:", fixed_exposure)))
      
      p <- plot_ly(df_filtered, x = ~Cycle, y = ~Signal, color = ~Sample, # Color by Sample
                   type = 'scatter', mode = 'lines+markers',
                   linetype = ~Array, # Use line type for array if multiple selected
                   hoverinfo = 'text',
                   text = ~paste("Sample:", Sample, "<br>Array:", Array, "<br>Cycle:", Cycle,
                                 "<br>Exposure:", fixed_exposure, "<br>Signal:", round(Signal, 2))) %>%
        layout(
          title = plot_title,
          annotations = list(x = 0.5, y = 1.05, text = paste("Fixed Exposure Time:", fixed_exposure),
                             showarrow = F, xref='paper', yref='paper', xanchor='center', yanchor='bottom'),
          xaxis = list(title = "Cycle"),
          yaxis = list(title = "Signal"),
          legend = list(title = list(text = 'Sample'))
        )
      return(p)
    }
  })
  
  # Plot: Kinetics vs Exposure Time
  output$plot_exposure_kinetics <- renderPlotly({
    plot_info <- kineticsPlotData()
    req(plot_info)
    
    df <- plot_info$data
    mode <- plot_info$mode
    peptide <- plot_info$peptide
    
    if(nrow(df) == 0) return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = "No data for selected filters"))
    
    plot_title <- paste("Kinetic vs Exposure Time:", peptide)
    
    if (mode == "single") {
      # Single Sample Mode: Plot Exposure vs Signal, color by Cycle
      p <- plot_ly(df, x = ~Exposure_time, y = ~Signal, color = ~factor(Cycle), # Color by Cycle
                   type = 'scatter', mode = 'lines+markers',
                   linetype = ~paste(Sample, Array, sep="-"), # Group by Sample/array combo
                   hoverinfo = 'text',
                   text = ~paste("Sample:", Sample, "<br>Array:", Array, "<br>Cycle:", Cycle,
                                 "<br>Exposure:", Exposure_time, "<br>Signal:", round(Signal, 2))) %>%
        layout(
          title = plot_title,
          xaxis = list(title = "Exposure Time"),
          yaxis = list(title = "Signal"),
          legend = list(title = list(text = 'Cycle'))
        )
      return(p)
      
    } else { # mode == "all"
      req(input$cycle) # Requires fixed cycle for comparison across Samples
      
      fixed_cycle <- as.numeric(input$cycle)
      df_filtered <- df %>% filter(Cycle == fixed_cycle)
      
      if(nrow(df_filtered) == 0) return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = paste("No data for Cycle:", fixed_cycle)))
      
      p <- plot_ly(df_filtered, x = ~Exposure_time, y = ~Signal, color = ~Sample, # Color by Sample
                   type = 'scatter', mode = 'lines+markers',
                   linetype = ~Array, # Line type for array
                   hoverinfo = 'text',
                   text = ~paste("Sample:", Sample, "<br>Array:", Array, "<br>Cycle:", fixed_cycle,
                                 "<br>Exposure:", Exposure_time, "<br>Signal:", round(Signal, 2))) %>%
        layout(
          title = plot_title,
          annotations = list(x = 0.5, y = 1.05, text = paste("Fixed Cycle:", fixed_cycle),
                             showarrow = F, xref='paper', yref='paper', xanchor='center', yanchor='bottom'),
          xaxis = list(title = "Exposure Time"),
          yaxis = list(title = "Signal"),
          legend = list(title = list(text = 'Sample'))
        )
      return(p)
    }
  })
  
  # --- Tab: Peptide Comparisons ---
  # Reactive for multi-peptide comparison data
  peptideComparisonData <- reactive({
    kdata_list <- kinaseData()
    # Require inputs needed regardless of which plot is shown
    req(kdata_list, input$peptide_multi, input$Sample, input$array)
    
    kdata <- kdata_list$data2
    selected_peptides <- input$peptide_multi
    peptide_cols_available <- kdata_list$peptide_choices
    
    # Check peptides exist
    actual_peptides_to_use <- intersect(selected_peptides, peptide_cols_available)
    missing_peptides <- setdiff(selected_peptides, actual_peptides_to_use)
    
    if(length(missing_peptides) > 0) {
      showNotification(paste("Warning: Ignoring selected peptides not found in data:", paste(missing_peptides, collapse=", ")), type="warning", duration = 5)
    }
    if(length(actual_peptides_to_use) == 0) {
      showNotification("No valid peptides selected or found.", type="error")
      return(NULL)
    }
    
    # Filter by Sample/array and select relevant columns
    filtered <- kdata %>%
      filter(Sample %in% input$Sample, Array %in% input$array) %>%
      select(Sample, Array, Cycle, Exposure_time, all_of(actual_peptides_to_use)) %>%
      mutate(across(c(Cycle, Exposure_time, all_of(actual_peptides_to_use)), ~ suppressWarnings(as.numeric(as.character(.)))))
    
    # Pivot for plotting
    pivoted <- filtered %>%
      pivot_longer(cols = all_of(actual_peptides_to_use), names_to = "Peptide", values_to = "Signal") %>%
      # Create unique ID for each original sample row for grouping lines
      mutate(SampleID = paste(Sample, Array, sep="_"))
    
    return(pivoted)
  })
  
  # Plot: Peptides vs Cycle
  output$plot_peptides_vs_cycle <- renderPlotly({
    plot_data <- peptideComparisonData()
    req(plot_data, input$exposure_time) # Need pivoted data and fixed exposure
    
    fixed_exposure <- as.numeric(input$exposure_time)
    
    df_filtered <- plot_data %>% filter(Exposure_time == fixed_exposure)
    
    if(nrow(df_filtered) == 0) return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = paste("No data for Exposure Time:", fixed_exposure)))
    
    # Plot: Cycle on X, Signal on Y, Color by Peptide, Lines grouped by SampleID
    plot_ly(df_filtered, x = ~Cycle, y = ~Signal, color = ~Peptide,
            type = 'scatter', mode = 'lines+markers',
            # Use linetype to distinguish different Sample/Array combinations
            linetype = ~SampleID,
            hoverinfo = 'text',
            text = ~paste("Peptide:", Peptide, "<br>Sample:", SampleID, "<br>Cycle:", Cycle, "<br>Signal:", round(Signal,2))) %>%
      layout(
        title = paste("Selected Peptides vs. Cycle"),
        annotations = list(x = 0.5, y = 1.05, text = paste("Fixed Exposure Time:", fixed_exposure),
                           showarrow = F, xref='paper', yref='paper', xanchor='center', yanchor='bottom'),
        xaxis = list(title = 'Cycle'),
        yaxis = list(title = 'Signal Value'),
        legend = list(title = list(text = 'Peptide'))
      )
  })
  
  # Plot: Peptides vs Exposure Time
  output$plot_peptides_vs_exposure <- renderPlotly({
    plot_data <- peptideComparisonData()
    req(plot_data, input$cycle) # Need pivoted data and fixed cycle
    
    fixed_cycle <- as.numeric(input$cycle)
    
    df_filtered <- plot_data %>% filter(Cycle == fixed_cycle)
    
    if(nrow(df_filtered) == 0) return(plotly_empty(type = "scatter", mode = "markers") %>% layout(title = paste("No data for Cycle:", fixed_cycle)))
    
    # Plot: Exposure Time on X, Signal on Y, Color by Peptide, Lines grouped by SampleID
    plot_ly(df_filtered, x = ~Exposure_time, y = ~Signal, color = ~Peptide,
            type = 'scatter', mode = 'lines+markers',
            linetype = ~SampleID,
            hoverinfo = 'text',
            text = ~paste("Peptide:", Peptide, "<br>Sample:", SampleID, "<br>Exposure:", Exposure_time, "<br>Signal:", round(Signal,2))) %>%
      layout(
        title = paste("Selected Peptides vs. Exposure Time"),
        annotations = list(x = 0.5, y = 1.05, text = paste("Fixed Cycle:", fixed_cycle),
                           showarrow = F, xref='paper', yref='paper', xanchor='center', yanchor='bottom'),
        xaxis = list(title = 'Exposure Time'),
        yaxis = list(title = 'Signal Value'),
        legend = list(title = list(text = 'Peptide'))
      )
  })
  
  
  # --- Tab: Heatmap ---
  output$plot_heatmap <- renderPlotly({
    kdata_list <- kinaseData()
    # Require data and the fixed cycle/exposure/array filters
    req(kdata_list, input$exposure_time, input$cycle, input$array)
    
    kdata <- kdata_list$data2
    fixed_exposure <- as.numeric(input$exposure_time)
    fixed_cycle <- as.numeric(input$cycle)
    selected_arrays <- input$array # Use the selected arrays
    peptide_cols <- kdata_list$peptide_choices
    
    # Filter data by the selected Exposure_time, Cycle, and selected Array(s)
    filtered_data <- kdata %>%
      filter(
        suppressWarnings(as.numeric(Exposure_time)) == fixed_exposure,
        suppressWarnings(as.numeric(Cycle)) == fixed_cycle,
        Array %in% selected_arrays # Filter by selected arrays
      ) %>%
      # Create a unique Sample ID for global clustering
      mutate(SampleID = paste(Sample, Array, sep="_")) %>%
      # Convert peptide columns to numeric *after* filtering
      mutate(across(all_of(peptide_cols), ~ suppressWarnings(as.numeric(as.character(.)))))
    
    if(nrow(filtered_data) == 0 || length(peptide_cols) == 0) {
      return(plotly_empty() %>% layout(title = "No data available for the selected Cycle, Exposure Time, and Array(s)"))
    }
    
    # --- Step 1: Prepare GLOBAL Matrix for Clustering ---
    # Matrix: Samples (SampleID = Sample_Array) as rows, Peptides as columns
    global_matrix_data <- filtered_data %>%
      select(SampleID, all_of(peptide_cols)) %>%
      # Handle potential duplicate SampleIDs (if raw data had them) by averaging
      group_by(SampleID) %>%
      summarise(across(all_of(peptide_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      tibble::column_to_rownames("SampleID")
    
    # Check if enough samples/peptides exist globally
    if(nrow(global_matrix_data) < 2 || ncol(global_matrix_data) < 2) {
      return(plotly_empty() %>% layout(title = "Fewer than 2 samples or peptides found globally for clustering"))
    }
    
    # Replace NA signal values with 0 (or consider other imputation)
    global_matrix_data[is.na(global_matrix_data)] <- 0
    
    # Remove peptides (columns) that are all zero globally
    non_zero_peptides <- colnames(global_matrix_data)[colSums(global_matrix_data != 0) > 0]
    if(length(non_zero_peptides) < 2){
      return(plotly_empty() %>% layout(title = "Fewer than 2 peptides with non-zero signal remain globally for clustering"))
    }
    global_matrix_data <- global_matrix_data[, non_zero_peptides, drop = FALSE]
    
    # Remove samples (rows) that are all zero globally (can happen after NA->0)
    non_zero_samples <- rownames(global_matrix_data)[rowSums(global_matrix_data != 0) > 0]
    if(length(non_zero_samples) < 2){
      return(plotly_empty() %>% layout(title = "Fewer than 2 samples with non-zero signal remain globally for clustering"))
    }
    global_matrix_data <- global_matrix_data[non_zero_samples, , drop = FALSE]
    
    # --- Step 2: Global Normalization & Clustering ---
    # Normalize globally (0-100 per peptide column)
    global_normalized <- percentize(global_matrix_data)
    
    # Cluster Peptides Globally (columns)
    ordered_peptides <- tryCatch({
      peptide_dist <- dist(t(global_normalized))
      peptide_clust <- hclust(peptide_dist, method = "ward.D2")
      colnames(global_normalized)[peptide_clust$order]
    }, error = function(e){
      warning("Peptide clustering failed: ", e$message)
      colnames(global_normalized) # Fallback to original order if clustering fails
    })
    
    # --- REMOVE: Clustering of samples based on Sample is eliminated ---
    # Use the original sample order (Sample_Array) without additional hierarchical clustering
    ordered_samples <- rownames(global_normalized)
    
    # --- Step 3: Create Subplots using GLOBAL Order ---
    unique_arrays_in_data <- sort(intersect(selected_arrays, unique(filtered_data$Array)))
    heatmap_plots <- list()
    
    for(arr in unique_arrays_in_data) {
      # Identify which globally ordered samples belong to this array
      # Use endsWith to check SampleID format "Sample_Array"
      ordered_samples_for_arr <- ordered_samples[endsWith(ordered_samples, paste0("_", arr))]
      
      # Skip if no samples for this array remain after global filtering/ordering
      if(length(ordered_samples_for_arr) == 0) next
      
      # Extract data for this array, maintaining GLOBAL row/col order
      peptides_to_plot <- intersect(ordered_peptides, colnames(global_normalized))
      if(length(peptides_to_plot) < 2) next # Need at least 2 peptides to plot
      
      # Subset the globally normalized matrix
      final_matrix_for_arr <- global_normalized[ordered_samples_for_arr, peptides_to_plot, drop = FALSE]
      
      # Skip if matrix dimensions are invalid after subsetting
      if(nrow(final_matrix_for_arr) < 1 || ncol(final_matrix_for_arr) < 1) next
      
      # Determine y-axis label dynamically based on selected array
      current_ylab <- paste("Samples (", arr, ")", sep = "")
      
      # Generate heatmap for this array using heatmaply, WITHOUT internal clustering
      hmap <- tryCatch({
        heatmaply(
          final_matrix_for_arr, # Pass pre-ordered matrix
          main = paste("Array:", arr),  # Simpler title per subplot
          xlab = "", # No xlab per subplot, use overall
          ylab = current_ylab,   # Dynamic y-axis label based on array
          
          # --- CRITICAL: Disable internal clustering, use pre-calculated global order ---
          Rowv = FALSE,
          Colv = FALSE,
          seriate = "none", # Ensure no reordering
          
          # Appearance
          colors = colorRampPalette(c("blue", "white", "red"))(100),
          scale = "none", # Data is already normalized (0-100)
          limits = c(0, 100), # Consistent color scale
          
          # Dendrograms - DO NOT SHOW (as clustering is external)
          show_dendrogram = c(FALSE, FALSE),
          
          # Labels and Font Size
          fontsize_row = 8,
          fontsize_col = max(4, 10 - length(peptides_to_plot) %/% 20), # Dynamic font size
          margins = c(60, 100, 30, 30) # Adjusted margins 
        )
      }, error = function(e) {
        message("Error generating heatmap for Array ", arr, ": ", e$message)
        plotly_empty() %>% layout(title = paste("Error - Array:", arr))
      })
      
      heatmap_plots[[as.character(arr)]] <- hmap # Use array name as list index
    }
    
    # --- Step 4: Combine Plots ---
    if(length(heatmap_plots) == 0) {
      return(plotly_empty() %>% layout(title = "No valid heatmaps could be generated based on global clustering"))
    } else if (length(heatmap_plots) == 1) {
      # If only one heatmap, add overall title 
      return(heatmap_plots[[1]] %>%
               layout(title = paste("Heatmap (Peptides Globally Clustered) - Cycle:", fixed_cycle, ", Exposure:", fixed_exposure, "- Array:", names(heatmap_plots)[1]),
                      xaxis = list(title = "Peptides"))) # Add overall X axis title here
    } else {
      # Use subplot for multiple arrays; plots will be arranged vertically
      plots_in_order <- heatmap_plots[as.character(unique_arrays_in_data)]
      combined_plot <- subplot(
        plots_in_order,
        nrows = length(plots_in_order), # Stack vertically
        shareX = TRUE, # Share Peptide axis (order is global)
        shareY = FALSE, # Do not share Sample axis
        titleX = TRUE,  # Show shared X axis title at bottom
        titleY = TRUE,
        margin = 0.04
      ) %>% layout(
        title = paste("Peptide Signal Heatmaps by Array (Peptides Globally Clustered)\nCycle:",
                      fixed_cycle, ", Exposure:", fixed_exposure),
        xaxis = list(title = "Peptides"), # Overall X axis title
        height = max(600, 300 * length(plots_in_order)) # Dynamic height, slightly reduced multiplier
      )
      return(combined_plot)
    }
  })
  
  
  # --- Tab: Dimensionality Reduction ---
  # Reactive expression for calculating PCA/UMAP
  clustering_results <- eventReactive(input$btn_clustering, {
    kdata_list <- kinaseData()
    # Require data and the fixed cycle/exposure filters from main sidebar
    req(kdata_list, input$cycle, input$exposure_time)
    
    kdata <- kdata_list$data2
    fixed_cycle <- as.numeric(input$cycle)
    fixed_exposure <- as.numeric(input$exposure_time)
    peptide_cols <- kdata_list$peptide_choices
    method <- input$clustering_method
    
    showNotification(paste("Running", method, "..."), type="message", duration=2)
    
    # Filter data
    data_filtered <- kdata %>%
      filter(
        suppressWarnings(as.numeric(Cycle)) == fixed_cycle,
        suppressWarnings(as.numeric(Exposure_time)) == fixed_exposure
      )
    
    if(nrow(data_filtered) < 2) {
      showNotification("Need at least 2 samples (rows) for dimensionality reduction.", type="warning")
      return(NULL)
    }
    
    # Prepare numeric matrix: Samples as rows, Peptides as columns
    # Create unique sample ID for rownames
    data_filtered <- data_filtered %>% mutate(SampleID = make.unique(paste(Sample, Array, sep="_")))
    
    numeric_data <- data_filtered %>%
      select(SampleID, all_of(peptide_cols)) %>%
      tibble::column_to_rownames("SampleID") %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
    
    
    # Impute NAs (e.g., with 0) before scaling/reduction
    numeric_data[is.na(numeric_data)] <- 0
    
    # Remove columns with zero variance
    col_vars <- apply(numeric_data, 2, var, na.rm = TRUE)
    # Use a small tolerance for variance check
    cols_to_keep <- !is.na(col_vars) & col_vars > 1e-8
    if(sum(!cols_to_keep) > 0) {
      warning("Removing ", sum(!cols_to_keep) , " peptide columns with zero/low variance or all NAs.")
      numeric_data <- numeric_data[, cols_to_keep, drop = FALSE]
    }
    
    if(ncol(numeric_data) < 2) {
      showNotification("Need at least 2 peptides with variance.", type="warning")
      return(NULL)
    }
    
    
    # Scale data (center and scale each peptide column)
    # Scaling is important for PCA/UMAP unless features are already comparable
    scaled_data <- scale(numeric_data, center = TRUE, scale = TRUE)
    # Check for issues after scaling (e.g., all NaN if input was bad)
    if(all(is.na(scaled_data))) {
      showNotification("Data scaling resulted in NA/NaN values. Check input data.", type="error")
      return(NULL)
    }
    
    # Ensure rownames are preserved after scaling
    rownames(scaled_data) <- rownames(numeric_data)
    
    
    # Add sample identifiers back (Sample, Array) - ensure order matches scaled_data rows
    # Retrieve original Sample/Array using the SampleID rownames
    sample_info <- data_filtered %>% 
      select(SampleID, Sample, Array) %>%
      filter(SampleID %in% rownames(scaled_data)) %>%
      # Arrange sample_info to match the order of rownames in scaled_data
      arrange(match(SampleID, rownames(scaled_data)))
    
    
    # Perform dimensionality reduction
    tryCatch({
      if (method == "PCA") {
        # Run PCA on scaled data
        # center=FALSE, scale.=FALSE because data is already scaled
        pca_res <- prcomp(scaled_data, center=FALSE, scale.=FALSE)
        res_df <- data.frame(
          Dim1 = pca_res$x[,1],
          Dim2 = pca_res$x[,2]
        )
        # Get variance explained
        pc_variance <- summary(pca_res)$importance[2, 1:2] * 100
        axis_labels <- c(paste0("PC1 (", round(pc_variance[1], 1), "%)"),
                         paste0("PC2 (", round(pc_variance[2], 1), "%)"))
      } else { # UMAP
        # Run UMAP on scaled data
        # Add umap.pars potentially if needed (n_neighbors, min_dist etc.)
        umap_res <- umap(scaled_data)
        res_df <- data.frame(
          Dim1 = umap_res$layout[,1],
          Dim2 = umap_res$layout[,2]
        )
        axis_labels <- c("UMAP 1", "UMAP 2")
      }
      
      # Combine results with sample info
      # Ensure row counts match after filtering and arrangement
      if(nrow(sample_info) != nrow(res_df)){
        stop(paste("Mismatch between sample info (", nrow(sample_info), 
                   ") and reduction results (", nrow(res_df), ") dimensions. Rownames might be lost or mismatched." ,
                   "Sample Rownames:", paste(head(rownames(scaled_data)), collapse=", "),
                   "Info SampleIDs:", paste(head(sample_info$SampleID), collapse=", ")
        ))
      }
      results <- bind_cols(sample_info %>% select(-SampleID), res_df) # Drop SampleID after matching
      showNotification(paste(method, "calculation complete."), type="message", duration=3)
      
      return(list(method = method, data = results, labels = axis_labels))
      
    }, error = function(e) {
      showNotification(paste("Error during", method, ":", e$message), type="error", duration=10)
      return(NULL)
    })
  })
  
  # Plot: Dimensionality Reduction
  output$plot_clustering <- renderPlotly({
    results <- clustering_results()
    req(results) # Need results from the eventReactive
    
    df <- results$data
    method <- results$method
    axis_labels <- results$labels
    
    # Define color and symbol mapping based on 'Array'
    # Ensure 'Array' is treated as a factor/character for mapping
    df <- df %>% mutate(Array = as.character(Array))
    unique_arrays <- sort(unique(df$Array))
    
    # Colors: Use RColorBrewer Set1, expanding if necessary
    num_colors_needed <- length(unique_arrays)
    if (num_colors_needed <= 9) {
      # Handle case with only 1 or 2 arrays for brewer.pal minimum requirement
      array_colors <- brewer.pal(max(3, num_colors_needed), "Set1")[1:num_colors_needed]
    } else {
      # Expand palette if more than 9 arrays
      color_palette_func <- colorRampPalette(brewer.pal(9, "Set1"))
      array_colors <- color_palette_func(num_colors_needed)
    }
    array_color_map <- setNames(array_colors, unique_arrays)
    
    # Symbols: Define a set of symbols
    available_symbols <- c("circle", "square", "diamond", "cross", "x",
                           "triangle-up", "triangle-down", "pentagon", "star")
    num_symbols_needed <- length(unique_arrays)
    if (num_symbols_needed <= length(available_symbols)) {
      array_symbols <- available_symbols[1:num_symbols_needed]
    } else {
      # Repeat symbols if more arrays than symbols
      array_symbols <- rep(available_symbols, length.out = num_symbols_needed)
    }
    array_symbol_map <- setNames(array_symbols, unique_arrays)
    
    # Create the plot using mapped values
    p <- plot_ly(df, x = ~Dim1, y = ~Dim2,
                 type = 'scatter', mode = 'markers',
                 # Map aesthetics directly using the named vectors/lists
                 color = ~Array, colors = array_color_map,
                 symbol = ~Array, symbols = array_symbol_map,
                 # Define hover text
                 hoverinfo = 'text',
                 text = ~paste("Sample:", Sample, "<br>Array:", Array),
                 marker = list(size = 10, opacity=0.8) # Adjust marker size/opacity
    ) %>%
      layout(
        title = paste(method, "Analysis (Cycle:", input$cycle, ", Exposure:", input$exposure_time, ")"),
        xaxis = list(title = axis_labels[1]),
        yaxis = list(title = axis_labels[2]),
        legend = list(title = list(text = 'Array'))
      )
    
    return(p)
  })
  
  
  # --- Tab: Network Analysis ---
  # Reactive expression for generating network layout and plot
  network_results <- eventReactive(input$btn_layout, {
    kdata_list <- kinaseData()
    req(kdata_list) # Need main data
    
    kdata <- kdata_list$data2
    top_n <- as.numeric(input$network_top_peptides)
    
    showNotification("Generating network...", type="message", duration=2)
    
    # Use the generate_network_plot function
    # Pass the full data and let the function handle peptide selection and calculation
    tryCatch({
      net_output <- generate_network_plot(
        data = kdata,
        layout_type = input$layout_type,
        threshold = input$threshold,
        seed = input$seed,
        top_n_peptides = top_n # Pass N for top peptides
      )
      showNotification("Network generation complete.", type="message", duration=3)
      return(net_output)
    }, error = function(e) {
      showNotification(paste("Network Error:", e$message), type = "error", duration = 10)
      # Return NULL or a specific structure indicating error
      return(list(error = e$message))
    })
  })
  
  # Plot: Network
  output$network_plot <- renderPlotly({
    results <- network_results()
    # Check for error structure first
    if(!is.null(results$error)){
      return(plotly_empty() %>% layout(title=paste("Network Error:", results$error)))
    }
    # If no error, require the plot element
    req(results$plot)
    results$plot # Return the plotly object from the list
  })
  
  # Download: Network Layout
  output$download_layout <- downloadHandler(
    filename = function() {
      paste0("network_layout_", input$layout_type, "_thresh", input$threshold, ".csv")
    },
    content = function(file) {
      results <- network_results()
      # Ensure results are valid and not an error structure before downloading
      req(results, !is.null(results$layout), is.data.frame(results$layout)) # Make sure layout is a dataframe
      write.csv(results$layout, file, row.names = TRUE)
    }
  )
  
  
  # --- Tab: Data Summary ---
  output$text_data_summary <- renderText({
    kdata_list <- kinaseData()
    req(kdata_list)
    kdata <- kdata_list$data2 # Use data2 for counts
    
    # Safely get unique lengths, handle potential NULLs or empty choices
    n_Samples <- length(kdata_list$Sample_choices)
    n_arrays <- length(kdata_list$array_choices)
    n_exposures <- length(kdata_list$exposure_choices)
    n_cycles <- length(kdata_list$cycle_choices)
    n_peptides <- length(kdata_list$peptide_choices)
    n_rows <- nrow(kdata)
    
    paste(
      "--- Dataset Dimensions ---",
      paste("Total Rows:", format(n_rows, big.mark=",")),
      paste("Peptide Columns:", format(n_peptides, big.mark=",")),
      "\n--- Unique Metadata Values ---",
      paste("Samples:", n_Samples),
      paste("Arrays:", n_arrays),
      paste("Exposure Times:", n_exposures),
      paste("Cycles:", n_cycles),
      sep = "\n"
    )
  })
  
  # Plot: Top 10 Peptides by Global Mean, grouped by Array
  output$plot_top10_peptides_by_array <- renderPlotly({
    kdata_list <- kinaseData()
    req(kdata_list)
    kdata <- kdata_list$data2
    peptide_cols <- kdata_list$peptide_choices
    
    if(length(peptide_cols) == 0) return(plotly_empty() %>% layout(title="No peptide columns found in data"))
    
    # Calculate global mean for each peptide
    # Ensure columns are numeric before calculating mean
    numeric_peptide_data <- kdata[, peptide_cols, drop=FALSE] %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))
    
    peptide_means <- colMeans(numeric_peptide_data, na.rm = TRUE)
    
    # Identify top 10 peptides (handle cases with < 10 peptides)
    peptide_means <- peptide_means[!is.na(peptide_means)] # Remove NAs resulting from all-NA columns
    num_to_show <- min(10, length(peptide_means))
    if(num_to_show == 0) return(plotly_empty() %>% layout(title="No peptides with valid average signal"))
    
    top_peptides_names <- names(sort(peptide_means, decreasing = TRUE))[1:num_to_show]
    
    # Aggregate mean intensity by Array for these top peptides
    agg_data <- kdata %>%
      select(Array, all_of(top_peptides_names)) %>%
      # Convert to numeric within the aggregation pipeline
      mutate(across(all_of(top_peptides_names), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
      group_by(Array) %>%
      # Calculate mean, ensuring na.rm=TRUE
      summarize(across(all_of(top_peptides_names), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      pivot_longer(cols = all_of(top_peptides_names), names_to = "Peptide", values_to = "Mean_Intensity") %>%
      # Ensure Peptide order matches original top N order for plotting
      mutate(Peptide = factor(Peptide, levels = top_peptides_names))
    
    if(nrow(agg_data) == 0) return(plotly_empty() %>% layout(title="No aggregated data to plot"))
    
    # Create grouped bar chart
    plot_ly(agg_data,
            x = ~Peptide,
            y = ~Mean_Intensity,
            color = ~as.factor(Array), # Ensure Array is factor for discrete colors
            type = 'bar',
            hoverinfo = 'text',
            text = ~paste("Array:", Array, "<br>Peptide:", Peptide, "<br>Mean Signal:", round(Mean_Intensity, 2))) %>%
      layout(
        title = paste("Top", num_to_show, "Peptides: Average Intensity by Array"),
        barmode = 'group', # Group bars for the same peptide together
        xaxis = list(title = "Peptide", tickangle = -45),
        yaxis = list(title = "Mean Intensity"),
        legend = list(title = list(text = 'Array'))
      )
  })
  
} # End Server function

# Run the application
shinyApp(ui = ui, server = server)
