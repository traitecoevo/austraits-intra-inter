plot_australia_map <- function(data,
                               facet_by_trait = TRUE,
                               trait_to_colour = NULL,
                               point_size = 2,
                               alpha = 0.6) {
  # Create a map(s) of Australia with data points overlaid
  #
  # Parameters:
  # - data: data frame containing decimal_longitude, decimal_latitude, and trait_name columns
  # - facet_by_trait: if TRUE, create faceted plots for each focal trait (default: TRUE)
  # - trait_to_colour: single trait name to plot (only used if facet_by_trait = FALSE)
  # - point_size: size of points on map
  # - alpha: transparency of points
  #
  # Returns: ggplot object of Australia map(s) with data points
  
  # Get Australia coastline using rnaturalearth
  australia <- rnaturalearth::ne_countries(country = "Australia", scale = 50, returnclass = "sf")
  
  # Remove na growth form
  data <- data %>% drop_na(woody)
  # Define trait labels using text notation that renders reliably
  trait_labels <- c(
    huber_value = "SA:LA",
    leaf_area = "Leaf area",
    leaf_delta13C = "Ci:Ca",
    leaf_mass_per_area = "LMA",
    leaf_N_per_area = "Leaf N per area",
    wood_density = "Wood density",
    seed_dry_mass = "Seed mass"
  )
  
  # Create labeller
  trait_labeller <- ggplot2::labeller(trait_name = ggplot2::as_labeller(trait_labels))
  
  # Prepare data - remove any NA coordinates
  data_map <- data |>
    drop_na(decimal_longitude, decimal_latitude)
  
  # Define focal traits
  focal_traits <- c("leaf_delta13C", "leaf_N_per_area", "leaf_mass_per_area",
                    "leaf_area", "huber_value", "wood_density", "seed_dry_mass")
  
  if (facet_by_trait) {
    # Filter to focal traits only
    data_map <- data_map |>
      filter(trait_name %in% focal_traits) |>
      mutate(trait_name = factor(trait_name, levels = focal_traits)) |>
      distinct(trait_name, decimal_longitude, decimal_latitude, woody, .keep_all = FALSE)
    # Compute per-trait, per-growth-form observation and taxa counts
    label_data <- data |>
      drop_na(decimal_longitude, decimal_latitude) |>
      filter(trait_name %in% focal_traits) |>
      group_by(trait_name, woody) |>
      summarise(
        n_obs = n(),
        n_taxa = n_distinct(taxon_name),
        .groups = "drop"
      ) |>
      mutate(
        trait_name = factor(trait_name, levels = focal_traits),
        label = paste0("n = ", n_obs, " (", n_taxa, ")")
      ) |>
      arrange(trait_name, desc(woody)) |>  # woody first (bottom), non-woody above
      group_by(trait_name) |>
      mutate(y_pos = -43.5 + (row_number() - 1) * 2.5) |>
      ungroup()
    
    # Create faceted map
    p <- ggplot() +
      geom_sf(data = australia, fill = "light grey", colour = "black", linewidth = 0.3) +
      geom_point(
        data = data_map,
        aes(x = decimal_longitude, y = decimal_latitude, colour = woody),
        size = point_size,
        alpha = alpha
      ) +
      geom_text(
        data = label_data,
        aes(x = 113, y = y_pos, label = label, colour = woody),
        hjust = 0,
        vjust = 0,
        size = 6,
        fontface = "bold",
        show.legend = FALSE
      ) +
      coord_sf(xlim = c(112, 154), ylim = c(-44, -10), expand = FALSE) +
      scale_colour_okabeito(order = c(3, 6), labels = c("non-woody" = "Non-woody", "woody" = "Woody")) +
      facet_wrap(~trait_name, ncol = 3, labeller = trait_labeller) +
      theme_minimal() +
      labs(
        x = "Longitude",
        y = "Latitude",
        colour = "Growth form"
      ) +
      theme(
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17)
      )
    
  } else {
    # Single trait plot
    if (!is.null(trait_to_colour)) {
      data_map <- data_map |>
        filter(trait_name == trait_to_colour)
    }
    
    # Get unique lat/lon combinations
    data_map <- data_map |>
      distinct(decimal_longitude, decimal_latitude, woody, .keep_all = FALSE)
    
    # Get trait label if specified
    trait_label_text <- if (!is.null(trait_to_colour) && trait_to_colour %in% names(trait_labels)) {
      as.character(trait_labels[[trait_to_colour]])
    } else if (!is.null(trait_to_colour)) {
      trait_to_colour
    } else {
    }
    
    # Create single map
    p <- ggplot() +
      geom_sf(data = australia, fill = "light grey", colour = "black", linewidth = 0.3) +
      geom_point(
        data = data_map,
        aes(x = decimal_longitude, y = decimal_latitude, colour = woody),
        size = point_size,
        alpha = alpha
      ) +
      coord_sf(xlim = c(112, 154), ylim = c(-44, -10), expand = FALSE) +
      scale_colour_okabeito(order = c(3, 6), labels = c("non-woody" = "Non-woody", "woody" = "Woody")) +
      theme_minimal() +
      labs(
        x = "Longitude",
        y = "Latitude",
        colour = "Growth form"
      ) +
      theme(
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)
      )
  }
  
  return(p)
}

