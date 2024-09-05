#################################################
###Author: Cássio Rampinelli & Saulo Aires
##################ANA-COMUC######################




############################################################################
###########################FIRST ATTEMPT TO GEV_MML#########################
############################################################################

#############################################################################
#################################XAVIER######################################
#############################################################################

#Clean R memory
rm(list=ls(all=TRUE))


# Load necessary libraries
library(sf)
library(dplyr)
library(ggplot2)
library(ggpattern)

##############################
#########FIGURA 21a a 25a#####
##############################




# Load necessary libraries
library(sf)
library(dplyr)
library(ggplot2)
library(ggpattern)

#Reading polygon shapefile
polygon_shapefile <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/xavier/RS_GEV_MML.shp")

#Set wd
setwd("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro")

#Reading points shapefile
points_shapefile <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GEV_MML_HIDRO.shp")


# Function to ensure CRS is assigned
assign_crs_if_missing <- function(sf_object, crs) {
  if (is.na(st_crs(sf_object))) {
    st_crs(sf_object) <- crs
  }
  return(sf_object)
}

# Define a CRS (for example, WGS 84)
crs_wgs84 <- 4326


# Assign CRS if missing
polygon_shapefile <- assign_crs_if_missing(polygon_shapefile, crs_wgs84)
points_shapefile <- assign_crs_if_missing(points_shapefile, crs_wgs84)

# Ensure the coordinate reference systems match
if (st_crs(polygon_shapefile) != st_crs(points_shapefile)) {
  points_shapefile <- st_transform(points_shapefile, st_crs(polygon_shapefile))
}

#Provide id references for points and polygons
points_shapefile$id<-seq(1:nrow(points_shapefile))
polygon_shapefile$id<-seq(1:nrow(polygon_shapefile))

# Check if each point is within or on the boundary of any polygon and assign polygon ID
points_shapefile$polygon_id <- sapply(1:nrow(points_shapefile), function(i) {
  point <- points_shapefile[i, ]
  containing_polygons <- which(st_intersects(point, polygon_shapefile, sparse = FALSE)[1, ])
  if (length(containing_polygons) > 0) {
    return(containing_polygons[1])
  } else {
    return(NA)
  }
})

# Identify which polygons are matched by at least one point
matched_polygons <- unique(points_shapefile$polygon_id[!is.na(points_shapefile$polygon_id)])
unmatched_polygons <- polygon_shapefile$id[!polygon_shapefile$id %in% matched_polygons]

# Create a new polygon layer with unmatched polygons
unmatched_polygons_layer <- polygon_shapefile %>%
  filter(id %in% unmatched_polygons)

# Create a new point layer with unmatched points
unmatched_points_layer <- points_shapefile %>%
  filter(is.na(polygon_id))

# Assign nearest polygon to unmatched points
# Compute distances between unmatched points and all polygons
distances <- st_distance(unmatched_points_layer, polygon_shapefile)

# Find the nearest polygon for each unmatched point
nearest_polygon_ids <- apply(distances, 1, function(row) {
  which.min(row)
})

# Update unmatched points with the nearest polygon ID
points_shapefile$polygon_id[is.na(points_shapefile$polygon_id)] <- nearest_polygon_ids

# Recalculate matched polygons and unmatched points
matched_polygons <- unique(points_shapefile$polygon_id[!is.na(points_shapefile$polygon_id)])
unmatched_polygons <- polygon_shapefile$id[!polygon_shapefile$id %in% matched_polygons]

# Create a new polygon layer with unmatched polygons
unmatched_polygons_layer <- polygon_shapefile %>%
  filter(id %in% unmatched_polygons)

# Create a new point layer with unmatched points
unmatched_points_layer <- points_shapefile %>%
  filter(is.na(polygon_id))


#Selecting the matched polygons for Xavier Grids
xavier.matched.polygons<-polygon_shapefile[matched_polygons,]



#######- Computing Return Periods by using IDFs equations

##-Hidro- Station

# Define the function to calculate IDF for a rainfall of 24 hs
calculate_rainfall <- function(idf_a, idf_b, idf_c,idf_d, Tr) {
  
  (idf_a*(Tr)^idf_b)/(24*60+idf_c)^idf_d
}

# User-specified values for Tr
Tr_values <- c(2,5,25,100,500)


# Apply the function for each value of Tr and store results
stations_results <- lapply(Tr_values, function(Tr) {
  points_shapefile %>%
    mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c,idf_d, Tr)) %>%
    mutate(Tr_values = as.character(Tr)) %>% 
    dplyr::select(codigo,lon, lat, rainfall, Tr_values)
})


# Combine results into a single dataframe
stations_results_df <- bind_rows(stations_results)

stations_results_df$Tr_values <- factor(stations_results_df$Tr_values, levels = as.character(Tr_values))  # Convert x_value to factor with proper levels


#Reading polygon shapefile
boundary_sf <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/GIS/RS_UF.shp")

# Assign CRS if missing
boundary_sf <- assign_crs_if_missing(boundary_sf, crs_wgs84)

min(stations_results_df$rainfall)
max(stations_results_df$rainfall)

# Define color and size scale based on the entire dataset
rainfall_range <- range(stations_results_df$rainfall, na.rm = TRUE)


#############################
##-Xavier Matched - Polygons#
#############################



# Apply the function for each value of Tr and store results
  xavier.matched.polygons_results <- lapply(Tr_values, function(Tr) {
    xavier.matched.polygons %>%
    mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c,idf_d, Tr)) %>%
    mutate(Tr_values = as.character(Tr)) %>% 
    dplyr::select(codigo,lon, lat, rainfall, Tr_values)
})


# Combine results into a single dataframe
xavier.matched.polygons_results_df <- bind_rows(xavier.matched.polygons_results)
xavier.matched.polygons_results_df$Tr_values <- factor(xavier.matched.polygons_results_df$Tr_values, levels = as.character(Tr_values))  # Convert x_value to factor with proper levels


xavier.matched.polygons_results_df
stations_results_df

library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)


# Create a base plot with polygons
p <- ggplot() +
  geom_sf(data = xavier.matched.polygons_results_df, aes(fill = rainfall), color = NA) +
  scale_fill_viridis_c(name = "Rainfall") +
  geom_sf(data = stations_results_df, color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Polygons and Stations",
       subtitle = "Polygons with rainfall data and stations",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "bottom")

# Print the plot
print(p)



# Calculate distances between each polygon and all stations
distances <- st_distance(xavier.matched.polygons_results_df, stations_results_df)

# Find the index of the closest station for each polygon
nearest_stations_indices <- apply(distances, 1, which.min)

# Extract the nearest stations based on the indices
nearest_stations <- stations_results_df[nearest_stations_indices, ]

nearest_stations
stations_results_df

# Combine the polygons with their nearest station data
matched_data <- xavier.matched.polygons_results_df %>%
  mutate(
    nearest_station = nearest_stations$codigo,
    nearest_rainfall = nearest_stations$rainfall
  )


# Assuming CRS is the same and distances have been computed correctly

# Convert the sf object to a data frame
matched_data_df <- st_drop_geometry(matched_data)

# Compute the percentage difference
matched_data_df <- matched_data_df %>%
  mutate(
    rainfall_diff =  nearest_rainfall-rainfall,
    percentage_diff = (rainfall_diff / nearest_rainfall) * 100
  )

# View the first few rows of the data frame
head(matched_data_df[c("codigo", "nearest_station", "rainfall", "nearest_rainfall", "rainfall_diff", "percentage_diff")])



# Create a list of unique Tr_values for plotting
unique_Tr_values <- unique(matched_data_df$Tr_values)






library(ggplot2)
library(dplyr)

# Create a directory to save the plots if it doesn't exist
dir.create("bubble_plots", showWarnings = FALSE)

# Loop through each unique Tr value
for (Tr in unique(matched_data_df$Tr_values)) {
  
  # Filter data for the current Tr value
  plot_data <- matched_data_df %>%
    filter(Tr_values == Tr)
  
  # Create the bubble plot
  plot <- ggplot(data = plot_data) +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
    geom_point(aes(x = lon, y = lat, size = abs(percentage_diff), color = percentage_diff), alpha = 0.7) +
    scale_size_continuous(name = "Diferença %", range = c(2, 30)) +
    scale_color_gradient(low = "lightblue", high = "red", name = "Diferença %") +
    labs(title = "XAVIER - Diferença percentual entre grade e estação pluviométrica",
         subtitle = paste("Tempo de Retorno =", Tr, "anos"),
         x = "Longitude",
         y = "Latitude") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  
  print(plot)
  
  # Save the plot
  #ggsave(filename = paste0("bubble_plots/Percentage_Difference_Tr=", Tr, ".png"), plot = plot, width = 10, height = 8, dpi = 500)
}




library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)

# Create a directory to save the plots if it doesn't exist
dir.create("bubble_plots", showWarnings = FALSE)

# Initialize a counter for filenames
counter <- 21

setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_21_25")

# Loop through each unique Tr value
for (Tr in unique(matched_data_df$Tr_values)) {
  
  # Filter data for the current Tr value
  plot_data <- matched_data_df %>%
    filter(Tr_values == Tr)
  
  # Create the bubble plot
  plot <- ggplot(data = plot_data) +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
    geom_point(aes(x = lon, y = lat, size = abs(percentage_diff), color = percentage_diff), alpha = 0.7) +
    scale_size_continuous(name = "Diferença %", range = c(1, 10), guide = "none") +  # Reduced bubble size and removed size legend
    scale_color_gradient(low = "red", high = "lightblue", name = "Diferença %") +  # Inverted color pattern
    labs(title = "XAVIER - Diferença percentual entre grade e estação pluviométrica",
         subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno =", Tr, "anos"),
         x = "Longitude",
         y = "Latitude") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  
  # Print the plot
  print(plot)
  
  # Save the plot with a sequential filename
  ggsave(filename = paste0("Figura", counter, "a.png"), plot = plot, width = 10, height = 8, dpi = 500)
  
  # Increment the counter
  counter <- counter + 1
}



#############################################################################
#################################CHIRPS######################################
#############################################################################

#Clean R memory
rm(list=ls(all=TRUE))


# Load necessary libraries
library(sf)
library(dplyr)
library(ggplot2)
library(ggpattern)

##############################
#########FIGURA 21b a 25b#####
##############################

#Reading polygon shapefile
polygon_shapefile <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/chirps/RS_GEV_MML.shp")

#Set wd
setwd("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro")

#Reading points shapefile
points_shapefile <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GEV_MML_HIDRO.shp")


# Function to ensure CRS is assigned
assign_crs_if_missing <- function(sf_object, crs) {
  if (is.na(st_crs(sf_object))) {
    st_crs(sf_object) <- crs
  }
  return(sf_object)
}

# Define a CRS (for example, WGS 84)
crs_wgs84 <- 4326


# Assign CRS if missing
polygon_shapefile <- assign_crs_if_missing(polygon_shapefile, crs_wgs84)
points_shapefile <- assign_crs_if_missing(points_shapefile, crs_wgs84)

# Ensure the coordinate reference systems match
if (st_crs(polygon_shapefile) != st_crs(points_shapefile)) {
  points_shapefile <- st_transform(points_shapefile, st_crs(polygon_shapefile))
}

#Provide id references for points and polygons
points_shapefile$id<-seq(1:nrow(points_shapefile))
polygon_shapefile$id<-seq(1:nrow(polygon_shapefile))

# Check if each point is within or on the boundary of any polygon and assign polygon ID
points_shapefile$polygon_id <- sapply(1:nrow(points_shapefile), function(i) {
  point <- points_shapefile[i, ]
  containing_polygons <- which(st_intersects(point, polygon_shapefile, sparse = FALSE)[1, ])
  if (length(containing_polygons) > 0) {
    return(containing_polygons[1])
  } else {
    return(NA)
  }
})

# Identify which polygons are matched by at least one point
matched_polygons <- unique(points_shapefile$polygon_id[!is.na(points_shapefile$polygon_id)])
unmatched_polygons <- polygon_shapefile$id[!polygon_shapefile$id %in% matched_polygons]

# Create a new polygon layer with unmatched polygons
unmatched_polygons_layer <- polygon_shapefile %>%
  filter(id %in% unmatched_polygons)

# Create a new point layer with unmatched points
unmatched_points_layer <- points_shapefile %>%
  filter(is.na(polygon_id))

# Assign nearest polygon to unmatched points
# Compute distances between unmatched points and all polygons
distances <- st_distance(unmatched_points_layer, polygon_shapefile)

# Find the nearest polygon for each unmatched point
nearest_polygon_ids <- apply(distances, 1, function(row) {
  which.min(row)
})

# Update unmatched points with the nearest polygon ID
points_shapefile$polygon_id[is.na(points_shapefile$polygon_id)] <- nearest_polygon_ids

# Recalculate matched polygons and unmatched points
matched_polygons <- unique(points_shapefile$polygon_id[!is.na(points_shapefile$polygon_id)])
unmatched_polygons <- polygon_shapefile$id[!polygon_shapefile$id %in% matched_polygons]

# Create a new polygon layer with unmatched polygons
unmatched_polygons_layer <- polygon_shapefile %>%
  filter(id %in% unmatched_polygons)

# Create a new point layer with unmatched points
unmatched_points_layer <- points_shapefile %>%
  filter(is.na(polygon_id))


#Selecting the matched polygons for Xavier Grids
chirps.matched.polygons<-polygon_shapefile[matched_polygons,]



#######- Computing Return Periods by using IDFs equations

##-Hidro- Station

# Define the function to calculate IDF for a rainfall of 24 hs
calculate_rainfall <- function(idf_a, idf_b, idf_c,idf_d, Tr) {
  
  (idf_a*(Tr)^idf_b)/(24*60+idf_c)^idf_d
}

# User-specified values for Tr
Tr_values <- c(2,5,25,100,500)


# Apply the function for each value of Tr and store results
stations_results <- lapply(Tr_values, function(Tr) {
  points_shapefile %>%
    mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c,idf_d, Tr)) %>%
    mutate(Tr_values = as.character(Tr)) %>% 
    dplyr::select(codigo,lon, lat, rainfall, Tr_values)
})


# Combine results into a single dataframe
stations_results_df <- bind_rows(stations_results)

stations_results_df$Tr_values <- factor(stations_results_df$Tr_values, levels = as.character(Tr_values))  # Convert x_value to factor with proper levels


#Reading polygon shapefile
boundary_sf <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/GIS/RS_UF.shp")

# Assign CRS if missing
boundary_sf <- assign_crs_if_missing(boundary_sf, crs_wgs84)

min(stations_results_df$rainfall)
max(stations_results_df$rainfall)

# Define color and size scale based on the entire dataset
rainfall_range <- range(stations_results_df$rainfall, na.rm = TRUE)


#############################
##-Chirps Matched - Polygons#
#############################



# Apply the function for each value of Tr and store results
chirps.matched.polygons_results <- lapply(Tr_values, function(Tr) {
  chirps.matched.polygons %>%
    mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c,idf_d, Tr)) %>%
    mutate(Tr_values = as.character(Tr)) %>% 
    dplyr::select(codigo,lon, lat, rainfall, Tr_values)
})


# Combine results into a single dataframe
chirps.matched.polygons_results_df <- bind_rows(chirps.matched.polygons_results)
chirps.matched.polygons_results_df$Tr_values <- factor(chirps.matched.polygons_results_df$Tr_values, levels = as.character(Tr_values))  # Convert x_value to factor with proper levels


chirps.matched.polygons_results_df
stations_results_df

library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)


# Create a base plot with polygons
p <- ggplot() +
  geom_sf(data = chirps.matched.polygons_results_df, aes(fill = rainfall), color = NA) +
  scale_fill_viridis_c(name = "Rainfall") +
  geom_sf(data = stations_results_df, color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Polygons and Stations",
       subtitle = "Polygons with rainfall data and stations",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "bottom")

# Print the plot
print(p)



# Calculate distances between each polygon and all stations
distances <- st_distance(chirps.matched.polygons_results_df, stations_results_df)

# Find the index of the closest station for each polygon
nearest_stations_indices <- apply(distances, 1, which.min)

# Extract the nearest stations based on the indices
nearest_stations <- stations_results_df[nearest_stations_indices, ]

nearest_stations
stations_results_df

# Combine the polygons with their nearest station data
matched_data <- chirps.matched.polygons_results_df %>%
  mutate(
    nearest_station = nearest_stations$codigo,
    nearest_rainfall = nearest_stations$rainfall
  )


# Assuming CRS is the same and distances have been computed correctly

# Convert the sf object to a data frame
matched_data_df <- st_drop_geometry(matched_data)

# Compute the percentage difference
matched_data_df <- matched_data_df %>%
  mutate(
    rainfall_diff =  nearest_rainfall-rainfall,
    percentage_diff = (rainfall_diff / nearest_rainfall) * 100
  )

# View the first few rows of the data frame
head(matched_data_df[c("codigo", "nearest_station", "rainfall", "nearest_rainfall", "rainfall_diff", "percentage_diff")])



# Create a list of unique Tr_values for plotting
unique_Tr_values <- unique(matched_data_df$Tr_values)

library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)

# Initialize a counter for filenames
counter <- 21

setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_21_25")

# Loop through each unique Tr value
for (Tr in unique(matched_data_df$Tr_values)) {
  
  # Filter data for the current Tr value
  plot_data <- matched_data_df %>%
    filter(Tr_values == Tr)
  
  # Create the bubble plot
  plot <- ggplot(data = plot_data) +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
    geom_point(aes(x = lon, y = lat, size = abs(percentage_diff), color = percentage_diff), alpha = 0.7) +
    scale_size_continuous(name = "Diferença %", range = c(1, 10), guide = "none") +  # Reduced bubble size and removed size legend
    scale_color_gradient(low = "red", high = "lightblue", name = "Diferença %") +  # Inverted color pattern
    labs(title = "CHIRPS - Diferença percentual entre grade e estação pluviométrica",
         subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno =", Tr, "anos"),
         x = "Longitude",
         y = "Latitude") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_line(color = "grey90")
    )
  
  # Print the plot
  print(plot)
  
  # Save the plot with a sequential filename
  ggsave(filename = paste0("Figura", counter, "b.png"), plot = plot, width = 10, height = 8, dpi = 500)
  
  # Increment the counter
  counter <- counter + 1
}

