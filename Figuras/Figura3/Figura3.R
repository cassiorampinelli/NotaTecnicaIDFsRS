#################################################
###Author: Cássio Rampinelli & Saulo Aires
##################ANA-COMUC######################

#Clean R memory
rm(list=ls(all=TRUE))


# Load necessary libraries
library(sf)
library(dplyr)
library(ggplot2)
library(ggpattern)


############################################################################
###########################FIRST ATTEMPT TO GEV_MML#########################
############################################################################

############################################################################
####################################FIGURAS 3 A 7###########################
############################################################################

#Clean R memory
rm(list=ls(all=TRUE))

#Set wd
setwd("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/chirps")

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



#Provide id references for points and polygons
points_shapefile$id<-seq(1:nrow(points_shapefile))
polygon_shapefile$id<-seq(1:nrow(polygon_shapefile))


# Assign CRS if missing
polygon_shapefile <- assign_crs_if_missing(polygon_shapefile, crs_wgs84)
points_shapefile <- assign_crs_if_missing(points_shapefile, crs_wgs84)

# Ensure the coordinate reference systems match
if (st_crs(polygon_shapefile) != st_crs(points_shapefile)) {
  points_shapefile <- st_transform(points_shapefile, st_crs(polygon_shapefile))
}

# Provide id references for points and polygons
points_shapefile$id <- seq(1:nrow(points_shapefile))
polygon_shapefile$id <- seq(1:nrow(polygon_shapefile))


###########
#HIDRO#####FIGURA 3a, 4a, 5a, 6a e 7a
###########




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

setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura3")



# Create and save plots for each Tr_value
for (Tr in Tr_values) {
  # Filter data for current Tr_value
  plot_data <- stations_results_df %>% filter(Tr_values ==  as.character(Tr))
  
  p <- ggplot(data = plot_data) +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
    geom_sf(aes(size = rainfall, color = rainfall), show.legend = TRUE) +
    scale_size_continuous(name = "Intensidade (mm/h)",guide="none") +
    scale_color_gradient(low = "lightblue", high = "red", name = "Intensidade (mm/h)") +
    labs(
      title = paste("HIDRO. Chuvas com duração de 24 hs. Tempo de Retorno =", Tr,"anos"),
      x = "Longitude",
      y = "Latitude"
    ) +
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
  
  print(p)
  # Save plot
  ggsave(filename = paste0("HIDRO_24hs_Tr=", Tr, ".png"), plot = p, width = 10, height = 8, dpi = 500)
}


###########
#XAVIER####FIGURA 3b, 4b, 5b,6b,7b
###########



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


#Bublle Plot

#######- Computing Return Periods by using IDFs equations


##-XAVIER SELECTED GRID ELEMENTS

#Selecting the matched grids elements
matched_elements_xavier<-polygon_shapefile$id %in% matched_polygons
matched_polygons_xavier<-polygon_shapefile[matched_elements_xavier, ]


# Apply the function for each value of Tr and store results
xavier_results <- lapply(Tr_values, function(Tr) {
  matched_polygons_xavier %>%
    mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c,idf_d, Tr)) %>%
    mutate(Tr_values = as.character(Tr)) %>% 
    dplyr::select(codigo,lon, lat, rainfall, Tr_values)
})


# Combine results into a single dataframe
xavier_results_df <- bind_rows(xavier_results)

xavier_results_df$Tr_values <- factor(xavier_results_df$Tr_values, levels = as.character(Tr_values))  # Convert x_value to factor with proper levels


# Create and save plots for each Tr_value
for (Tr in Tr_values) {
  # Filter data for current Tr_value
  plot_data <- xavier_results_df %>% filter(Tr_values ==  as.character(Tr))
  
  # Create plot
  p <- ggplot() +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed",linewidth=1) +  # Plot the boundary
    geom_point(data = plot_data, aes(x = lon, y = lat, size = rainfall,  color=rainfall), alpha = 0.7) +
    scale_size_continuous(name = "Intensidade (mm/h)",guide="none") +
    scale_color_gradient(low = "lightblue", high = "red", name = "Intensidade (mm/h)") +
    labs(
      title = paste("XAVIER. Chuvas com duração de 24 hs.Tempo de Retorno =", Tr,"anos"),
      x = "Longitude",
      y = "Latitude"
      ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)  # Centralize the title
    )
  
  print(p)
  # Save plot
  ggsave(filename = paste0("XAVIER_Grid_24hs_Tr=", Tr, ".png"), plot = p, width = 10, height = 8, dpi = 500)
}


###########
#CHIRPS####FIGURA 3c, 4c, 5c,6c,7c
###########



#######- Identifying the polygons that match with the gauging stations



#Set wd
setwd("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/chirps")

#Reading polygon shapefile
polygon_shapefile <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/chirps/RS_GEV_MML.shp")


#Set wd
setwd("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro")

#Reading polygon shapefile
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



#######- Computing Return Periods by using IDFs equations

##-CHIRPS SELECTED GRID ELEMENTS

#Selecting the matched grids elements
matched_elements_chirps<-polygon_shapefile$id %in% matched_polygons
matched_polygons_chirps<-polygon_shapefile[matched_elements_chirps, ]


# Apply the function for each value of Tr and store results
chirps_results <- lapply(Tr_values, function(Tr) {
  matched_polygons_chirps %>%
    mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c,idf_d, Tr)) %>%
    mutate(Tr_values = as.character(Tr)) %>% 
    dplyr::select(lon, lat, rainfall, Tr_values)
})


# Combine results into a single dataframe
chirps_results_df <- bind_rows(chirps_results)

chirps_results_df$Tr_values <- factor(chirps_results_df$Tr_values, levels = as.character(Tr_values))  # Convert x_value to factor with proper levels

min(chirps_results_df$rainfall)
max(chirps_results_df$rainfall)

setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura3")

# Create and save plots for each Tr_value
for (Tr in Tr_values) {
  # Filter data for current Tr_value
  plot_data <- chirps_results_df %>% filter(Tr_values ==  as.character(Tr))
  
  # Create plot
  p <- ggplot() +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed",linewidth=1) +  # Plot the boundary
    geom_point(data = plot_data, aes(x = lon, y = lat, size = rainfall,  color=rainfall), alpha = 0.7) +
    scale_size_continuous(name = "Intensidade (mm/h)",guide="none") +
    scale_color_gradient(low = "lightblue", high = "red", name = "Intensidade (mm/h)") +
    
    
    labs(
      title = paste("CHIRPS. Chuvas com duração de 24 hs.Tempo de Retorno =", Tr,"anos"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)  # Centralize the title
    )
  
  print(p)
  # Save plot
  ggsave(filename = paste0("CHIRPS_Grid_24hs_Tr=", Tr, ".png"), plot = p, width = 10, height = 8, dpi = 500)
}



