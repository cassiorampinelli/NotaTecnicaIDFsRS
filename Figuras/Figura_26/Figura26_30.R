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


###########################################################################
#######################################TR=2#################################
############################################################################
####################################FIGURA 26###############################


####################
#########FIGURA 10a#
####################

###IDW



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


library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)




############################################################
## Creating Grid Size with the same size of XAVIER resolution
############################################################


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height),crs=crs(polygon_shapefile))

#Convert raster to points
#raster_sites_points <- rasterToPoints(rasterized_precip_4326, spatial = TRUE)
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# 2. Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(polygon_shapefile))

plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(2))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

###############################################################################
###############################################################################
#####################################IDW#######################################
###############################################################################
###############################################################################

# Interpolate the grid cells using a power value of 2 (idp=2.5)
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
P.idw <- st_transform(P.idw, crs = st_crs(boundary_sf_transformed))

# Crop the IDW results using the boundary
P.idw_cropped <- st_intersection(P.idw, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
P.idw_cropped <- P.idw_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot<-ggplot() +
  # Tile layer for IDW interpolation
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for IDW interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade IDW") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação IDW - Chuva com duração de 24 h. Tempo de Retorno = 2 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura10a.png"), plot = plot, width = 10, height = 8, dpi = 500)



#RMSE e MAPE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Perform IDW interpolation
  idw_result <- idw(rainfall ~ 1, train_data, test_data, idp = 2.5)
  
  # Store the prediction
  IDW.out[i] <- idw_result$var1.pred
}

# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape



metrics.data.frame<-cbind(points_shapefile,rmse)
metrics.data.frame<-cbind(metrics.data.frame,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame


# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "IDW - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 2 anos. MAPE médio: 8.471%"),
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

counter=26
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "a.png"), plot = plot, width = 10, height = 8, dpi = 500)

################################################################################
################################################################################
###################################KRIGING######################################
################################################################################
################################################################################

# Check and transform CRS if necessary
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Verify CRS consistency
if (st_crs(plot_data_sf) != st_crs(sf_grid_xavier)) {
  stop("CRS mismatch between plot_data_sf and sf_grid_xavier")
}

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define variogram models with specific parameters
v_linear_model <- vgm(psill = 0.8, model = "Lin", range = 500000, nugget = 0.12)
v_sph_model <- vgm(psill = 0.8, model = "Sph", range = 500000, nugget = 0.12)
v_exp_model <- vgm(psill = 0.8, model = "Exp", range = 500000, nugget = 0.12)
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit different variogram models
fit_linear <- fit.variogram(v_empirical, model = v_linear_model)
fit_sph <- fit.variogram(v_empirical, model = v_sph_model)
fit_exp <- fit.variogram(v_empirical, model = v_exp_model)
fit_gauss <- fit.variogram(v_empirical, model = v_gau_model)

# Plot the empirical variogram and the fitted models
plot(v_empirical, model = fit_linear, main = "Fitted Linear Variogram")
plot(v_empirical, model = fit_sph, main = "Fitted Spherical Variogram")
plot(v_empirical, model = fit_exp, main = "Fitted Exponential Variogram")
plot(v_empirical, model = fit_gauss, main = "Fitted Gaussian Variogram")


#Performing Kriging with the the gaussian variogram


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height), crs = crs(polygon_shapefile))

# Convert raster to points
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(raster_grid_xavier))

plot_data <- stations_results_df %>% filter(Tr_values == as.character(2))

# Transform stations points to the same CRS as sf_grid_xavier
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, plot_data_sf)

# Define the Gaussian variogram model
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit the Gaussian variogram model to the empirical variogram
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)

# Perform Ordinary Kriging interpolation
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)

# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))

# Crop the Kriging results using the boundary
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for Kriging interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade Kriging") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação OK - Chuva com duração de 24 h. Tempo de Retorno = 2 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)





#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

###########################
####LEAVE ONE OUT#########
##########################

# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Prepare empirical variogram using training data
  v_empirical <- variogram(rainfall ~ 1, train_data)
  
  # Define the Gaussian variogram model
  v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
  
  # Fit the Gaussian variogram model to the empirical variogram
  fit_gau <- fit.variogram(v_empirical, model = v_gau_model)
  
  # Perform Ordinary Kriging interpolation
  kriging_result <- krige(rainfall ~ 1, train_data, newdata = test_data, model = fit_gau)
  
  # Store the prediction
  kriging_out[i] <- kriging_result$var1.pred
}



# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape



metrics.data.frame.kriging<-cbind(points_shapefile,rmse)
metrics.data.frame.kriging<-cbind(metrics.data.frame.kriging,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.kriging


# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "OK - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 2 anos. MAPE Médio: 8.214%"),
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

counter=26
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "b.png"), plot = plot, width = 10, height = 8, dpi = 500)








###########################################################################
#######################################TR=5#################################
############################################################################
####################################FIGURA 27###############################

counter=27
####################
#########FIGURA 10a#
####################

###IDW



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


library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)




############################################################
## Creating Grid Size with the same size of XAVIER resolution
############################################################


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height),crs=crs(polygon_shapefile))

#Convert raster to points
#raster_sites_points <- rasterToPoints(rasterized_precip_4326, spatial = TRUE)
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# 2. Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(polygon_shapefile))

plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(5))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

###############################################################################
###############################################################################
#####################################IDW#######################################
###############################################################################
###############################################################################

# Interpolate the grid cells using a power value of 5 (idp=2.5)
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
P.idw <- st_transform(P.idw, crs = st_crs(boundary_sf_transformed))

# Crop the IDW results using the boundary
P.idw_cropped <- st_intersection(P.idw, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
P.idw_cropped <- P.idw_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot<-ggplot() +
  # Tile layer for IDW interpolation
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for IDW interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade IDW") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação IDW - Chuva com duração de 24 h. Tempo de Retorno = 5 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")



#RMSE e MAPE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Perform IDW interpolation
  idw_result <- idw(rainfall ~ 1, train_data, test_data, idp = 2.5)
  
  # Store the prediction
  IDW.out[i] <- idw_result$var1.pred
}

# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame<-0
metrics.data.frame<-cbind(plot_data,rmse)
metrics.data.frame<-cbind(metrics.data.frame,mape)



# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "IDW - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 5 anos. MAPE Médio: 7.265%"),
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

counter=27
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "a.png"), plot = plot, width = 10, height = 8, dpi = 500)

################################################################################
################################################################################
###################################KRIGING######################################
################################################################################
################################################################################

# Check and transform CRS if necessary
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Verify CRS consistency
if (st_crs(plot_data_sf) != st_crs(sf_grid_xavier)) {
  stop("CRS mismatch between plot_data_sf and sf_grid_xavier")
}

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define variogram models with specific parameters
v_linear_model <- vgm(psill = 0.8, model = "Lin", range = 500000, nugget = 0.12)
v_sph_model <- vgm(psill = 0.8, model = "Sph", range = 500000, nugget = 0.12)
v_exp_model <- vgm(psill = 0.8, model = "Exp", range = 500000, nugget = 0.12)
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit different variogram models
fit_linear <- fit.variogram(v_empirical, model = v_linear_model)
fit_sph <- fit.variogram(v_empirical, model = v_sph_model)
fit_exp <- fit.variogram(v_empirical, model = v_exp_model)
fit_gauss <- fit.variogram(v_empirical, model = v_gau_model)

# Plot the empirical variogram and the fitted models
plot(v_empirical, model = fit_linear, main = "Fitted Linear Variogram")
plot(v_empirical, model = fit_sph, main = "Fitted Spherical Variogram")
plot(v_empirical, model = fit_exp, main = "Fitted Exponential Variogram")
plot(v_empirical, model = fit_gauss, main = "Fitted Gaussian Variogram")


#Performing Kriging with the the gaussian variogram


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height), crs = crs(polygon_shapefile))

# Convert raster to points
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(raster_grid_xavier))

plot_data <- stations_results_df %>% filter(Tr_values == as.character(5))

# Transform stations points to the same CRS as sf_grid_xavier
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, plot_data_sf)

# Define the Gaussian variogram model
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit the Gaussian variogram model to the empirical variogram
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)

# Perform Ordinary Kriging interpolation
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)

# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))

# Crop the Kriging results using the boundary
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for Kriging interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade Kriging") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação OK - Chuva com duração de 24 h. Tempo de Retorno = 5 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)





#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

###########################
####LEAVE ONE OUT#########
##########################

# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Prepare empirical variogram using training data
  v_empirical <- variogram(rainfall ~ 1, train_data)
  
  # Define the Gaussian variogram model
  v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
  
  # Fit the Gaussian variogram model to the empirical variogram
  fit_gau <- fit.variogram(v_empirical, model = v_gau_model)
  
  # Perform Ordinary Kriging interpolation
  kriging_result <- krige(rainfall ~ 1, train_data, newdata = test_data, model = fit_gau)
  
  # Store the prediction
  kriging_out[i] <- kriging_result$var1.pred
}



# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape



metrics.data.frame.kriging<-cbind(plot_data,rmse)
metrics.data.frame.kriging<-cbind(metrics.data.frame.kriging,mape)




# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "OK - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 5 anos. MAPE Médio: 7.123%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "b.png"), plot = plot, width = 10, height = 8, dpi = 500)



###########################################################################
#######################################TR=25################################
############################################################################
####################################FIGURA 28###############################

counter=28
####################
#########FIGURA 10a#
####################

###IDW



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


library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)




############################################################
## Creating Grid Size with the same size of XAVIER resolution
############################################################


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height),crs=crs(polygon_shapefile))

#Convert raster to points
#raster_sites_points <- rasterToPoints(rasterized_precip_4326, spatial = TRUE)
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# 2. Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(polygon_shapefile))

plot_data<-0
plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(25))


#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

###############################################################################
###############################################################################
#####################################IDW#######################################
###############################################################################
###############################################################################

# Interpolate the grid cells using a power value of 5 (idp=2.5)
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
P.idw <- st_transform(P.idw, crs = st_crs(boundary_sf_transformed))

# Crop the IDW results using the boundary
P.idw_cropped <- st_intersection(P.idw, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
P.idw_cropped <- P.idw_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot<-ggplot() +
  # Tile layer for IDW interpolation
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for IDW interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade IDW") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação IDW - Chuva com duração de 24 h. Tempo de Retorno = 25 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")



#RMSE e MAPE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Perform IDW interpolation
  idw_result <- idw(rainfall ~ 1, train_data, test_data, idp = 2)
  
  # Store the prediction
  IDW.out[i] <- idw_result$var1.pred
}

# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame.idw<-0
metrics.data.frame.idw<-cbind(plot_data,rmse)
metrics.data.frame.idw<-cbind(metrics.data.frame.idw,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.idw

# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "IDW - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 25 anos. MAPE Média: 9.779%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "a.png"), plot = plot, width = 10, height = 8, dpi = 500)

################################################################################
################################################################################
###################################KRIGING######################################
################################################################################
################################################################################

# Check and transform CRS if necessary
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Verify CRS consistency
if (st_crs(plot_data_sf) != st_crs(sf_grid_xavier)) {
  stop("CRS mismatch between plot_data_sf and sf_grid_xavier")
}

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define variogram models with specific parameters
v_linear_model <- vgm(psill = 0.8, model = "Lin", range = 500000, nugget = 0.12)
v_sph_model <- vgm(psill = 0.8, model = "Sph", range = 500000, nugget = 0.12)
v_exp_model <- vgm(psill = 0.8, model = "Exp", range = 500000, nugget = 0.12)
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit different variogram models
fit_linear <- fit.variogram(v_empirical, model = v_linear_model)
fit_sph <- fit.variogram(v_empirical, model = v_sph_model)
fit_exp <- fit.variogram(v_empirical, model = v_exp_model)
fit_gauss <- fit.variogram(v_empirical, model = v_gau_model)

# Plot the empirical variogram and the fitted models
plot(v_empirical, model = fit_linear, main = "Fitted Linear Variogram")
plot(v_empirical, model = fit_sph, main = "Fitted Spherical Variogram")
plot(v_empirical, model = fit_exp, main = "Fitted Exponential Variogram")
plot(v_empirical, model = fit_gauss, main = "Fitted Gaussian Variogram")


#Performing Kriging with the the gaussian variogram


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height), crs = crs(polygon_shapefile))

# Convert raster to points
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(raster_grid_xavier))

plot_data <- stations_results_df %>% filter(Tr_values == as.character(25))

# Transform stations points to the same CRS as sf_grid_xavier
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, plot_data_sf)

# Define the Gaussian variogram model
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit the Gaussian variogram model to the empirical variogram
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)

# Perform Ordinary Kriging interpolation
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)

# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))

# Crop the Kriging results using the boundary
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for Kriging interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade Kriging") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação OK - Chuva com duração de 24 h. Tempo de Retorno = 25 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)





#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

###########################
####LEAVE ONE OUT#########
##########################

# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Prepare empirical variogram using training data
  v_empirical <- variogram(rainfall ~ 1, train_data)
  
  # Define the Gaussian variogram model
  v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
  
  # Fit the Gaussian variogram model to the empirical variogram
  fit_gau <- fit.variogram(v_empirical, model = v_gau_model)
  
  # Perform Ordinary Kriging interpolation
  kriging_result <- krige(rainfall ~ 1, train_data, newdata = test_data, model = fit_gau)
  
  # Store the prediction
  kriging_out[i] <- kriging_result$var1.pred
}



# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame.kriging<-0
metrics.data.frame.kriging<-cbind(plot_data,rmse)
metrics.data.frame.kriging<-cbind(metrics.data.frame.kriging,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.kriging


# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "OK - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 25 anos. MAPE Média: 9.821%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "b.png"), plot = plot, width = 10, height = 8, dpi = 500)



###########################################################################
#######################################TR=100################################
############################################################################
####################################FIGURA 29###############################

counter=29
####################
#########FIGURA 10a#
####################

###IDW



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


library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)




############################################################
## Creating Grid Size with the same size of XAVIER resolution
############################################################


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height),crs=crs(polygon_shapefile))

#Convert raster to points
#raster_sites_points <- rasterToPoints(rasterized_precip_4326, spatial = TRUE)
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# 2. Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(polygon_shapefile))

plot_data<-0
plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(100))


#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

###############################################################################
###############################################################################
#####################################IDW#######################################
###############################################################################
###############################################################################

# Interpolate the grid cells using a power value of 5 (idp=1)
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 1)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
P.idw <- st_transform(P.idw, crs = st_crs(boundary_sf_transformed))

# Crop the IDW results using the boundary
P.idw_cropped <- st_intersection(P.idw, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
P.idw_cropped <- P.idw_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot<-ggplot() +
  # Tile layer for IDW interpolation
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for IDW interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade IDW") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação IDW - Chuva com duração de 24 h. Tempo de Retorno = 100 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")



#RMSE e MAPE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Perform IDW interpolation
  idw_result <- idw(rainfall ~ 1, train_data, test_data, idp = 1)
  
  # Store the prediction
  IDW.out[i] <- idw_result$var1.pred
}

# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame.idw<-0
metrics.data.frame.idw<-cbind(plot_data,rmse)
metrics.data.frame.idw<-cbind(metrics.data.frame.idw,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.idw

# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "IDW - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 100 anos. MAPE Média: 14.754%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "a.png"), plot = plot, width = 10, height = 8, dpi = 500)

################################################################################
################################################################################
###################################KRIGING######################################
################################################################################
################################################################################

# Check and transform CRS if necessary
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Verify CRS consistency
if (st_crs(plot_data_sf) != st_crs(sf_grid_xavier)) {
  stop("CRS mismatch between plot_data_sf and sf_grid_xavier")
}

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define variogram models with specific parameters
v_linear_model <- vgm(psill = 0.8, model = "Lin", range = 500000, nugget = 0.12)
v_sph_model <- vgm(psill = 0.8, model = "Sph", range = 500000, nugget = 0.12)
v_exp_model <- vgm(psill = 0.8, model = "Exp", range = 500000, nugget = 0.12)
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit different variogram models
fit_linear <- fit.variogram(v_empirical, model = v_linear_model)
fit_sph <- fit.variogram(v_empirical, model = v_sph_model)
fit_exp <- fit.variogram(v_empirical, model = v_exp_model)
fit_gauss <- fit.variogram(v_empirical, model = v_gau_model)

# Plot the empirical variogram and the fitted models
plot(v_empirical, model = fit_linear, main = "Fitted Linear Variogram")
plot(v_empirical, model = fit_sph, main = "Fitted Spherical Variogram")
plot(v_empirical, model = fit_exp, main = "Fitted Exponential Variogram")
plot(v_empirical, model = fit_gauss, main = "Fitted Gaussian Variogram")


#Performing Kriging with the the gaussian variogram


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height), crs = crs(polygon_shapefile))

# Convert raster to points
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(raster_grid_xavier))

plot_data <-0
plot_data <- stations_results_df %>% filter(Tr_values == as.character(100))

# Transform stations points to the same CRS as sf_grid_xavier
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, plot_data_sf)

# Define the Gaussian variogram model
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit the Gaussian variogram model to the empirical variogram
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)

# Perform Ordinary Kriging interpolation
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)

# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))

# Crop the Kriging results using the boundary
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for Kriging interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade Kriging") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação OK - Chuva com duração de 24 h. Tempo de Retorno = 100 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)





#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

###########################
####LEAVE ONE OUT#########
##########################

# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Prepare empirical variogram using training data
  v_empirical <- variogram(rainfall ~ 1, train_data)
  
  # Define the Gaussian variogram model
  v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
  
  # Fit the Gaussian variogram model to the empirical variogram
  fit_gau <- fit.variogram(v_empirical, model = v_gau_model)
  
  # Perform Ordinary Kriging interpolation
  kriging_result <- krige(rainfall ~ 1, train_data, newdata = test_data, model = fit_gau)
  
  # Store the prediction
  kriging_out[i] <- kriging_result$var1.pred
}



# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame.kriging<-0
metrics.data.frame.kriging<-cbind(plot_data,rmse)
metrics.data.frame.kriging<-cbind(metrics.data.frame.kriging,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.kriging

# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "OK - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 100 anos. MAPE Média: 15.467%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "b.png"), plot = plot, width = 10, height = 8, dpi = 500)




###########################################################################
#######################################TR=500################################
############################################################################
####################################FIGURA 30###############################

counter=30
####################
#########FIGURA 10a#
####################

###IDW



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


library(sf)
library(stars)
library(gstat)
library(automap)
library(leaflet)
library(leafem)
library(terra)
library(raster)




############################################################
## Creating Grid Size with the same size of XAVIER resolution
############################################################


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height),crs=crs(polygon_shapefile))

#Convert raster to points
#raster_sites_points <- rasterToPoints(rasterized_precip_4326, spatial = TRUE)
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# 2. Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(polygon_shapefile))

plot_data<-0
plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(500))


#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

###############################################################################
###############################################################################
#####################################IDW#######################################
###############################################################################
###############################################################################

# Interpolate the grid cells using a power value of 5 (idp=1)
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 0.5)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
P.idw <- st_transform(P.idw, crs = st_crs(boundary_sf_transformed))

# Crop the IDW results using the boundary
P.idw_cropped <- st_intersection(P.idw, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
P.idw_cropped <- P.idw_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot<-ggplot() +
  # Tile layer for IDW interpolation
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for IDW interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade IDW") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação IDW - Chuva com duração de 24 h. Tempo de Retorno = 500 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")



#RMSE e MAPE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Perform IDW interpolation
  idw_result <- idw(rainfall ~ 1, train_data, test_data, idp = 0.5)
  
  # Store the prediction
  IDW.out[i] <- idw_result$var1.pred
}

# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame.idw<-0
metrics.data.frame.idw<-cbind(plot_data,rmse)
metrics.data.frame.idw<-cbind(metrics.data.frame.idw,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.idw

# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "IDW - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 500 anos. MAPE Média: 21.805%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "a.png"), plot = plot, width = 10, height = 8, dpi = 500)

################################################################################
################################################################################
###################################KRIGING######################################
################################################################################
################################################################################

# Check and transform CRS if necessary
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Verify CRS consistency
if (st_crs(plot_data_sf) != st_crs(sf_grid_xavier)) {
  stop("CRS mismatch between plot_data_sf and sf_grid_xavier")
}

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define variogram models with specific parameters
v_linear_model <- vgm(psill = 0.8, model = "Lin", range = 500000, nugget = 0.12)
v_sph_model <- vgm(psill = 0.8, model = "Sph", range = 500000, nugget = 0.12)
v_exp_model <- vgm(psill = 0.8, model = "Exp", range = 500000, nugget = 0.12)
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit different variogram models
fit_linear <- fit.variogram(v_empirical, model = v_linear_model)
fit_sph <- fit.variogram(v_empirical, model = v_sph_model)
fit_exp <- fit.variogram(v_empirical, model = v_exp_model)
fit_gauss <- fit.variogram(v_empirical, model = v_gau_model)

# Plot the empirical variogram and the fitted models
plot(v_empirical, model = fit_linear, main = "Fitted Linear Variogram")
plot(v_empirical, model = fit_sph, main = "Fitted Spherical Variogram")
plot(v_empirical, model = fit_exp, main = "Fitted Exponential Variogram")
plot(v_empirical, model = fit_gauss, main = "Fitted Gaussian Variogram")


#Performing Kriging with the the gaussian variogram


# Calculate the resolution based on the size of the squared polygons
# Assuming the polygons are regular and have the same size
poly_bbox <- st_bbox(polygon_shapefile)
poly_width <- diff(poly_bbox[c("xmin", "xmax")]) / sqrt(nrow(polygon_shapefile))
poly_height <- diff(poly_bbox[c("ymin", "ymax")]) / sqrt(nrow(polygon_shapefile))

# Create an empty raster with the calculated resolution
raster_grid_xavier <- raster(extent(poly_bbox), res = c(poly_width, poly_height), crs = crs(polygon_shapefile))

# Convert raster to points
df_grid_xavier <- as.data.frame(raster_grid_xavier, xy = TRUE)

# Convert to sf object
sf_grid_xavier <- st_as_sf(df_grid_xavier, coords = c("x", "y"), crs = crs(raster_grid_xavier))

plot_data <-0
plot_data <- stations_results_df %>% filter(Tr_values == as.character(500))

# Transform stations points to the same CRS as sf_grid_xavier
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, plot_data_sf)

# Define the Gaussian variogram model
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)

# Fit the Gaussian variogram model to the empirical variogram
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)

# Perform Ordinary Kriging interpolation
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)

# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Ensure CRS consistency
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))

# Crop the Kriging results using the boundary
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao, color = after_stat(level)), size = 0.5) +
  # Fill gradient for Kriging interpolation
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Intensidade Kriging") +
  # Color scale for contour lines
  scale_color_viridis_c(name = "Contours") +
  # Custom color scale for points reflecting rainfall values with multiple transitions
  scale_color_gradientn(
    colors = c("lightblue", "lightgreen", "yellow", "orange", "red"),
    name = "Intensidade - Estações(mm/h)",
    guide = guide_colorbar(title = "Intensidade - Estações(mm/h)")
  ) +
  # Add boundary
  geom_sf(data = boundary_sf_transformed, fill = NA, color = "black", size = 1) +
  # Coordinate system
  coord_sf() +
  # Labels and theme
  labs(title = "Interpolação OK - Chuva com duração de 24 h. Tempo de Retorno = 100 anos",
       fill = "Intensidade - Interpolação(mm/h)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(plot)





#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

###########################
####LEAVE ONE OUT#########
##########################

# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")

# Initialize vector for predictions
kriging_out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  # Prepare empirical variogram using training data
  v_empirical <- variogram(rainfall ~ 1, train_data)
  
  # Define the Gaussian variogram model
  v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
  
  # Fit the Gaussian variogram model to the empirical variogram
  fit_gau <- fit.variogram(v_empirical, model = v_gau_model)
  
  # Perform Ordinary Kriging interpolation
  kriging_result <- krige(rainfall ~ 1, train_data, newdata = test_data, model = fit_gau)
  
  # Store the prediction
  kriging_out[i] <- kriging_result$var1.pred
}



# Extract observed values
observed_values <- plot_data_sp$rainfall

# Create a data frame for plotting
validation_df <- data.frame(
  Observed = observed_values,
  Predicted = IDW.out
)

rmse<-sqrt((validation_df$Predicted-validation_df$Observed)^2)
rmse

mape <- (abs((validation_df$Predicted - validation_df$Observed) / validation_df$Observed)) * 100
mape


metrics.data.frame.kriging<-0
metrics.data.frame.kriging<-cbind(plot_data,rmse)
metrics.data.frame.kriging<-cbind(metrics.data.frame.kriging,mape)



# Filter data for the current Tr value
plot_data <- metrics.data.frame.kriging

# Create the bubble plot
plot <- ggplot(data = plot_data) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed", linewidth = 1) +  # Plot the boundary
  geom_point(aes(x = lon, y = lat, size = mape, color = mape), alpha = 0.7) +
  scale_size_continuous(name = "MAPE %", guide = "none") +  # Reduced bubble size and removed size legend
  scale_color_gradient(low = "lightblue", high = "red", name = "MAPE %") +  # Inverted color pattern
  labs(title = "OK - LOOCV - MAPE para cada estação",
       subtitle = paste("Chuva de Duração de 24 h. Tempo de Retorno = 500 anos. MAPE Média: 22.276%"),
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


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_26")
# Save the plot with a sequential filename
ggsave(filename = paste0("Figura", counter, "b.png"), plot = plot, width = 10, height = 8, dpi = 500)


