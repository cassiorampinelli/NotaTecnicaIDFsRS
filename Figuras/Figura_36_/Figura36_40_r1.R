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
####################################FIGURA 36###############################


####################
#########FIGURA 36a#
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

#Reading polygon shapefile
sedes.municipais_sf <- st_read("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\IDF_Project\\Rio_Grande_do_Sul\\GIS\\Sedes_Municipais_RS.shp")


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
sedes.municipais_sf <-assign_crs_if_missing(sedes.municipais_sf , crs_wgs84)


# Ensure the coordinate reference systems match
if (st_crs(polygon_shapefile) != st_crs(points_shapefile)) {
  points_shapefile <- st_transform(points_shapefile, st_crs(polygon_shapefile))
}



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
    dplyr::select(codigo,lon, lat, rainfall, Tr_values,idf_a, idf_b, idf_c,idf_d)
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

plot_data <- stations_results_df %>% filter(Tr_values == as.character(2))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(polygon_shapefile ))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(polygon_shapefile ))

# Interpolate the grid cells using a power value of 2 (idp=2.5) for rainfall
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sedes.municipais_sf, idp = 2.5)

# Interpolate idf_a
idf_a_idw <- idw(idf_a ~ 1, plot_data_sf, newdata=sedes.municipais_sf, idp = 2.5)

# Interpolate idf_b
idf_b_idw <- idw(idf_b ~ 1, plot_data_sf, newdata=sedes.municipais_sf, idp = 2.5)

# Interpolate idf_c
idf_c_idw <- idw(idf_c ~ 1, plot_data_sf, newdata=sedes.municipais_sf, idp = 2.5)

# Interpolate idf_d
idf_d_idw <- idw(idf_d ~ 1, plot_data_sf, newdata=sedes.municipais_sf, idp = 2.5)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)



Tr=2

P.idw2<-calculate_rainfall(idf_a_idw$var1.pred, idf_b_idw$var1.pred, idf_c_idw$var1.pred,idf_d_idw$var1.pred, Tr) 


P.idw$precipitacao2<-P.idw2


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
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "IDW com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 2 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura31a.png"), plot = plot, width = 10, height = 8, dpi = 500)




library(gstat)

# Function to perform IDW interpolation with a specific idp value and maximum number of neighbors
interpolate_idw <- function(variable, idp_value, plot_data_sf, newdata, nmax) {
  # Create the gstat object for the variable of interest
  gstat_obj <- gstat(formula = as.formula(paste(variable, "~ 1")),
                     locations = plot_data_sf,
                     nmax = nmax,   # Number of nearest neighbors to consider
                     set = list(idp = idp_value))
  
  # Perform IDW interpolation
  interpolated <- predict(gstat_obj, newdata = newdata)
  
  return(interpolated)
}

# Interpolate the grid cells for each variable with a power value of 2.5 and maximum 15 neighbors
nmax_neighbors <- 15
idp_value <- 2.5

# Interpolate rainfall
P.idw <- interpolate_idw("rainfall", idp_value, plot_data_sf, sedes.municipais_sf, nmax_neighbors)

# Interpolate idf_a
idf_a_idw <- interpolate_idw("idf_a", idp_value, plot_data_sf, sedes.municipais_sf, nmax_neighbors)

# Interpolate idf_b
idf_b_idw <- interpolate_idw("idf_b", idp_value, plot_data_sf, sedes.municipais_sf, nmax_neighbors)

# Interpolate idf_c
idf_c_idw <- interpolate_idw("idf_c", idp_value, plot_data_sf, sedes.municipais_sf, nmax_neighbors)

# Interpolate idf_d
idf_d_idw <- interpolate_idw("idf_d", idp_value, plot_data_sf, sedes.municipais_sf, nmax_neighbors)


sedes.municipais_sf$idf_a_idw<-idf_a_idw $var1.pred
sedes.municipais_sf$idf_b_idw<-idf_b_idw $var1.pred
sedes.municipais_sf$idf_c_idw<-idf_c_idw $var1.pred
sedes.municipais_sf$idf_d_idw<-idf_d_idw $var1.pred

head(sedes.municipais_sf)

output_file<-"C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_36_\\sedes.municipais_sf.shp"
# Save the sf object as a shapefile
st_write(sedes.municipais_sf, output_file)



###############################
#############################
##########################












####################
#########FIGURA 10b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)


# Ensure CRS consistency between sf_grid_xavier and plot_data_sf
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define and fit the Gaussian variogram model to the empirical variogram
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)


# Perform Kriging interpolation for rainfall
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Perform Kriging interpolation for idf_a, idf_b, idf_c, idf_d
idf_a_krig <- krige(idf_a ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_b_krig <- krige(idf_b ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_c_krig <- krige(idf_c ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_d_krig <- krige(idf_d ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Calculate rainfall using interpolated IDF parameters
Tr <- 2
kriging_result$precipitacao2 <- calculate_rainfall(idf_a_krig$var1.pred, idf_b_krig$var1.pred, idf_c_krig$var1.pred, idf_d_krig$var1.pred, Tr)

# Ensure CRS consistency and crop the Kriging results using the boundary
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "OK com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 2 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura31b.png"), plot = plot, width = 10, height = 8, dpi = 500)




###########################################################################
#######################################TR=5#################################
############################################################################
####################################FIGURA 32###############################


####################
#########FIGURA 32a#
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
    dplyr::select(codigo,lon, lat, rainfall, Tr_values,idf_a, idf_b, idf_c,idf_d)
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

plot_data <- stations_results_df %>% filter(Tr_values == as.character(5))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))



#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Interpolate the grid cells using a power value of 2 (idp=2.5) for rainfall
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)

# Interpolate idf_a
idf_a_idw <- idw(idf_a ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)

# Interpolate idf_b
idf_b_idw <- idw(idf_b ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)

# Interpolate idf_c
idf_c_idw <- idw(idf_c ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)

# Interpolate idf_d
idf_d_idw <- idw(idf_d ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2.5)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)



Tr=5

P.idw2<-calculate_rainfall(idf_a_idw$var1.pred, idf_b_idw$var1.pred, idf_c_idw$var1.pred,idf_d_idw$var1.pred, Tr) 


P.idw$precipitacao2<-P.idw2


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
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "IDW com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 5 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura32a.png"), plot = plot, width = 10, height = 8, dpi = 500)




####################
#########FIGURA 32b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)


# Ensure CRS consistency between sf_grid_xavier and plot_data_sf
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define and fit the Gaussian variogram model to the empirical variogram
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)


# Perform Kriging interpolation for rainfall
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Perform Kriging interpolation for idf_a, idf_b, idf_c, idf_d
idf_a_krig <- krige(idf_a ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_b_krig <- krige(idf_b ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_c_krig <- krige(idf_c ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_d_krig <- krige(idf_d ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Calculate rainfall using interpolated IDF parameters
Tr <- 5
kriging_result$precipitacao2 <- calculate_rainfall(idf_a_krig$var1.pred, idf_b_krig$var1.pred, idf_c_krig$var1.pred, idf_d_krig$var1.pred, Tr)

# Ensure CRS consistency and crop the Kriging results using the boundary
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "OK com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 5 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura32b.png"), plot = plot, width = 10, height = 8, dpi = 500)




###########################################################################
#######################################TR=25#################################
############################################################################
####################################FIGURA 33###############################


####################
#########FIGURA 33a#
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
    dplyr::select(codigo,lon, lat, rainfall, Tr_values,idf_a, idf_b, idf_c,idf_d)
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

plot_data <- stations_results_df %>% filter(Tr_values == as.character(25))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))



#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Interpolate the grid cells using a power value of 2 (idp=2) for rainfall
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2)

# Interpolate idf_a
idf_a_idw <- idw(idf_a ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2)

# Interpolate idf_b
idf_b_idw <- idw(idf_b ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2)

# Interpolate idf_c
idf_c_idw <- idw(idf_c ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2)

# Interpolate idf_d
idf_d_idw <- idw(idf_d ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 2)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)



Tr=25

P.idw2<-calculate_rainfall(idf_a_idw$var1.pred, idf_b_idw$var1.pred, idf_c_idw$var1.pred,idf_d_idw$var1.pred, Tr) 


P.idw$precipitacao2<-P.idw2


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
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "IDW com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 25 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura33a.png"), plot = plot, width = 10, height = 8, dpi = 500)




####################
#########FIGURA 33b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)


# Ensure CRS consistency between sf_grid_xavier and plot_data_sf
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define and fit the Gaussian variogram model to the empirical variogram
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)


# Perform Kriging interpolation for rainfall
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Perform Kriging interpolation for idf_a, idf_b, idf_c, idf_d
idf_a_krig <- krige(idf_a ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_b_krig <- krige(idf_b ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_c_krig <- krige(idf_c ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_d_krig <- krige(idf_d ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Calculate rainfall using interpolated IDF parameters
Tr <- 25
kriging_result$precipitacao2 <- calculate_rainfall(idf_a_krig$var1.pred, idf_b_krig$var1.pred, idf_c_krig$var1.pred, idf_d_krig$var1.pred, Tr)

# Ensure CRS consistency and crop the Kriging results using the boundary
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "OK com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 25 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura33b.png"), plot = plot, width = 10, height = 8, dpi = 500)



###########################################################################
#######################################TR=100#################################
############################################################################
####################################FIGURA 34###############################


####################
#########FIGURA 33a#
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
    dplyr::select(codigo,lon, lat, rainfall, Tr_values,idf_a, idf_b, idf_c,idf_d)
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

plot_data <- stations_results_df %>% filter(Tr_values == as.character(100))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))



#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Interpolate the grid cells using a power value of 1 (idp=1) for rainfall
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 1)

# Interpolate idf_a
idf_a_idw <- idw(idf_a ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 1)

# Interpolate idf_b
idf_b_idw <- idw(idf_b ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 1)

# Interpolate idf_c
idf_c_idw <- idw(idf_c ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 1)

# Interpolate idf_d
idf_d_idw <- idw(idf_d ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 1)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)



Tr=100

P.idw2<-calculate_rainfall(idf_a_idw$var1.pred, idf_b_idw$var1.pred, idf_c_idw$var1.pred,idf_d_idw$var1.pred, Tr) 


P.idw$precipitacao2<-P.idw2


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
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "IDW com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 100 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura34a.png"), plot = plot, width = 10, height = 8, dpi = 500)




####################
#########FIGURA 34b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)


# Ensure CRS consistency between sf_grid_xavier and plot_data_sf
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define and fit the Gaussian variogram model to the empirical variogram
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)


# Perform Kriging interpolation for rainfall
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Perform Kriging interpolation for idf_a, idf_b, idf_c, idf_d
idf_a_krig <- krige(idf_a ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_b_krig <- krige(idf_b ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_c_krig <- krige(idf_c ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_d_krig <- krige(idf_d ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Calculate rainfall using interpolated IDF parameters
Tr <- 100
kriging_result$precipitacao2 <- calculate_rainfall(idf_a_krig$var1.pred, idf_b_krig$var1.pred, idf_c_krig$var1.pred, idf_d_krig$var1.pred, Tr)

# Ensure CRS consistency and crop the Kriging results using the boundary
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "OK com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 100 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura34b.png"), plot = plot, width = 10, height = 8, dpi = 500)




###########################################################################
#######################################TR=500#################################
############################################################################
####################################FIGURA 35###############################


####################
#########FIGURA 34a#
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
    dplyr::select(codigo,lon, lat, rainfall, Tr_values,idf_a, idf_b, idf_c,idf_d)
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

plot_data <- stations_results_df %>% filter(Tr_values == as.character(500))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))



#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Interpolate the grid cells using a power value of 1 (idp=1) for rainfall
P.idw <- idw(rainfall ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 0.5)

# Interpolate idf_a
idf_a_idw <- idw(idf_a ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 0.5)

# Interpolate idf_b
idf_b_idw <- idw(idf_b ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 0.5)

# Interpolate idf_c
idf_c_idw <- idw(idf_c ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 0.5)

# Interpolate idf_d
idf_d_idw <- idw(idf_d ~ 1, plot_data_sf, newdata=sf_grid_xavier, idp = 0.5)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(P.idw)) {
  # If not correctly identified, you may need to explicitly assign it
  st_geometry(P.idw) <- "geometry"
}

# Create a new column for precipitation values
P.idw <- P.idw %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)



Tr=500

P.idw2<-calculate_rainfall(idf_a_idw$var1.pred, idf_b_idw$var1.pred, idf_c_idw$var1.pred,idf_d_idw$var1.pred, Tr) 


P.idw$precipitacao2<-P.idw2


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
  geom_tile(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for IDW interpolation
  geom_contour(data = P.idw_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "IDW com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 100 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura35a.png"), plot = plot, width = 10, height = 8, dpi = 500)




####################
#########FIGURA 35b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)


# Ensure CRS consistency between sf_grid_xavier and plot_data_sf
if (st_crs(sf_grid_xavier) != st_crs(plot_data_sf)) {
  sf_grid_xavier <- st_transform(sf_grid_xavier, crs = st_crs(plot_data_sf))
}

# Transform boundary to the same CRS as sf_grid_xavier
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))

# Define the correct CRS based on your data location
crs_utm_21 <- "+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs"

# Prepare empirical variogram
v_empirical <- variogram(rainfall ~ 1, st_transform(plot_data_sf, crs = crs_utm_21))

# Define and fit the Gaussian variogram model to the empirical variogram
v_gau_model <- vgm(psill = 0.8, model = "Gau", range = 500000, nugget = 0.12)
fit_gau <- fit.variogram(v_empirical, model = v_gau_model)


# Perform Kriging interpolation for rainfall
kriging_result <- krige(rainfall ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Perform Kriging interpolation for idf_a, idf_b, idf_c, idf_d
idf_a_krig <- krige(idf_a ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_b_krig <- krige(idf_b ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_c_krig <- krige(idf_c ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)
idf_d_krig <- krige(idf_d ~ 1, plot_data_sf, newdata = sf_grid_xavier, model = fit_gau)


# Ensure the geometry column is correctly identified
if (!"geometry" %in% names(kriging_result)) {
  st_geometry(kriging_result) <- "geometry"
}

# Create a new column for precipitation values
kriging_result <- kriging_result %>%
  mutate(precipitacao = var1.pred) %>%
  dplyr::select(geometry, precipitacao)

# Calculate rainfall using interpolated IDF parameters
Tr <- 500
kriging_result$precipitacao2 <- calculate_rainfall(idf_a_krig$var1.pred, idf_b_krig$var1.pred, idf_c_krig$var1.pred, idf_d_krig$var1.pred, Tr)

# Ensure CRS consistency and crop the Kriging results using the boundary
kriging_result <- st_transform(kriging_result, crs = st_crs(boundary_sf_transformed))
kriging_result_cropped <- st_intersection(kriging_result, boundary_sf_transformed)

# Remove any duplicate levels in the precipitation data
kriging_result_cropped <- kriging_result_cropped %>%
  distinct(geometry, .keep_all = TRUE)

# Generate the plot
plot <- ggplot() +
  # Tile layer for Kriging interpolation
  geom_tile(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], fill = precipitacao2)) +
  # Points layer for stations with color reflecting the magnitude of rainfall
  geom_point(data = plot_data, aes(x = lon, y = lat, color = rainfall), size = 2) +
  # Contour lines for Kriging interpolation
  geom_contour(data = kriging_result_cropped, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], z = precipitacao2, color = after_stat(level)), size = 0.5) +
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
  labs(title = "OK com parâmetros da IDF - Chuva com duração de 24 h. Tempo de Retorno = 500 anos",
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
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figura_31")
ggsave(filename = paste0("Figura35b.png"), plot = plot, width = 10, height = 8, dpi = 500)




