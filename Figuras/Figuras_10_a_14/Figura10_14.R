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

matrix.results.metrics.IDW<-matrix(NA,nrow=4,ncol=5)
matrix.results.metrics.KRIGING<-matrix(NA,nrow=4,ncol=5)

rownames(matrix.results.metrics.IDW)=c("Tr","RMSE","MAPE","MODELO")
rownames(matrix.results.metrics.KRIGING)=c("Tr","RMSE","MAPE","MODELO")

matrix.results.metrics.IDW=as.data.frame(matrix.results.metrics.IDW)
matrix.results.metrics.KRIGING=as.data.frame(matrix.results.metrics.KRIGING)

matrix.results.metrics.IDW[1,]=c(2,5,25,100,500)
matrix.results.metrics.IDW[4,]=rep("IDW",5)

matrix.results.metrics.KRIGING[1,]=c(2,5,25,100,500)
matrix.results.metrics.KRIGING[4,]=rep("KRIGING",5)
############################################################################
###########################FIRST ATTEMPT TO GEV_MML#########################
############################################################################


###########################################################################
#######################################TR=2#################################
############################################################################
####################################FIGURA 10###############################


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


####################
#########FIGURA 10c#
####################



#RMSE-LEAVE ONE OUT


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
    
    
    # Compute RMSE
    rmse <- sqrt(sum((IDW.out - observed_values)^2) / length(observed_values))
    matrix.results.metrics.IDW[2,1]=rmse
    # Compute MAPE
    mape <- mean(abs((IDW.out - observed_values) / observed_values)) * 100
    matrix.results.metrics.IDW[3,1]=mape
    
   
    # Plot observed vs. predicted values
    grafico<-ggplot(validation_df, aes(x = Observed, y = Predicted)) +
      geom_point(color = "black", alpha = 0.5) +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      labs(x = "Observado", y = "Previsto", title =paste("IDW. Leave-One-Out - Cross Validation: Observado X Previsto. Tr=2 anos. RMSE=",round(rmse,3),"; MAPE=",round(mape,3),"%" )) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90")
      )
    
print(grafico)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura10c.png"), plot = grafico, width = 10, height = 8, dpi = 500)



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

#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura10b.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 10d#
####################



#RMSE-LEAVE ONE OUT


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
  Predicted = kriging_out
)


# Compute RMSE
rmse <- sqrt(sum((kriging_out - observed_values)^2) / length(observed_values))
matrix.results.metrics.KRIGING[2,1]=rmse

# Compute MAPE
mape <- mean(abs((kriging_out - observed_values) / observed_values)) * 100
matrix.results.metrics.KRIGING[3,1]=mape


# Plot observed vs. predicted values
grafico <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title = paste("OK. Leave-One-Out - Cross Validation: Observado X Previsto. Tr= 2 anos. RMSE=", round(rmse, 3), "; MAPE=", round(mape, 3), "%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)


#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura10d.png"), plot = grafico, width = 10, height = 8, dpi = 500)














###########################################################################
#######################################TR=5#################################
############################################################################
####################################FIGURA 11###############################



####################
#########FIGURA 11a#
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
ggsave(filename = paste0("Figura11a.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 11c#
####################



#RMSE-LEAVE ONE OUT


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


# Compute RMSE
rmse <- sqrt(sum((IDW.out - observed_values)^2) / length(observed_values))
matrix.results.metrics.IDW[2,2]=rmse
# Compute MAPE
mape <- mean(abs((IDW.out - observed_values) / observed_values)) * 100
matrix.results.metrics.IDW[3,2]=mape


# Plot observed vs. predicted values
grafico<-ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title =paste("IDW. Leave-One-Out - Cross Validation: Observado X Previsto. Tr= 5 anos. RMSE=",round(rmse,3),"; MAPE=",round(mape,3),"%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura11c.png"), plot = grafico, width = 10, height = 8, dpi = 500)



####################
#########FIGURA 11b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)

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

#TROCAR TR AQUI!!!!!!!!!!!
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

#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura11b.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 11d#
####################



#RMSE-LEAVE ONE OUT


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
  Predicted = kriging_out
)


# Compute RMSE
rmse <- sqrt(sum((kriging_out - observed_values)^2) / length(observed_values))
matrix.results.metrics.KRIGING[2,2]=rmse

# Compute MAPE
mape <- mean(abs((kriging_out - observed_values) / observed_values)) * 100
matrix.results.metrics.KRIGING[3,2]=mape


# Plot observed vs. predicted values
grafico <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title = paste("OK. Leave-One-Out - Cross Validation: Observado X Previsto.Tr= 5 anos. RMSE=", round(rmse, 3), "; MAPE=", round(mape, 3), "%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)


#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura11d.png"), plot = grafico, width = 10, height = 8, dpi = 500)


























###########################################################################
#######################################TR=25#################################
############################################################################
####################################FIGURA 12###############################



####################
#########FIGURA 12a#
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



#TROCAR O TR AQUI!!!!!!!!!!!!!!!!
plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(25))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))


#TROCAR O VALOR DE p AQUI!!!!!!!!!!!
# Interpolate the grid cells using a power value of 2 (idp=2)
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
ggsave(filename = paste0("Figura12a.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 12c#
####################



#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  #TROCAR VALOR DE P AQUI!!!!!!!!!!!!!!
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


# Compute RMSE
rmse <- sqrt(sum((IDW.out - observed_values)^2) / length(observed_values))
matrix.results.metrics.IDW[2,3]=rmse
# Compute MAPE
mape <- mean(abs((IDW.out - observed_values) / observed_values)) * 100
matrix.results.metrics.IDW[3,3]=mape


# Plot observed vs. predicted values
grafico<-ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title =paste("IDW. Leave-One-Out - Cross Validation: Observado X Previsto. Tr= 25 anos. RMSE=",round(rmse,3),"; MAPE=",round(mape,3),"%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura12c.png"), plot = grafico, width = 10, height = 8, dpi = 500)



####################
#########FIGURA 12b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)

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

#TROCAR TR AQUI!!!!!!!!!!
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

#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura12b.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 12d#
####################



#RMSE-LEAVE ONE OUT


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
  Predicted = kriging_out
)


# Compute RMSE
rmse <- sqrt(sum((kriging_out - observed_values)^2) / length(observed_values))
matrix.results.metrics.KRIGING[2,3]=rmse

# Compute MAPE
mape <- mean(abs((kriging_out - observed_values) / observed_values)) * 100
matrix.results.metrics.KRIGING[3,3]=mape


# Plot observed vs. predicted values
grafico <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title = paste("OK. Leave-One-Out - Cross Validation: Observado X Previsto.Tr= 25 anos. RMSE=", round(rmse, 3), "; MAPE=", round(mape, 3), "%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)


#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura12d.png"), plot = grafico, width = 10, height = 8, dpi = 500)




















###########################################################################
#######################################TR=100#################################
############################################################################
####################################FIGURA 13###############################



####################
#########FIGURA 13a#
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



#TROCAR O TR AQUI!!!!!!!!!!!!!!!!
plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(100))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))


#TROCAR O VALOR DE p AQUI!!!!!!!!!!!
# Interpolate the grid cells using a power value of 2 (idp=2)
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
ggsave(filename = paste0("Figura13a.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 13c#
####################



#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  #TROCAR VALOR DE P AQUI!!!!!!!!!!!!!!
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


# Compute RMSE
rmse <- sqrt(sum((IDW.out - observed_values)^2) / length(observed_values))
matrix.results.metrics.IDW[2,4]=rmse
# Compute MAPE
mape <- mean(abs((IDW.out - observed_values) / observed_values)) * 100
matrix.results.metrics.IDW[3,4]=mape


# Plot observed vs. predicted values
grafico<-ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title =paste("IDW. Leave-One-Out - Cross Validation: Observado X Previsto. Tr= 100 anos. RMSE=",round(rmse,3),"; MAPE=",round(mape,3),"%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura13c.png"), plot = grafico, width = 10, height = 8, dpi = 500)



####################
#########FIGURA 13b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)

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

#TROCAR TR AQUI!!!!!!!!!!!!!!!!
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

#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura13b.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 13d#
####################



#RMSE-LEAVE ONE OUT


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
  Predicted = kriging_out
)


# Compute RMSE
rmse <- sqrt(sum((kriging_out - observed_values)^2) / length(observed_values))
matrix.results.metrics.KRIGING[2,4]=rmse

# Compute MAPE
mape <- mean(abs((kriging_out - observed_values) / observed_values)) * 100
matrix.results.metrics.KRIGING[3,4]=mape


# Plot observed vs. predicted values
grafico <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title = paste("OK. Leave-One-Out - Cross Validation: Observado X Previsto.Tr= 100 anos. RMSE=", round(rmse, 3), "; MAPE=", round(mape, 3), "%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)


#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura13d.png"), plot = grafico, width = 10, height = 8, dpi = 500)


















###########################################################################
#######################################TR=500#################################
############################################################################
####################################FIGURA 14###############################



####################
#########FIGURA 14a#
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



#TROCAR O TR AQUI!!!!!!!!!!!!!!!!
plot_data<-plot_data <- stations_results_df %>% filter(Tr_values == as.character(500))

#Transform hidro stations points to the same CRS as sf_grid_chirps
plot_data_sf <- st_transform(plot_data, st_crs(sf_grid_xavier))

#Transform sf_hidro_stations_points to the same CRS as sf_grid_chirps
boundary_sf_transformed <- st_transform(boundary_sf, st_crs(sf_grid_xavier))


#TROCAR O VALOR DE p AQUI!!!!!!!!!!!
# Interpolate the grid cells using a power value of 2 (idp=2)
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
ggsave(filename = paste0("Figura14a.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 14c#
####################



#RMSE-LEAVE ONE OUT


# Convert to Spatial object for gstat functions
plot_data_sp <- as(plot_data, "Spatial")
# Initialize vector for predictions
IDW.out <- numeric(length = nrow(plot_data_sp))

for (i in 1:nrow(plot_data_sp)) {
  # Exclude the ith point for training
  train_data <- plot_data_sp[-i, ]
  test_data <- plot_data_sp[i, ]
  
  #TROCAR VALOR DE P AQUI!!!!!!!!!!!!!!
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


# Compute RMSE
rmse <- sqrt(sum((IDW.out - observed_values)^2) / length(observed_values))
matrix.results.metrics.IDW[2,5]=rmse
# Compute MAPE
mape <- mean(abs((IDW.out - observed_values) / observed_values)) * 100
matrix.results.metrics.IDW[3,5]=mape


# Plot observed vs. predicted values
grafico<-ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title =paste("IDW. Leave-One-Out - Cross Validation: Observado X Previsto. Tr= 500 anos. RMSE=",round(rmse,3),"; MAPE=",round(mape,3),"%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)
#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura14c.png"), plot = grafico, width = 10, height = 8, dpi = 500)



####################
#########FIGURA 14b#
####################


###############
###KRIGING####
##############

# Load required libraries
library(gstat)
library(sf)
library(ggplot2)
library(dplyr)

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
  labs(title = "Interpolação OK - Chuva com duração de 24 h. Tempo de Retorno = 500 anos",
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
ggsave(filename = paste0("Figura14b.png"), plot = plot, width = 10, height = 8, dpi = 500)


####################
#########FIGURA 14d#
####################



#RMSE-LEAVE ONE OUT


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
  Predicted = kriging_out
)


# Compute RMSE
rmse <- sqrt(sum((kriging_out - observed_values)^2) / length(observed_values))
matrix.results.metrics.KRIGING[2,5]=rmse

# Compute MAPE
mape <- mean(abs((kriging_out - observed_values) / observed_values)) * 100
matrix.results.metrics.KRIGING[3,5]=mape


# Plot observed vs. predicted values
grafico <- ggplot(validation_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Observado", y = "Previsto", title = paste("OK. Leave-One-Out - Cross Validation: Observado X Previsto.Tr= 500 anos. RMSE=", round(rmse, 3), "; MAPE=", round(mape, 3), "%" )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

print(grafico)


#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura14d.png"), plot = grafico, width = 10, height = 8, dpi = 500)




final.comparison.data.frame<-t(matrix.results.metrics.IDW)

final.comparison.data.frame<-rbind(final.comparison.data.frame,t(matrix.results.metrics.KRIGING))

final.comparison.data.frame<-as.data.frame(final.comparison.data.frame)

final.comparison.data.frame$RMSE<-as.numeric(final.comparison.data.frame$RMSE)
final.comparison.data.frame$MAPE<-as.numeric(final.comparison.data.frame$MAPE)
final.comparison.data.frame$Tr<-as.numeric(final.comparison.data.frame$Tr)


# Create the plot
plot <- ggplot(final.comparison.data.frame, aes(x = Tr, y = RMSE, color = MODELO, shape = MODELO, linetype = MODELO)) +
  geom_line(size = 1) +  # Draw lines
  geom_point(size = 3) +  # Add points
  scale_x_log10() +  # Use log scale for x-axis
  scale_y_continuous() +  # Linear scale for y-axis
  labs(title = "Tempo de Retorno X RMSE",
       x = "Tempo de Retorno (anos)",
       y = "RMSE",
       color = "Modelo",
       shape = "Modelo",
       linetype = "Modelo") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

# Print the plot
print(plot)

#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura15a.png"), plot = plot, width = 10, height = 8, dpi = 500)




# Create the plot
plot <- ggplot(final.comparison.data.frame, aes(x = Tr, y = MAPE, color = MODELO, shape = MODELO, linetype = MODELO)) +
  geom_line(size = 1) +  # Draw lines
  geom_point(size = 3) +  # Add points
  scale_x_log10() +  # Use log scale for x-axis
  scale_y_continuous() +  # Linear scale for y-axis
  labs(title = "Tempo de Retorno X MAPE",
       x = "Tempo de Retorno (anos)",
       y = "MAPE (%)",
       color = "Modelo",
       shape = "Modelo",
       linetype = "Modelo") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

# Print the plot
print(plot)


#Set working directory
setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\gitHUB\\idf_RS_Technical_Note2\\Figuras\\Figuras_10_a_14")
ggsave(filename = paste0("Figura15b.png"), plot = plot, width = 10, height = 8, dpi = 500)


