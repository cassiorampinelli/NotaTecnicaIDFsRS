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

###########
#XAVIER####FIGURA 2a
###########

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

# Plotting
p <- ggplot() +
  # Plot all polygons
  geom_sf(data = polygon_shapefile, aes(fill = ifelse(id %in% matched_polygons, "Selecionado", "Não Selecionado")), alpha = 0.5, color = "black") +
  # Highlight unmatched polygons with light blue and dashed borders
  geom_sf(data = unmatched_polygons_layer, fill = "lightblue", alpha = 0.5, color = "black", linetype = "dashed") +
  # Plot all points
  geom_sf(data = points_shapefile, aes(color = ifelse(is.na(polygon_id), "Não Selecionado", "Selecionado")), size = 1.5) +
  # Highlight unmatched points
  geom_sf(data = unmatched_points_layer, color = "darkblue", size = 1.5, shape = 17) +
  # Set colors for matched and unmatched features
  scale_fill_manual(values = c("Selecionado" = "lightcoral", "Não Selecionado" = "lightblue")) +
  scale_color_manual(values = c("Selecionado" = "darkred", "Não Selecionado" = "darkblue"), labels = c("Estações Pluviométricas", "Não Selecionado")) +
  theme_minimal() +
  labs(fill = "Grid", color = "") +
  labs(title = "XAVIER - Grids selecionados",
       x = expression("Longitude (" * degree * ")"),
       y = expression("Latitude (" * degree * ")"))

p

# Save the plot with high quality
setwd("C:/Users/cassi/OneDrive/Documents/ANA/gitHUB/idf_RS_Technical_Note2/Figuras")
ggsave(filename = "Figura2a.png", plot = p, width = 10, height = 8, dpi = 500)





###########
#CHIRPS####FIGURA 2b
###########

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

# Plotting
p <- ggplot() +
  # Plot all polygons
  geom_sf(data = polygon_shapefile, aes(fill = ifelse(id %in% matched_polygons, "Selecionado", "Não Selecionado")), alpha = 0.5, color = "black") +
  # Highlight unmatched polygons with light blue and dashed borders
  geom_sf(data = unmatched_polygons_layer, fill = "lightblue", alpha = 0.5, color = "black", linetype = "dashed") +
  # Plot all points
  geom_sf(data = points_shapefile, aes(color = ifelse(is.na(polygon_id), "Não Selecionado", "Selecionado")), size = 3) +
  # Highlight unmatched points
  geom_sf(data = unmatched_points_layer, color = "darkblue", size = 3, shape = 17) +
  # Set colors for matched and unmatched features
  scale_fill_manual(values = c("Selecionado" = "lightcoral", "Não Selecionado" = "lightblue")) +
  scale_color_manual(values = c("Selecionado" = "darkred", "Não Selecionado" = "darkblue"), labels = c("Estações Pluviométricas", "Não Selecionado")) +
  theme_minimal() +
  labs(fill = "Grid", color = "") +
  labs(title = "CHIRPS - Grids selecionados",
       x = expression("Longitude (" * degree * ")"),
       y = expression("Latitude (" * degree * ")"))

p


# Save the plot with high quality
setwd("C:/Users/cassi/OneDrive/Documents/ANA/gitHUB/idf_RS_Technical_Note2/Figuras")
ggsave(filename = "Figura2b.png", plot = p, width = 10, height = 8, dpi = 500)


nrow(polygon_shapefile)

length(matched_polygon)
