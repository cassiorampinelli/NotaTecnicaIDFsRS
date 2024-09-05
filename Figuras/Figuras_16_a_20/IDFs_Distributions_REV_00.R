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

#Set wd
setwd("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro")



##########################################################
#READING DATA AT SITES
##########################################################

#########################
#####MML#################
#########################


###GEV
#Reading polygon shapefile
points_shapefile.MML.GEV <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GEV_MML_HIDRO.shp")

###GAM
#Reading polygon shapefile
points_shapefile.MML.GAM <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GAM_MML_HIDRO.shp")


###GUM
#Reading polygon shapefile
points_shapefile.MML.GUM <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GUM_MML_HIDRO.shp")



#########################
#####MOM#################
#########################


###GEV
#Reading polygon shapefile
points_shapefile.MOM.GEV <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GEV_MOM_HIDRO.shp")

###GAM
#Reading polygon shapefile
points_shapefile.MOM.GAM <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GAM_MOM_HIDRO.shp")


###GUM
#Reading polygon shapefile
points_shapefile.MOM.GUM <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/enviados_Cassio/hidro/RS_GUM_MOM_HIDRO.shp")


points_shapefile<-list(points_shapefile.MML.GEV,points_shapefile.MML.GAM,points_shapefile.MML.GUM,points_shapefile.MOM.GEV,points_shapefile.MOM.GAM,points_shapefile.MOM.GUM)


# Define a CRS (for example, WGS 84)
crs_wgs84 <- 4326

##################################################
#Assigining CRS
##################################################

# Function to ensure CRS is assigned
assign_crs_if_missing <- function(sf_object, crs) {
  if (is.na(st_crs(sf_object))) {
    st_crs(sf_object) <- crs
  }
  return(sf_object)
}


i=1
st_crs(points_shapefile[[1]])

for(i in 1:length(points_shapefile)){
  
  
  points_shapefile[[i]] <- assign_crs_if_missing(points_shapefile[[i]], crs_wgs84)
  
}

i=1
st_crs(points_shapefile[[1]])




#Provide id references for points

for(i in 1:length(points_shapefile)){
  
  
  points_shapefile[[i]]$id<-seq(1:nrow(points_shapefile[[i]]))
  
}


#######- Computing Return Periods by using IDFs equations

##-Hidro- Station

# Define the function to calculate IDF for a rainfall of 24 hs
calculate_rainfall <- function(idf_a, idf_b, idf_c,idf_d, Tr) {
  
  (idf_a*(Tr)^idf_b)/(24*60+idf_c)^idf_d
}

# User-specified values for Tr
Tr_values <- c(2,5,25,100,500)

#Record data
# Initialize an empty list to store results
stations_results_df <- list()
stations_results_df.temp<-0

# Define the categories
categories <- c("MOM.GEV", "MOM.GAM", "MOM.GUM", "MML.GEV", "MML.GAM", "MML.GUM")



for (i in 1:length(points_shapefile)) {
  # Apply the function for each value of Tr and store results
  stations_results.temp <- lapply(Tr_values, function(Tr) {
    points_shapefile[[i]] %>%
      mutate(rainfall = calculate_rainfall(idf_a, idf_b, idf_c, idf_d, Tr)) %>%
      mutate(Tr_values = as.character(Tr)) %>%
      dplyr::select(codigo, lon, lat, rainfall, Tr_values)
  })
  
  # Combine the results into a single data frame
  stations_results_df.temp <- bind_rows(stations_results.temp) %>%
    mutate(category = categories[i]) # Add the category column
  
  # Store the combined result in the list
  stations_results_df[[i]] <- stations_results_df.temp
}

# Combine all data frames in the list into a single data frame
final_results_df <- bind_rows(stations_results_df)

# Print the combined data frame
print(final_results_df)




# Calculate the maximum, minimum, and difference between rainfall values for each codigo and Tr_values
rainfall_stats_df <- final_results_df %>%
  group_by(codigo, Tr_values) %>%
  summarise(
    max_rainfall = max(rainfall),
    min_rainfall = min(rainfall),
    rainfall_diff = max_rainfall - min_rainfall,
    percentage_diff = (rainfall_diff / min_rainfall) * 100
  ) %>%
  ungroup()


# Assuming final_results_df is already available with percentage_diff calculated

# Convert to sf object if not already
rainfall_stats_sf <- st_as_sf(rainfall_stats_df, coords = c("lon", "lat"), crs = st_crs(final_results_df))


#Reading polygon shapefile
boundary_sf <- st_read("C:/Users/cassi/OneDrive/Documents/ANA/IDF_Project/Rio_Grande_do_Sul/GIS/RS_UF.shp")


setwd("C:\\Users\\cassi\\OneDrive\\Documents\\ANA\\IDF_Project\\Rio_Grande_do_Sul\\enviados_Cassio\\Figuras\\all_distributions")


# Assign CRS if missing
boundary_sf <- assign_crs_if_missing(boundary_sf, crs_wgs84)

# Create a list of unique Tr_values for plotting
unique_Tr_values <- unique(rainfall_stats_sf$Tr_values)

# Plot for each Tr_value
for (Tr in unique_Tr_values) {
  # Filter data for the current Tr_value
  plot_data <- rainfall_stats_sf %>% filter(Tr_values == Tr)
  
  # Generate the plot
  plot <- ggplot(data = plot_data) +
    geom_sf(data = boundary_sf, fill = NA, color = "black", linetype = "dashed",linewidth=1) +  # Plot the boundary
    geom_sf(aes(size = percentage_diff, color = percentage_diff), show.legend = TRUE) +
    scale_size_continuous(name = "Diferença %",guide="none") +
    scale_color_gradient(low = "lightblue", high = "red", name = "Diferença %") +
    labs(title = "Diferença Percentual entre os Maiores e Menores Valores das Distribuições",
         subtitle = paste("Tempo de Retorno =", Tr,"anos"),
         color = "Diferença %") +
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
  
  # Save the plot
  ggsave(filename = paste0("Percentage_Difference_Tr=", Tr, ".png"), plot = plot, width = 10, height = 8, dpi = 500)
}


