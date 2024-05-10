# Filtering ROS to arctic land points only and exporting plots
# 14-01-24

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(sp)
library(raster) 
library(vegan)
library(geosphere)
library(s2dv)

# Read the shapefile for the Arctic
arctic_shape <- st_read("Arctic_Climate_Project/data/evidence-map-scope/evidence-map-scope.shp")

# Clean the geometry of the Arctic shapefile
arctic_shape <- st_make_valid(arctic_shape)
arctic_shape <- st_union(arctic_shape)

# Create a bounding box and subtract the arctic_shape to get the complement
bounding_box <- st_bbox(arctic_shape) %>%
  st_as_sfc() %>%
  st_sf(geometry = .)

arctic_shape <- st_difference(bounding_box, arctic_shape)

# Importing ocean shapefile 
ocean_polygons <- read_sf("Arctic_Climate_Project/data/ne_10m_ocean/ne_10m_ocean.shp")

coastline <- read_sf("Arctic_Climate_Project/data/ne_10m_coastline/ne_10m_coastline.shp")

# Define bounding box
bounding_box <- st_bbox(c(xmin = -180, xmax = 180, ymin = 50, ymax = 90), crs = st_crs(coastline))
# Crop coastline data using bounding box
coastline <- st_crop(coastline, bounding_box)

# Function to calculate average ROS events per year
yearly_gridded_count <- function(ROS_file) {
  
  # Extracting year
  filename_without_path <- basename(ROS_file)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- as.integer(filename_parts[length(filename_parts)])
  
  print(paste("processing", year, sep = " "))
    
  # Read the CSV file with lon-lat points
  ROS_data <- read.csv(ROS_file)
    
  # Convert lon values to -180 to 180 range
  ROS_data$lon <- ifelse(ROS_data$lon > 180, ROS_data$lon - 360, ROS_data$lon)
  
  ROS_data <- ROS_data[ROS_data$lat > 50, ]
    
  total_rf <- ROS_data %>%
    group_by(lat, lon) %>%
    summarise(total_rf = sum(rf))
    
    # Add year information to the result
  total_rf$year <- year
    
  return(total_rf)
  
}

# Specify the path to the folder containing CSV files
ROS_folder <- "Arctic_Climate_Project/output/ROS_10mm_1t2m_sd"

# Get a list of all CSV files in the folder
ROS_files <- list.files(path = ROS_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
ROS_data_list <- lapply(ROS_files, yearly_gridded_count)

# Combine the results into a single data frame
ROS_data <- bind_rows(ROS_data_list)

ROS_data$total_rf <- ROS_data$total_rf * 1000

## Filtering out land and out of arctic data for baseline
filter_fun <- function(global_data, ocean_polygons, arctic_shape) {
  
  coords_df <- global_data[, c("lon", "lat")]
  coords_df <- distinct(coords_df, lat, lon)
  
  coords_df <- st_as_sf(coords_df, coords=c("lon","lat"))
  
  st_crs(coords_df) <- st_crs(ocean_polygons)
  sf_use_s2(FALSE)
  
  ##find where out points intersect with the ocean
  tmp <- sapply(st_intersects(coords_df, ocean_polygons), function(z) if (length(z)==0) NA_integer_ else z[1])
  
  if (sum(!is.na(tmp))>0) {
    coords_df<-data.frame(st_coordinates(coords_df[is.na(tmp),]))} else {
      coords_df<-data.frame(st_coordinates(coords_df))}
  
  colnames(coords_df) <- c("lon","lat")
  
  coords_df <- st_as_sf(coords_df, coords=c("lon","lat"))
  
  st_crs(coords_df) <- st_crs(arctic_shape)
  sf_use_s2(FALSE)
  
  ##find where out points intersect with the arctic area
  tmp <- sapply(st_intersects(coords_df, arctic_shape), function(z) if (length(z)==0) NA_integer_ else z[1])
  
  if (sum(!is.na(tmp))>0) {
    coords_df <- data.frame(st_coordinates(coords_df[is.na(tmp),]))} else {
      coords_df <- data.frame(st_coordinates(coords_df))}
  
  colnames(coords_df) <- c("lon","lat")
  
  filtered_data <- global_data %>%
    inner_join(coords_df, by = c("lat", "lon"))  
  
  return(filtered_data)

}

filtered_ROS <- filter_fun(ROS_data, ocean_polygons, arctic_shape)
rm(ROS_data)

filtered_baseline <- filtered_ROS %>%
  subset(year > 1950 & year <= 1980) %>%
  group_by(lon, lat) %>%
  summarise(ROS_baseline = sum(total_rf)/30)

filtered_current <- filtered_ROS %>%
  subset(year > 1990 & year <= 2020) %>%
  group_by(lon, lat) %>%
  summarise(ROS_current = sum(total_rf)/30)

ROS_absolute_change <- full_join(filtered_baseline, filtered_current, by = c("lon", "lat"))
rm(filtered_baseline, filtered_current)

ROS_absolute_change[is.na(ROS_absolute_change)] <- 0

ROS_absolute_change <- ROS_absolute_change %>%
  mutate(Absolute_diff = ROS_current - ROS_baseline)

legend_limits <- c(0, 100)

gg_baseline <- ggplot() +
  geom_tile(data = ROS_absolute_change, aes(x = lon, y = lat, fill = ROS_baseline), width = 1, height = 1) +
  geom_sf(data = coastline, fill = "transparent", color = "black") +
  geom_sf(data = arctic_shape, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(mako(50)), limits = legend_limits, trans = "sqrt") +
  geom_tile(data = subset(ROS_absolute_change, ROS_baseline > 100), 
            aes(x = lon, y = lat), 
            fill = "black", width = 1, height = 1) +
  theme_minimal() +
  #labs(title = "Average mm of ROS Between 1950-1980 (Square Root Transformation)",
   #    fill = "") + 
  theme(
    plot.title = element_blank(),
    legend.key.size = unit(1, "cm"),  # Adjust the size of legend color key
    legend.text = element_text(size = 12),  # Adjust the size of legend text
    legend.title = element_blank(),  # Adjust the size of legend title
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.spacing.y = unit(0.1, "cm")  # Adjust vertical spacing in the legend
  ) 

#+
#  coord_sf(crs = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m +no_defs")

plot(gg_baseline)

gg_current <- ggplot() +
  geom_tile(data = ROS_absolute_change, aes(x = lon, y = lat, fill = ROS_current), width = 1, height = 1) +
  geom_sf(data = coastline, fill = "transparent", color = "black") +  # Add coastline
  geom_sf(data = arctic_shape, fill = "transparent", color = "black") +
  geom_tile(data = subset(ROS_absolute_change, ROS_current > 100), 
            aes(x = lon, y = lat), 
            fill = "black", width = 1, height = 1) +
  scale_fill_gradientn(colors = rev(mako(50)), limits = legend_limits, trans = "sqrt") +
  theme_minimal() +
 # labs(title = "Average mm of ROS Between 1990-2020 (Square Root Transformation)",
  #     fill = "") + 
  theme(
    plot.title = element_blank(),
    legend.key.size = unit(1, "cm"),  # Adjust the size of legend color key
    legend.text = element_text(size = 12),  # Adjust the size of legend text
    legend.title = element_text(size = 14),  # Adjust the size of legend title
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.spacing.y = unit(0.1, "cm")  # Adjust vertical spacing in the legend
  )

plot(gg_current)

# Specify the RdBu color palette
rdBu_palette <- colorRampPalette(brewer.pal(11, "RdBu"))
legend_limits_change <- c(-5, 5)

study_locations <- read.csv("Arctic_Climate_Project/data/study_locations.csv")

gg_change <- ggplot() +
  geom_tile(data = ROS_absolute_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  geom_sf(data = coastline, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change, name = "Absolute difference (mm)") +
  geom_tile(data = subset(ROS_absolute_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(ROS_absolute_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_sf(data = arctic_shape, fill = "transparent", color = "black") +
  geom_point(data = subset(ROS_absolute_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 5, shape = "*") +
  geom_point(data = subset(study_locations), 
             aes(x = Lon, y = Lat), color = "black", size =2, shape = 17) +
  geom_rect(aes(xmin = -168, xmax = -141, ymin = 66, ymax = 72),
            fill = "#CCCCCC50", color = "black", linetype = "dashed")  + # Alaska
  geom_rect(aes(xmin = 16, xmax = 32, ymin = 67, ymax = 71),
            fill = "#CCCCCC50", color = "black", linetype = "dashed") + # Northern Fennoscandia
   geom_rect(aes(xmin = 10, xmax = 35, ymin = 76, ymax = 82),
           fill = "#CCCCCC50", color = "black", linetype = "dashed") + # Svalbard
  geom_rect(aes(xmin = 64, xmax = 82, ymin = 64, ymax = 74),
           fill = "#CCCCCC50", color = "black", linetype = "dashed") + # Yamal
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.spacing.y = unit(0.1, "cm")
    )

plot(gg_change)

bounding_box_alaska <- st_bbox(c(xmin = -168, xmax = -141, ymin = 66, ymax = 72), crs = st_crs(coastline))
alaska_coastline <- st_crop(coastline, bounding_box_alaska)

alaska_change <- subset(ROS_absolute_change, 
                        lat >= 66 & lat <= 72 & 
                          lon >= -168 & lon <= -141)

bounding_box_northern_fennoscandia <- st_bbox(c(xmin = 16, xmax = 32, ymin = 67, ymax = 71), crs = st_crs(coastline))
northern_fennoscandia_coastline <- st_crop(coastline, bounding_box_northern_fennoscandia)

northern_fennoscandia_change <- subset(ROS_absolute_change, 
                         lat >= 67 & lat <= 71 & 
                           lon >= 16 & lon <= 32)

bounding_box_svalbard <- st_bbox(c(xmin = 10, xmax = 35, ymin = 76, ymax = 82), crs = st_crs(coastline))
svalbard_coastline <- st_crop(coastline, bounding_box_svalbard)

svalbard_change <- subset(ROS_absolute_change, 
                          lat >= 76 & lat <= 82 & 
                            lon >= 10 & lon <= 35)

bounding_box_yamal <- st_bbox(c(xmin = 64, xmax = 82, ymin = 64, ymax = 74), crs = st_crs(coastline))
yamal_coastline <- st_crop(coastline, bounding_box_yamal)

yamal_change <- subset(ROS_absolute_change, 
                       lat >= 64 & lat <= 73 & 
                         lon >= 64 & lon <= 82)

gg_change_alaska <- ggplot() +
  geom_tile(data = alaska_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  geom_sf(data = alaska_coastline, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(alaska_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_point(data = subset(alaska_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  ) 

plot(gg_change_alaska)

gg_change_northern_fennoscandia <- ggplot() +
  geom_tile(data = northern_fennoscandia_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  geom_sf(data = northern_fennoscandia_coastline, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_point(data = subset(northern_fennoscandia_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  ) 

plot(gg_change_northern_fennoscandia)

gg_change_svalbard <- ggplot() +
  geom_tile(data = svalbard_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  geom_sf(data = svalbard_coastline, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(svalbard_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(svalbard_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_point(data = subset(svalbard_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  )

plot(gg_change_svalbard)

gg_change_yamal <- ggplot() +
  geom_tile(data = yamal_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  geom_sf(data = yamal_coastline, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(yamal_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(yamal_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_point(data = subset(yamal_change, ROS_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 15, shape = "*") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Add dashed grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.title = element_blank(),
    axis.text.x = element_blank(),  # Rotate x-axis labels
    axis.text.y = element_blank(),  # Adjust y-axis labels alignment
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove axis ticks
    plot.title = element_blank(),  # Remove the plot title
    legend.position = "none"  # Remove the legend
  )

gg_change_yamal

plot_name <- "Arctic_Climate_Project/output/iceland.jpg"
ggsave(plot_name, plot = gg_change_iceland, width = 20, height = 7, units = "in", dpi = 300)

write.csv(filtered_ROS, "Arctic_Climate_Project/output/mantel_test/yearly_total_ROS_df.csv", row.names = FALSE)

## MANTEL TEST

coords <- ROS_absolute_change[c("lon", "lat")]
spatial_dist <- distm(coords)

# Compute distance matrix for ROS and time
ROS_absolute_dist <- dist(ROS_absolute_change$Absolute_diff)

mantel_test_result <- mantel(ROS_absolute_dist, spatial_dist, method = "pearson", permutations = 40)

print(mantel_test_result)

correlation_value <- mantel_test_result$statistic
space_dist_vector <- as.vector(spatial_dist[lower.tri(spatial_dist)])

densities <- densCols(space_dist_vector, ROS_absolute_dist, colramp = colorRampPalette(c("black", "white")))

# Extract densities
dens_values <- col2rgb(densities)[1,] + 1L

# Create a dataframe with spatial distance, exceedance distance, and densities
df <- data.frame(space_dist = space_dist_vector,
                 ROS_absolute_dist = ROS_absolute_dist,
                 dens = dens_values)

# Plot it
mantel_plot <- ggplot(df, aes(x = space_dist, y = ROS_absolute_dist, color = dens)) +
  geom_point(size = 2) +
  scale_color_viridis(option = "inferno") +  # Using viridis color palette
  labs(
    x = "Spatial Distance",
    y = "Rain-On-Snow Absolute Change Distance"
    ) +
  theme_minimal()

plot(mantel_plot)

## PAIRED T-TEST

result <- t.test(ROS_absolute_change$ROS_current, ROS_absolute_change$ROS_baseline, paired = TRUE)
result_alaska <- t.test(alaska_change$ROS_current, alaska_change$ROS_baseline, paired = TRUE)
result_northern_fennoscandia <- t.test(northern_fennoscandia_change$ROS_current, northern_fennoscandia_change$ROS_baseline, paired = TRUE)
result_svalbard <- t.test(svalbard_change$ROS_current, svalbard_change$ROS_baseline, paired = TRUE)
result_yamal <- t.test(yamal_change$ROS_current, yamal_change$ROS_baseline, paired = TRUE)

print(result)
print(result_alaska)
print(result_northern_fennoscandia)
print(result_svalbard)
print(result_yamal)



