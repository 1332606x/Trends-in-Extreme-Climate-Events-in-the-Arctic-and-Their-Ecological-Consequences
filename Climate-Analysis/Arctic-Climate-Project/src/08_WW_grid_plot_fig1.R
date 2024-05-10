# Filtering WW to arctic land points only and exporting plots
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

# Function to calculate average WW events per year
yearly_gridded_count <- function(WW_file) {
  
  # Extracting year
  filename_without_path <- basename(WW_file)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- as.integer(filename_parts[length(filename_parts)])
  
  print(paste("processing", year, sep = " "))
  
  # Read the CSV file with lon-lat points
  WW_data <- read.csv(WW_file)
  
  # Convert lon values to -180 to 180 range
  WW_data$lon <- ifelse(WW_data$lon > 180, WW_data$lon - 360, WW_data$lon)
  
  WW_data <- WW_data[WW_data$lat > 50, ]
  
  total_WW <- WW_data %>%
    group_by(lat, lon) %>%
    summarise(total_exceedance = sum(exceedance))
  
  # Add year information to the result
  total_WW$year <- year
  
  return(total_WW)
  
}

# Specify the path to the folder containing CSV files
WW_folder <- "Arctic_Climate_Project/output/WW_gridded_t2m"

# Get a list of all CSV files in the folder
WW_files <- list.files(path = WW_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
WW_data_list <- lapply(WW_files, yearly_gridded_count)

# Combine the results into a single data frame
WW_data <- bind_rows(WW_data_list)

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

filtered_WW <- filter_fun(WW_data, ocean_polygons, arctic_shape)
rm(WW_data)

filtered_baseline <- filtered_WW %>%
  subset(year > 1950 & year <= 1980) %>%
  group_by(lon, lat) %>%
  summarise(WW_baseline = sum(total_exceedance)/30)

filtered_current <- filtered_WW %>%
  subset(year > 1990 & year <= 2020) %>%
  group_by(lon, lat) %>%
  summarise(WW_current = sum(total_exceedance)/30)

WW_absolute_change <- full_join(filtered_baseline, filtered_current, by = c("lon", "lat"))
rm(filtered_baseline, filtered_current)

WW_absolute_change[is.na(WW_absolute_change)] <- 0

WW_absolute_change <- WW_absolute_change %>%
  mutate(Absolute_diff = WW_current - WW_baseline)

legend_limits <- c(0, 25)

gg_baseline <- ggplot() +
  geom_tile(data = WW_absolute_change, aes(x = lon, y = lat, fill = WW_baseline), width = 1, height = 1) +
  geom_sf(data = coastline, fill = "transparent", color = "black") +
  geom_sf(data = arctic_shape, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(inferno(50)), limits = legend_limits) + 
  geom_tile(data = subset(WW_absolute_change, WW_baseline > 25), 
            aes(x = lon, y = lat), 
            fill = "black", width = 1, height = 1) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    legend.key.size = unit(1, "cm"),  # Adjust the size of legend color key
    legend.text = element_text(size = 12),  # Adjust the size of legend text
    legend.title = element_text(size = 14),  # Adjust the size of legend title
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.spacing.y = unit(0.1, "cm")  # Adjust vertical spacing in the legend
  )

plot(gg_baseline)

gg_current <- ggplot() +
  geom_tile(data = WW_absolute_change, aes(x = lon, y = lat, fill = WW_current), width = 1, height = 1) +
  geom_sf(data = coastline, fill = "transparent", color = "black") +  # Add coastline
  geom_sf(data = arctic_shape, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(inferno(50)), limits = legend_limits) + 
  geom_tile(data = subset(WW_absolute_change, WW_current > 25), 
            aes(x = lon, y = lat), 
            fill = "black", width = 1, height = 1) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    legend.key.size = unit(1, "cm"),  # Adjust the size of legend color key
    legend.text = element_text(size = 12),  # Adjust the size of legend text
    legend.title = element_blank(),  # Adjust the size of legend title
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.spacing.y = unit(0.1, "cm")  # Adjust vertical spacing in the legend
  )

plot(gg_current)

# Specify the RdBu color palette
rdBu_palette <- colorRampPalette(brewer.pal(11, "RdBu"))
legend_limits_change <- c(-10, 10)

gg_change <- ggplot() +
  geom_tile(data = WW_absolute_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  geom_sf(data = coastline, fill = "transparent", color = "black") +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(WW_absolute_change, Absolute_diff < -10), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(WW_absolute_change, Absolute_diff > 10), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_sf(data = arctic_shape, fill = "transparent", color = "black") +
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

alaska_change <- subset(WW_absolute_change, 
                        lat >= 66 & lat <= 72 & 
                          lon >= -168 & lon <= -141)

bounding_box_northern_fennoscandia <- st_bbox(c(xmin = 16, xmax = 32, ymin = 67, ymax = 71), crs = st_crs(coastline))
northern_fennoscandia_coastline <- st_crop(coastline, bounding_box_northern_fennoscandia)

northern_fennoscandia_change <- subset(WW_absolute_change, 
                                       lat >= 67 & lat <= 71 & 
                                         lon >= 16 & lon <= 32)

bounding_box_svalbard <- st_bbox(c(xmin = 10, xmax = 35, ymin = 76, ymax = 82), crs = st_crs(coastline))
svalbard_coastline <- st_crop(coastline, bounding_box_svalbard)

svalbard_change <- subset(WW_absolute_change, 
                          lat >= 76 & lat <= 82 & 
                            lon >= 10 & lon <= 35)

bounding_box_yamal <- st_bbox(c(xmin = 64, xmax = 82, ymin = 64, ymax = 74), crs = st_crs(coastline))
yamal_coastline <- st_crop(coastline, bounding_box_yamal)

yamal_change <- subset(WW_absolute_change, 
                       lat >= 64 & lat <= 73 & 
                         lon >= 64 & lon <= 82)

gg_change_alaska <- ggplot() +
  geom_tile(data = alaska_change, aes(x = lon, y = lat, fill = Absolute_diff), width = 1, height = 1) +
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(alaska_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(alaska_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_sf(data = alaska_coastline, fill = "transparent", color = "black") +
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
  scale_fill_gradientn(colors = rev(rdBu_palette(100)), limits = legend_limits_change) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff < -5), 
            aes(x = lon, y = lat), 
            fill = "#1A1F43", width = 1, height = 1) +
  geom_tile(data = subset(northern_fennoscandia_change, Absolute_diff > 5), 
            aes(x = lon, y = lat), 
            fill = "firebrick", width = 1, height = 1) +
  geom_point(data = subset(northern_fennoscandia_change, WW_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 5, shape = "*") +
  geom_sf(data = northern_fennoscandia_coastline, fill = "transparent", color = "black") +
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
  geom_point(data = subset(svalbard_change, WW_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 5, shape = "*") +
  geom_sf(data = svalbard_coastline, fill = "transparent", color = "black") +
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
  geom_point(data = subset(yamal_change, WW_baseline == 0), 
             aes(x = lon, y = lat), color = "darkred", size = 5, shape = "*") +
  geom_sf(data = yamal_coastline, fill = "transparent", color = "black") +
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

write.csv(filtered_WW, "Arctic_Climate_Project/output/mantel_test/yearly_total_WW_df.csv", row.names = FALSE)

## MANTEL TEST

coords <- WW_absolute_change[c("lon", "lat")]
spatial_dist <- distm(coords)

# Compute distance matrix for WW and time
WW_absolute_dist <- dist(WW_absolute_change$Absolute_diff)

mantel_test_result <- mantel(WW_absolute_dist, spatial_dist, method = "pearson", permutations = 99)

print(mantel_test_result)

# Calculate distance matrix based on latitude only
latitude_matrix <- dist(WW_absolute_change$lat)

mantel_test_result_lat <- mantel(WW_absolute_dist, latitude_matrix, method = "pearson", permutations = 20)

# View the Mantel test results
print(mantel_test_result_lat)

correlation_value <- mantel_test_result$statistic
space_dist_vector <- as.vector(spatial_dist[lower.tri(spatial_dist)])

densities <- densCols(space_dist_vector, WW_absolute_dist, colramp = colorRampPalette(c("black", "white")))

# Extract densities
dens_values <- col2rgb(densities)[1,] + 1L

# Create a dataframe with spatial distance, exceedance distance, and densities
df <- data.frame(space_dist = space_dist_vector,
                 WW_absolute_dist = WW_absolute_dist,
                 dens = dens_values)

# Plot it
mantel_plot <- ggplot(df, aes(x = space_dist, y = WW_absolute_dist, color = dens)) +
  geom_point(size = 2) +
  scale_color_viridis(option = "inferno") +  # Using viridis color palette
  labs(
    x = "Spatial Distance",
    y = "Cumulative Winter Warming Absolute Change Distance"
  ) +
  theme_minimal()

plot(mantel_plot)

plot_name <- "Arctic_Climate_Project/output/Mantel_test_WW.jpg"
ggsave(plot_name, plot = mantel_plot, width = 20, height = 7, units = "in", dpi = 300)


## PAIRED T-TEST

result <- t.test(WW_absolute_change$WW_current, WW_absolute_change$WW_baseline, paired = TRUE)
result_alaska <- t.test(alaska_change$WW_current, alaska_change$WW_baseline, paired = TRUE)
result_northern_fennoscandia <- t.test(northern_fennoscandia_change$WW_current, northern_fennoscandia_change$WW_baseline, paired = TRUE)
result_svalbard <- t.test(svalbard_change$WW_current, svalbard_change$WW_baseline, paired = TRUE)
result_yamal <- t.test(yamal_change$WW_current, yamal_change$WW_baseline, paired = TRUE)

print(result)
print(result_alaska)
print(result_northern_fennoscandia)
print(result_svalbard)
print(result_yamal)

