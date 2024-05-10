## Calculating the baselines for observational paper (16) in meta-analysis
## 28-01-24

library(dplyr)
library(sf)
library(lubridate)

# Specify the path to the folder containing CSV files
ROS_folder <- "Arctic_Climate_Project/data/t2m_rf"

# Get a list of all CSV files in the folder
ROS_files <- list.files(path = ROS_folder, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty dataframe to store the merged data
baseline_data <- data.frame(lon = numeric(0), lat = numeric(0), 
                            total_days_rf = numeric(0))

# Read and merge data from each file
for (ROS_file in ROS_files) {
  
  filename_without_path <- basename(ROS_file)
  
  # Remove the ".csv" extension
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the format "era5_daily_2m_temperature_YEAR.csv")
  year <- as.numeric(filename_parts[length(filename_parts)])
  
  print(paste("processing", year, sep = " "))
  
  if (year >= 1950 && year <= 1980) {
    
    ROS_data <- read.csv(ROS_file) %>%
      mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
      filter(
        t2m > 274.15,
        rf > 0.001,
        lon > 13,
        lon < 17,
        lat > 76,
        lat < 80,
        month(date) %in% c(12, 1, 2, 3)
      ) %>%
      group_by(lon, lat) %>%
      summarise(total_days_rf = n()) %>%
      select(lon, lat, total_days_rf)
    
    # Merge the summarized data into baseline_data
    baseline_data <- merge(baseline_data, ROS_data, all = TRUE)
  }
  
}

# Remove data which intersects with the ocean
ocean_polygons <- read_sf("Arctic_Climate_Project/data/ne_10m_ocean/ne_10m_ocean.shp")

results_coords <- st_as_sf(baseline_data, coords = c("lon", "lat"), crs = st_crs(ocean_polygons))

sf_use_s2(FALSE)

# Filter out any grid points which are not within the ocean
tmp <- sapply(st_intersects(results_coords, ocean_polygons), function(z) if (length(z) == 0) NA_integer_ else z[1])

# Remove points that don't intersect with the ocean and convert back to a table of coordinates
if (sum(!is.na(tmp)) > 0) {
  results_coords <- data.frame(st_coordinates(results_coords[is.na(tmp), ]))
} else {
  results_coords <- data.frame(st_coordinates(results_coords))
}

colnames(results_coords) <- c("lon", "lat")

print_this <- results_coords

results_coords <- left_join(results_coords, baseline_data, by = c("lon", "lat"))

rm(baseline_data)

print("calculating 90th percentile")

# Calculate the 90th percentile
result_summary <- results_coords %>%
  group_by(lon, lat) %>%
  summarise(total_days_rf = quantile(total_days_rf, 0.9, na.rm = TRUE))

print(sum(result_summary$total_days_rf) / unique(results_coords[, c("lat", "lon")]))

print("writing")

csv_fname <- ("Arctic_Climate_Project/output/1950_1980_gridded_rf_days.csv")

# Write the filtered dataframe to a CSV file
write.table(print_this, file = csv_fname, row.names = FALSE, sep = ",")
