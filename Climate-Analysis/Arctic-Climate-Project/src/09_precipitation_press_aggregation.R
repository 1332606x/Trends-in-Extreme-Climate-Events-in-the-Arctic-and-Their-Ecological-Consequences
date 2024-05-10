## Aggregating daily data to winter season data for climate press (precipitation)
## 20-04-24

library(dplyr)
library(lubridate)
library(purrr)

# Function to aggregate temperatures for winter season
aggregate_winter_season <- function(precip_file_pair) {
  
  file1 <- precip_file_pair[[1]]
  file2 <- precip_file_pair[[2]]
  
  print(file1)
  print(file2)
  
  filename_without_path <- basename(file1)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- as.integer(filename_parts[length(filename_parts)])
  print(paste("season being processed:", year, sep = " "))
  
  # Read data from file1 and file2
  data1 <- read.csv(file1)
  data2 <- read.csv(file2)
  
  data1$date <- as.Date(data1$time)
  data2$date <- as.Date(data2$time)
  
  print("filtering")
  
  # Filter data for months 10, 11, 12 from file1 and 1, 2, 3 from file2
  winter_data <- rbind(
    data1 %>% filter(month(date) %in% c(10, 11, 12)),
    data2 %>% filter(month(date) %in% c(1, 2, 3))
  )
  
  print("aggregating")
  
  # Aggregate temperatures for winter season
  aggregated_data <- winter_data %>%
    group_by(lat, lon) %>%
    summarise(sf = sum(sf),
              rf = sum(tp) - sum(sf)) %>%
    select(lat, lon, rf, sf)
  
  # Specify the filename for the output CSV file
  csv_fname <- paste("Arctic_Climate_Project/output/winter_season_precipitation_press/era5_gridded_winter_season_precipitation_", year, ".csv", sep="")
  
  print("writing")
  
  # Write the aggregated dataframe to a CSV file
  write.csv(aggregated_data, file = csv_fname, row.names = FALSE)
}

# Specify the folder containing the files
precip_folder <- "Arctic_Climate_Project/data/precipitation"

# Get a list of all files in the folder
precip_files <- list.files(precip_folder, full.names = TRUE)

# Create pairs of consecutive files
precip_file_pairs <- lapply(1:(length(precip_files)-1), function(i) list(precip_files[i], precip_files[i + 1]))

# Process each pair of files
mapply(FUN = aggregate_winter_season, precip_file_pairs, SIMPLIFY = FALSE)
