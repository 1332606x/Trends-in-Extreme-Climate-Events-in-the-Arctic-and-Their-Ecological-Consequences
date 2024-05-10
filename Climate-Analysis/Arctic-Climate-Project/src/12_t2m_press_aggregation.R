## Aggregating daily data to winter season data for climate press (t2m)
## 20-04-24

library(dplyr)
library(lubridate)
library(purrr)

# Function to aggregate temperatures for winter season
aggregate_winter_season <- function(t2m_file_pair) {
  
  file1 <- t2m_file_pair[[1]]
  file2 <- t2m_file_pair[[2]]
  
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
  
  data1$date <- as.Date(data1$date)
  data2$date <- as.Date(data2$date)
  
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
    summarise(t2m = mean(t2m, na.rm = TRUE)) %>%
    select(lat, lon, t2m)
  
  # Specify the filename for the output CSV file
  csv_fname <- paste("Arctic_Climate_Project/output/mean_max_winter_season_t2m_press/era5_gridded_winter_season_t2m_", year, ".csv", sep="")
  
  print("writing")
  
  # Write the aggregated dataframe to a CSV file
  write.csv(aggregated_data, file = csv_fname, row.names = FALSE)
}


# Specify the folder containing the files
t2m_folder <- "Arctic_Climate_Project/data/2m_temperature"

# Get a list of all files in the folder
t2m_files <- list.files(t2m_folder, full.names = TRUE)

# Create pairs of consecutive files
t2m_file_pairs <- lapply(1:(length(t2m_files)-1), function(i) list(t2m_files[i], t2m_files[i + 1]))

# Process each pair of files
mapply(FUN = aggregate_winter_season, t2m_file_pairs, SIMPLIFY = FALSE)
