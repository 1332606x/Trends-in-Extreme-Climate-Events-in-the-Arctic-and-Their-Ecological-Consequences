# Calculating the 99th percentile for temperature from baseline years for each
# grid box
# 06-03-24

library(dplyr)

# Specify the path to the folder containing CSV files
t2m_folder <- "Arctic_Climate_Project/data/2m_temperature"

# Get a list of all CSV files in the folder
t2m_files <- list.files(path = t2m_folder, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty dataframe to store the merged data
baseline_data <- data.frame(lon = numeric(0), lat = numeric(0), t2m = numeric(0))

# Read and merge data from each file
for (t2m_file in t2m_files) {
  
  filename_without_path <- basename(t2m_file)
  
  # Remove the ".csv" extension
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the format "era5_daily_2m_temperature_YEAR.csv")
  year <- as.numeric(filename_parts[length(filename_parts)])
  
  print(paste("processing", year, sep = " "))
  
  if(year > 1950 && year <= 1980) {
    t2m_data <- read.csv(t2m_file)
    baseline_data <- merge(baseline_data, t2m_data[, c("lon", "lat", "t2m")], all = TRUE)
    rm(t2m_data)
  }
 
}

print("calculating 99th percentile")

# Calculate the 99th percentile for each latitude and longitude
result_summary <- baseline_data %>%
  group_by(lon, lat) %>%
  summarise(p99_temperature = quantile(t2m, 0.99, na.rm = TRUE))

rm(baseline_data)

print("writing")

csv_fname <- ("Arctic_Climate_Project/output/1951_1980_gridded_0.99_percentile.csv")

# Write the filtered dataframe to a CSV file
write.table(result_summary, file = csv_fname, row.names = FALSE, sep = ",")
