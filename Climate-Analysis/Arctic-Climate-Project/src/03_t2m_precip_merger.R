## Merging temperature, rainfall and snowdepth in 1 file for each year
## 30-01-24

library(dplyr)

# Specify the path to the folder containing CSV files
precip_folder <- "Arctic_Climate_Project/data/precipitation"

# Get a list of all CSV files in the folder
precip_files <- list.files(path = precip_folder, pattern = "\\.csv$", full.names = TRUE)

# Read and merge data from each file
for (precip_file in precip_files) {
  
  filename_without_path <- basename(precip_file)
  
  # Remove the ".csv" extension
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the format "era5_daily_2m_temperature_YEAR.csv")
  year <- as.numeric(filename_parts[length(filename_parts)])
  
  if(year >= 1950 & year <= 2021) {
    
    print(paste("processing", year, sep = " "))
    
    precip_data <- read.csv(precip_file) %>% 
      mutate(date = as.Date(time, format = "%Y-%m-%d"),
             rf = tp - sf) %>% 
      select(lon, lat, date, rf, sd)
    
    t2m_file <- paste("Arctic_Climate_Project/data/2m_temperature/era5_gridded_daily_max_t2m_", year, ".csv", sep = "")
    t2m_data <- read.csv(t2m_file) %>% 
      mutate(date = as.Date(date, format = "%Y-%m-%d"))
    
    t2m_data <- left_join(t2m_data, precip_data[, c("lon", "lat", "date", "rf", "sd")], by = c("lon", "lat", "date"))
    
    print("writing")
    
    csv_fname <- paste("Arctic_Climate_Project/data/t2m_rf_sd/era5_gridded_daily_t2m_rf_sd_", year, ".csv", sep="")
    
    # Write the filtered dataframe to a CSV file
    write.table(t2m_data, file = csv_fname, row.names = FALSE, sep = ",")
  }
  
}
