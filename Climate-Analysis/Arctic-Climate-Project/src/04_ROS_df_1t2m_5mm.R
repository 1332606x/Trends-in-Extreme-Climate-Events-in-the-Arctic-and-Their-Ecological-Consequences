# Producing a data frame with the ROS events for t2m > 274.15 and rf > 10mm
# Scale - 1 by 1 degree grid boxes for entire arctic region
# 30-01-24

library(dplyr)

ROS_count <- function(ROS_file) {

  result_data <- data.frame(lon = integer(), lat = integer(), date = character(), rf = numeric()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 )
  
  # Importing data
  year_data <- read.csv(ROS_file) %>%
    na.omit()  # Remove NAs
  
  # Extracting year
  filename_without_path <- basename(ROS_file)
  
  # Remove the ".nc" extension
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the same format
  year <- filename_parts[length(filename_parts)]
  
  print(paste("processing", year, sep = " "))
  
  # Filter data based on conditions for ROS event
  result_data <- year_data %>%
    filter(rf > 0.005, t2m > 274.15, sd > 0.003) %>%
    select(lon, lat, date, rf)
  
  if (nrow(result_data) > 0) {
    csv_fname <- paste("Arctic_Climate_Project/output/ROS_5mm_1t2m_sd/ROS_gridded_5mm_1t2m_sd_", year, ".csv", sep="")
    
    print("writing")
    
    # Write the filtered dataframe to a CSV file
    write.table(result_data, file = csv_fname, row.names = FALSE, sep = ",")
    
  }

}

# Specify the path to the folder containing NetCDF files
ROS_folder <- "Arctic_Climate_Project/data/t2m_rf_sd"

# Get a list of all NetCDF files in the folder
ROS_files <- list.files(path = ROS_folder, pattern = "\\.csv$", full.names = TRUE)

mapply(FUN = ROS_count, ROS_files)

