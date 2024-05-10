## Converting nc to csv for ROS
## 12-01-24

library(ncdf4)
library(tidyr)
library(lubridate)
library(dplyr)

# function

convert_to_csv <- function(sf_file) {
  
  # Extracting year
  # Specify the output CSV file with the year from the NetCDF file name
  filename_without_path <- basename(sf_file)
  
  # Remove the ".nc" extension
  filename_without_extension <- sub("\\.nc$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the format "era5_daily_2m_temperature_YEAR.nc")
  year <- filename_parts[length(filename_parts)]
  
  print(paste("processing", year, sep = " "))
  
  sd_file <- paste("Arctic_Climate_Project/data/precipitation/snow_depth/era5_daily_snow_depth_", year, ".nc", sep = "")
  tp_file <- paste("Arctic_Climate_Project/data/precipitation/total_precipitation/era5_daily_total_precipitation_", year, ".nc", sep = "")
  
  #Extracting data for snowfall
  
  # Read ERA5 NetCDF data
  nc_ds <- nc_open(sf_file)
  
  print(sf_file)
  print("sf file open")
  
  dim_lon <- ncvar_get(nc_ds, "longitude")
  dim_lat <- ncvar_get(nc_ds, "latitude")
  dim_time <- ncvar_get(nc_ds, "time")
  units <- ncatt_get(nc_ds, "time", "units")
  ntt <- dim(dim_time)
  
  # Find the indices corresponding to the desired latitude range
  lat_indices <- which(dim_lat >= 50)
  
  # Extract the latitude subset
  lat_subset <- dim_lat[lat_indices]
  
  # cropping tp 
  sf <- ncvar_get(nc_ds, "sf", start = c(1, lat_indices[1], 1), 
                   count = c(length(dim_lon), length(lat_indices), length(dim_time)), 
                   collapse_degen = FALSE)
  
  coords <- as.matrix(expand.grid(dim_lon, lat_subset, dim_time))
  
  print("creating sf dataframe")
  
  # Create sf dataframe
  nc_df <- data.frame(cbind(coords, sf))
  names(nc_df) <- c("lon", "lat",  "time", "sf")
  
  print("filtering months")
  
  nc_df <- nc_df %>%
    mutate(time = as.POSIXct(time * 3600, origin = '1900-01-01 00:00')) %>%
    filter(month(time) %in% c(10, 11, 12, 01, 02, 03))
  
  #aggregating sf
  print("aggregating sf")
  
  nc_df <- nc_df %>%
    group_by(round(lon), round(lat), time) %>%
    summarise(sf = mean(sf)) %>%
    rename(lon = `round(lon)`, lat = `round(lat)`)
  
  # Read ERA5 NetCDF data
  nc_ds <- nc_open(sd_file)
  
  print(sd_file)
  print("sd file open")
  
  dim_lon <- ncvar_get(nc_ds, "longitude")
  dim_lat <- ncvar_get(nc_ds, "latitude")
  dim_time <- ncvar_get(nc_ds, "time")
  units <- ncatt_get(nc_ds, "time", "units")
  ntt <- dim(dim_time)
  
  # Find the indices corresponding to the desired latitude range
  lat_indices <- which(dim_lat >= 50)
  
  # Extract the latitude subset
  lat_subset <- dim_lat[lat_indices]
  
  # cropping sd
  sd <- ncvar_get(nc_ds, "sd", start = c(1, lat_indices[1], 1), 
                  count = c(length(dim_lon), length(lat_indices), length(dim_time)), 
                  collapse_degen = FALSE)
  
  coords <- as.matrix(expand.grid(dim_lon, lat_subset, dim_time))
  
  print("creating dataframe")
  
  # Create a dataframe
  nc_df_sd <- data.frame(cbind(coords, sd))
  names(nc_df_sd) <- c("lon", "lat",  "time", "sd")
  
  print("filtering months")
  
  nc_df_sd <- nc_df_sd %>%
    mutate(time = as.POSIXct(time * 3600, origin = '1900-01-01 00:00')) %>%
    filter(month(time) %in% c(10, 11, 12, 01, 02, 03))
  
  #aggregating sd
  print("aggregating sd")
  
  nc_df_sd <- nc_df_sd %>%
    group_by(round(lon), round(lat), time) %>%
    summarise(sd = mean(sd)) %>%
    rename(lon = `round(lon)`, lat = `round(lat)`)
  
  nc_df <- merge(nc_df, nc_df_sd, by = c("lon", "lat", "time"), all.x = TRUE)
  rm(nc_df_sd)
  
  # Read ERA5 NetCDF data
  nc_ds <- nc_open(tp_file)
  
  print(tp_file)
  print("tp file open")
  
  dim_lon <- ncvar_get(nc_ds, "longitude")
  dim_lat <- ncvar_get(nc_ds, "latitude")
  dim_time <- ncvar_get(nc_ds, "time")
  units <- ncatt_get(nc_ds, "time", "units")
  ntt <- dim(dim_time)
  
  # Find the indices corresponding to the desired latitude range
  lat_indices <- which(dim_lat >= 50)
  
  # Extract the latitude subset
  lat_subset <- dim_lat[lat_indices]
  
  # cropping tp 
  tp <- ncvar_get(nc_ds, "tp", start = c(1, lat_indices[1], 1), 
                  count = c(length(dim_lon), length(lat_indices), length(dim_time)), 
                  collapse_degen = FALSE)
  
  coords <- as.matrix(expand.grid(dim_lon, lat_subset, dim_time))
  
  print("creating dataframe")
  
  # Create a dataframe
  nc_df_tp <- data.frame(cbind(coords, tp))
  names(nc_df_tp) <- c("lon", "lat",  "time", "tp")
  
  print("filtering months")
  
  nc_df_tp <- nc_df_tp %>%
    mutate(time = as.POSIXct(time * 3600, origin = '1900-01-01 00:00')) %>%
    filter(month(time) %in% c(10, 11, 12, 01, 02, 03))
  
  #aggregating tp
  print("aggregating tp")
  
  nc_df_tp <- nc_df_tp %>%
    group_by(round(lon), round(lat), time) %>%
    summarise(tp = mean(tp)) %>%
    rename(lon = `round(lon)`, lat = `round(lat)`)
  
  nc_df <- merge(nc_df, nc_df_tp, by = c("lon", "lat", "time"), all.x = TRUE)
  rm(nc_df_tp)
  
  csv_fname <- paste("Arctic_Climate_Project/data/precipitation/era5_gridded_daily_precipitation_", year, ".csv", sep="")
  
  print("writing")
  
  # Write the filtered dataframe to a CSV file
  write.table(nc_df, file = csv_fname, row.names = FALSE, sep = ",")
  
  # Close the NetCDF file
  nc_close(nc_ds) 
  
}

# Specify the path to the folder containing NetCDF files
sf_folder <- "Arctic_Climate_Project/data/precipitation/snowfall"

# Get a list of all NetCDF files in the folder
sf_files <- list.files(path = sf_folder, pattern = "\\.nc$", full.names = TRUE)

mapply(FUN = convert_to_csv, sf_files)
