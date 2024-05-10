## Converting nc to csv for t2m
## 10-01-24

library(ncdf4)
library(tidyr)
library(lubridate)
library(dplyr)
library(sf)

# function
convert_to_csv <- function(nc_file) {
  
  print(nc_file)
  
  # Read ERA5 NetCDF data
  nc_ds <- nc_open(nc_file)
  
  print("file open")
  
  dim_lon <- ncvar_get(nc_ds, "longitude")
  dim_lat <- ncvar_get(nc_ds, "latitude")
  dim_time <- ncvar_get(nc_ds, "time")
  units <- ncatt_get(nc_ds, "time", "units")
  ntt <- dim(dim_time)
  
  # Find the indices corresponding to the desired latitude range
  lat_indices <- which(dim_lat > 50)
  
  # Extract the latitude subset
  lat_subset <- dim_lat[lat_indices]
  
  # cropping t2m 
  t2m <- ncvar_get(nc_ds, "t2m", start = c(1, lat_indices[1], 1), 
                   count = c(length(dim_lon), length(lat_indices), length(dim_time)), 
                   collapse_degen = FALSE)
  
  coords <- as.matrix(expand.grid(dim_lon, lat_subset, dim_time))
  
  print("creating dataframe")
  
  # Create a dataframe
  nc_df <- data.frame(cbind(coords, t2m))
  
  names(nc_df) <- c("lon", "lat",  "time", "t2m")
  
  print("filtering months")
  
  nc_df <- nc_df %>%
    mutate(time = as.POSIXct(time * 3600, origin = '1900-01-01 00:00')) %>%
    filter(month(time) %in% c(10, 11, 12, 01, 02, 03))
  
  print("filtering out ocean")
  
  # Convert lon values to -180 to 180 range
  nc_df$lon <- ifelse(nc_df$lon > 180, nc_df$lon - 360, nc_df$lon)
  
  # Convert lon-lat points to sf object
  count_coords <- st_as_sf(nc_df, coords = c("lon", "lat"), crs = st_crs(ocean_polygons))
  
  sf_use_s2(FALSE)
  
  # Filter out any grid points which are not within the Arctic boundary
  tmp <- sapply(st_intersects(count_coords, ocean_polygons), function(z) if (length(z) == 0) NA_integer_ else z[1])
  
  # Remove points that don't intersect with the Arctic shapefile and convert back to a table of coordinates
  if (sum(!is.na(tmp)) > 0) {
    count_coords <- data.frame(st_coordinates(count_coords[is.na(tmp), ]))
  } else {
    count_coords <- data.frame(st_coordinates(count_coords))
  }
  colnames(count_coords) <- c("lon", "lat")
  
  nc_df <- left_join(count_coords, nc_df, by = c("lon", "lat"))
  
  print("aggregating")
  
  # aggregating data into 1 degree by 1 degree grid
  nc_df <- nc_df %>%
    group_by(round(lon), round(lat), date(time)) %>%
    summarise(t2m = max(t2m)) %>%
    rename(lon = `round(lon)`, lat = `round(lat)`, date = `date(time)`)
  
  # Specify the output CSV file with the year from the NetCDF file name
  filename_without_path <- basename(nc_file)
  
  # Remove the ".nc" extension
  filename_without_extension <- sub("\\.nc$", "", filename_without_path)
  
  # Split the filename using "_" as a separator
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  
  # Extract the year (assuming it's always in the format "era5_daily_2m_temperature_YEAR.nc")
  year <- filename_parts[length(filename_parts)]
  
  csv_fname <- paste("Arctic_Climate_Project/data/2m_temperature/era5_gridded_daily_max_t2m_", year, ".csv", sep="")
  
  print("writing")
  
  # Write the filtered dataframe to a CSV file
  write.table(nc_df, file = csv_fname, row.names = FALSE, sep = ",")
  
  # Close the NetCDF file
  nc_close(nc_ds)
  print(paste(year, "complete", sep = " "))
  
}

ocean_polygons <- read_sf("Arctic_Climate_Project/data/ne_110m_ocean/ne_110m_ocean.shp")

# Specify the path to the folder containing NetCDF files
nc_folder <- "/soge-home/data/analysis/era5/0.28125x0.28125/daily/2m_temperature/nc"

# Get a list of all NetCDF files in the folder
nc_files <- list.files(path = nc_folder, pattern = "\\.nc$", full.names = TRUE)

mapply(FUN = convert_to_csv, nc_files)
