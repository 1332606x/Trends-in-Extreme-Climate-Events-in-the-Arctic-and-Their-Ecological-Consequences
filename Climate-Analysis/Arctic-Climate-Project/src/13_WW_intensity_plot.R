# Grouping into indivdual events and plotting intensity for WW
# 15-03-24

library(ggplot2)
library(viridis)
library(sf)
library(dplyr)
library(gridExtra)
library(lubridate)
library(ggpubr)
library(ggpmisc)
library(lmtest)
library(car)
library(tseries)
library(trend)

ocean_polygons <- read_sf("Arctic_Climate_Project/data/ne_10m_ocean/ne_10m_ocean.shp")

area_filter <- function(data_file) {
  
  # Extracting year
  filename_without_path <- basename(data_file)
  filename_without_extension <- sub("\\.csv$", "", filename_without_path)
  filename_parts <- strsplit(filename_without_extension, "_")[[1]]
  year <- as.integer(filename_parts[length(filename_parts)])
  
  print(paste("Processing", year))
  
  # Read the CSV file with lon-lat points
  current_data <- read.csv(data_file)
  
  current_data$year <- year
  
  # Convert lon values to -180 to 180 range
  current_data$lon <- ifelse(current_data$lon > 180, current_data$lon - 360, current_data$lon)
  
  # Convert lon-lat points to sf object
  data_coords <- st_as_sf(current_data, coords = c("lon", "lat"), crs = st_crs(ocean_polygons))
  
  sf_use_s2(FALSE)
  
  # Find where our points intersect with the ocean
  tmp <- sapply(st_intersects(data_coords, ocean_polygons), function(z) if (length(z) == 0) NA_integer_ else z[1])
  
  # Remove points that intersect with the ocean and convert back to a table of coordinates
  if (sum(!is.na(tmp)) > 0) {
    data_coords <- data.frame(st_coordinates(data_coords[is.na(tmp), ]))
  } else {
    data_coords <- data.frame(st_coordinates(data_coords))
  }
  colnames(data_coords) <- c("lon", "lat")
  
  # Extract unique lon-lat combinations from data_coords
  unique_coords <- data_coords %>% 
    distinct(lon, lat)
  
  # Filter current_data based on unique lon-lat combinations for each region
  alaska_subset <- current_data %>%
    filter(lat >= 66 & lat <= 72 & lon >= -168 & lon <= -141)
  
  northern_fennoscandia_subset <- current_data %>%
    filter(lat >= 67 & lat <= 71 & lon >= 16 & lon <= 32)
  
  svalbard_subset <- current_data %>%
    filter(lat >= 76 & lat <= 82 & lon >= 10 & lon <= 35)
  
  yamal_subset <- current_data %>%
    filter(lat >= 64 & lat <= 74 & lon >= 64 & lon <= 82)
  
  return(
    list(
      alaska = alaska_subset,
      northern_fennoscandia = northern_fennoscandia_subset,
      svalbard = svalbard_subset,
      yamal = yamal_subset
    )
  )
}

# Specify the path to the folder containing CSV files
data_folder <- "Arctic_Climate_Project/output/WW_gridded_t2m"

# Get a list of all CSV files in the folder
data_files <- list.files(path = data_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
data_list <- lapply(data_files, area_filter)

# Access the dataframes for each area and each year separately
# Filter the data for each year based on the specified conditions
alaska_data_list <- lapply(data_list, function(data) data$alaska)
northern_fennoscandia_data_list <- lapply(data_list, function(data) data$northern_fennoscandia)
svalbard_data_list <- lapply(data_list, function(data) data$svalbard)
yamal_data_list <- lapply(data_list, function(data) data$yamal)

# Combine the lists into dataframes
alaska_data_WW <- bind_rows(alaska_data_list)
northern_fennoscandia_data_WW <- bind_rows(northern_fennoscandia_data_list)
svalbard_data_WW <- bind_rows(svalbard_data_list)
yamal_data_WW <- bind_rows(yamal_data_list)

group_dates <- function(data) {
  data %>%
    mutate(date = as.Date(date)) %>%
    arrange(lon, lat, date) %>%
    group_by(lat, lon, grp = cumsum(c(0, diff(date) > 1))) %>%
    summarise(
      start_date = min(date),  # Add a check for non-empty data
      end_date = max(date),  # Add a check for non-empty data
      event_length = n(),
      intensity = ifelse(n() > 0, sum(exceedance * seq_along(exceedance)), NA),  # Add a check for non-empty data
      max_t2m = ifelse(n() > 0, max(t2m) - 273.15, NA)  # Add a check for non-empty data
    ) %>%
    ungroup() %>%
    dplyr::select(-grp)
}

alaska_data_WW <- group_dates(alaska_data_WW)
northern_fennoscandia_data_WW <- group_dates(northern_fennoscandia_data_WW)
svalbard_data_WW <- group_dates(svalbard_data_WW)
yamal_data_WW <- group_dates(yamal_data_WW)

alaska_data_WW$region <- "Northern Alaska"
northern_fennoscandia_data_WW$region <- "Northern Fennoscandia"
svalbard_data_WW$region <- "Svalbard"
yamal_data_WW$region <- "Yamal"

season_assignment <- function(WW_data) {
  WW_data <- WW_data %>%
    mutate(
      month = month(start_date),  # Extract month from the start_date
      season = if_else(month %in% 1:3, year(start_date) - 1, year(start_date))  # Assign season value based on month
    ) %>%
    dplyr::select(-month)  # Remove the month column
}

alaska_data_WW <- season_assignment(alaska_data_WW)
northern_fennoscandia_data_WW <- season_assignment(northern_fennoscandia_data_WW)
svalbard_data_WW <- season_assignment(svalbard_data_WW)
yamal_data_WW <- season_assignment(yamal_data_WW)

combined_data <- rbind(alaska_data_WW, northern_fennoscandia_data_WW, svalbard_data_WW, yamal_data_WW)

## PRESS

# Specify the path to the folder containing CSV files
data_folder <- "Arctic_Climate_Project/output/mean_max_winter_season_t2m_press"

# Get a list of all CSV files in the folder
data_files <- list.files(path = data_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
data_list <- lapply(data_files, area_filter)

# Access the dataframes for each area and each year separately
# Filter the data for each year based on the specified conditions
alaska_data_list <- lapply(data_list, function(data) data$alaska)
northern_fennoscandia_data_list <- lapply(data_list, function(data) data$northern_fennoscandia)
svalbard_data_list <- lapply(data_list, function(data) data$svalbard)
yamal_data_list <- lapply(data_list, function(data) data$yamal)

# Combine the lists into dataframes
alaska_data_press <- bind_rows(alaska_data_list)
northern_fennoscandia_data_press <- bind_rows(northern_fennoscandia_data_list)
svalbard_data_press <- bind_rows(svalbard_data_list)
yamal_data_press <- bind_rows(yamal_data_list)

alaska_data_press <- alaska_data_press %>%
  group_by(year) %>%
  summarise(t2m = mean(t2m))

northern_fennoscandia_data_press <- northern_fennoscandia_data_press %>%
  group_by(year) %>%
  summarise(t2m = mean(t2m))

svalbard_data_press <- svalbard_data_press %>%
  group_by(year) %>%
  summarise(t2m = mean(t2m))

yamal_data_press <- yamal_data_press %>%
  group_by(year) %>%
  summarise(t2m = mean(t2m))

alaska_data_press$region <- "Northern Alaska"
northern_fennoscandia_data_press$region <- "Northern Fennoscandia"
svalbard_data_press$region <- "Svalbard"
yamal_data_press$region <- "Yamal"

WW_combined_data_press <- rbind(alaska_data_press, northern_fennoscandia_data_press, svalbard_data_press, yamal_data_press)
WW_combined_data_press$t2m <- WW_combined_data_press$t2m - 273.15

WW_combined_data_press <- WW_combined_data_press %>%
  rename(season = year)

WW_combined_data <- full_join(combined_data, WW_combined_data_press, by = c("season", "region"))

# Intensity plot
ww_intensity_combined <- ggplot(WW_combined_data, aes(x = start_date, y = max_t2m, size = event_length, color = intensity)) +
  geom_point() +
  scale_size_continuous(range = c(1, 5), breaks = c(1, 10, 20, 29), labels = c("1 day", "10 days", "20 days", "30 days")) +
  scale_color_viridis_c(option = "rocket", direction = -1, trans = "sqrt") +  # Square root transformation
  labs(
    x = "Date",
    y = "Max Temperature (°C)",  # Add units to y-axis label
    size = "Event Length",
    color = "Intensity"
  ) +
  theme_minimal() +
  facet_grid(region ~ ., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 12),  # Increase size of axis text
    axis.title = element_text(size = 14, angle = 0),  # Increase size of axis text
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 12, angle = 90)  # Increase size of facet labels
  ) +
  coord_cartesian(xlim = as.Date(c("1950-10-01", "2021-04-01")))  +
  geom_vline(xintercept = seq(as.Date("1950-01-01"), as.Date("2021-01-01"), by = "5 years"), linetype = "dashed", color = "grey")  # Add lines every 5 years

plot(ww_intensity_combined)

# Climate press plot
gg_press <- ggplot(data = WW_combined_data_press, aes(x = season, y = t2m)) +
  geom_line(color = "darkred") +
  #  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) +  # Adding a linear regression line
  labs(
    x = "Winter Season",
    y = "Mean Maximum Daily Winter Temperature (°C)"  # Add units to y-axis label
  ) +
  theme_minimal() +
  facet_grid(region ~., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # Center the title
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 12),  # Increase size of axis text
    axis.title = element_text(size = 14, angle = 0),  # Increase size of axis text
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 12, angle = 90)  # Increase size of facet labels
  ) +
  # stat_poly_eq(
  #   aes(
  #     label = paste(
  #       after_stat(eq.label), "*\" with \"*", 
  #       after_stat(rr.label), "*\", \"*", 
  #       after_stat(p.value.label), "*\".\"",
  #       sep = ""
  #     )
  #   ),
  #   formula = y ~ x, 
  #   size = 3
  # ) +
  geom_vline(xintercept = seq(min(WW_combined_data_press$season), max(WW_combined_data_press$season), by = 5), color = "gray", linetype = "dashed")

 plot(gg_press)

plot_name <- "Arctic_Climate_Project/output/Plots/WW_Intensity_plots.jpg"
ggsave(plot_name, plot = ww_intensity_combined, width = 20, height = 9, units = "in", dpi = 300)

plot_name <- "Arctic_Climate_Project/output/Plots/Temperature_press_plot.jpg"
ggsave(plot_name, plot = gg_press, width = 20, height = 9, units = "in", dpi = 300)

## Visualising trends

climate_extremes_season <- WW_combined_data %>%
  group_by(lat, lon, season, t2m, region) %>%
  summarise(
    max_t2m = mean(max_t2m, na.rm = TRUE),
    sum_intensity = mean(intensity, na.rm = TRUE),
    sum_event_length = mean(event_length, na.rm = TRUE),
    num_events = n()
  ) %>% mutate(decade = floor(season/10)*10) %>% 
  group_by(decade, region) %>% 
  summarise(
    max_t2m = mean(max_t2m, na.rm = TRUE),
    sum_intensity = mean(sum_intensity, na.rm = TRUE),
    sum_event_length = mean(sum_event_length, na.rm = TRUE),
    num_events = sum(num_events, na.rm = TRUE),
    mean_t2m = mean(t2m, na.rm = TRUE)
  ) %>%
  ungroup()  # Remove grouping

climate_extremes_season$decade_change <- (climate_extremes_season$decade - 1950) / 10

gg_press_time <- ggplot(data = climate_extremes_season, aes(x = decade_change, y = mean_t2m)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Season",
    y = "Mean Maximum Daily Temperature"  # Add units to y-axis label
  ) +
  theme_minimal() +
  facet_grid(region ~., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # Center the title
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black")  # Add axis lines
  ) +
  stat_poly_eq(
    aes(
      label = paste(
        after_stat(eq.label), "*\" with \"*", 
        after_stat(rr.label), "*\", \"*", 
        after_stat(p.value.label), "*\".\"",
        sep = ""
      )
    ),
    formula = y ~ x,  
    size = 3, 
    parse = TRUE, 
    coef.digits = 4,
    coef.keep.zeros = TRUE,
    rr.digits = 4,
    f.digits = 4,
    p.digits = 4
  ) 

plot(gg_press_time)

gg_intensity_time <- ggplot(data = climate_extremes_season, aes(x = decade_change, y = sum_intensity)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Winter Season",
    y = "Cumulative Intensity"  # Add units to y-axis label
  ) +
  theme_minimal() +
  facet_grid(region ~., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # Center the title
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black")  # Add axis lines
  ) +
  stat_poly_eq(
    aes(
      label = paste(
        after_stat(eq.label), "*\" with \"*", 
        after_stat(rr.label), "*\", \"*", 
        after_stat(p.value.label), "*\".\"",
        sep = ""
      )
    ),
    formula = y ~ x,
    size = 3,
    parse = TRUE, 
    coef.digits = 4,
    coef.keep.zeros = TRUE,
    rr.digits = 4,
    f.digits = 4,
    p.digits = 4
  ) 

plot(gg_intensity_time)

gg_length_time <- ggplot(data = climate_extremes_season, aes(x = decade_change, y = sum_event_length)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Winter Season",
    y = "Event length"  # Add units to y-axis label
  ) +
  theme_minimal() +
  facet_grid(region ~., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # Center the title
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black")  # Add axis lines
  ) +
  stat_poly_eq(
    aes(
      label = paste(
        after_stat(eq.label), "*\" with \"*", 
        after_stat(rr.label), "*\", \"*", 
        after_stat(p.value.label), "*\".\"",
        sep = ""
      )
    ),
    formula = y ~ x,
    size = 3,
    parse = TRUE, 
    coef.digits = 4,
    coef.keep.zeros = TRUE,
    rr.digits = 4,
    f.digits = 4,
    p.digits = 4
  ) 

plot(gg_length_time)

gg_frequency_time <- ggplot(data = climate_extremes_season, aes(x = mean_t2m, y = sum_event_length)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Winter Season",
    y = "Cumulative Intensity"  # Add units to y-axis label
  ) +
  theme_minimal() +
  facet_grid(region ~., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # Center the title
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black")  # Add axis lines
  ) +
  stat_poly_eq(
    aes(
      label = paste(
        after_stat(eq.label), "*\" with \"*", 
        after_stat(rr.label), "*\", \"*", 
        after_stat(p.value.label), "*\".\"",
        sep = ""
      )
    ),
    formula = y ~ x,
    size = 3
  ) 

plot(gg_frequency_time)
