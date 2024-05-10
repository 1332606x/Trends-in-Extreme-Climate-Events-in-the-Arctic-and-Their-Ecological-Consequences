# Grouping into indivdual events and plotting intensity for ROS
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
data_folder <- "Arctic_Climate_Project/output/ROS_5mm_1t2m_sd"

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
alaska_data <- bind_rows(alaska_data_list)
northern_fennoscandia_data <- bind_rows(northern_fennoscandia_data_list)
svalbard_data <- bind_rows(svalbard_data_list)
yamal_data <- bind_rows(yamal_data_list)

group_dates <- function(data) {
  data %>%
    mutate(date = as.Date(date)) %>%
    arrange(lon, lat, date) %>%
    group_by(lat, lon, grp = cumsum(c(0, diff(date) > 1))) %>%
    summarise(
      start_date = min(date), 
      end_date = max(date), 
      event_length = n(),
      intensity = sum(rf * seq_along(rf)),
      total_rf = sum(rf) 
    ) %>%
    ungroup() %>%
    dplyr::select(-grp)
}

alaska_data <- group_dates(alaska_data)
northern_fennoscandia_data <- group_dates(northern_fennoscandia_data)
svalbard_data <- group_dates(svalbard_data)
yamal_data <- group_dates(yamal_data)

alaska_data$region <- "Northern Alaska"
northern_fennoscandia_data$region <- "Northern Fennoscandia"
svalbard_data$region <- "Svalbard"
yamal_data$region <- "Yamal"

season_assignment <- function(WW_data) {
  WW_data <- WW_data %>%
    mutate(
      month = month(start_date),  # Extract month from the start_date
      season = if_else(month %in% 1:3, year(start_date) - 1, year(start_date))  # Assign season value based on month
    ) %>%
    dplyr::select(-month)  # Remove the month column
}

alaska_data <- season_assignment(alaska_data)
northern_fennoscandia_data <- season_assignment(northern_fennoscandia_data)
svalbard_data <- season_assignment(svalbard_data)
yamal_data <- season_assignment(yamal_data)

combined_data <- rbind(alaska_data, northern_fennoscandia_data, svalbard_data, yamal_data)

combined_data$total_rf <- combined_data$total_rf * 1000

## PRESS TESTING

# Specify the path to the folder containing CSV files
data_folder <- "Arctic_Climate_Project/output/winter_season_precipitation_press"

# Get a list of all CSV files in the folder
data_files <- list.files(path = data_folder, pattern = "\\.csv$", full.names = TRUE)

# Apply the function to all files
data_list <- lapply(data_files, area_filter)

# Access the dataframes for each area and each year separately
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
  summarise(sf = mean(sf), 
            rf = mean(rf)) %>%
  mutate(ratio = rf/sf)

northern_fennoscandia_data_press <- northern_fennoscandia_data_press %>%
  group_by(year) %>%
  summarise(sf = mean(sf), 
            rf = mean(rf)) %>%
  mutate(ratio = rf/sf)

svalbard_data_press <- svalbard_data_press %>%
  group_by(year) %>%
  summarise(sf = mean(sf), 
            rf = mean(rf)) %>%
  mutate(ratio = rf/sf)

yamal_data_press <- yamal_data_press %>%
  group_by(year) %>%
  summarise(sf = mean(sf), 
            rf = mean(rf)) %>%
  mutate(ratio = rf/sf)

alaska_data_press$region <- "Northern Alaska"
northern_fennoscandia_data_press$region <- "Northern Fennoscandia"
svalbard_data_press$region <- "Svalbard"
yamal_data_press$region <- "Yamal"

combined_data_press <- rbind(alaska_data_press, northern_fennoscandia_data_press, svalbard_data_press, yamal_data_press)

combined_data_press <- combined_data_press %>%
  rename(season = year) %>%
  filter(season >= 1950 & season <= 2020)

combined_data_ROS <- full_join(combined_data, combined_data_press, c("season", "region"))

# Plotting
ROS_intensity_combined <- ggplot(combined_data_ROS, aes(x = start_date, y = total_rf, size = event_length, color = intensity)) +
  geom_point() +
  scale_size_continuous(range = c(1, 5), breaks = c(1, 5, 10, 15), labels = c("1 day", "5 days", "10 days", "15 days")) +
  scale_color_viridis_c(option = "mako", direction = -1, trans = "sqrt") +  # Using "mako" color palette
  labs(
    x = "Date",
    y = "Total Rainfall (mm)",
    size = "Event Length",
    color = "Intensity"
  ) +
  theme_minimal() +
  facet_grid(region ~ ., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # Center the title
    panel.border = element_rect(color = "black", fill = NA),  # Add border around each panel
    panel.spacing = unit(1, "lines"),  # Add space between panels
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 12),  # Increase size of axis text
    axis.title = element_text(size = 14, angle = 0),  # Increase size of axis text
    strip.text = element_text(size = 12),  # Increase size of facet labels
    legend.text = element_text(size = 12)
    ) +
  coord_cartesian(xlim = as.Date(c("1950-10-01", "2021-04-01"))) +
  geom_vline(xintercept = seq(as.Date("1950-01-01"), as.Date("2021-01-01"), by = "5 years"), linetype = "dashed", color = "grey")  # Add lines every 5 years

plot(ROS_intensity_combined)

legend_colors <- c("darkblue", "lightblue")  # Snow, Rain
legend_title <- "Precipitation Type"

combined_data_press$sf <- combined_data_press$sf * 1000
combined_data_press$rf <- combined_data_press$rf * 1000
combined_data_press$ratio <- combined_data_press$rf / combined_data_press$sf


# Climate press plot
gg_press <- ggplot(data = combined_data_press, aes(x = season)) +
  geom_line(aes(y = sf, color = "Snowfall")) +
  geom_line(aes(y = rf, color = "Rainfall")) +
  labs(
    x = "Winter Season",
    y = "Total Precipitation (mm)"  # Add units to y-axis label
  ) +
  theme_minimal() +
  facet_grid(region ~ ., scales = "free_y") +  # Adjusting scales for y-axis
  theme(
    plot.title = element_blank(),  # No title
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
  scale_color_manual(values = legend_colors, name = legend_title) +  # Manually set colors, labels, and title for legend
  geom_vline(xintercept = seq(min(combined_data_ROS$season), max(combined_data_ROS$season), by = 5), color = "gray", linetype = "dashed")  # Adding vertical lines every 5 years

plot(gg_press)

ROS_extremes_yearly <- combined_data_ROS %>%
  group_by(season, rf, sf, region) %>%
  summarise(
    mean_ros = mean(total_rf, na.rm = TRUE),
    intensity = mean(intensity, na.rm = TRUE),
    event_length = mean(event_length, na.rm = TRUE),
    num_events = mean(n()),
    mean_ratio = mean(ratio, na.rm = TRUE)
  )

ROS_extremes_decades <- combined_data_ROS %>%
  group_by(season, rf, sf, region) %>%
  summarise(
    mean_ros = mean(total_rf, na.rm = TRUE),
    intensity = mean(intensity, na.rm = TRUE),
    event_length = mean(event_length, na.rm = TRUE),
    num_events = mean(n()),
    mean_ratio = mean(ratio, na.rm = TRUE)
  ) %>% mutate(decade = floor(season/10)*10) %>% 
  group_by(decade, region) %>% 
  summarise(
    mean_ros = mean(mean_ros, na.rm = TRUE),
    intensity = mean(intensity, na.rm = TRUE),
    event_length = mean(event_length, na.rm = TRUE),
    num_events = sum(num_events, na.rm = TRUE),
    mean_sf = mean(sf, na.rm = TRUE),
    mean_rf = mean(rf, na.rm = TRUE), 
    mean_ratio = mean(mean_ratio, na.rm = TRUE)
  )  %>% 
  ungroup()  # Remove grouping

ROS_extremes_decades$decade_change <- (ROS_extremes_decades$decade - 1940) / 10

# STATS

gg_rf_time <- ggplot(data = combined_data_press, aes(x = season, y = rf)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Season",
    y = "Mean rainfall"  # Add units to y-axis label
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

plot(gg_rf_time)

gg_sf_time <- ggplot(data = combined_data_press, aes(x = season, y = sf)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Season",
    y = "Mean snowfall"  # Add units to y-axis label
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

plot(gg_sf_time)

gg_ratio_time <- ggplot(data = combined_data_press, aes(x = season, y = ratio)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Season",
    y = "Mean ratio"  # Add units to y-axis label
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

plot(gg_ratio_time)

gg_intensity_time <- ggplot(data = ROS_extremes_decades, aes(x = decade_change, y = intensity)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) + 
  labs(
    x = "Decade",
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

gg_length_time <- ggplot(data = ROS_extremes_decades, aes(x = decade_change, y = event_length)) +
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




