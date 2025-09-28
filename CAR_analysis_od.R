# Load required libraries 
library(sp)
library(sf)
library(spdep)
library(CARBayes)
library(gstat)
library(raster)
library(dplyr)
library(purrr)

# Read brinjal volatility dataset with coordinates, convert types, and filter for years 2018–2024
df <- read.csv2("C:/Users/Abbin/Desktop/JNU/Paper-1/data/BRINJAL-CONDVOL-COORD(GPR).csv", sep = ',')
df$Date <- as.Date(df$Date, format = "%Y-%m-%d")
df$Latitude <- as.numeric(df$Latitude)
df$Longitude <- as.numeric(df$Longitude)
df$Volatility <- as.numeric(df$Volatility)
df$District <- toupper(df$District)
df_filter <- df %>% filter(format(Date, "%Y") %in% 2018:2024)

# Load boundary shapefile for India, subset to Odisha and create a raster template
india_map <- st_read("C:/Users/Abbinav Sankar/Desktop/JNU/Paper-1/gadm41_IND_shp/gadm41_IND_1.shp", quiet = TRUE)
odisha <- filter(india_map, NAME_1 == "Odisha")

r_template <- raster(extent(odisha), res = 0.01, crs = st_crs(odisha)$proj4string)

# Function to process a single date: fit CAR model, interpolate volatility with IDW, and return raster data
process_date <- function(date_val) {
  df_date <- df_filter[df_filter$Date == date_val, ]
  if (nrow(df_date) < 3) return(NULL)  # Skip if insufficient data points
  
  coordinates(df_date) <- ~Longitude + Latitude
  proj4string(df_date) <- CRS("+proj=longlat +datum=WGS84")
  
  coords <- coordinates(df_date)
  knn <- knearneigh(coords, k = min(2, nrow(df_date) - 1)) # Tried with k = 2, 3, 4, 5, 6
  nb <- knn2nb(knn)
  W <- nb2mat(nb, style = "B")
  W <- (W + t(W)) / 2
  
  df_car <- as.data.frame(df_date)
  df_car$region <- 1:nrow(df_car)

  model <- S.CARleroux(Volatility ~ 1, family = "gaussian", data = df_car,
                       W = W, burnin = 5000, n.sample = 20000)

  df_date@data$Fitted <- model$fitted.values
  idw_model <- gstat::gstat(formula = Fitted ~ 1, data = df_date, nmax = 3, set = list(idp = 2.0))
  raster_interp <- interpolate(r_template, idw_model)

  raster_df <- as.data.frame(raster_interp, xy = TRUE)
  colnames(raster_df)[3] <- "Volatility"
  raster_df$Date <- date_val
  return(raster_df)
}

# Apply the process_date function to all unique dates and collect successful outputs
unique_dates <- unique(df_filter$Date)
results <- map(unique_dates, safely(process_date))
successful_rasters <- transpose(results)$result %>% compact()
all_raster_data <- bind_rows(successful_rasters)

# Save the necessary files
write.csv(all_raster_data, "C:/Users/Abbin/Desktop/JNU/Paper-1/Export/all_raster_data_OD.csv", row.names = FALSE)

df_coords <- df %>% select(District, Latitude, Longitude) %>% distinct()
write.csv(df_coords[, c("District", "Longitude", "Latitude")], "C:/Users/Abbin/Desktop/JNU/Paper-1/Export/district_coords_OD.csv", row.names = FALSE)

st_write(odisha, "C:/Users/Abbin/Desktop/JNU/Paper-1/Export/odisha_boundary.geojson", delete_dsn = TRUE)