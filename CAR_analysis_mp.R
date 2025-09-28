# Load required libraries
library(sp)
library(sf)
library(spdep)
library(CARBayes)
library(gstat)
library(raster)
library(dplyr)
library(purrr)

# Read soybean volatility dataset with coordinates, convert types, and filter by year (2018–2024)
df <- read.csv2("data/SOYBEAN-CONDVOL-COORD(GPR).csv", sep = ',')
df$Date <- as.Date(df$Date, format = "%Y-%m-%d")
df$Latitude <- as.numeric(df$Latitude)
df$Longitude <- as.numeric(df$Longitude)
df$Volatility <- as.numeric(df$Volatility)
df$District <- toupper(df$District)
df_filter <- df %>% filter(format(Date, "%Y") %in% 2018:2024)

# Load boundary shapefile for India, subset to Madhya Pradesh and create a raster template
india_map <- st_read("gadm41_IND_shp/gadm41_IND_1.shp", quiet = TRUE)
madhya_pradesh <- filter(india_map, NAME_1 == "Madhya Pradesh")

r_template <- raster(extent(madhya_pradesh), res = 0.01, crs = st_crs(madhya_pradesh)$proj4string)

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
  idw_model <- gstat::gstat(formula = Fitted ~ 1, data = df_date, nmax = 7, set = list(idp = 2.0))
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
write.csv(all_raster_data, "Export/all_raster_data_mp.csv", row.names = FALSE)

df_coords <- df %>% select(District, Latitude, Longitude) %>% distinct()
write.csv(df_coords[, c("District", "Longitude", "Latitude")], "Export/district_coords_mp.csv", row.names = FALSE)

st_write(madhya_pradesh, "Export/mp_boundary.geojson", delete_dsn = TRUE)
