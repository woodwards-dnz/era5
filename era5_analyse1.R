
# downloading ERA5 climate data

# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
# global climate data reanalysis
# 1940 -
# hourly
# uncertainty estimate dervied from 10 member ensemble
# name = "ERA5 hourly data on single levels from 1940 to present"
# 
# Main variables
# 

# simple R example
# https://confluence.ecmwf.int/display/CUSF/Download+CDS+ERA5+data+using+R
# ERA-Interim vs ERA5 vs ERA5-Land 
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=74764925
# ERA5-Land vs AgERA5
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=239338353#:~:text=ERA5%2DLand%20is%20a%20reanalysis,interpolated%20to%200.1%20degree%20resolution.
# Ethiopia example in R
# https://rpubs.com/Ajeep05/era5#:~:text=Download%20ERA5%20data%20using%20R&text=Set%20a%20key%20to%20a,profile%20on%20the%20ERA5%20webpage.&text=Get%20a%20command%20line%20request%20to%20provide%20the%20required%20details.
# ERA5-Land official overview
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.e2161bac?tab=overview
# ERA5-Land variable list
# https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation
# query arguments
# https://cds.climate.copernicus.eu/cdsapp#!/software/app-c3s-daily-era5-statistics?tab=overview
# query examples
# https://cds.climate.copernicus.eu/toolbox/doc/how-to/1_how_to_retrieve_data/1_how_to_retrieve_data.html

library(tidyr)
library(janitor)
library(sf)
library(lubridate)
library(ncdf4) 
library(raster)
library(rgdal) 
library(ggplot2) 
library(dplyr)
library(stringr)
library(foreach)
library(ecmwfr)

# variables ####
variables <- c("2m_temperature" = "t2m", # Kelvin 
               "2m_dewpoint_temperature" = "d2m", # Kelvin, Rh = 100 - 5(T-Td), Bonshoms et al 2022
               "surface_solar_radiation_downwards" = "ssrd", # J/m2  
               "10m_u_component_of_wind" = "u10", # m/s
               "10m_v_component_of_wind" = "v10" # m/s
)

# north island data ####
suffix = "N_12"
kelvin_to_celsius <- 273.16
t2m <- readRDS(paste0("downloaded/t2m_", suffix, ".rds")) - kelvin_to_celsius
d2m <- readRDS(paste0("downloaded/d2m_", suffix, ".rds")) - kelvin_to_celsius
rh2m <- 100 - 5 * (t2m - d2m) # Bonshoms et al 2022
ssrd <- readRDS(paste0("downloaded/ssrd_", suffix, ".rds")) / 1000000 # J/m2/h to MJ/m2/h
u10 <- readRDS(paste0("downloaded/u10_", suffix, ".rds"))
v10 <- readRDS(paste0("downloaded/v10_", suffix, ".rds"))
w2m <- sqrt(u10 ^ 2 + v10 ^ 2) * 0.75 # https://www.researchgate.net/post/Is-that-possible-to-convert-wind-speed-measured-in-10-m-height-to-a-possible-2-m-height-wind-speed
# get meta data 
nc_data <- nc_open(paste0("downloaded/2m_temperature_", suffix, ".nc"))
lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")
dt <- (ymd_hms("1900-01-01 00:00:00") + hours(t)) %>% with_tz("Pacific/Auckland")
nc_close(nc_data)
data_N <- data.frame(dt = dt) %>% 
  cross_join(data.frame(lat = lat)) %>% 
  cross_join(data.frame(lon = lon)) 
data_N <- data_N %>% 
  mutate(
    t2m = as.vector(t2m), 
    rh2m = as.vector(rh2m), 
    ssrd = as.vector(ssrd), 
    w2m = as.vector(w2m) 
  ) %>% 
  drop_na()
when <- ymd_hm("2022-02-12 18:00", tz = "Pacific/Auckland") %>% with_tz("UTC")
data_N %>% 
  dplyr::filter(dt == when) %>% 
  ggplot() +
  geom_point(aes(x = lon, y = lat, colour = t2m)) # yep!

# south island data ####
suffix = "S_12"
kelvin_to_celsius <- 273.16
t2m <- readRDS(paste0("downloaded/t2m_", suffix, ".rds")) - kelvin_to_celsius
d2m <- readRDS(paste0("downloaded/d2m_", suffix, ".rds")) - kelvin_to_celsius
rh2m <- 100 - 5 * (t2m - d2m) # Bonshoms et al 2022
ssrd <- readRDS(paste0("downloaded/ssrd_", suffix, ".rds")) / 1000000 # J/m2/h to MJ/m2/h
u10 <- readRDS(paste0("downloaded/u10_", suffix, ".rds"))
v10 <- readRDS(paste0("downloaded/v10_", suffix, ".rds"))
w2m <- sqrt(u10 ^ 2 + v10 ^ 2) * 0.75 # https://www.researchgate.net/post/Is-that-possible-to-convert-wind-speed-measured-in-10-m-height-to-a-possible-2-m-height-wind-speed
# get meta data 
nc_data <- nc_open(paste0("downloaded/2m_temperature_", suffix, ".nc"))
lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")
dt <- (ymd_hms("1900-01-01 00:00:00") + hours(t)) %>% with_tz("Pacific/Auckland")
nc_close(nc_data)
data_S <- data.frame(dt = dt) %>% 
  cross_join(data.frame(lat = lat)) %>% 
  cross_join(data.frame(lon = lon)) 
data_S <- data_S %>% 
  mutate(
    t2m = as.vector(t2m), 
    rh2m = as.vector(rh2m), 
    ssrd = as.vector(ssrd), 
    w2m = as.vector(w2m) 
  ) %>% 
  drop_na()
when <- ymd_hm("2022-02-12 18:00", tz = "Pacific/Auckland") %>% with_tz("UTC")
data_S %>% 
  dplyr::filter(dt == when) %>% 
  ggplot() +
  geom_point(aes(x = lon, y = lat, colour = t2m)) # yep!

# combine ####
data <- bind_rows(data_N, data_S)

