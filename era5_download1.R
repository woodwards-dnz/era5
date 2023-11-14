
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

# https://confluence.ecmwf.int/display/CUSF/Download+CDS+ERA5+data+using+R
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=239338353#:~:text=ERA5%2DLand%20is%20a%20reanalysis,interpolated%20to%200.1%20degree%20resolution.
# https://confluence.ecmwf.int/pages/viewpage.action?pageId=74764925
# https://www.ecmwf.int/en/computing/software/ecmwf-web-api
# https://rpubs.com/Ajeep05/era5#:~:text=Download%20ERA5%20data%20using%20R&text=Set%20a%20key%20to%20a,profile%20on%20the%20ERA5%20webpage.&text=Get%20a%20command%20line%20request%20to%20provide%20the%20required%20details.
# list of variables
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.e2161bac?tab=overview

library(lubridate)
library(ncdf4) 
library(raster)
library(rgdal) 
library(ggplot2) 
library(dplyr)
library(stringr)
library(foreach)
library(ecmwfr)

cds <- readLines("cds_key.dat") # creds from account page on website
wf_set_key(user = cds[1], # "Insert_your_CDS_UID_here", 
           key = cds[2], # "Insert_your_CDS_API_KEY_here"
           service = "cds")
catalogue <- wf_datasets(user = cds[1], service = "cds") %>% 
  dplyr::filter(str_detect(name, "era5"))
# dataset_short_name <- "reanalysis-era5-single-levels" # 31km 1940ff tons of atmospheric stuff
dataset_short_name <- "reanalysis-era5-land" # 9km 1950ff focused on land surface
dataset_info <- wf_product_info(dataset = dataset_short_name, user = cds[1], service = "cds")
names(dataset_info)
abstract <- dataset_info$rich_abstract
variables <- c("2m_temperature")
when <- ymd_hm("2022-02-12 18:00", tz = "Pacific/Auckland") %>% with_tz("UTC")
years <- sprintf("%04d", year(when))
months <- sprintf("%02d", month(when))
days <- sprintf("%02d", day(when))
times <- c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00")
times <- sprintf("%02d:00", hour(when))
fname <- "download.nc"

request <- list(
  dataset_short_name = dataset_short_name,
  product_type   = "reanalysis",
  format = "netcdf",
  variable = variables[1],
  year = years[1],
  month = months[1],
  day = days[1],
  time = times,
  # area is specified as N, W, S, E
  # area = c(-34.3, 172.5, -38.8, 177.5), # `1_uni` = c(xmin = 172.5, ymin = -38.8, xmax = 177.5, ymax = -34.3)
  area = c(-34.12972, 166.42555, -48.06417, 178.55083), # `5_nz` = c(xmin = 166.42555, ymin = -48.06417, xmax = 178.55083, ymax = -34.12972)
  target = fname
)

file <- wf_request(user = cds[1], # "Insert_your_CDS_UID_here",
                   request = request,
                   transfer = TRUE,
                   path = "downloaded",
                   verbose = TRUE)

# https://rpubs.com/boyerag/297592

# get view meta data ####
nc_data <- nc_open(paste0("downloaded/", fname))
{
  sink(str_replace(paste0("downloaded/", fname), "\\.nc$", ".txt"))
  print(nc_data)
  sink()
}
lon <- ncvar_get(nc_data, "longitude")
length(lon)
lat <- ncvar_get(nc_data, "latitude", verbose = F)
length(lat)
t <- ncvar_get(nc_data, "time")
length(t)
t2m <- ncvar_get(nc_data, "t2m") # 3d matrix lon * lat * t
fillvalue <- ncatt_get(nc_data, "t2m", "_FillValue")
nc_close(nc_data)
t2m[t2m == fillvalue$value] <- NA
kelvin_to_celsius <- 273.16
slice <- t2m[,] - kelvin_to_celsius
max(slice, na.rm = TRUE)
min(slice, na.rm = TRUE)
crs_wgs <- 4326
r <- raster(t(slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
            crs=crs_wgs)
# r <- flip(r, direction='y')
plot(r)
