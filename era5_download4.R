
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

# access era5 ####
cds <- readLines("cds_key.dat") # creds from account page on website
wf_set_key(user = cds[1], # "Insert_your_CDS_UID_here", 
           key = cds[2], # "Insert_your_CDS_API_KEY_here"
           service = "cds")
catalogue <- wf_datasets(user = cds[1], service = "cds") %>% 
  dplyr::filter(str_detect(name, "era5"))

# connect era5-land ####
# dataset_short_name <- "reanalysis-era5-single-levels" # 31km 1940ff tons of atmospheric stuff
dataset_short_name <- "reanalysis-era5-land" # 9km 1950ff focused on land surface
dataset_info <- wf_product_info(dataset = dataset_short_name, user = cds[1], service = "cds")
names(dataset_info)
abstract <- dataset_info$rich_abstract

# get dyna_vcsn locations ####
crs_wgs <- 4326
dyna_vcsn <- read.csv("../deep_south/mockup/dyna_vcs_locations2.csv") %>% 
  clean_names() %>% 
  st_as_sf(coords = c("lon", "lat"), remove = FALSE) %>% 
  st_set_crs(crs_wgs) %>% 
  mutate(
    lat = round((lat - 0.025) * 20) / 20 + 0.025, # avoid rounding error
    lon = round((lon - 0.025) * 20) / 20 + 0.025, # avoid rounding error
    island = ifelse(lat < -40.0 & lon < 174.5, "S", "N") # separate islands to reduce unused space
  )
grid <- c(.2,.2) # rough grid
grid <- c(.5,.5) # rough grid
island <- "S"
# area is specified as N, W, S, E
area <- c(max(dyna_vcsn$lat[dyna_vcsn$island == island]), 
          min(dyna_vcsn$lon[dyna_vcsn$island == island]), 
          min(dyna_vcsn$lat[dyna_vcsn$island == island]), 
          max(dyna_vcsn$lon[dyna_vcsn$island == island])) # dyna only
area <- round(area / grid[1]) * grid[1] # round to grid

# set up query ####
variables <- c("2m_temperature" = "t2m", # Kelvin 
               "2m_dewpoint_temperature" = "d2m", # Kelvin, Rh = 100 - 5(T-Td), Bonshoms et al 2022
               "surface_solar_radiation_downwards" = "ssrd", # J/m2  
               "10m_u_component_of_wind" = "u10", # m/s
               "10m_v_component_of_wind" = "v10" # m/s
               )
hours <- sprintf("%02d:00", 0:23)
when <- ymd_hm("2022-02-12 18:00", tz = "Pacific/Auckland") %>% with_tz("UTC")
# when <- ymd_hm(paste("2022-02-12", hours), tz = "Pacific/Auckland") %>% with_tz("UTC")
years <- sprintf("%04d", year(when))
months <- sprintf("%02d", month(when))
months <- sprintf("%02d", 1:12) # full year
days <- sprintf("%02d", day(when))
days <- sprintf("%02d", 1:31) # full year
hours <- sprintf("%02d:00", hour(when))
hours <- sprintf("%02d:00", 0:23) # full year
nitems <- ceiling((area[4]-area[2])/grid[1]) * ceiling((area[1]-area[3])/grid[2]) * length(years) * length(months) * length(days) * length(hours)
# limit to 12000 items per query? recommend GRIB over netcdf  

# echo grid details ####
print(years)
print(months)
print(days)
print(hours)
print(island)
print(area)
print(grid)
print(variables)
print(nitems)

# loop through variables ####
variables <- variables[1:5] # subsetting
for (variable in names(variables)){
  
  suffix <- paste0(island, "_", max(months))
  fname <- paste0(variable, "_", suffix, ".nc")
  print(fname)
    
  request <- list(
    dataset_short_name = dataset_short_name,
    product_type   = "reanalysis",
    format = "netcdf",
    variable = variable,
    year = years,
    month = months,
    day = days,
    time = hours,
    # area is specified as N, W, S, E
    # area = c(-34.3, 172.5, -38.8, 177.5)
    area = area,
    grid = grid,
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
  dt <- (ymd_hms("1900-01-01 00:00:00") + hours(t)) %>% with_tz("Pacific/Auckland")
  length(t)
  vname <- unname(variables[variable])
  data <- ncvar_get(nc_data, vname) # 3d matrix lon * lat * t
  fillvalue <- ncatt_get(nc_data, vname, "_FillValue")
  data[data == fillvalue$value] <- NA
  print(paste("Assigning", vname))
  assign(vname, data)
  print(paste("Saving", vname))
  saveRDS(data, paste0("downloaded/", vname, "_", suffix, ".rds"))
  nc_close(nc_data)
  
}

kelvin_to_celsius <- 273.16
t2m <- readRDS(paste0("downloaded/t2m_", suffix, ".rds")) - kelvin_to_celsius
d2m <- readRDS(paste0("downloaded/d2m_", suffix, ".rds")) - kelvin_to_celsius
rh2m <- 100 - 5 * (t2m - d2m) # Bonshoms et al 2022
ssrd <- readRDS(paste0("downloaded/ssrd_", suffix, ".rds")) / 1000000 # J/m2/h to MJ/m2/h
u10 <- readRDS(paste0("downloaded/u10_", suffix, ".rds"))
v10 <- readRDS(paste0("downloaded/v10_", suffix, ".rds"))
w2m <- sqrt(u10 ^ 2 + v10 ^ 2) * 0.75 # https://www.researchgate.net/post/Is-that-possible-to-convert-wind-speed-measured-in-10-m-height-to-a-possible-2-m-height-wind-speed

vars <- c("t2m", "rh2m", "ssrd", "w2m")
for (var in vars){
  i <- which(dt == when)
  slice <- get(var)[,,i] # 2022-02-12 18:00
  crs_wgs <- 4326
  r <- raster(t(slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=crs_wgs)
  png(paste0("downloaded/", var, "_", suffix, ".png"))
  plot(r, main = var)
  dev.off()
}
