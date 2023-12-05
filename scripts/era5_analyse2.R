
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

library(GGally)
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

# function s####
myggsave <- function(fname, height = 210, width = 297, units = "mm", ...){
  print(paste("Plotting", fname))
  ggsave(fname, height = height, width = width, units = units, ...)
}
safemean <- function(x){
  mean(x, na.rm = !all(is.na(x)))
}
safemin <- function(x){
  min(x, na.rm = !all(is.na(x)))
}
safemax <- function(x){
  max(x, na.rm = !all(is.na(x)))
}
safesum <- function(x){
  sum(x, na.rm = !all(is.na(x)))
}
safeget <- function(x, i){ # get last element with i == TRUE, i is a boolean vector
  j <- i & !is.na(x)
  ifelse(length(j) == 0, NA_real_, tail(x[j], 1))
}

# variables ####
variables <- c("2m_temperature" = "t2m", # Kelvin 
               "2m_dewpoint_temperature" = "d2m", # Kelvin, Rh = 100 - 5(T-Td), Bonshoms et al 2022
               "surface_solar_radiation_downwards" = "ssrd", # J/m2 cumulative 
               "10m_u_component_of_wind" = "u10", # m/s
               "10m_v_component_of_wind" = "v10" # m/s
)

# get ERA5-Land data ####
# north island data ####
suffix = "N_12"
kelvin_to_celsius <- 273.16
t2m <- readRDS(paste0("downloaded/t2m_", suffix, ".rds")) - kelvin_to_celsius
d2m <- readRDS(paste0("downloaded/d2m_", suffix, ".rds")) - kelvin_to_celsius
rh2m <- 100 - 5 * (t2m - d2m) # Bonshoms et al 2022
ssrd <- readRDS(paste0("downloaded/ssrd_", suffix, ".rds")) / 1000000 # J/m2/h to MJ/m2/h (cumulative)
u10 <- readRDS(paste0("downloaded/u10_", suffix, ".rds"))
v10 <- readRDS(paste0("downloaded/v10_", suffix, ".rds"))
w2m <- sqrt(u10 ^ 2 + v10 ^ 2) * 0.75 # https://www.researchgate.net/post/Is-that-possible-to-convert-wind-speed-measured-in-10-m-height-to-a-possible-2-m-height-wind-speed
# get meta data 
nc_data <- nc_open(paste0("downloaded/2m_temperature_", suffix, ".nc"))
lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")
dt <- (ymd_hms("1900-01-01 00:00:00") + hours(t)) # %>% with_tz("Pacific/Auckland")
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
  drop_na() %>% 
  dplyr::filter(lat > -41.5 | lon > 174) # remove mistaken bit
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
ssrd <- readRDS(paste0("downloaded/ssrd_", suffix, ".rds")) / 1000000 # J/m2/h to MJ/m2/h (cumulative)
u10 <- readRDS(paste0("downloaded/u10_", suffix, ".rds"))
v10 <- readRDS(paste0("downloaded/v10_", suffix, ".rds"))
w2m <- sqrt(u10 ^ 2 + v10 ^ 2) * 0.75 # https://www.researchgate.net/post/Is-that-possible-to-convert-wind-speed-measured-in-10-m-height-to-a-possible-2-m-height-wind-speed
# get meta data 
nc_data <- nc_open(paste0("downloaded/2m_temperature_", suffix, ".nc"))
lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")
dt <- (ymd_hms("1900-01-01 00:00:00") + hours(t)) # %>% with_tz("Pacific/Auckland")
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
data <- bind_rows(data_N, data_S) %>% 
  arrange(lon, lat, dt) %>% 
  group_by(lon, lat) %>% 
  mutate(
    ssrd = ifelse(hour(dt) == 1, ssrd, ssrd - lag(ssrd))
  ) %>% 
  ungroup() %>% 
  mutate(
    dt = dt %>% with_tz("Pacific/Auckland")
  )

## max vs means (code from smaxdata) ####
weathlocal <- data %>% 
  transmute(
    STATION_ID = paste(lon, lat),
    DATE_TIME = dt,
    DATE = as.Date(dt),
    HOUR = hour(dt),
    MONTH = month(dt),
    AIR_C = t2m,
    HUM_PC = rh2m,
    SOLAR_MJH = ssrd,
    WIND_MPS = w2m,
    THI = (1.8 * AIR_C + 32) - 0.55 * (1 - HUM_PC / 100) * (1.8 * AIR_C - 26),
    GHLI = 61.78 + 4.21 * (AIR_C - 22.48) - 1.70 * (WIND_MPS - 7.05) + 5.89 * (SOLAR_MJH - 2.41),
    GHLI2 = 59.67 + 8.87 * (AIR_C - 23.5) + 2.78 * (HUM_PC - 58.34) - 3.15 * (WIND_MPS - 7.3) + 7.65 * (SOLAR_MJH - 2.14)
  ) %>% 
  # dplyr::filter(HOUR >= 11 & HOUR <= 16) %>%
  # dplyr::filter(MONTH == 12 | MONTH <= 4) %>%
  arrange(STATION_ID, DATE_TIME) 
n_distinct(weathlocal$STATION_ID)
weathlocal2 <- weathlocal %>% 
  group_by(STATION_ID, DATE) %>% 
  # dplyr::filter(n() == 24) %>% 
  summarise(
    # daily max
    AIR_MIN = safemin(AIR_C),
    AIR_MAX = safemax(AIR_C),
    AIR_HOUR = safeget(HOUR, AIR_C == safemax(AIR_C)),
    THI_MAX = safemax(THI),
    THI_HOUR = safeget(HOUR, THI == safemax(THI)), 
    GHLI_MAX = safemax(GHLI),
    GHLI_HOUR = safeget(HOUR, GHLI == safemax(GHLI)),
    GHLI2_MAX = safemax(GHLI2),
    GHLI2_HOUR = safeget(HOUR, GHLI2 == safemax(GHLI2)),
    # daily mean
    AIR_MEAN = safemean(AIR_C),
    HUM_MEAN = safemean(HUM_PC),
    SOL_MEAN = safemean(SOLAR_MJH),
    WIND_MEAN = safemean(WIND_MPS),
    # vars to model 
    HUM_MAX = safeget(HUM_PC, HOUR == AIR_HOUR),
    SOL_MAX = safeget(SOLAR_MJH, HOUR == AIR_HOUR),
    WIND_MAX = safeget(WIND_MPS, HOUR == AIR_HOUR),
    # cumulative
    THI_HRS = safesum(THI > 72),
    THI_HRS = ifelse(THI_HRS == 24, NA, THI_HRS),
    THI_ACC = safesum(pmax(0, THI - 72)),
    GHLI_HRS = safesum(GHLI > 70),
    GHLI_HRS = ifelse(GHLI_HRS == 24, NA, GHLI_HRS),
    GHLI_ACC = safesum(pmax(0, GHLI - 70)),
    GHLI2_HRS = safesum(GHLI2 > 70),
    GHLI2_HRS = ifelse(GHLI2_HRS == 24, NA, GHLI2_HRS),
    GHLI2_ACC = safesum(pmax(0, GHLI2 - 70)),
  ) %>% 
  ungroup() %>% 
  mutate(
    SOL_MEAN = SOL_MEAN * 24, # convert to MJ/m2/d
    SOL_MAX = ifelse(SOL_MAX > 0, SOL_MAX, NA), # MJ/m2/h
  ) %>% 
  drop_na(SOL_MAX)
# hist(weathlocal2$SOL_MAX) # plausible

# plot corr ####
if (FALSE){
  png("plots/weathmeanmaxcorr.png", width = 297, height = 210, units = "mm", res = 150)
  g <- ggpairs(weathlocal2, columns = which(str_detect(names(weathlocal2), "MAX|MEAN")))
  print(g)
  dev.off()
}

# phase shift ####
set.seed(123)
weathlocal2 %>% 
  dplyr::select(matches("_HOUR")) %>% 
  pivot_longer(-AIR_HOUR) %>% 
  mutate(
    name = str_extract(name, "^[^_]+"),
    name = factor(name, levels = c("THI", "GHLI", "GHLI2"))
    ) %>% 
  ggplot() +
  theme_bw() +
  labs(x = "Hour of Maximimum Air Temperature", y = "Hour of Maximum Index Value", colour = "Series") +
  geom_jitter(aes(x = AIR_HOUR, y = value, colour = "Heat Index"), alpha = 0.2) +
  geom_smooth(aes(x = AIR_HOUR, y = value, colour = "Regression"), method = "lm", se = FALSE) +
  geom_line(aes(x = AIR_HOUR, y = AIR_HOUR, colour = "1:1"), linewidth = 1) +
  # geom_abline(colour = "red") +
  # guides(colour = "none") +
  # coord_fixed() +
  scale_colour_manual(values = c("Heat Index" = "steelblue", "1:1" = "black", "Regression" = "firebrick")) +
  facet_wrap( ~ name, ncol = 3, scales = "fixed")
myggsave("plots/weathphase.png", height = 297/3)
  
# predict solar max ####
sol_mod <- lm(SOL_MAX ~ AIR_MAX + AIR_MIN + HUM_MEAN + WIND_MEAN + SOL_MEAN, data = weathlocal2)
summary(sol_mod)
sol_mod <- lm(SOL_MAX ~ 1 + SOL_MEAN, data = weathlocal2 %>% filter(SOL_MEAN>1))
summary(sol_mod)
sol_mod <- lm(SOL_MAX ~ 0 + SOL_MEAN, data = weathlocal2 %>% filter(SOL_MEAN>1))
summary(sol_mod)
weathlocal2$SOL_MAX2 <- predict(sol_mod, newdata = weathlocal2)
ggplot(weathlocal2) +
  geom_point(aes(x = SOL_MEAN, y = SOL_MAX)) +
  geom_abline(colour = "red") 
ggplot(weathlocal2) +
  labs(x = "Model for SOL_MAX") +
  geom_point(aes(x = SOL_MAX2, y = SOL_MAX)) +
  geom_abline(colour = "red") +
  coord_fixed()
myggsave("plots/weathsolmaxmodel.png")
temp1 <- weathlocal2 %>% 
  dplyr::select(MEAN = SOL_MEAN, MAX = SOL_MAX, MODEL = SOL_MAX2) %>% 
  mutate(var = "Solar Radiation (MJ/m2/h)")

# predict wind max ####
wind_mod <- lm(WIND_MAX ~ AIR_MAX + AIR_MIN + HUM_MEAN + WIND_MEAN + SOL_MEAN, data = weathlocal2)
summary(wind_mod)
wind_mod <- lm(WIND_MAX ~ 0 + WIND_MEAN, data = weathlocal2)
summary(wind_mod)
# wind_mod <- lm(WIND_MAX ~ 0 + WIND_MEAN + I(WIND_MEAN^2), data = weathlocal2)
# summary(wind_mod)
weathlocal2$WIND_MAX2 <- predict(wind_mod, newdata = weathlocal2)
ggplot(weathlocal2) +
  geom_point(aes(x = WIND_MEAN, y = WIND_MAX)) +
  geom_abline(colour = "red") 
ggplot(weathlocal2) +
  labs(x = "Model for WIND_MAX") +
  geom_point(aes(x = WIND_MAX2, y = WIND_MAX)) +
  geom_abline(colour = "red") +
  coord_fixed()
myggsave("plots/weathwindmaxmodel.png")
temp2 <- weathlocal2 %>% 
  dplyr::select(MEAN = WIND_MEAN, MAX = WIND_MAX, MODEL = WIND_MAX2) %>% 
  mutate(var = "Wind Speed (m/s)")

# predict humidity max ####
hum_mod <- lm(HUM_MAX ~ AIR_MAX + AIR_MIN + HUM_MEAN + WIND_MEAN + SOL_MEAN, data = weathlocal2)
summary(hum_mod)
hum_mod <- lm(I(100-HUM_MAX) ~ 0 + I(100-HUM_MEAN), data = weathlocal2)
summary(hum_mod)
weathlocal2$HUM_MAX2 <- 100-predict(hum_mod, newdata = weathlocal2)
ggplot(weathlocal2) +
  geom_point(aes(x = HUM_MEAN, y = HUM_MAX)) +
  geom_abline(colour = "red") 
ggplot(weathlocal2) +
  labs(x = "Model for HUM_MAX") +
  geom_point(aes(x = HUM_MAX2, y = HUM_MAX)) +
  geom_abline(colour = "red") +
  coord_fixed()
myggsave("plots/weathhummaxmodel.png")
temp3 <- weathlocal2 %>% 
  dplyr::select(MEAN = HUM_MEAN, MAX = HUM_MAX, MODEL = HUM_MAX2) %>% 
  dplyr::filter(MAX > 0, MODEL > 0) %>% 
  mutate(var = "Relative Humidity (%)")

# combined plot ####
set.seed(123)
bind_rows(temp1, temp2, temp3) %>% 
  # slice_sample(prop = 0.1) %>% 
  ggplot() +
  theme_bw() +
  labs(x = "Daily Mean Value", y = "Maximum Hourly Value", colour = "Series") +
  geom_point(aes(x = MEAN, y = MAX, colour = "Data")) +
  geom_line(aes(x = MEAN, y = MODEL, colour = "Model"), linewidth = 1) +
  # geom_abline(colour = "red") +
  # guides(colour = "none") +
  # coord_fixed() +
  scale_colour_manual(values = c("Data" = "steelblue", "Model" = "black")) +
  facet_wrap( ~ var, ncol = 3, scales = "free")
myggsave("plots/weathmodelall.png", height = 297/3)
bind_rows(temp1, temp2, temp3) %>% 
  rename(
    Data = MAX,
    Model = MODEL,
    Var = var
  ) %>% 
  group_by(Var) %>% 
  summarise(
    mean_value = mean(Data),
    mean_model = mean(Model),
    bias = mean(Model - Data),
    rmse = sqrt(mean((Model - Data)^2)),
    nse = 1 - sum((Model - Data)^2) / sum((mean(Data) - Data)^2),
    gain = as.numeric(coef(lm(Data ~ Model))[2]), # slope of regression
    rsq = summary(lm(Data ~ Model))$adj.r.squared
  ) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  print()

# validate model ####
temp <- weathlocal2 %>% 
  # dplyr::filter(AIR_MAX >= 20) %>% 
  mutate(
    THI_MOD = (1.8 * AIR_MAX + 32) - 0.55 * (1 - HUM_MAX2 / 100) * (1.8 * AIR_MAX - 26),
    GHLI_MOD = 61.78 + 4.21 * (AIR_MAX - 22.48) - 1.70 * (WIND_MAX2 - 7.05) + 5.89 * (SOL_MAX2 - 2.41),
    GHLI2_MOD = 59.67 + 8.87 * (AIR_MAX - 23.5) + 2.78 * (HUM_MAX2 - 58.34) - 3.15 * (WIND_MAX2 - 7.3) + 7.65 * (SOL_MAX2 - 2.14)
  ) %>%
  dplyr::select(THI_MAX, GHLI_MAX, GHLI2_MAX, THI_MOD, GHLI_MOD, GHLI2_MOD) %>% 
  pivot_longer(ends_with("_MAX"), names_to = "Index", values_to = "Data") %>% 
  pivot_longer(ends_with("_MOD"), names_to = "Dummy", values_to = "Model") %>% 
  mutate(
    Index = str_extract(Index, "^[^_]+"),
    Dummy = str_extract(Dummy, "^[^_]+")
  ) %>% 
  dplyr::filter(Index == Dummy) %>% 
  mutate(
    Index = factor(Index, levels = c("THI", "GHLI", "GHLI2")),
    Dummy = NULL,
    Model = unname(Model)
    )
temp %>% 
  ggplot() +
  theme_bw() +
  labs(x = "Model of Daily Maximum Value", y = "Maximum Hourly Value", colour = "Series") +
  geom_point(aes(x = Model, y = Data, colour = "Data")) +
  geom_blank(aes(y = Model, x = Data, colour = "Data")) +
  geom_smooth(aes(x = Model, y = Data, colour = "Regression"), method = "lm", se = FALSE) +
  geom_abline(aes(slope = 1, intercept = 0, colour = "1:1"), linewidth = 1, linetype = 2) +
  # geom_abline(colour = "red") +
  # guides(colour = "none") +
  # coord_fixed() +
  scale_colour_manual(values = c("Data" = "steelblue", "1:1" = "black", "Regression" = "steelblue4")) +
  facet_wrap( ~ Index, ncol = 3, scales = "free")
myggsave("plots/weathvalidmodel.png", height = 297/3)
temp %>% 
  group_by(Index) %>% 
  summarise(
    mean_value = mean(Data),
    mean_model = mean(Model),
    bias = mean(Model - Data),
    rmse = sqrt(mean((Model - Data)^2)),
    nse = 1 - sum((Model - Data)^2) / sum((mean(Data) - Data)^2),
    gain = as.numeric(coef(lm(Data ~ Model))[2]), # slope of regression
    rsq = summary(lm(Data ~ Model))$adj.r.squared
  ) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  print()

# validate hours ####
temp <- weathlocal2 %>% 
  # dplyr::filter(AIR_MAX >= 20) %>% 
  mutate(
    THI_MOD = (1.8 * AIR_MAX + 32) - 0.55 * (1 - HUM_MAX2 / 100) * (1.8 * AIR_MAX - 26),
    GHLI_MOD = 61.78 + 4.21 * (AIR_MAX - 22.48) - 1.70 * (WIND_MAX2 - 7.05) + 5.89 * (SOL_MAX2 - 2.41),
    GHLI2_MOD = 59.67 + 8.87 * (AIR_MAX - 23.5) + 2.78 * (HUM_MAX2 - 58.34) - 3.15 * (WIND_MAX2 - 7.3) + 7.65 * (SOL_MAX2 - 2.14),
    # THI_MOD = THI_MAX,
    # GHLI_MOD = GHLI_MAX,
    # GHLI2_MOD = GHLI2_MAX,
    THI_MOD = ifelse(THI_MOD > 72, THI_MOD, NA),
    GHLI_MOD = ifelse(GHLI_MOD > 70, GHLI_MOD, NA),
    GHLI2_MOD = ifelse(GHLI2_MOD > 70, GHLI2_MOD, NA),
  ) %>%
  dplyr::select(THI_HRS, GHLI_HRS, GHLI2_HRS, THI_MOD, GHLI_MOD, GHLI2_MOD) %>% 
  pivot_longer(ends_with("_HRS"), names_to = "Index", values_to = "Hours") %>% 
  pivot_longer(ends_with("_MOD"), names_to = "Dummy", values_to = "Model") %>% 
  mutate(
    Index = str_extract(Index, "^[^_]+"),
    Dummy = str_extract(Dummy, "^[^_]+")
  ) %>% 
  dplyr::filter(Index == Dummy) %>% 
  mutate(
    Index = factor(Index, levels = c("THI", "GHLI", "GHLI2")),
    Dummy = NULL,
    Model = unname(Model)
  )
temp %>% 
  ggplot() +
  theme_bw() +
  labs(x = "Model of Daily Maximum Value", y = "Hours Above Threshold", colour = "Series") +
  geom_point(aes(x = Model, y = Hours, colour = "Data")) +
  # geom_blank(aes(y = Model, x = Hours, colour = "Data")) +
  geom_smooth(aes(x = Model, y = Hours, colour = "Regression"), method = "lm", se = FALSE) +
  # geom_abline(aes(slope = 1, intercept = 0, colour = "1:1"), linewidth = 1, linetype = 2) +
  # geom_abline(colour = "red") +
  # guides(colour = "none") +
  # coord_fixed() +
  scale_y_continuous(limits = c(0,24)) +
  scale_colour_manual(values = c("Data" = "steelblue", "1:1" = "black", "Regression" = "steelblue4")) +
  facet_wrap( ~ Index, ncol = 3, scales = "free_x")
myggsave("plots/weathvalidhours.png", height = 297/3)

# validate accum ####
temp <- weathlocal2 %>% 
  # dplyr::filter(AIR_MAX >= 20) %>% 
  mutate(
    THI_MOD = (1.8 * AIR_MAX + 32) - 0.55 * (1 - HUM_MAX2 / 100) * (1.8 * AIR_MAX - 26),
    GHLI_MOD = 61.78 + 4.21 * (AIR_MAX - 22.48) - 1.70 * (WIND_MAX2 - 7.05) + 5.89 * (SOL_MAX2 - 2.41),
    GHLI2_MOD = 59.67 + 8.87 * (AIR_MAX - 23.5) + 2.78 * (HUM_MAX2 - 58.34) - 3.15 * (WIND_MAX2 - 7.3) + 7.65 * (SOL_MAX2 - 2.14),
    # THI_MOD = THI_MAX,
    # GHLI_MOD = GHLI_MAX,
    # GHLI2_MOD = GHLI2_MAX,
    THI_MOD = ifelse(THI_MOD > 72, THI_MOD, NA),
    GHLI_MOD = ifelse(GHLI_MOD > 70, GHLI_MOD, NA),
    GHLI2_MOD = ifelse(GHLI2_MOD > 70, GHLI2_MOD, NA),
  ) %>%
  dplyr::select(THI_ACC, GHLI_ACC, GHLI2_ACC, THI_MOD, GHLI_MOD, GHLI2_MOD) %>% 
  pivot_longer(ends_with("_ACC"), names_to = "Index", values_to = "Accum") %>% 
  pivot_longer(ends_with("_MOD"), names_to = "Dummy", values_to = "Model") %>% 
  mutate(
    Index = str_extract(Index, "^[^_]+"),
    Dummy = str_extract(Dummy, "^[^_]+")
  ) %>% 
  dplyr::filter(Index == Dummy) %>% 
  mutate(
    Index = factor(Index, levels = c("THI", "GHLI", "GHLI2")),
    Dummy = NULL,
    Model = unname(Model)
  )
temp %>% 
  ggplot() +
  theme_bw() +
  labs(x = "Model of Daily Maximum Value", y = "Accumulation Above Threshold", colour = "Series") +
  geom_point(aes(x = Model, y = Accum, colour = "Data")) +
  # geom_blank(aes(y = Model, x = Hours, colour = "Data")) +
  geom_smooth(aes(x = Model, y = Accum, colour = "Regression"), method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
  # geom_abline(aes(slope = 1, intercept = 0, colour = "1:1"), linewidth = 1, linetype = 2) +
  # geom_abline(colour = "red") +
  # guides(colour = "none") +
  # coord_fixed() +
  scale_colour_manual(values = c("Data" = "steelblue", "1:1" = "black", "Regression" = "steelblue4")) +
  facet_wrap( ~ Index, ncol = 3, scales = "free")
myggsave("plots/weathvalidaccum.png", height = 297/3)
