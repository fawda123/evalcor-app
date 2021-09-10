box::use(
  WtRegDO[ecometab, wtreg, meteval, met_day_fun, evalcor],
  here[...], 
  doParallel[registerDoParallel], 
  parallel[detectCores],
  dplyr[...],
  tidyr[...],
  SWMPr[import_local, qaqc, comb], 
  ggplot2[...], 
  tibble[enframe], 
  lubridate[...], 
  ggforce[facet_zoom],
  oce[...],
  foreach[...]
)

source(here('R/funcs.R'))

# apa observed with cordat floored sun angle ------------------------------

# original apaobs created in appalachicola.Rmd

data(apaobs)

apaobs <- apaobs %>%
  filter(year(DateTimeStamp) == 2012) %>% 
  select(-cordat)

# site metadata
locs <- SWMPr::stat_locs %>% 
  filter(station_code == 'apacp')
lat <- locs$latitude
long <- locs$longitude
tz <- attr(apaobs$DateTimeStamp, 'tzone')

# setup parallel backend
ncores <- detectCores()  
registerDoParallel(cores = ncores - 1)

# get evalcor results to check sun angle/tidal height correlations
cordatflr <- evalcorflr(apaobs, tz, lat, long, progress = T, plot = F)
apaobs$cordatflr <- cordatflr

# get regular evalcor results
cordat <- evalcor(apaobs, tz, lat, long, progress = T, plot = F)
apaobs$cordat <- cordat
apaevlcr <- apaobs

save(apaevlcr, file = here('data/apaevlcr.RData'))

