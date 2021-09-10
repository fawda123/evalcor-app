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

# apa cat point 2012 with evalcor and floor sun angle evalcor -------------

# import apacp data, clean up with qaqc
apacpwq <- import_local(station_code = 'apacpwq', path = here('data/222454.zip'))
apacpwq <- qaqc(apacpwq, qaqc_keep = c('0', '1', '2', '3', '4', '5'))
apaebmet <- import_local(station_code = 'apaebmet', path = here('data/222454.zip'))
apaebmet <- qaqc(apaebmet, qaqc_keep = c('0', '1', '2', '3', '4', '5'))

# combine wx and wq, select/rename relevant columns
apa <- comb(apacpwq, apaebmet, timestep = 60, method = 'union') %>% 
  select(
    DateTimeStamp = datetimestamp, 
    Temp = temp, 
    Sal = sal, 
    DO_obs = do_mgl,
    WSpd = wspd,
    BP = bp,
    Tide = depth,
    PAR = totpar
  )

# site metadata
locs <- SWMPr::stat_locs %>% 
  filter(station_code == 'apacp')
lat <- locs$latitude
long <- locs$longitude
tz <- attr(apa$DateTimeStamp, 'tzone')

# setup parallel backend
ncores <- detectCores()  
registerDoParallel(cores = ncores - 1)

# get evalcor results to check sun angle/tidal height correlations
corstd <- evalcor(apa, tz, lat, long, progress = T, plot = F)
apa$corstd <- corstd

# get evalcor results to check sun angle/tidal height correlations
corflr <- evalcorflr(apa, tz, lat, long, progress = T, plot = F)
apa$corflr <- corflr

# filter na tide
apa <- apa %>% 
  filter(!is.na(Tide))

save(apa, file = here('data/apa.RData'))

