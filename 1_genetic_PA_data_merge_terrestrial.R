# 1. Prepare data

# Libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)

source('functions.R')


# Load data ------
gdata <- read.csv('clean_scripts/squeaky_clean_scripts/genetdata_terr.csv', header = TRUE)
gsf <- st_as_sf(gdata, coords = c('lon', 'lat'), crs = 4326)

allPA <- read_sf('global_pa/global_pa.shp') # WDPA data
pas <- filter(allPA, MARINE != 'marine')
rm(allPA)

# Put in same crs
gsf.t <- st_transform(gsf, crs = st_crs(pas))

# Spatial join genetic to PA data -----
site_in_pas <- st_intersects(gsf.t, pas)

sip <- sites_to_PAs(gdata$pop, site_in_pas)

# Binary PA variable (in/out) -----
binn <- sip[[1]]
binn$PA_bin <- ifelse(binn$PA_num > 0, 1, 0)

gsf.t2 <- merge(gsf.t, binn %>% select(pop, PA_bin), by = 'pop', all = TRUE)

# Distance to nearest PA -----
outpa <- filter(gsf.t2, PA_bin==0) # filter for sites that aren't in PA; inside PA set distance to 0

# Find the nearest PA
nearest.pa <- st_nearest_feature(outpa, pas)

# Distance between site and nearest PA
# Distance (m) to nearest PA only (by.element = TRUE)
padist <- st_distance(outpa, pas[nearest.pa,], by_element = TRUE)

# To reattach to data, get index of sites with PA_bin=0
notinpa <- which(gsf.t2$PA_bin==0)
gsf.t2$PA_dist <- 0                                                # set all distances to 0
gsf.t2$PA_dist[notinpa] <- as.numeric(padist)                      # replace 0s with distances


# PA metadata ------

# IDs of PAs sites are in:
pa_index <- sip[[2]]$PA
pa_index <- na.omit(pa_index)
pa_index <- unique(pa_index)

pas_for_conn <- pas[pa_index,]

# make the index an ID to match with genetic data:
pas_for_conn$index <- pa_index

pa_merge <- merge(sip[[2]], pas_for_conn[,c('index', 'WDPA_PI', 'AREA_KM', 'IUCN_CA', 'DESIG_E')], by.x='PA', by.y = 'index', all.x = TRUE, all.y = FALSE)

## PA area and IUCN designation ---------
# If a site is in multiple PAs, take one with bigger area (can be in multiple if nested PAs or on the border)

areacatdesig <- pa_merge %>% 
  group_by(pop) %>% 
  slice_max(AREA_KM) %>% 
  summarise(WDPA_PI = WDPA_PI,
            AREA_KM = AREA_KM,
            IUCN_mode = Mode(IUCN_CA)) %>%
  distinct(AREA_KM, DESIG_E_mode, IUCN_mode, .keep_all = TRUE)

## Proportion of populations in PA ---------
gsf.t2 <- gsf.t2 %>% 
  group_by(species) %>% 
  mutate(n_pops = n(),
         PA_prop = sum(PA_bin)/n_pops) %>% 
  ungroup()

# Add PA metrics to data -----
gdata_temp <- merge(gsf.t2, areacatdesig, all.x = TRUE, all.y = FALSE, by = 'pop')

# Check:
nrow(filter(gdata_temp, PA_bin == 0, PA_dist == 0)) == 0
nrow(filter(gdata_temp, PA_bin == 1, PA_dist != 0)) == 0

# attach body size ------
adult_mass_pantheria <- read.csv('body_mass_pantheria.csv', h= T)
gdata_temp2 <- merge(gdata_temp, adult_mass_pantheria, all = TRUE)

gdata_temp2 <- gdata_temp2 %>% 
              st_transform(crs = 4326)

# Convert to dataframe
gdata <- gdata_temp2 %>% st_drop_geometry()
gdata$lon <- st_coordinates(gdata_temp2)[,1]
gdata$lat <- st_coordinates(gdata_temp2)[,2]



# Export -----
write.csv(gdata, 'genetdata_PA_terr.csv', quote = F, row.names = F)


