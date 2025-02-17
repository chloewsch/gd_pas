# 1. Prepare terrestrial data

# Libraries ----
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(sp)

source('functions.R')


# Load data ------
gdata <- read.csv('genetdata_terr.csv', header = TRUE)
gsf <- st_as_sf(gdata, coords = c('lon', 'lat'), crs = 4326)

allPA <- read_sf('global_pa.shp') # WDPA data
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

# Distance to nearest PA [edge] -----
# Find the nearest PA
nearest.pa <- st_nearest_feature(gsf.t2, pas)

# Distance between site and nearest PA
# Distance (m) to nearest PA only (by.element = TRUE)
padist <- st_distance(gsf.t2, pas[nearest.pa,], by_element = TRUE)

gsf.t2$PA_dist <- as.numeric(padist)

# Distance to nearest edge for all sites ----
near_pas <- pas[nearest.pa,] # Find nearest PAs
near_pas_line <- near_pas %>% 
  st_geometry() %>% 
  st_cast(to='MULTILINESTRING')

# Distance (m) to nearest PA edge
padist2 <- st_distance(gsf.t2, near_pas_line, by_element = TRUE)

# Get index of sites with PA_bin=1 to reattach to data
gsf.t2$PA_dist_neg <- as.numeric(padist2)
inpa <- which(gsf.t2$PA_bin==1)                            # identify sites inside PAs
gsf.t2$PA_dist_neg[inpa] <- gsf.t2$PA_dist_neg[inpa]*-1    # make distances negative if inside PA


# Distance to nearest PA [centroid] -----
# Distance between site and centroid of nearest PA: all sites, not only ones outside PA.
nearest_pa_all <- st_nearest_feature(gsf.t2, pas)
near_pas <- pas[nearest_pa_all,]
near_pas_cent <- st_centroid(near_pas)

padistcent <- st_distance(gsf.t2, near_pas_cent, by_element = TRUE)

gsf.t2$PA_dist_cent <- as.numeric(padistcent)

# Number of PAs within buffers ------
# Create buffers of 10, 50, 100, and 150 km
gsf.moll <- st_transform(gsf.t2, crs='ESRI:54009') # Mollweide units = m
buffs_km <- as.list(c(10, 50, 100, 150)*1000)      # buffer distances in m
sitebuff <- lapply(buffs_km, function(d) st_buffer(gsf.moll, dist = d))

site_buff.t <- lapply(sitebuff, function(x) st_transform(x, crs = st_crs(pas)))

# Count number of PAs within each buffer size
pas_in_buffer <- lapply(site_buff.t, function(x) st_intersects(x, pas))

pa_in_buff_count <- lapply(pas_in_buffer, function(x) sites_to_PAs(gdata$pop, x))

# Count PAs inside buffer
pa_count <- lapply(pa_in_buff_count, function(x) x[[1]])
names(pa_count) <- c('PA_count_10km', 'PA_count_50km', 'PA_count_100km', 'PA_count_150km')

gsf.t2$PA_count_10km <- pa_count$PA_count_10km$PA_num
gsf.t2$PA_count_50km <- pa_count$PA_count_50km$PA_num
gsf.t2$PA_count_100km <- pa_count$PA_count_100km$PA_num
gsf.t2$PA_count_150km <- pa_count$PA_count_150km$PA_num

# PA metadata ------

# IDs of PAs sites are in:
pa_index <- sip[[2]]$PA
pa_index <- na.omit(pa_index)
pa_index <- unique(pa_index)

pas_for_conn <- pas[pa_index,]

# make the index an ID to match with genetic data:
pas_for_conn$index <- pa_index

pa_merge <- merge(sip[[2]], pas_for_conn[,c('index', 'WDPA_PI', 'STATUS_', 'AREA_KM', 'IUCN_CA', 'DESIG_E')], by.x='PA', by.y = 'index', all.x = TRUE, all.y = FALSE)

## PA area ---------
# If a site is in multiple PAs, take one with bigger area (can be in multiple if nested PAs or on the border)

areacatdesig <- pa_merge %>% 
  group_by(pop) %>% 
  slice_max(AREA_KM) %>% 
  summarise(WDPA_PI = WDPA_PI,
            AREA_KM = AREA_KM) %>%
  distinct(AREA_KM, .keep_all = TRUE)

## PA IUCN category -------
# If a site is in multiple PAs, record the most restrictive IUCN category
IUCN_cat_data = as.data.frame(table(pa_merge$IUCN_CA)/sum(table(pa_merge$IUCN_CA))*100)
sum(IUCN_cat_data[c(1:5,9,10),2]) # 59% of PAs have IUCN category data (without removing populations in multiple PAs)

iucncatdat <- pa_merge %>% 
  filter(PA_num !=0) %>% 
  filter(!IUCN_CA %in% c('Not Applicable', 'Not Assigned', 'Not Reported')) %>% 
  mutate(IUCN_CA_num = case_when(IUCN_CA == 'Ia' ~ 1,
                                 IUCN_CA == 'Ib' ~ 1,
                                 IUCN_CA=='II' ~ 2,
                                 IUCN_CA=='III' ~ 3,
                                 IUCN_CA=='IV' ~ 4,
                                 IUCN_CA=='V' ~ 5,
                                 IUCN_CA=='VI' ~ 6)) %>% 
  group_by(pop) %>% 
  slice_min(IUCN_CA_num) %>% 
  summarise(IUCN_CA_num = Mode(IUCN_CA_num))

## Oldest year of PA ------
agedat <- pa_merge %>% 
  filter(PA_num !=0) %>% 
  group_by(pop) %>% 
  slice_min(STATUS_) %>% 
  rename(year = STATUS_) %>% 
  summarise(year = year) %>%
  distinct(year, .keep_all = TRUE)

## Proportion of populations in PA ---------
gsf.t2 <- gsf.t2 %>% 
  group_by(species) %>% 
  mutate(n_pops = n(),
         PA_prop = sum(PA_bin)/n_pops) %>% 
  ungroup()

# Add PA metrics to data -----
gdata_temp <- Reduce(function(x, y) merge(x, y, all.x = TRUE, all.y = FALSE, by = 'pop'), list(gsf.t2, areacatdesig, agedat, iucncatdat))

# Check:
nrow(filter(gdata_temp, PA_bin == 0, PA_dist == 0)) == 0
nrow(filter(gdata_temp, PA_bin == 1, PA_dist_neg >= 0)) == 0

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
#write.csv(gdata, 'data/genetdata_PA_terr.csv', quote = F, row.names = F)
