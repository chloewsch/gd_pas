# Marine data ----------------
#dir.create('marine')

# 2. Compile protected areas predictors

library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(rfishbase)
library(rredlist)

source('functions.R')

# Load data ------
gdata <- read.csv('marine/genetdata_mar.csv', header = TRUE)

gdata <- gdata %>% 
          drop_na(lon, lat)

gsf <- st_as_sf(gdata, coords = c('lon', 'lat'), crs = 4326)

allPA <- read_sf('global_pa/global_pa.shp') 
pas <- filter(allPA, MARINE == 'marine')
rm(allPA)

# Put in same crs
gsf.t <- st_transform(gsf, crs = st_crs(pas))

# Spatial join genetic to PA data -----
site_in_pas <- st_intersects(gsf.t, pas)

sip <- sites_to_PAs(gdata$pop, site_in_pas)

# Binary PA variable (in/out) -----
binn <- sip[[1]]
binn$PA_bin <- ifelse(binn$PA_num > 0, 1, 0)

gr_temp <- merge(gsf.t %>% select(pop), binn %>% select(pop, PA_bin), by = 'pop', all = TRUE)

gsf.t$PA_bin <- gr_temp$PA_bin

# Distance to nearest PA -----
outpa <- filter(gsf.t, PA_bin==0) # filter for sites that aren't in PA; inside PA set distance to 0

# Find the nearest PA
nearest.pa <- st_nearest_feature(outpa, pas)

# Distance (m) between site and nearest PA
padist <- st_distance(outpa, pas[nearest.pa,], by_element = TRUE)

# To reattach to data, get index of sites with PA_bin=0
notinpa <- which(gsf.t$PA_bin==0)
gsf.t$PA_dist <- 0                                                # set all distances to 0
gsf.t$PA_dist[notinpa] <- as.numeric(padist)                      # replace 0s with distances

# PA metadata ------
# IDs of PAs sites are in:
pa_index <- sip[[2]]$PA
pa_index <- na.omit(pa_index)
pa_index <- unique(pa_index)

pas_for_conn <- pas[pa_index,]

## Attach connectivity values to data:
# make the index an ID to match with genetic data:
pas_for_conn$index <- pa_index

pa_merge <- merge(sip[[2]], pas_for_conn[,c('index', 'WDPA_PI', 'AREA_KM', 'IUCN_CA', 'DESIG_E')], by.x='PA', by.y = 'index', all.x = TRUE, all.y = FALSE)

# If a site is in multiple PAs take the maximum value (can be in multiple if nested PAs or on the border)
areacatdesig <- pa_merge %>% 
  group_by(pop) %>% 
  slice_max(AREA_KM) %>% 
  summarise(WDPA_PI = WDPA_PI,
            AREA_KM = AREA_KM,
            IUCN_mode = Mode(IUCN_CA)) %>%
  distinct(AREA_KM, IUCN_mode, .keep_all = TRUE)

# Add new columns to data ------
gdata$PA_bin <- gsf.t$PA_bin
gdata$PA_dist <- gsf.t$PA_dist

# Check:
nrow(filter(gdata, PA_bin == 0, PA_dist == 0)) == 0
nrow(filter(gdata, PA_bin == 1, PA_dist != 0)) == 0

## Proportion of populations in PA ---------
gdata <- gdata %>% 
  group_by(species) %>% 
  mutate(n_pops = n(),
         PA_prop = sum(PA_bin)/n_pops) %>% 
  ungroup()

# Add PA metrics to data -----
gdata_temp <- merge(gdata, areacatdesig, all.x = TRUE, all.y = FALSE, by = 'pop')

gdata <- gdata_temp

# Body size -------------
library(rfishbase)

fish <- unique(gdata$species)

fishbs <- fb_tbl("species") %>% 
  mutate(sci_name = paste0(Genus, '_', Species)) %>%
  filter(sci_name %in% fish) %>% 
  select(sci_name, Length, Weight)

setdiff(fish, fishbs$sci_name)

# Gobiusculus flavescens: 6cm
fishbs <- rbind(fishbs, data.frame(sci_name = 'Gobiusculus_flavescens', Length = 6, Weight = NA))

grdat_envbs <- merge(grdat_env, fishbs %>% select(sci_name, Length), by.x = 'species',
                     by.y = 'sci_name', all = TRUE)

grdat_envbs <- grdat_envbs %>% 
  rename(length_cm = Length) %>% 
  select(names(grdat_env), length_cm)


# Export -----
#write.csv(grdat_envbs, 'marine/genetdata_PA_mar.csv', quote = F, row.names = F)



