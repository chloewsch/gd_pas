# 1. Prepare marine data

library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(gdistance)
library(rfishbase)
library(rredlist)

source('functions.R')

# Load data ------
gdata <- read.csv('genetdata_mar.csv', header = TRUE)

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

gsf.t2 <- merge(gsf.t %>% dplyr::select(pop), binn %>% dplyr::select(pop, PA_bin), by = 'pop', all = TRUE)

gsf.t$PA_bin <- gsf.t2$PA_bin

# Distance to nearest PA -----

## Find the nearest PA
nearest.pa <- st_nearest_feature(gsf.t2, pas)
nearest.point <- st_nearest_points(gsf.t2, pas[nearest.pa,], pairwise=TRUE)
near_pas <- pas[nearest.pa,] # Nearest PAs only to reduce computation

world <- rnaturalearth::ne_countries(scale = 'small', returnclass = 'sf')
world.t <- st_transform(world, crs = st_crs(pas))

# Identify sites with clear overland paths:
gsf.t2$index <- 1:nrow(gsf.t2)

ggplot() +
  geom_sf(data=world.t, fill='gray85', color=NA) +
  geom_sf(data=gsf.t2, color ='red') +
  geom_sf_text(data = gsf.t2, aes(label=index), size=3) +
  geom_sf(data=near_pas, color = 'blue') +
  geom_sf(data=nearest.point, color='red') +
  coord_sf(crs=4326, xlim=c(-100, -80), ylim=c(10,25))
# Points to fix: genetic sites 3, 8, 10, 12, 401, 905, 1027

# Distance to nearest edge for all sites ----
# Distance between site and nearest PA

near_pas_line <- near_pas %>% 
  st_geometry() %>% 
  st_cast(to='MULTILINESTRING') # Convert to linestring to get distances from sites inside PAs

# Distance (m) to nearest PA (by.element = TRUE)
padist2 <- st_distance(gsf.t2, near_pas_line, by_element = TRUE)

# Get index of sites with PA_bin=1 to reattach to data
gsf.t2$PA_dist <- as.numeric(padist2)
gsf.t2$PA_dist_neg <- as.numeric(padist2)
inpa <- which(gsf.t2$PA_bin==1)                            # identify sites inside PAs
gsf.t2$PA_dist_neg[inpa] <- gsf.t2$PA_dist_neg[inpa]*-1    # make distances negative if inside PA

## Marine distances -----
# Points to fix: genetic sites 3, 8, 10, 12, 401, 905, 1027
# Site 8 is closest to the same nearest feature of 10
# Site 905 is now closest to the same nearest feature of 85
# Subset df to sites that need fixing:
to_be_fixed <- which(gsf.t2$pop %in% c('Almojil1_Oman', 'Almojil1_Yemen', 'Almojil2_Oman', 'Almojil2_Yemen',
                                       'Kousteni_MYR', 'Sellas_Mexico', 'Vignaud_Djibouti'))
ptfx <- gsf.t2[to_be_fixed,]
ptfx_nearest <- near_pas[to_be_fixed,]
ptfx_nearest[c(2,6),] <- near_pas[c(10,85),]

near_pas_pt <- ptfx_nearest %>% 
  group_by(WDPAID) %>% 
  st_cast(to='MULTILINESTRING') %>% 
  sf::st_segmentize(10) %>%  
  st_cast(to='MULTIPOINT')

# Check PAs polygon to points:
ggplot() + 
  geom_sf(data=near_pas_pt) #+
#coord_sf(crs=4326, xlim=c(22, 23), ylim=c(36,37))
#coord_sf(crs=4326, xlim=c(-90, -80), ylim=c(20,25))
#coord_sf(crs=4326, xlim=c(57.5, 58.5), ylim=c(23.5,24))
#coord_sf(crs=4326, xlim=c(41, 41.2), ylim=c(18.9,19))

near_pas_pt <- near_pas_pt[-c(6:10),]

nearest.points.fix <- st_nearest_points(ptfx, near_pas_pt, pairwise=TRUE)

np4mar <- nearest.points.fix %>% 
  st_cast(to='POINT')  
np4mar <- np4mar[c(2,4,6,8,10,12,14)]

np4marsf <- st_sf(pop = ptfx$pop, geometry = st_geometry(np4mar))

genet_sites <- st_transform(ptfx, crs = st_crs(world))
near_pa_points <- st_transform(np4marsf, crs = st_crs(world))

r <- rast()
worldtvec <- vect(world)
r <- terra::rasterize(worldtvec, r, 1)

ptfv <- vect(genet_sites)
nppas <- vect(near_pa_points)

# Ensure all genetic sites and near PA points fall on water cells 
masking_points <- rbind(ptfv, nppas)

rm <- mask(r, masking_points, inverse=TRUE) # mask raster so the covered points become NA

# make all sea (NA) = -999
rm[is.na(rm)] <- -999
# turn all landmass to NA
rm[rm>-999] <- NA
# assign unit cost to all grid cells in water
rm[rm==-999] <- 1

# plot
ggplot() + 
  tidyterra::geom_spatraster(data=rm) 

## Check all points are on water cells:
extract(rm, genet_sites)
extract(rm, near_pa_points)

rm <- raster(rm) # convert to raster

# Create transition layer
tr <- transition(rm, mean, directions = 8)
tr <- geoCorrection(tr, "c")

# Find least cost path distance between sites and PAs
marine_distances <- list()
linez <- for(i in 1:nrow(genet_sites)){
  marine_distances[[i]] <- shortestPath(tr, st_coordinates(genet_sites)[i,], st_coordinates(near_pa_points)[i,], output = "SpatialLines")
}

# Compile results into dataframe
marine_distances <- lapply(marine_distances, function(x) geosphere::lengthLine(x))
marine_distances_m <- do.call('rbind.data.frame', marine_distances)
marine_distances_m <- data.frame(pop = genet_sites$pop, dist_m = marine_distances_m[,1])

# Replace distance values with new marine distances (all are outside PAs)
gsf.t2$PA_dist_neg[to_be_fixed] <- marine_distances_m$dist_m
gsf.t2$PA_dist[to_be_fixed] <- marine_distances_m$dist_m

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

## Attach connectivity values to data:
# make the index an ID to match with genetic data:
pas_for_conn$index <- pa_index

pa_merge <- merge(sip[[2]], pas_for_conn[,c('index', 'WDPA_PI', 'STATUS_', 'AREA_KM', 'IUCN_CA', 'DESIG_E')], by.x='PA', by.y = 'index', all.x = TRUE, all.y = FALSE)

# If a site is in multiple PAs, take one with bigger area (can be in multiple if nested PAs or on the border)
areacatdesig <- pa_merge %>% 
  group_by(pop) %>% 
  slice_max(AREA_KM) %>% 
  summarise(WDPA_PI = WDPA_PI,
            AREA_KM = AREA_KM) %>%
  distinct(AREA_KM, .keep_all = TRUE)

# PA IUCN category ----
# If a site is in multiple PAs, record the most restrictive IUCN category
areacatdesig <- pa_merge %>% 
  group_by(pop) %>% 
  slice_max(AREA_KM) %>% 
  summarise(WDPA_PI = WDPA_PI,
            AREA_KM = AREA_KM) %>%
  distinct(AREA_KM, .keep_all = TRUE)

IUCN_cat_data_marine = as.data.frame(table(pa_merge$IUCN_CA)/sum(table(pa_merge$IUCN_CA))*100)
sum(IUCN_cat_data_marine[c(1:4,8,9),2]) # 52% of PAs have IUCN category data (without removing populations in multiple PAs)

# Remove pops in multiple PAs:
iucncatdat <- pa_merge %>% 
  filter(PA_num !=0) %>% 
  filter(!IUCN_CA %in% c('Not Applicable', 'Not Assigned', 'Not Reported')) %>% 
  mutate(IUCN_CA_num = case_when(IUCN_CA == 'Ia' ~ 1,
                                 IUCN_CA == 'Ib' ~ 1,
                                 IUCN_CA=='II' ~ 2,
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


# Add PA metrics to data -----
gdata_temp <- Reduce(function(x, y) merge(x, y, all.x = TRUE, all.y = FALSE, by = 'pop'), list(gdata, areacatdesig, agedat, iucncatdat))


gsf.t2_4merge <- gsf.t2 %>% 
  dplyr::select(pop, PA_bin, PA_dist, PA_dist_neg, PA_count_10km, PA_count_50km, PA_count_100km, PA_count_150km) %>% 
  st_drop_geometry()

gdat <- merge(gdata_temp, gsf.t2_4merge, by = 'pop')

# Check:
nrow(filter(gdat, PA_bin == 0, PA_dist == 0)) == 0
nrow(filter(gdat, PA_bin == 1, PA_dist_neg >= 0)) == 0

## Proportion of populations in PA ---------
gdat <- gdat %>% 
  group_by(species) %>% 
  mutate(n_pops = n(),
         PA_prop = sum(PA_bin)/n_pops) %>% 
  ungroup()

# Body size -------------
library(rfishbase)

fish <- unique(gdat$species)

fishbs <- fb_tbl("species") %>% 
  mutate(sci_name = paste0(Genus, '_', Species)) %>%
  filter(sci_name %in% fish) %>% 
  dplyr::select(sci_name, Length, Weight)

setdiff(fish, fishbs$sci_name)

# Gobiusculus flavescens: 6cm
fishbs <- rbind(fishbs, data.frame(sci_name = 'Gobiusculus_flavescens', Length = 6, Weight = NA))

grdat_envbs <- merge(gdat, fishbs %>% dplyr::select(sci_name, Length), by.x = 'species',
                     by.y = 'sci_name', all = TRUE)

grdat_envbs <- grdat_envbs %>%
  rename(length_cm=Length )

# Export -----
#write.csv(grdat_envbs, 'data/genetdata_PA_mar.csv', quote = F, row.names = F)



