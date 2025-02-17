# 4. Body size and PA attribute analyses -- marine

# Libraries ------
library(tidyr)
library(dplyr)

library(sf)
library(INLA)
library(adespatial)
library(spdep)

source('functions.R')

# Data -----
grdat_env <- read.csv('genetdata_PA_mar.csv', h=T) 
frdat_both <- read.csv('model_output/marine/frdat_both_mar.csv', h=T)
frdat_allbuf <- read.csv('model_output/marine/frdat_allbufm.csv', h=T)

var_list <- list('gene diversity', 'allelic richness', 'allelic richness2', 'population-specific FST', 'effective population size')
frdat_allbuf$response_var <- factor(frdat_allbuf$response_var, levels = var_list[c(5,1,2,4)])

# Body size -----
frdat <- frdat_both %>% filter(predictor == 'scale_PA_dist_neg')
frdatb <- frdat_both %>% filter(predictor == 'PA_bin')

binbeta_BS <- merge(frdatb %>% select(species, response_var, sp_fx, sd, num_sites),
                    grdat_env %>% select(species, length_cm) %>% distinct(), by = 'species', all = TRUE)

binbeta_BS <- binbeta_BS %>% 
  mutate(predict_var = 'binary')

distbeta_BS <- merge(frdat %>% select(species, response_var, sp_fx, sd, num_sites),
                     grdat_env %>% select(species, length_cm) %>% distinct(), by = 'species', all = TRUE)

distbeta_BS <- distbeta_BS %>% 
  mutate(predict_var = 'distance')

beta_bs <- rbind(binbeta_BS, distbeta_BS)

# Buffers
bufbeta_BS <- merge(frdat_allbuf %>% select(species, predictor, response_var, sp_fx, sd, num_sites),
                    grdat_env %>% select(species, length_cm) %>% distinct(), by = 'species', all = TRUE)
bufbeta_BS <- bufbeta_BS %>% 
  rename(predict_var = predictor)

## Plots (Fig. S3) ------
merm_trunk <- '#0e1756'

bsbin <- ggplot(beta_bs %>%  filter(predict_var == 'binary'), aes(y = sp_fx, x = log10(length_cm))) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_point(color = "#55e4e7", alpha = 1) +
  geom_smooth(color = merm_trunk, fill = merm_trunk, method = 'lm', alpha = 0.25)+
  labs(y = 'species effect size', x = expression(log[10]~adult~length), title = 'A') +
  guides(color = 'none', fill = 'none') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium")) +
  facet_wrap(~response_var)

bsdist <- ggplot(beta_bs %>%  filter(predict_var == 'distance'), aes(y = sp_fx, x = log10(length_cm), color = response_var, fill = response_var)) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_point(color = "#252099") +
  geom_smooth(color = merm_trunk, fill = merm_trunk, method = 'lm') +
  labs(y = 'species effect size', x = expression(log[10]~adult~length), title = 'B') +
  guides(color = 'none', fill = 'none') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium")) +
  facet_wrap(~response_var)

bsbuff <- 
  ggplot(bufbeta_BS %>%  filter(predict_var == 'scale_log_PA10'), aes(y = sp_fx, x = log10(length_cm), color = response_var, fill = response_var)) +
  geom_hline(yintercept = 0, lty = 'dashed') +
  geom_point(color = '#9f6eec') +
  geom_smooth(color = merm_trunk, fill = merm_trunk, method = 'lm') +
  labs(y = 'species effect size', x = expression(log[10]~adult~length), title = 'C') +
  guides(color = 'none', fill = 'none') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium")) +
  facet_wrap(~response_var)

beta_bodysize_marine <- bsbin + bsdist + bsbuff
#ggsave(filename = 'figures/FigS3.pdf', plot = beta_bodysize_marine, device = cairo_pdf,
#       width = 11, height = 5, units = 'in')
#


## Models -----
# Binary
binbeta_BS_divgd <- filter(binbeta_BS, response_var == 'gene diversity')
binbeta_BS_divar <- filter(binbeta_BS, response_var == 'allelic richness')
binbeta_BS_divne <- filter(binbeta_BS, response_var == 'effective population size')
binbeta_BS_divfst <- filter(binbeta_BS, response_var == 'population-specific FST')
lapply(list(binbeta_BS_divgd, binbeta_BS_divar, binbeta_BS_divne, binbeta_BS_divfst), nrow)

bs_mod_gd <- lm(sp_fx ~ scale(log10(binbeta_BS_divgd$length_cm)), data = binbeta_BS_divgd)
bs_mod_ar <- lm(sp_fx ~ scale(log10(binbeta_BS_divar$length_cm)), data = binbeta_BS_divar)
bs_mod_ne <- lm(sp_fx ~ scale(log10(binbeta_BS_divne$length_cm)), data = binbeta_BS_divne)
bs_mod_fst <- lm(sp_fx ~ scale(log10(binbeta_BS_divfst$length_cm)), data = binbeta_BS_divfst)

lapply(list(bs_mod_gd, bs_mod_ar, bs_mod_fst, bs_mod_ne), summary)
lapply(list(bs_mod_gd, bs_mod_ar, bs_mod_fst, bs_mod_ne), coef)
lapply(list(bs_mod_gd, bs_mod_ar, bs_mod_fst, bs_mod_ne), confint)

# Distance
distbeta_BS_divgd <-  filter(distbeta_BS, response_var == 'gene diversity')
distbeta_BS_divar <-  filter(distbeta_BS, response_var == 'allelic richness')
distbeta_BS_divne <-  filter(distbeta_BS, response_var == 'effective population size')
distbeta_BS_divfst <- filter(distbeta_BS, response_var == 'population-specific FST')


bs_modd_gd <-  lm(sp_fx ~ scale(log10(distbeta_BS_divgd$length_cm)), data =  distbeta_BS_divgd)
bs_modd_ar <-  lm(sp_fx ~ scale(log10(distbeta_BS_divar$length_cm)), data =  distbeta_BS_divar)
bs_modd_ne <-  lm(sp_fx ~ scale(log10(distbeta_BS_divne$length_cm)), data =  distbeta_BS_divne)
bs_modd_fst <- lm(sp_fx ~ scale(log10(distbeta_BS_divfst$length_cm)), data = distbeta_BS_divfst)

lapply(list(bs_modd_gd, bs_modd_ar, bs_modd_fst, bs_modd_ne), summary)
lapply(list(bs_modd_gd, bs_modd_ar, bs_modd_fst, bs_modd_ne), coef)
lapply(list(bs_modd_gd, bs_modd_ar, bs_modd_fst, bs_modd_ne), confint)

# Buffer (10)
bufbeta_BS_divgd <-  filter(bufbeta_BS, predict_var == 'scale_log_PA10', response_var == 'gene diversity')
bufbeta_BS_divar <-  filter(bufbeta_BS, predict_var == 'scale_log_PA10', response_var == 'allelic richness')
bufbeta_BS_divne <-  filter(bufbeta_BS, predict_var == 'scale_log_PA10', response_var == 'effective population size')
bufbeta_BS_divfst <- filter(bufbeta_BS, predict_var == 'scale_log_PA10', response_var == 'population-specific FST')
lapply(list(bufbeta_BS_divgd, bufbeta_BS_divar, bufbeta_BS_divne, bufbeta_BS_divfst), nrow)

bs_modbu_gd <-  lm(sp_fx ~ scale(log10(bufbeta_BS_divgd$length_cm)), data =  bufbeta_BS_divgd)
bs_modbu_ar <-  lm(sp_fx ~ scale(log10(bufbeta_BS_divar$length_cm)), data =  bufbeta_BS_divar)
bs_modbu_ne <-  lm(sp_fx ~ scale(log10(bufbeta_BS_divne$length_cm)), data =  bufbeta_BS_divne)
bs_modbu_fst <- lm(sp_fx ~ scale(log10(bufbeta_BS_divfst$length_cm)), data = bufbeta_BS_divfst)

lapply(list(bs_modbu_gd, bs_modbu_ar, bs_modbu_fst, bs_modbu_ne), summary)
lapply(list(bs_modbu_gd, bs_modbu_ar, bs_modbu_fst, bs_modbu_ne), coef)
lapply(list(bs_modbu_gd, bs_modbu_ar, bs_modbu_fst, bs_modbu_ne), confint)


# Area, IUCN cat, year -----
# Filter populations inside PAs
popsinpa <- grdat_env %>% 
  filter(PA_bin == 1)

# Area dataset:
popsinpa_area <- popsinpa %>% 
  drop_na(AREA_KM) %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_gd = scale(gene_diversity),
         scale_ar = scale(allelic_richness),
         scale_arsp = scale(allelic_richness_rbysp))

popsinpa_area_fst <- popsinpa %>% 
  drop_na(AREA_KM, global_fst) %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_fst = scale(global_fst))

# Test without species with strong effect size
popsinpa_area_fst_N <- popsinpa %>% 
  drop_na(AREA_KM, global_fst) %>% 
  filter(species != 'Kryptolebias_marmoratus') %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_fst = scale(global_fst))

popsinpa_area_ne <- popsinpa %>% 
  drop_na(AREA_KM, Ne) %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_ne = scale(log10(Ne)))

lapply(list(popsinpa_area, popsinpa_area_fst, popsinpa_area_ne), nrow)

# IUCN designation dataset:
popsinpa_iucn <- popsinpa %>% 
  drop_na(IUCN_CA_num) %>% 
  mutate(PA_IUCN_bin = ifelse(IUCN_CA_num <5, 0, 1)) #0 = non multiple use, 1 = multiple use

ggplot(popsinpa_iucn, aes(y=gene_diversity, x = as.factor(PA_IUCN_bin), color=as.factor(PA_IUCN_bin))) +
  geom_jitter(width=0.25) +
  theme_classic() +
  facet_wrap(~species)

popsinpa_iucn$scale_gd <- scale(popsinpa_iucn$gene_diversity)
popsinpa_iucn$scale_ar <- scale(popsinpa_iucn$allelic_richness)
popsinpa_iucn$scale_arsp <- scale(popsinpa_iucn$allelic_richness_rbysp)

popsinpa_iucn_fst <- popsinpa_iucn %>% 
  drop_na(global_fst) %>% 
  mutate(scale_fst = scale(global_fst))

popsinpa_iucn_ne <- popsinpa_iucn %>% 
  drop_na(Ne) %>% 
  mutate(scale_ne = scale(log10(Ne)))

lapply(list(popsinpa_iucn, popsinpa_iucn_fst, popsinpa_iucn_ne), nrow)


# Year dataset:
popsinpa_year <- popsinpa %>% 
  drop_na(year)

popsinpa_year$scale_gd <- scale(popsinpa_year$gene_diversity)
popsinpa_year$scale_ar <- scale(popsinpa_year$allelic_richness)
popsinpa_year$scale_arsp <- scale(popsinpa_year$allelic_richness_rbysp)
popsinpa_year$scale_year <- scale(popsinpa_year$year)

popsinpa_year_fst <- popsinpa_year %>% 
  drop_na(global_fst) %>% 
  mutate(scale_fst = scale(global_fst),
         scale_year = scale(year)) 

popsinpa_year_ne <- popsinpa_year %>% 
  drop_na(Ne) %>% 
  mutate(scale_ne = scale(log10(Ne)),
         scale_year=scale_year)

lapply(list(popsinpa_year, popsinpa_year_fst, popsinpa_year_ne), nrow)

ggplot(popsinpa_year) + 
  geom_histogram(aes(x=year)) +
  theme_classic()

ggplot(popsinpa_year, aes(y = gene_diversity, x = year)) + 
  geom_point() +
  theme_classic() +
  facet_wrap(~species)

# Create adjacency matrices to test for spatial autocorrelation -----

## Area----
xygd_a = coord_jitter(popsinpa_area)
xygd_f = coord_jitter(popsinpa_area_fst)
xygd_n = coord_jitter(popsinpa_area_ne)

nb2knn = chooseCN(xygd_a, type = 6, result.type = "nb", k = 3, edit.nb = F) # GD and AR
nb2knnf = chooseCN(xygd_f, type = 6, result.type = "nb", k = 3, edit.nb = F)
nb2knnn = chooseCN(xygd_n, type = 6, result.type = "nb", k = 3, edit.nb = F)

nb_list <- list(nb2knn, nb2knn, nb2knn, nb2knnf, nb2knnn)

# INLA lattice
#nb2INLA('networks/marine/Lattice_knn_k3_area_gdar1.graph', nb2knn)

# without Killifish
xygd_fnk = coord_jitter(popsinpa_area_fst_N)
nb2knnfnk = chooseCN(xygd_fnk, type = 6, result.type = "nb", k = 3, edit.nb = F)

## IUCN ------
xygd_i = coord_jitter(popsinpa_iucn)
xygd_if = coord_jitter(popsinpa_iucn_fst)
xygd_in = coord_jitter(popsinpa_iucn_ne)

nb2knni6 = chooseCN(xygd_i, type = 6, result.type = "nb", k = 6, edit.nb = F) # GD and AR
nb2knni8 = chooseCN(xygd_i, type = 6, result.type = "nb", k = 8, edit.nb = F) # GD and AR
nb2knnif = chooseCN(xygd_if, type = 6, result.type = "nb", k = 5, edit.nb = F)
nb2knnin = chooseCN(xygd_in, type = 6, result.type = "nb", k = 3, edit.nb = F)

nb_listi <- list(nb2knni, nb2knni, nb2knni, nb2knnif, nb2knnin)

#nb2INLA('clean_scripts/squeaky_clean_scripts/revision1/marine_networksr1/Lattice_knn_k6_iucn_gdar1.graph', nb2knn)
#nb2INLA('clean_scripts/squeaky_clean_scripts/revision1/marine_networksr1/Lattice_knn_k8_iucn_gdar1.graph', nb2knn)
#nb2INLA('clean_scripts/squeaky_clean_scripts/revision1/marine_networksr1/Lattice_knn_k5_iucn_fst1.graph', nb2knnif)


## Year -----
xygd_y = coord_jitter(popsinpa_year)
xygd_yf = coord_jitter(popsinpa_year_fst)
xygd_yn = coord_jitter(popsinpa_year_ne)

nb2knny4 = chooseCN(xygd_y, type = 6, result.type = "nb", k = 4, edit.nb = F) # GD and AR
nb2knny5 = chooseCN(xygd_y, type = 6, result.type = "nb", k = 5, edit.nb = F) # GD and AR
nb2knnyf = chooseCN(xygd_yf, type = 6, result.type = "nb", k = 5, edit.nb = F)
nb2knnyn = chooseCN(xygd_yn, type = 6, result.type = "nb", k = 5, edit.nb = F)

nb_listy <- list(nb2knny, nb2knny, nb2knny, nb2knnyf, nb2knnyn)

# INLA lattice
#nb2INLA("networks/marine/Lattice_knn_k4_year_gdar.graph", nb2knny4)
#nb2INLA("networks/marine/Lattice_knn_k5_year_gdar.graph", nb2knny5)
#nb2INLA("networks/marine/Lattice_knn_k4_year_fst.graph", nb2knnyf)


# INLA setup: Area -----------
inla.setOption(scale.model.default = F)

## create an index for sites
popsinpa_area$site_ID = 1:nrow(popsinpa_area)
popsinpa_area_fst$site_ID = 1:nrow(popsinpa_area_fst)
popsinpa_area_fst_N$site_ID = 1:nrow(popsinpa_area_fst_N)
popsinpa_area_ne$site_ID = 1:nrow(popsinpa_area_ne)

# Index for species
popsinpa_area$sp_id <- as.numeric(as.factor(popsinpa_area$species))     # intercept
popsinpa_area$sp_id_s <- popsinpa_area$sp_id + max(popsinpa_area$sp_id) # slope

popsinpa_area_fst$sp_id <- as.numeric(as.factor(popsinpa_area_fst$species))         # intercept
popsinpa_area_fst$sp_id_s <- popsinpa_area_fst$sp_id + max(popsinpa_area_fst$sp_id) # slope

popsinpa_area_fst_N$sp_id <- as.numeric(as.factor(popsinpa_area_fst_N$species))           # intercept
popsinpa_area_fst_N$sp_id_s <- popsinpa_area_fst_N$sp_id + max(popsinpa_area_fst_N$sp_id) # slope

popsinpa_area_ne$sp_id <- as.numeric(as.factor(popsinpa_area_ne$species))        # intercept
popsinpa_area_ne$sp_id_s <- popsinpa_area_ne$sp_id + max(popsinpa_area_ne$sp_id) # slope

## priors
pcprior <- list(prec = list(prior = "pc.prec",param = c(1, .1)))
hyper_besag <- list(prec = list(prior = "pc.prec", param = c(1, .1)))

# data lists
response_all <- c('scale_gd', 'scale_ar', 'scale_arsp', 'scale_fst', 'scale_ne')
data_list <- list(popsinpa_area, popsinpa_area, popsinpa_area, popsinpa_area_fst, popsinpa_area_ne)

# Area models ------

## spatial --------
# only AR is spatially autocorrelated
lattice_listA <- list(NA, 'networks/marine/Lattice_knn_k3_area_gdar1.graph', 
                      'networks/marine/Lattice_knn_k3_area_gdar1.graph', 
                      NA, NA)

areamod <- Map(function(dat, r, lati, spat) run_INLA_model_r1(data = dat, 
                                                              predictor = 'scale_area', 
                                                              response = r,
                                                              spatial = spat,
                                                              Lattice.adj = lati),
               dat = data_list,
               r = response_all, 
               spat = list(FALSE, TRUE, TRUE, FALSE, FALSE),
               lati = lattice_listA)

moranisp <- Map(function(dat, resp, model, nb){
  obs = dat[,resp]
  res = obs - model$summary.fitted.values$mean
  listw_temp = nb2listw(nb, style = 'W')
  moran.test(res, listw_temp)
}, dat = data_list, resp = response_all, model = areamod, nb = nb_list)
moranisp

lapply(areamod, summary)

# Save model output
#saveRDS(areamod, 'model_output/marine/area_modelsm.rds')

# FST model without mangrove killifish:

mk_fst_mod <- run_INLA_model(data = popsinpa_area_fst_N, 
                             predictor = 'scale_area', 
                             response = response_all[4],
                             spatial = FALSE)

summary(mk_fst_mod)
#saveRDS(mk_fst_mod, 'model_output/marine/area_fst_nokilli.rds')

obs = popsinpa_area_fst_N$scale_fst
res = obs - mk_fst_mod$summary.fitted.values$mean
listw_temp = nb2listw(nb2knnfnk, style = 'W')
moran.test(res, listw_temp)

# INLA setup: IUCN -----------
## create an index for sites
popsinpa_iucn$site_ID = 1:nrow(popsinpa_iucn)
popsinpa_iucn_fst$site_ID = 1:nrow(popsinpa_iucn_fst)
popsinpa_iucn_ne$site_ID = 1:nrow(popsinpa_iucn_ne)

# Index for species
popsinpa_iucn$sp_id <- as.numeric(as.factor(popsinpa_iucn$species))     # intercept
popsinpa_iucn$sp_id_s <- popsinpa_iucn$sp_id + max(popsinpa_iucn$sp_id) # slope

popsinpa_iucn_fst$sp_id <- as.numeric(as.factor(popsinpa_iucn_fst$species))         # intercept
popsinpa_iucn_fst$sp_id_s <- popsinpa_iucn_fst$sp_id + max(popsinpa_iucn_fst$sp_id) # slope

popsinpa_iucn_ne$sp_id <- as.numeric(as.factor(popsinpa_iucn_ne$species))        # intercept
popsinpa_iucn_ne$sp_id_s <- popsinpa_iucn_ne$sp_id + max(popsinpa_iucn_ne$sp_id) # slope

# Number of random effect levels (species)
n.species.ga <- max(popsinpa_iucn$sp_id)
n.species.fst <- max(popsinpa_iucn_fst$sp_id)
n.species.ne <- max(popsinpa_iucn_ne$sp_id)

# data lists
data_listi <- list(popsinpa_iucn, popsinpa_iucn, popsinpa_iucn, popsinpa_iucn_fst, popsinpa_iucn_ne)

# IUCN models ------

## non spatial ------
iucnmod_nosp <- Map(function(dat, r, lati) run_INLA_model_r1(data = dat, predictor = 'PA_IUCN_bin', response = r, spatial = FALSE),
                    dat = data_listi,
                    r = response_all)

moranii <- Map(function(dat, resp, model, nb){
  obs = dat[,resp]
  res = obs - model$summary.fitted.values$mean
  listw_temp = nb2listw(nb, style = 'W')
  moran.test(res, listw_temp)
}, dat = data_listi, resp = response_all, model = iucnmod_nosp, nb = nb_listi)
moranii # autocorr for 1,2,3,4

lapply(iucnmod_nosp, summary)

lattice_listi <- c(rep('networks/marine/Lattice_knn_k6_iucn_gdar1.graph', 2),
                   'networks/marine/Lattice_knn_k8_iucn_gdar1.graph',
                   'networks/marine/Lattice_knn_k5_iucn_fst1.graph',
                   NA)
iucnmod <- Map(function(dat, r, lati, spat) run_INLA_model_r1(data = dat, 
                                                              predictor = 'PA_IUCN_bin', 
                                                              response = r,
                                                              spatial = spat,
                                                              Lattice.adj = lati),
               dat = data_listi,
               r = response_all, 
               spat = list(TRUE, TRUE, TRUE, TRUE, FALSE),
               lati = lattice_listi) 
moranii <- Map(function(dat, resp, model, nb){
  obs = dat[,resp]
  res = obs - model$summary.fitted.values$mean
  listw_temp = nb2listw(nb, style = 'W')
  moran.test(res, listw_temp)
}, dat = data_listi, resp = response_all, model = iucnmod, nb = nb_listi)
moranii

lapply(iucnmod, summary)

# Save model output
#saveRDS(iucnmod, 'model_output/marine/pa_IUCN_modelsm.rds')

# Year models ------
## create an index for sites
popsinpa_year$site_ID = 1:nrow(popsinpa_year)
popsinpa_year_fst$site_ID = 1:nrow(popsinpa_year_fst)
popsinpa_year_ne$site_ID = 1:nrow(popsinpa_year_ne)

# Index for species
popsinpa_year$sp_id <- as.numeric(as.factor(popsinpa_year$species))     # intercept
popsinpa_year$sp_id_s <- popsinpa_year$sp_id + max(popsinpa_year$sp_id) # slope

popsinpa_year_fst$sp_id <- as.numeric(as.factor(popsinpa_year_fst$species))         # intercept
popsinpa_year_fst$sp_id_s <- popsinpa_year_fst$sp_id + max(popsinpa_year_fst$sp_id) # slope

popsinpa_year_ne$sp_id <- as.numeric(as.factor(popsinpa_year_ne$species))        # intercept
popsinpa_year_ne$sp_id_s <- popsinpa_year_ne$sp_id + max(popsinpa_year_ne$sp_id) # slope

# Number of random effect levels (species)
n.species.ga <- max(popsinpa_year$sp_id)
n.species.fst <- max(popsinpa_year_fst$sp_id)
n.species.ne <- max(popsinpa_year_ne$sp_id)

data_listy <- list(popsinpa_year, popsinpa_year, popsinpa_year, popsinpa_year_fst, popsinpa_year_ne)

## Spatial ------
lattice_listY <- as.list(c(rep('networks/marine/Lattice_knn_k5_year_gdar.graph', 2),
                           'networks/marine/Lattice_knn_k4_year_gdar.graph',
                           'networks/marine/Lattice_knn_k5_year_fst.graph',
                           NA))
yearmod <- Map(function(dat, r, lati, spat) run_INLA_model_r1(data = dat, predictor = 'scale_year', 
                                                              response = r, spatial = spat,
                                                              Lattice.adj = lati),
               dat = data_listy,
               r = response_all,
               spat= list(TRUE,TRUE,TRUE,TRUE,FALSE),
               lati = lattice_listY)

moranii <- Map(function(dat, resp, model, nb){
  obs = dat[,resp]
  res = obs - model$summary.fitted.values$mean
  listw_temp = nb2listw(nb, style = 'W')
  moran.test(res, listw_temp)
}, dat = data_listy, resp = response_all, model = yearmod, nb = nb_listy)
moranii # autocorr in 1:4 removed


lapply(yearmod, summary)


# Save model output
#saveRDS(yearmod, 'model_output/marine/pa_year_models_m.rds')


# New plot ------
var_list <- list('gene diversity', 'allelic richness', 'population-specific FST', 'effective population size')

# Species names and sample sizes from data:
# Area
sshgd <- popsinpa_area %>% 
  group_by(species, sp_id_s) %>% 
  mutate(n = n()) %>% 
  summarise(num_sites = n())

sshfst <- popsinpa_area_fst %>% 
  group_by(species, sp_id_s) %>% 
  mutate(n = n()) %>% 
  summarise(num_sites = n())

sshne <- popsinpa_area_ne %>% 
  group_by(species, sp_id_s) %>% 
  mutate(n = n()) %>% 
  summarise(num_sites = n())

ss_list <- list(sshgd, sshgd, sshfst, sshne)

models_without_arsp <- c(1,2,4,5)

# Summarize fixed effects
fixefa <- Map(summary_fixed, model = areamod[models_without_arsp], resp_var = response_all[models_without_arsp], pred_var = 'scale_area')
fixef_dfa <- do.call('rbind', fixefa)
rownames(fixef_dfa) <- NULL

fixef_dfa$response_var <- unlist(var_list)

fixefi <- Map(summary_fixed, model = iucnmod[models_without_arsp], resp_var = response_all[models_without_arsp], pred_var = 'PA_IUCN_bin')
fixef_dfi <- do.call('rbind', fixefi)
rownames(fixef_dfi) <- NULL

fixef_dfi$response_var <- unlist(var_list)

fixefy <- Map(summary_fixed, model = yearmod[models_without_arsp], resp_var = response_all[models_without_arsp], pred_var = 'scale_year')
fixef_dfy <- do.call('rbind', fixefy)
rownames(fixef_dfy) <- NULL

fixef_dfy$response_var <- unlist(var_list)

## Summarize species random slopes
ranefa <- Map(summary_random, model = areamod[models_without_arsp], 
              resp_var = var_list, pred_var='scale_area',
              summary_fixed = fixefa, dat1=data_list[models_without_arsp])


frdata <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefa, ss=ss_list)

frdata <- do.call('rbind', frdata)

ranefi <- Map(summary_random, model = iucnmod[models_without_arsp], 
              resp_var = var_list, pred_var='PA_IUCN_bin',
              summary_fixed = fixefi, dat1=data_list[models_without_arsp])


frdati <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefi, ss=ss_list)

frdati <- do.call('rbind', frdati)

ranefy <- Map(summary_random, model = yearmod[models_without_arsp], 
              resp_var = var_list, pred_var='scale_year',
              summary_fixed = fixefy, dat1=data_list[models_without_arsp])


frdaty <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefy, ss=ss_list)

frdaty <- do.call('rbind', frdaty)

# Reorder factor levels for plot
frdata$response_var <- factor(frdata$response_var, levels = var_list[c(3,2,1,4)])
frdati$response_var <- factor(frdati$response_var, levels = var_list[c(3,2,1,4)])
frdaty$response_var <- factor(frdaty$response_var, levels = var_list[c(3,2,1,4)])


# Merge species effect sizes for predictors
fixef_both <- rbind(fixef_dfa, fixef_dfi, fixef_dfy)
frdat_both <- rbind(frdata, frdati, frdaty)


# Reorder factor levels for plot
fixef_both$response_var <- factor(fixef_both$response_var, levels = var_list[c(3,2,1,4)])

frdat_both$response_var <- factor(frdat_both$response_var, levels = var_list[c(3,2,1,4)])

mar_area_IUCN <- 
  ggplot() +
  # Horizontal axis lines between groups
  geom_hline(yintercept=seq(1.5, length(unique(fixef_both$response_var))-0.5, 1),
             lwd=0.5, colour="gray90") +
  # Zero line
  geom_vline(xintercept = 0, color = '#000000', lty = 'dashed') +
  # Random effect densities
  ggdist::stat_slab(data=frdat_both, 
                    aes(x = sp_fx, y = response_var, 
                        fill = predictor), 
                    alpha = 0.4, 
                    position = 'dodge') +
  # 90% CIs on overall effect size
  geom_linerange(data = fixef_both, aes(y = response_var, xmin = lo_90CI, xmax = up_90CI, color = predictor), 
                 lwd = 2.5, position = position_dodge(width = 1)) +
  
  # 95% CIs + overall coefficient
  geom_pointrange(data = fixef_both, aes(y = response_var, x = mean, xmin = lo_95CI,
                                         xmax = up_95CI, color = predictor),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  scale_fill_manual(values = c('#023373', '#0477BF', '#04B2D9'), 
                    guide = 'none') +
  scale_color_manual(values = c('#023373', '#0477BF', '#04B2D9'), 
                     labels = c('IUCN category', 'area', 'year'),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = c(expression(population-specific~F[ST]),
                              'allelic richness', 'gene diversity', 
                              'effective population size')) +
  labs(y= "", x = "model coefficients", title = "", color = "", fill = "") +
  theme_classic(base_size = 14) +
  theme(text=element_text(family="Roboto Medium"),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "gray90"))


## write files -----
#ggsave(filename = 'figures/FigS5.pdf', plot = mar_area_IUCN, device = cairo_pdf,
#       width = 9, height = 7, units = 'in')
