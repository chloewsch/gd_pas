# Libraries ------
library(tidyr)
library(dplyr)

library(sf)
library(INLA)
library(adegenet)
library(spdep)

source('functions.R')

# Data -----
grdat_env <- read.csv('genetdata_PA_terr.csv', h=T)

# Filter populations inside PAs
popsinpa <- grdat_env %>% 
  filter(PA_bin == 1)

# Area dataset:
popsinpa_area <- popsinpa %>% 
  drop_na(AREA_KM) %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_gd = scale(gene_diversity),
         scale_ar = scale(allelic_richness))

popsinpa_area_fst <- popsinpa %>% 
  drop_na(AREA_KM, global_fst) %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_fst = scale(global_fst))

popsinpa_area_ne <- popsinpa %>% 
  drop_na(AREA_KM, Ne) %>% 
  mutate(scale_area = scale(log10(AREA_KM)),
         scale_ne = scale(log10(Ne)))

lapply(list(popsinpa_area, popsinpa_area_fst, popsinpa_area_ne), nrow)

# IUCN designation dataset:
popsinpa_iucn <- popsinpa %>% 
  drop_na(IUCN_mode) %>% 
  filter(!IUCN_mode %in% c('Not Applicable', 'Not Assigned',
                           'Not Reported')) %>% 
  mutate(PA_IUCN_bin = ifelse(IUCN_mode %in% c('Ia', 'Ib', 'II', 'III', 'IV'), 0, 1)) #0 = non multiple use, 1 = multiple use

popsinpa_iucn$scale_gd <- scale(popsinpa_iucn$gene_diversity)
popsinpa_iucn$scale_ar <- scale(popsinpa_iucn$allelic_richness)

popsinpa_iucn_fst <- popsinpa_iucn %>% 
                  drop_na(global_fst) %>% 
                  mutate(scale_fst = scale(global_fst))

popsinpa_iucn_ne <- popsinpa_iucn %>% 
                  drop_na(Ne) %>% 
                  mutate(scale_ne = scale(log10(Ne)))

lapply(list(popsinpa_iucn, popsinpa_iucn_fst, popsinpa_iucn_ne), nrow)

# INLA setup: Area -----------
inla.setOption(scale.model.default = F)

## create an index for sites
popsinpa_area$site_ID = 1:nrow(popsinpa_area)
popsinpa_area_fst$site_ID = 1:nrow(popsinpa_area_fst)
popsinpa_area_ne$site_ID = 1:nrow(popsinpa_area_ne)

# Index for species
popsinpa_area$sp_id <- as.numeric(as.factor(popsinpa_area$species)) # intercept
popsinpa_area$sp_id_s <- popsinpa_area$sp_id + max(popsinpa_area$sp_id) # slope

popsinpa_area_fst$sp_id <- as.numeric(as.factor(popsinpa_area_fst$species)) # intercept
popsinpa_area_fst$sp_id_s <- popsinpa_area_fst$sp_id + max(popsinpa_area_fst$sp_id) # slope

popsinpa_area_ne$sp_id <- as.numeric(as.factor(popsinpa_area_ne$species)) # intercept
popsinpa_area_ne$sp_id_s <- popsinpa_area_ne$sp_id + max(popsinpa_area_ne$sp_id) # slope

# Number of random effect levels (species)
n.species.ga <- max(popsinpa_area$sp_id)
n.species.fst <- max(popsinpa_area_fst$sp_id)
n.species.ne <- max(popsinpa_area_ne$sp_id)

## priors
prior.prec <- list(prior = 'normal', param = c(0, 0.1)) # mean and *precision*
prior <- list(prec = prior.prec)

# data lists
response_all <- c('scale_gd', 'scale_ar', 'scale_fst', 'scale_ne')
data_list <- list(popsinpa_area, popsinpa_area, popsinpa_area_fst, popsinpa_area_ne)

# Area models ------
lattice_listA <- as.list(rep('networks/Lattice_knn_k3_area_gd.graph', 2))

areamod <- Map(function(dat, r, lati, spat) run_INLA_model(data = dat, 
                                                     predictor = 'scale_area', 
                                                     response = r,
                                                     spatial = spat,
                                                     Lattice.adj = lati),
              dat = data_list,
              r = response_all, 
              spat = list(TRUE, TRUE, FALSE, FALSE),
              lati = lattice_listA)

moranisp <- Map(function(dat, resp, model, nb){
  obs = dat[,resp]
  res = obs - model$summary.fitted.values$mean
  listw_temp = nb2listw(nb, style = 'W')
  moran.test(res, listw_temp)
}, dat = data_list, resp = response_all, model = areamod, nb = nb_list) # no autocorrelation

lapply(areamod, summary)

# Save model output
#saveRDS(areamod, 'model_output/area_models.rds')


# INLA setup: IUCN -----------

## create an index for sites
popsinpa_iucn$site_ID = 1:nrow(popsinpa_iucn)
popsinpa_iucn_fst$site_ID = 1:nrow(popsinpa_iucn_fst)
popsinpa_iucn_ne$site_ID = 1:nrow(popsinpa_iucn_ne)

# Index for species
popsinpa_iucn$sp_id <- as.numeric(as.factor(popsinpa_iucn$species)) # intercept
popsinpa_iucn$sp_id_s <- popsinpa_iucn$sp_id + max(popsinpa_iucn$sp_id) # slope

popsinpa_iucn_fst$sp_id <- as.numeric(as.factor(popsinpa_iucn_fst$species)) # intercept
popsinpa_iucn_fst$sp_id_s <- popsinpa_iucn_fst$sp_id + max(popsinpa_iucn_fst$sp_id) # slope

popsinpa_iucn_ne$sp_id <- as.numeric(as.factor(popsinpa_iucn_ne$species)) # intercept
popsinpa_iucn_ne$sp_id_s <- popsinpa_iucn_ne$sp_id + max(popsinpa_iucn_ne$sp_id) # slope

# Number of random effect levels (species)
n.species.ga <- max(popsinpa_iucn$sp_id)
n.species.fst <- max(popsinpa_iucn_fst$sp_id)
n.species.ne <- max(popsinpa_iucn_ne$sp_id)

# data lists
data_listi <- list(popsinpa_iucn, popsinpa_iucn, popsinpa_iucn_fst, popsinpa_iucn_ne)

# IUCN models ------

## non spatial ------
iucnmod_nosp <- Map(function(dat, r, lati) run_INLA_model(data = dat, predictor = 'PA_IUCN_bin', response = r, spatial = FALSE),
                    dat = data_listi,
                    r = response_all)

moranii <- Map(function(dat, resp, model, nb){
  obs = dat[,resp]
  res = obs - model$summary.fitted.values$mean
  listw_temp = nb2listw(nb, style = 'W')
  moran.test(res, listw_temp)
}, dat = data_listi, resp = response_all, model = iucnmod_nosp, nb = nb_listi)
moranii # no autocorrelation

lapply(iucnmod_nosp, summary)

# Save model output
# saveRDS(iucnmod_nosp, 'model_output/pa_iucn_models.rds')


# Plots ------
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


# Summarize fixed effects
fixefa <- Map(summary_fixed, model = areamod, resp_var = response_all, pred_var = 'scale_area')
fixef_dfa <- do.call('rbind', fixefa)
rownames(fixef_dfa) <- NULL

fixef_dfa$response_var <- unlist(var_list)

fixefi <- Map(summary_fixed, model = iucnmod_nosp, resp_var = response_all, pred_var = 'PA_IUCN_bin')
fixef_dfi <- do.call('rbind', fixefi)
rownames(fixef_dfi) <- NULL

fixef_dfi$response_var <- unlist(var_list)

## Summarize species random slopes
ranefa <- Map(summary_random, model = areamod, 
              resp_var = var_list, pred_var='scale_area',
              summary_fixed = fixefa, dat1=data_list)


frdata <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefa, ss=ss_list)

frdata <- do.call('rbind', frdata)

ranefi <- Map(summary_random, model = iucnmod_nosp, 
              resp_var = var_list, pred_var='PA_IUCN_bin',
              summary_fixed = fixefi, dat1=data_list)


frdati <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefi, ss=ss_list)

frdati <- do.call('rbind', frdati)

# Reorder factor levels for plot
frdata$response_var <- factor(frdata$response_var, levels = var_list[c(3,2,1,4)])
frdati$response_var <- factor(frdati$response_var, levels = var_list[c(3,2,1,4)])


# Merge species effect sizes for predictors and write file
fixef_both <- rbind(fixef_dfa, fixef_dfi)
frdat_both <- rbind(frdata, frdati)


# Reorder factor levels for plot
fixef_both$response_var <- factor(fixef_both$response_var, levels = var_list[c(3,2,1,4)])

frdat_both$response_var <- factor(frdat_both$response_var, levels = var_list[c(3,2,1,4)])

terr_area_IUCN <- 
  ggplot() +
  # Horizontal axis lines between groups
  geom_hline(yintercept=seq(1.5, length(unique(fixef_both$response_var))-0.5, 1),
             lwd=0.5, colour="gray90") +
  # Zero line
  geom_vline(xintercept = 0, color = '#0A360C', lty = 'dashed') +
  # Random effect densities
  ggdist::stat_slab(data=frdat_both, 
                    aes(x = sp_fx, y = response_var, 
                        fill = predictor), 
                    alpha = 0.9, 
                    position = 'dodge') +
  # 90% CIs on overall effect size
  geom_linerange(data = fixef_both, aes(y = response_var, xmin = lo_90CI, xmax = up_90CI, color = predictor), 
                 lwd = 2.5, position = position_dodge(width = 1)) +
  
  # 95% CIs + overall coefficient
  geom_pointrange(data = fixef_both, aes(y = response_var, x = mean, xmin = lo_95CI,
                                         xmax = up_95CI, color = predictor),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  scale_fill_manual(values = c('#a8c693', '#8dd8a6'), 
                    guide = 'none') +
  scale_color_manual(values = c('#1a3607', '#168039'), 
                     labels = c('IUCN category', 'area'),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = c(expression(population-specific~F[ST]),
                              'allelic richness', 'gene diversity', 
                              'effective population size')) +
  labs(y= "", x = "model coefficients", title = "", color = "", fill = "") +
  theme_classic(base_size = 14) +
  theme(text=element_text(family="Roboto Medium"),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "gray90"))

### write files -----
#ggsave(filename = 'figures/SI_PAarea_models.pdf', plot = terr_area_IUCN, device = cairo_pdf,
#       width = 9, height = 7, units = 'in')
#
#png(filename = 'figures/SI_PAarea_models.png', width = 9, height = 7, units = 'in', res=300)
#terr_area_IUCN
#dev.off()




