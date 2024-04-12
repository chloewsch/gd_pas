# Libraries ------
library(tidyr)
library(dplyr)

library(sf)
library(INLA)
library(adegenet)
library(spdep)

source('functions.R')

# Data -----
grdat_env <- read.csv('marine/genetdata_PA_mar.csv', h=T)

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

lapply(list(popsinpa_iucn, popsinpa_iucn_fst, popsinpa_iucn_fst), nrow)

# INLA setup: Area -----------
inla.setOption(scale.model.default = F)

## create an index for sites
popsinpa_area$site_ID = 1:nrow(popsinpa_area)
popsinpa_area_fst$site_ID = 1:nrow(popsinpa_area_fst)
popsinpa_area_fst_N$site_ID = 1:nrow(popsinpa_area_fst_N)
popsinpa_area_ne$site_ID = 1:nrow(popsinpa_area_ne)

# Index for species
popsinpa_area$sp_id <- as.numeric(as.factor(popsinpa_area$species)) # intercept
popsinpa_area$sp_id_s <- popsinpa_area$sp_id + max(popsinpa_area$sp_id) # slope

popsinpa_area_fst$sp_id <- as.numeric(as.factor(popsinpa_area_fst$species)) # intercept
popsinpa_area_fst$sp_id_s <- popsinpa_area_fst$sp_id + max(popsinpa_area_fst$sp_id) # slope

popsinpa_area_fst_N$sp_id <- as.numeric(as.factor(popsinpa_area_fst_N$species)) # intercept
popsinpa_area_fst_N$sp_id_s <- popsinpa_area_fst_N$sp_id + max(popsinpa_area_fst_N$sp_id) # slope

popsinpa_area_ne$sp_id <- as.numeric(as.factor(popsinpa_area_ne$species)) # intercept
popsinpa_area_ne$sp_id_s <- popsinpa_area_ne$sp_id + max(popsinpa_area_ne$sp_id) # slope

## priors
prior.prec <- list(prior = 'normal', param = c(0, 0.1)) # mean and precision
prior <- list(prec = prior.prec)

# data lists
response_all <- c('scale_gd', 'scale_ar', 'scale_fst', 'scale_ne')
data_list <- list(popsinpa_area, popsinpa_area, popsinpa_area_fst, popsinpa_area_ne)

# Area models ------
## non spatial ------

## spatial --------
# only AR is spatially autocorrelated
lattice_listA <- list(NA, 'marine/networks/Lattice_knn_k5_area_gd.graph', NA, NA)

areamod <- Map(function(dat, r, lati, spat) run_INLA_model(data = dat, 
                                                     predictor = 'scale_area', 
                                                     response = r,
                                                     spatial = spat,
                                                     Lattice.adj = lati),
              dat = data_list,
              r = response_all, 
              spat = list(FALSE, TRUE, FALSE, FALSE),
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
#saveRDS(areamod, 'marine/model_output/area_models.rds')

# FST model without mangrove killifish:

mk_fst_mod <- run_INLA_model(data = popsinpa_area_fst_N, 
               predictor = 'scale_area', 
               response = response_all[3],
               spatial = FALSE)

summary(mk_fst_mod)
#saveRDS(mk_fst_mod, 'marine/model_output/area_fst_nokilli.rds')


# Plots ------
var_list <- list('gene diversity', 'allelic richness', 'population-specific FST', 'effective population size')

# Species names and sample sizes from data:
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
fixefb <- Map(summary_fixed, model = areamod, resp_var = response_all, pred_var = 'scale_area')
fixef_dfb <- do.call('rbind', fixefb)
rownames(fixef_dfb) <- NULL

fixef_dfb$response_var <- unlist(var_list)

## Summarize species random slopes
ranefb <- Map(summary_random, model = areamod, 
              resp_var = var_list, pred_var='scale_area',
              summary_fixed = fixefb, dat1=data_list)


frdatb <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefb, ss=ss_list)

frdatb <- do.call('rbind', frdatb)

# Reorder factor levels for plot
frdatb$response_var <- factor(frdatb$response_var, levels = var_list[c(3,2,1,4)])


marine_area <- 
  ggplot() +
  # Horizontal axis lines between groups
  geom_hline(yintercept=seq(1.5, length(unique(fixef_dfb$response_var))-0.5, 1),
             lwd=0.5, colour="gray90") +
  # Zero line
  geom_vline(xintercept = 0, color = '#252099', lty = 'dashed') +
  # Random effect densities
  ggdist::stat_slab(data=frdatb, 
                    aes(x = sp_fx, y = response_var),
                    fill = '#89a5ef',
                    #color='#0F4BE8',
                    alpha = 0.9, 
                    position = 'dodge') +
  # 90% CIs on overall effect size
  geom_linerange(data = fixef_dfb, aes(y = response_var, xmin = lo_90CI, xmax = up_90CI), 
                 color='#0F4BE8', lwd = 2.5, position = position_dodge(width = 1)) +
  
  # 95% CIs + overall coefficient
  geom_pointrange(data = fixef_dfb, aes(y = response_var, x = mean, xmin = lo_95CI,
                                         xmax = up_95CI),
                  lwd = 1, position = position_dodge(width = 1), color='#0F4BE8',
                  shape = 21, fill = "white", stroke = 3) +
  scale_y_discrete(labels = c(expression(population-specific~F[ST]),
                              'allelic richness', 'gene diversity', 
                              'effective population size')) +
  labs(y= "", x = "model coefficients", title = "", color = "", fill = "") +
  theme_classic(base_size = 14) +
  theme(text=element_text(family="Roboto Medium"),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "gray90"))

## write files -----
#ggsave(filename = 'marine/figures/SI_marine_PAarea_models.pdf', plot = marine_area, device = cairo_pdf,
#       width = 9, height = 5, units = 'in')
#
#png(filename = 'marine/figures/SI_marine_PAarea_models.png', width = 9, height = 5, units = 'in', res=300)
#marine_area
#dev.off()





