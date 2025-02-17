# 3. Protected area models -- marine
# Binary (in/out)
# Distance to PA
# Number of nearby PAs

# Libraries -----
library(tidyr)
library(dplyr)
library(sf)

library(INLA)

library(ggplot2)
library(extrafont)

source('functions.R')

# Data -----
grdat_env <- read.csv('data/genetdata_PA_mar.csv', h=T)
grdat_env$PA_bin <- as.factor(grdat_env$PA_bin)

# Data prep for INLA ------
grdat_env_fst <- grdat_env[!(is.na(grdat_env$global_fst)),]
grdat_env_ne <- grdat_env[!(is.na(grdat_env$Ne)),]

# Scale & center all variables:
# Data subsets
grdat_env$scale_gd <- scale(grdat_env$gene_diversity)
grdat_env$scale_ar <- scale(grdat_env$allelic_richness)
grdat_env$scale_arsp <- scale(grdat_env$allelic_richness_rbysp)
grdat_env$scale_PA_dist <- scale(log10(grdat_env$PA_dist+1))
grdat_env$scale_PA_dist_neg <- scale(grdat_env$PA_dist_neg)
grdat_env$scale_log_PA10 <- scale(log10(grdat_env$PA_count_10km+1))
grdat_env$scale_log_PA50 <- scale(log10(grdat_env$PA_count_50km+1))
grdat_env$scale_log_PA100 <- scale(log10(grdat_env$PA_count_100km+1))
grdat_env$scale_log_PA150 <- scale(log10(grdat_env$PA_count_150km+1))

# FST
grdat_env_fst$scale_fst <- scale(grdat_env_fst$global_fst)
grdat_env_fst$scale_PA_dist <- scale(log10(grdat_env_fst$PA_dist+1))
grdat_env_fst$scale_PA_dist_neg <- scale(grdat_env_fst$PA_dist_neg)
grdat_env_fst$scale_log_PA10 <-  scale(log10(grdat_env_fst$PA_count_10km+1))
grdat_env_fst$scale_log_PA50 <-  scale(log10(grdat_env_fst$PA_count_50km+1))
grdat_env_fst$scale_log_PA100 <- scale(log10(grdat_env_fst$PA_count_100km+1))
grdat_env_fst$scale_log_PA150 <- scale(log10(grdat_env_fst$PA_count_150km+1))

## Ne
grdat_env_ne$scale_ne <- scale(log10(grdat_env_ne$Ne))
grdat_env_ne$scale_PA_dist <- scale(log10(grdat_env_ne$PA_dist+1))
grdat_env_ne$scale_PA_dist_neg <- scale(grdat_env_ne$PA_dist_neg)
grdat_env_ne$scale_log_PA10  <- scale(log10(grdat_env_ne$PA_count_10km+1))
grdat_env_ne$scale_log_PA50  <- scale(log10(grdat_env_ne$PA_count_50km+1))
grdat_env_ne$scale_log_PA100 <- scale(log10(grdat_env_ne$PA_count_100km+1))
grdat_env_ne$scale_log_PA150 <- scale(log10(grdat_env_ne$PA_count_150km+1))

## Index for spatial effect
grdat_env$site_ID = 1:nrow(grdat_env)
grdat_env_fst$site_ID = 1:nrow(grdat_env_fst)
grdat_env_ne$site_ID = 1:nrow(grdat_env_ne)

## Index for species random effect
grdat_env$sp_id <- as.numeric(as.factor(grdat_env$species))  # intercept
grdat_env$sp_id_s <- grdat_env$sp_id + max(grdat_env$sp_id)  # slope
grdat_env$sp_id_s2 <- grdat_env$sp_id + max(grdat_env$sp_id) # slope (for buffer)

grdat_env_fst$sp_id <- as.numeric(as.factor(grdat_env_fst$species))      # intercept
grdat_env_fst$sp_id_s <- grdat_env_fst$sp_id + max(grdat_env_fst$sp_id)  # slope
grdat_env_fst$sp_id_s2 <- grdat_env_fst$sp_id + max(grdat_env_fst$sp_id) # slope (for buffer)

grdat_env_ne$sp_id <- as.numeric(as.factor(grdat_env_ne$species))     # intercept
grdat_env_ne$sp_id_s <- grdat_env_ne$sp_id + max(grdat_env_ne$sp_id)  # slope
grdat_env_ne$sp_id_s2 <- grdat_env_ne$sp_id + max(grdat_env_ne$sp_id) # slope (for buffer)

# Number of species
n.species.ga <- max(grdat_env$sp_id)
n.species.fst <- max(grdat_env_fst$sp_id)
n.species.ne <- max(grdat_env_ne$sp_id)

## priors ----
# fixed effect prior
hyper_besag <- list(prec = list(prior = "pc.prec", param = c(1, .1)))
pcprior <- list(prec = list(prior = "pc.prec",param = c(1, .1)))

# Run models ------
response_all <- c('scale_gd', 'scale_ar', 'scale_arsp', 'scale_fst', 'scale_ne')
data_list <- list(grdat_env, grdat_env, grdat_env, grdat_env_fst, grdat_env_ne)
marine_latti.dir <- 'networks/marine/'

marine_lattiB <- list('Lattice_knn_k10.graph', 'Lattice_knn_k10.graph', 'Lattice_knn_k8.graph', 'Lattice_knn_k10_FST.graph', 'Lattice_knn_k9_NE.graph')
lattice_listB <- lapply(marine_lattiB, function(f) paste0(marine_latti.dir, f))

marine_lattiD <- list('Lattice_knn_k10.graph', 'Lattice_knn_k10.graph', 'Lattice_knn_k9.graph', 'Lattice_knn_k9_FST.graph', 'Lattice_knn_k9_NE.graph')
lattice_listD <- lapply(marine_lattiD, function(f) paste0(marine_latti.dir, f))

## Binary PA ------
# directory to store model results:
#dir.create('model_output/marine')
mod_output <- 'model_output/marine'

binmod <- Map(function(dat, r, lati) run_INLA_model(data = dat, predictor = 'PA_bin', response = r, Lattice.adj = lati),
              dat = data_list,
              r = response_all, 
              lati = lattice_listB)
#saveRDS(binmod, paste0(mod_output, '/bin_models_marine_besag.rds'))


## Distance to PA neg ------
distmodneg <- Map(function(dat, r, lati) run_INLA_model(data = dat, predictor = 'scale_PA_dist_neg', response = r, Lattice.adj = lati),
                  dat = data_list,
                  r = response_all, 
                  lati = lattice_listD)
#saveRDS(distmodneg, paste0(mod_output, '/dist_models_marine_besag.rds'))

# Plots ------
var_list <- list('gene diversity', 'allelic richness', 'allelic richness2', 'population-specific FST', 'effective population size')


# Species names and sample sizes from data:
sshgd <- grdat_env %>% 
  group_by(species, sp_id_s) %>% 
  mutate(n = n()) %>% 
  summarise(num_sites = n())

sshfst <- grdat_env_fst %>% 
  group_by(species, sp_id_s) %>% 
  mutate(n = n()) %>% 
  summarise(num_sites = n())

sshne <- grdat_env_ne %>% 
  group_by(species, sp_id_s) %>% 
  mutate(n = n()) %>% 
  summarise(num_sites = n())

ss_list <- list(sshgd, sshgd, sshgd, sshfst, sshne)

## BINARY ----
# Summarize fixed effects
fixefb <- Map(summary_fixed, model = binmod, resp_var = response_all, pred_var = 'PA_bin')
fixef_dfb <- do.call('rbind', fixefb)
rownames(fixef_dfb) <- NULL

fixef_dfb$response_var <- unlist(var_list)

## Summarize species random slopes
ranefb <- Map(summary_random, model = binmod, 
              resp_var = var_list, pred_var='PA_bin',
              summary_fixed = fixefb, dat1=data_list)

frdatb <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
              r = ranefb, ss=ss_list)

frdatb <- do.call('rbind', frdatb)

## DISTANCE ----
# Summarize fixed effects
fixef <- Map(summary_fixed, model = distmodneg, resp_var = response_all, pred_var = 'scale_PA_dist_neg')
fixef_df <- do.call('rbind', fixef)
rownames(fixef_df) <- NULL

fixef_df$response_var <- unlist(var_list)

## Summarize species random slopes
ranef <- Map(summary_random, model = distmodneg, 
             resp_var = var_list, pred_var='scale_PA_dist_neg',
             summary_fixed = fixef, dat1=data_list)


frdat <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
             r = ranef, ss=ss_list)

frdat <- do.call('rbind', frdat)


# Merge species effect sizes for predictors
mfixef_both0 <- rbind(fixef_dfb, fixef_df)
mfrdat_both0 <- rbind(frdatb, frdat)

mfixef_both <- mfixef_both0 %>% filter(response_var != 'allelic richness2')
mfrdat_both <- mfrdat_both0 %>% filter(response_var != 'allelic richness2')

# Reorder factor levels for plot
mfixef_both$response_var <- factor(mfixef_both$response_var, levels = var_list[c(4,2,1,5)])
mfrdat_both$response_var <- factor(mfrdat_both$response_var, levels = var_list[c(4,2,1,5)])
#write.csv(mfrdat_both, paste0(mod_output, 'frdat_both_mar.csv'), quote = F, row.names = F)

## Fig. 2b -----
fig2b <- 
  ggplot() +
  # Horizontal axis lines between groups
  geom_hline(yintercept=seq(1.5, length(unique(mfixef_both$response_var))-0.5, 1),
             lwd=0.5, colour="gray90") +
  # Zero line
  geom_vline(xintercept = 0, color = '#252099', lty = 'dashed') +
  # Random effect densities
  ggdist::stat_slab(data=mfrdat_both, 
                    aes(x = sp_fx, y = response_var, 
                        fill = predictor), 
                    alpha = 0.9, 
                    position = 'dodge') +
  # 90% CIs on overall effect size
  geom_linerange(data = mfixef_both, aes(y = response_var, xmin = lo_90CI, xmax = up_90CI, color = predictor), 
                 lwd = 2.5, position = position_dodge(width = 1)) +
  
  # 95% CIs + overall coefficient
  geom_pointrange(data = mfixef_both, aes(y = response_var, x = mean, xmin = lo_95CI,
                                          xmax = up_95CI, color = predictor),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  scale_fill_manual(values = c('#89a5ef', '#96cdec'), 
                    guide = 'none') +
  scale_color_manual(values = c('#0F4BE8', '#0F9AE8'), 
                     labels = c('status', 'distance'),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = c(expression(population-specific~F[ST]),
                              'allelic richness', 'gene diversity', 
                              'effective population size')) +
  labs(y= "", x = "model coefficients", title = "", color = "", fill = "") +
  theme_classic(base_size = 14) +
  theme(text=element_text(family="Roboto Medium"),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "gray90"))


# Buffer models -----
buf10_binmodm <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA10', 
                                                                   response = r, Lattice.adj = lati)}, 
                     dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf10_binmodm, paste0(mod_output, '/buf10_binmodm.rds'))

buf50_binmodm <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA50', 
                                                                   response = r, Lattice.adj = lati)}, 
                     dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf50_binmodm, paste0(mod_output, '/buf50_binmodm.rds'))

buf100_binmodm <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA100', 
                                                                    response = r, Lattice.adj = lati)}, 
                      dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf100_binmodm, paste0(mod_output, '/buf100_binmodm.rds'))

buf150_binmodm <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA150', 
                                                                    response = r, Lattice.adj = lati)}, 
                      dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf150_binmodm, paste0(mod_output, '/buf150_binmodm.rds'))

## Fig. 3b ------

# Summarize fixed effects
fixef10m <- Map(summary_fixed, model = buf10_binmodm, resp_var = response_all, pred_var = 'scale_log_PA10')
fixef_df10m <- do.call('rbind', fixef10m)
fixef_df10m$predictor <- rep(c('PA_bin_10', 'scale_log_PA10'),5)
fixef_df10m$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df10m$scale <- rep(10,nrow(fixef_df10m))
fixef_df10m$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df10m) <- NULL

fixef50m <- Map(summary_fixed, model = buf50_binmodm, resp_var = response_all, pred_var = 'scale_log_PA50')
fixef_df50m <- do.call('rbind', fixef50m)
fixef_df50m$predictor <- rep(c('PA_bin_50', 'scale_log_PA50'),5)
fixef_df50m$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df50m$scale <- rep(50,nrow(fixef_df50m))
fixef_df50m$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df50m) <- NULL

fixef100m <- Map(summary_fixed, model = buf100_binmodm, resp_var = response_all, pred_var = 'scale_log_PA100')
fixef_df100m <- do.call('rbind', fixef100m)
fixef_df100m$predictor <- rep(c('PA_bin_100', 'scale_log_PA100'),5)
fixef_df100m$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df100m$scale <- rep(100,nrow(fixef_df100m))
fixef_df100m$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df100m) <- NULL

fixef150m <- Map(summary_fixed, model = buf150_binmodm, resp_var = response_all, pred_var = 'scale_log_PA150')
fixef_df150m <- do.call('rbind', fixef150m)
fixef_df150m$predictor <- rep(c('PA_bin_150', 'scale_log_PA150'),5)
fixef_df150m$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df150m$scale <- rep(150,nrow(fixef_df150m))
fixef_df150m$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df150m) <- NULL


fixef_allm <- rbind(fixef_df10m, fixef_df50m, fixef_df100m, fixef_df150m)

## Summarize species random slopes
# 10km
ranef10m <- Map(summary_random_buf, model = buf10_binmodm, 
                resp_var = var_list, pred_var='scale_log_PA10',
                summary_fixed = fixef10m, dat1=data_list)  

frdat10m <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
                r = ranef10m, ss=ss_list)

frdat10m <- do.call('rbind', frdat10m)

# 50 km
ranef50m <- Map(summary_random_buf, model = buf50_binmodm, 
                resp_var = var_list, pred_var='scale_log_PA50',
                summary_fixed = fixef50m, dat1=data_list)  

frdat50m <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
                r = ranef50m, ss=ss_list)

frdat50m <- do.call('rbind', frdat50m)

# 100 km
ranef100m <- Map(summary_random_buf, model = buf100_binmodm, 
                 resp_var = var_list, pred_var='scale_log_PA100',
                 summary_fixed = fixef100m, dat1=data_list)  

frdat100m <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
                 r = ranef100m, ss=ss_list)

frdat100m <- do.call('rbind', frdat100m)

# 150 km
ranef150m <- Map(summary_random_buf, model = buf150_binmodm, 
                 resp_var = var_list, pred_var='scale_log_PA150',
                 summary_fixed = fixef150m, dat1=data_list)  

frdat150m <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
                 r = ranef150m, ss=ss_list)

frdat150m <- do.call('rbind', frdat150m)


# Merge species effect sizes for predictors
fixef_allbuf0m <- rbind(fixef_df10m, fixef_df50m, fixef_df100m, fixef_df150m)
frdat_allbuf0m <- rbind(frdat10m, frdat50m, frdat100m, frdat150m)

fixef_allbufm <- fixef_allbuf0m %>% filter(response_var != 'allelic richness2')
frdat_allbufm <- frdat_allbuf0m %>% filter(response_var != 'allelic richness2')

# Reorder factor levels for plot
fixef_allbufm$response_var <- factor(fixef_allbufm$response_var, levels = var_list[c(4,2,1,5)])

frdat_allbufm$response_var <- factor(frdat_allbufm$response_var, levels = var_list[c(4,2,1,5)])

# Reorder factor levels for plot
fixef_countsm <- fixef_allbufm %>% filter(predictor_type=='PA count')
fixef_countsm$predictor <- factor(fixef_countsm$predictor, levels = c("scale_log_PA10", "scale_log_PA50", "scale_log_PA100", "scale_log_PA150"))

frdat_allbufm$predictor <- factor(frdat_allbufm$predictor, levels = c("scale_log_PA10", "scale_log_PA50", "scale_log_PA100", "scale_log_PA150"))
#write.csv(frdat_allbufm, paste0(mod_output, 'frdat_allbufm.csv'), quote = F, row.names = F)

fig3b <- 
  ggplot() +
  # Horizontal axis lines between groups
  geom_hline(yintercept=seq(1.5, length(unique(fixef_countsm$predictor))-0.5, 1),
             lwd=0.5, colour="gray90") +
  # Zero line
  geom_vline(xintercept = 0, lty = 'dashed') +
  
  # Random effect densities
  ggdist::stat_slab(data=frdat_allbufm, 
                    aes(x = sp_fx, y = predictor, 
                        fill = response_var), 
                    alpha = 0.4, 
                    position = 'dodge') +
  
  # 90% CIs on overall effect size
  geom_linerange(data = fixef_countsm, aes(y = predictor, xmin = lo_90CI, xmax = up_90CI, color = response_var), 
                 lwd = 2.5, position = position_dodge(width = 1)) +
  
  # 95% CIs + overall coefficient
  geom_pointrange(data = fixef_countsm, aes(y = predictor, x = mean, xmin = lo_95CI,
                                            xmax = up_95CI, color = response_var),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  scale_fill_manual(values = c('#5219E8', '#2929FF','#436cdf','#1C93FF'), 
                    guide = 'none') +
  scale_color_manual(values = c('#5219E8', '#2929FF','#436cdf','#1C93FF'), 
                     labels = c(expression(population-specific~F[ST]),
                                'allelic richness', 'gene diversity', 
                                'effective population size'),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = c('150 km', '100 km', '50 km', '10 km'), limits=rev) +
  labs(y= "", x = "model coefficients", title = "", color = "", fill = "") +
  theme_classic(base_size = 14) +
  theme(text=element_text(family="Roboto Medium"),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "gray90"))


# Proportion of alleles inside PAs per species ------
prop_alleles_df <- grdat_env %>% 
  filter(!prop_alleles_inPA %in% c(0,1))

prop_alleles_dfsp <- grdat_env %>% 
  filter(!prop_alleles_inPA %in% c(0,1)) %>% 
  group_by(species) %>% 
  summarise(mean_inpa = mean(prop_alleles_inPA),
            mean_gd = mean(gene_diversity))

summary(prop_alleles_dfsp$mean_inpa)
