# 3. Protected area models -- terrestrial
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
grdat_env <- read.csv('data/genetdata_PA_terr.csv', h=T)
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
lattice_listB <- list('networks/Lattice_knn_k9.graph', 'networks/Lattice_knn_k10.graph', 'networks/Lattice_knn_k9.graph', 'networks/Lattice_knn_k9_FST.graph', 'networks/Lattice_knn_k9_NE.graph')
lattice_listD_neg <- list('networks/Lattice_knn_k9.graph', 'networks/Lattice_knn_k10.graph', 'networks/Lattice_knn_k9.graph', 'networks/Lattice_knn_k9_FST.graph', 'networks/Lattice_knn_k9_NE.graph')

## Binary PA ------
# directory to store model results:
#dir.create('model_output')
mod_output <- 'model_output'

binmod <- Map(function(dat, r, lati) run_INLA_model_r1(data = dat, predictor = 'PA_bin', response = r, Lattice.adj = lati),
              dat = data_list,
              r = response_all, 
              lati = lattice_listB)
#saveRDS(binmod, paste0(mod_output, '/bin_models_besag.rds'))

## Distance to PA neg ------
distmodneg <- Map(function(dat, r, lati) run_INLA_model_r1(data = dat, predictor = 'scale_PA_dist_neg', response = r, Lattice.adj = lati),
                  dat = data_list,
                  r = response_all, 
                  lati = lattice_listD_neg)
#saveRDS(distmodneg, paste0(mod_output, '/dist_neg_models_besag.rds'))

# Plots ------
#dir.create('figures')
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


## Fig. 2 ------
# Merge species effect sizes for predictors
fixef_both0 <- rbind(fixef_dfb, fixef_df)
frdat_both0 <- rbind(frdatb, frdat)

fixef_both <- fixef_both0 %>% filter(response_var != 'allelic richness2')
frdat_both <- frdat_both0 %>% filter(response_var != 'allelic richness2')
#write.csv(frdat_both, paste0(mod_output, 'frdat_both_terr.csv'), quote = F, row.names = F)

# Reorder factor levels for plot
fixef_both$response_var <- factor(fixef_both$response_var, levels = var_list[c(4,2,1,5)])

frdat_both$response_var <- factor(frdat_both$response_var, levels = var_list[c(4,2,1,5)])

fig2a <- 
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
                     labels = c('status', 'distance'),
                     guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = c(expression(population-specific~F[ST]),
                              'allelic richness', 'gene diversity', 
                              'effective population size')) +
  labs(y= "", x = "model coefficients", title = "", color = "", fill = "") +
  theme_classic(base_size = 14) +
  theme(text=element_text(family="Roboto Medium"),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "gray90")) +
  expand_limits(y = 4.7)

### write file -----
Fig2 <- fig2a / fig2b + plot_annotation(tag_levels = 'A')

#ggsave(filename = 'figures/figure2.pdf', plot = Fig2, device = cairo_pdf,
#      width = 11, height = 10, units = 'in')


# Buffer models -----
buf10_binmod <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA10', 
                                                                  response = r, Lattice.adj = lati)}, 
                    dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf10_binmod, paste0(mod_output, '/buf10_binmod.rds'))

buf50_binmod <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA50', 
                                                                  response = r, Lattice.adj = lati)}, 
                    dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf50_binmod, paste0(mod_output, '/buf50_binmod.rds'))

buf100_binmod <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA100', 
                                                                   response = r, Lattice.adj = lati)}, 
                     dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf100_binmod, paste0(mod_output, '/buf100_binmod.rds'))

buf150_binmod <- Map(function(dat, r, lati) {run_INLA_model_buf_r1(data = dat, predictor = 'scale_log_PA150', 
                                                                   response = r, Lattice.adj = lati)}, 
                     dat = data_list, r = response_all, lati = lattice_listB)
#saveRDS(buf150_binmod, paste0(mod_output, '/buf150_binmod.rds'))

## Fig. 3 ------

# Summarize fixed effects
fixef10 <- Map(summary_fixed, model = buf10_binmod, resp_var = response_all, pred_var = 'scale_log_PA10')
fixef_df10 <- do.call('rbind', fixef10)
fixef_df10$predictor <- rep(c('PA_bin_10', 'scale_log_PA10'),5)
fixef_df10$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df10$scale <- rep(10,nrow(fixef_df10))
fixef_df10$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df10) <- NULL

fixef50 <- Map(summary_fixed, model = buf50_binmod, resp_var = response_all, pred_var = 'scale_log_PA50')
fixef_df50 <- do.call('rbind', fixef50)
fixef_df50$predictor <- rep(c('PA_bin_50', 'scale_log_PA50'),5)
fixef_df50$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df50$scale <- rep(50,nrow(fixef_df50))
fixef_df50$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df50) <- NULL

fixef100 <- Map(summary_fixed, model = buf100_binmod, resp_var = response_all, pred_var = 'scale_log_PA100')
fixef_df100 <- do.call('rbind', fixef100)
fixef_df100$predictor <- rep(c('PA_bin_100', 'scale_log_PA100'),5)
fixef_df100$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df100$scale <- rep(100,nrow(fixef_df100))
fixef_df100$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df100) <- NULL

fixef150 <- Map(summary_fixed, model = buf150_binmod, resp_var = response_all, pred_var = 'scale_log_PA150')
fixef_df150 <- do.call('rbind', fixef150)
fixef_df150$predictor <- rep(c('PA_bin_150', 'scale_log_PA150'),5)
fixef_df150$predictor_type <- rep(c('status', 'PA count'),5)
fixef_df150$scale <- rep(150,nrow(fixef_df150))
fixef_df150$response_var <- unlist(rep(var_list, each=2))
rownames(fixef_df150) <- NULL


fixef_all <- rbind(fixef_df10, fixef_df50, fixef_df100, fixef_df150)

## Summarize species random slopes
# 10km
ranef10 <- Map(summary_random_buf, model = buf10_binmod, 
               resp_var = var_list, pred_var='scale_log_PA10',
               summary_fixed = fixef10, dat1=data_list)  

frdat10 <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
               r = ranef10, ss=ss_list)

frdat10 <- do.call('rbind', frdat10)

# 50 km
ranef50 <- Map(summary_random_buf, model = buf50_binmod, 
               resp_var = var_list, pred_var='scale_log_PA50',
               summary_fixed = fixef50, dat1=data_list)  

frdat50 <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
               r = ranef50, ss=ss_list)

frdat50 <- do.call('rbind', frdat50)

# 100 km
ranef100 <- Map(summary_random_buf, model = buf100_binmod, 
                resp_var = var_list, pred_var='scale_log_PA100',
                summary_fixed = fixef100, dat1=data_list)  

frdat100 <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
                r = ranef100, ss=ss_list)

frdat100 <- do.call('rbind', frdat100)

# 150 km
ranef150 <- Map(summary_random_buf, model = buf150_binmod, 
                resp_var = var_list, pred_var='scale_log_PA150',
                summary_fixed = fixef150, dat1=data_list)  

frdat150 <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
                r = ranef150, ss=ss_list)

frdat150 <- do.call('rbind', frdat150)


# Merge species effect sizes for predictors
fixef_allbuf0 <- rbind(fixef_df10, fixef_df50, fixef_df100, fixef_df150)
frdat_allbuf0 <- rbind(frdat10, frdat50, frdat100, frdat150)

fixef_allbuf <- fixef_allbuf0 %>% filter(response_var != 'allelic richness2')
frdat_allbuf <- frdat_allbuf0 %>% filter(response_var != 'allelic richness2')

# Reorder factor levels for plot
fixef_allbuf$response_var <- factor(fixef_allbuf$response_var, levels = var_list[c(4,2,1,5)])

frdat_allbuf$response_var <- factor(frdat_allbuf$response_var, levels = var_list[c(4,2,1,5)])

# Reorder factor levels for plot
fixef_counts <- fixef_allbuf %>% filter(predictor_type=='PA count')
fixef_counts$predictor <- factor(fixef_counts$predictor, levels = c("scale_log_PA10", "scale_log_PA50", "scale_log_PA100", "scale_log_PA150"))

frdat_allbuf$predictor <- factor(frdat_allbuf$predictor, levels = c("scale_log_PA10", "scale_log_PA50", "scale_log_PA100", "scale_log_PA150"))
#write.csv(frdat_allbuf, paste0(mod_output, 'frdat_allbuf_terr.csv'), quote = F, row.names = F)

fig3a <- ggplot() +
  # Horizontal axis lines between groups
  geom_hline(yintercept=seq(1.5, length(unique(fixef_counts$predictor))-0.5, 1),
             lwd=0.5, colour="gray90") +
  # Zero line
  geom_vline(xintercept = 0, lty = 'dashed') +
  
  # Random effect densities
  ggdist::stat_slab(data=frdat_allbuf, 
                    aes(x = sp_fx, y = predictor, 
                        fill = response_var), 
                    alpha = 0.4, 
                    position = 'dodge') +
  
  # 90% CIs on overall effect size
  geom_linerange(data = fixef_counts, aes(y = predictor, xmin = lo_90CI, xmax = up_90CI, color = response_var), 
                 lwd = 2.5, position = position_dodge(width = 1)) +
  
  # 95% CIs + overall coefficient
  geom_pointrange(data = fixef_counts, aes(y = predictor, x = mean, xmin = lo_95CI,
                                           xmax = up_95CI, color = response_var),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3) +
  scale_fill_manual(values = c('#14260C', '#214013','#678C30','#99BF73'), 
                    guide = 'none') +
  scale_color_manual(values = c('#14260C', '#214013', '#678C30', '#8baf69'), 
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


Fig3 <- fig3a / fig3b + plot_annotation(tag_levels = 'A')
#ggsave(filename = 'figures/figure3.pdf', plot = Fig3, device = cairo_pdf,
#       width = 11, height = 10, units = 'in')


# Proportion of alleles inside PAs per species ------
prop_alleles_df <- grdat_env %>% 
  filter(!prop_alleles_inPA %in% c(0,1))

prop_alleles_dfsp <- grdat_env %>% 
  filter(!prop_alleles_inPA %in% c(0,1)) %>% 
  group_by(species) %>% 
  summarise(mean_inpa = mean(prop_alleles_inPA),
            mean_gd = mean(gene_diversity))

summary(prop_alleles_dfsp$mean_inpa)
