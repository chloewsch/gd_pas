# 2. Protected area models
# Binary (in/out)
# Distance to PA


# Libraries -----
library(tidyr)
library(dplyr)
library(sf)

#library(adespatial) # build neighborhood network
library(INLA)

library(ggplot2)
library(extrafont)

source('functions.R')

# Data -----
grdat_env <- read.csv('genetdata_PA_terr.csv', h=T)
grdat_env$PA_bin <- as.factor(grdat_env$PA_bin)

# Data prep for INLA ------
grdat_env_fst <- grdat_env[!(is.na(grdat_env$global_fst)),]
grdat_env_ne <- grdat_env[!(is.na(grdat_env$Ne)),]

# Scale & center all variables:
# Data subsets
grdat_env$scale_gd <- scale(grdat_env$gene_diversity)
grdat_env$scale_ar <- scale(grdat_env$allelic_richness)
grdat_env$scale_PA_dist <- scale(log10(grdat_env$PA_dist+1))


# FST
grdat_env_fst$scale_fst <- scale(grdat_env_fst$global_fst)
grdat_env_fst$scale_PA_dist <- scale(log10(grdat_env_fst$PA_dist+1))


## Ne
grdat_env_ne$scale_ne <- scale(log10(grdat_env_ne$Ne))
grdat_env_ne$scale_PA_dist <- scale(log10(grdat_env_ne$PA_dist+1))


## Index for spatial effect
grdat_env$site_ID = 1:nrow(grdat_env)
grdat_env_fst$site_ID = 1:nrow(grdat_env_fst)
grdat_env_ne$site_ID = 1:nrow(grdat_env_ne)

## Index for species random effect
grdat_env$sp_id <- as.numeric(as.factor(grdat_env$species)) # intercept
grdat_env$sp_id_s <- grdat_env$sp_id + max(grdat_env$sp_id) # slope

grdat_env_fst$sp_id <- as.numeric(as.factor(grdat_env_fst$species)) # intercept
grdat_env_fst$sp_id_s <- grdat_env_fst$sp_id + max(grdat_env_fst$sp_id) # slope

grdat_env_ne$sp_id <- as.numeric(as.factor(grdat_env_ne$species)) # intercept
grdat_env_ne$sp_id_s <- grdat_env_ne$sp_id + max(grdat_env_ne$sp_id) # slope

# Number of species
n.species.ga <- max(grdat_env$sp_id)
n.species.fst <- max(grdat_env_fst$sp_id)
n.species.ne <- max(grdat_env_ne$sp_id)

## priors ----
# fixed effect prior
prior.prec <- list(prior = 'normal', param = c(0, 0.1)) # mean and precision
prior <- list(prec = prior.prec)


# Run models ------
response_all <- c('scale_gd', 'scale_ar', 'scale_fst', 'scale_ne')
data_list <- list(grdat_env, grdat_env, grdat_env_fst, grdat_env_ne)
lattice_listB <- list('networks_26Feb/Lattice_knn_k8.graph', 'networks_26Feb/Lattice_knn_k7.graph','networks_26Feb/Lattice_knn_k8_FST.graph', 'networks_26Feb/Lattice_knn_k4_NE.graph')
lattice_listD <- list('networks_26Feb/Lattice_knn_k5.graph', 'networks_26Feb/Lattice_knn_k8.graph','networks_26Feb/Lattice_knn_k6_FST.graph', 'networks_26Feb/Lattice_knn_k6_NE.graph')

## Binary PA ------
# directory to store model results:
#dir.create('model_output')

binmod <- Map(function(dat, r, lati) run_INLA_model(data = dat, predictor = 'PA_bin', response = r, Lattice.adj = lati),
              dat = data_list,
              r = response_all, 
              lati = lattice_listB)
#saveRDS(binmod, 'model_output/bin_models.rds')

## Distance to PA ------
distmod <- Map(function(dat, r, lati) run_INLA_model(data = dat, predictor = 'scale_PA_dist', response = r, Lattice.adj = lati),
               dat = data_list,
               r = response_all, 
               lati = lattice_listD)
#saveRDS(distmod, 'model_output/dist_models.rds')

# Plots ------
#model_list <- list(gd_mod, ar_mod, fst_mod, ne_mod)
var_list <- list('gene diversity', 'allelic richness', 'population-specific FST', 'effective population size')

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

ss_list <- list(sshgd, sshgd, sshfst, sshne)

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
fixef <- Map(summary_fixed, model = distmod, resp_var = response_all, pred_var = 'scale_PA_dist')
fixef_df <- do.call('rbind', fixef)
rownames(fixef_df) <- NULL

fixef_df$response_var <- unlist(var_list)

## Summarize species random slopes
ranef <- Map(summary_random, model = distmod, 
             resp_var = var_list, pred_var='scale_PA_dist',
             summary_fixed = fixef, dat1=data_list)


frdat <- Map(function(r, ss) merge(r, ss, by = 'sp_id_s'),
             r = ranef, ss=ss_list)

frdat <- do.call('rbind', frdat)


## Fig 2a ------
# Merge species effect sizes for predictors
fixef_both <- rbind(fixef_dfb, fixef_df)
frdat_both <- rbind(frdatb, frdat)

# Reorder factor levels for plot
fixef_both$response_var <- factor(fixef_both$response_var, levels = var_list[c(3,2,1,4)])

frdat_both$response_var <- factor(frdat_both$response_var, levels = var_list[c(3,2,1,4)])

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
        axis.line = element_line(size = 0.5, colour = "gray90"))

### write files -----
#Fig2 <- fig2a / fig2b + plot_annotation(tag_levels = 'A')
#ggsave(filename = 'figures/figure2.pdf', plot = Fig2, device = cairo_pdf,
#       width = 11, height = 10, units = 'in')
#
#png(filename = 'figures/figure2.png', width = 11, height = 10, units = 'in', res=300)
#Fig2
#dev.off()
