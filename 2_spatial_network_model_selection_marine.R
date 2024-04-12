# 2. Select spatial network to address autocorrelation in main models (in/out PA and distance to PA)

# Libraries

library(tidyr)
library(dplyr)
library(spdep)
library(adespatial)
library(INLA)
library(sf)
library(ggplot2)

source('functions.R')

# Data -------------
grdat_env <- read.csv('marine/genetdata_PA_mar.csv', h=T)
grdat_env$PA_bin <- as.factor(grdat_env$PA_bin)

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

# Adjacency matrices with different settings -------------
# For gene diversity, allelic richness; FST, Ne done the same way separately

# Jitter points at the same location:
xygd = coord_jitter(grdat_env)

# FST
xyfst = coord_jitter(grdat_env_fst)

# Ne
xygd_ne = coord_jitter(grdat_env_ne)

## Make directory to store graphs
#dir.create('marine/networks')

##K-nearest neighbor to remove the cross-continent links
#with varying settings for the number of adjacencies per sites
for (i in 1:8) {
  #gene diversity, allelic richness
  nb2knn = chooseCN(xygd, type = 6, result.type = "nb", k = i, edit.nb = F)
  nb2INLA(paste0("marine/networks/Lattice_knn_k", i, ".graph"), nb2knn)
  
  #FST
  nb2knn_FST = chooseCN(xyfst, type = 6, result.type = "nb", k = i, edit.nb = F)
  nb2INLA(paste0("marine/networks/Lattice_knn_k", i, "_FST.graph"), nb2knn_FST)
  
  #NE
  nb2knn_NE = chooseCN(xygd_ne, type = 6, result.type = "nb", k = i, edit.nb = F)
  nb2INLA(paste0("marine/networks/Lattice_knn_k", i, "_NE.graph"), nb2knn_NE)
}



## Run models for each response variable and each adjacency matrix
# INLA -------
# Path to graph
inla.setOption(scale.model.default = F)


# Match species to create index for the phylogeny

## create an index for sites because it doesn't take factors
grdat_env$site_ID = 1:nrow(grdat_env)
grdat_env_fst$site_ID = 1:nrow(grdat_env_fst)
grdat_env_ne$site_ID = 1:nrow(grdat_env_ne)

# Index for species
grdat_env$sp_id <- as.numeric(as.factor(grdat_env$species))
grdat_env$sp_id_s <- grdat_env$sp_id + max(grdat_env$sp_id) # slope

grdat_env_fst$sp_id <- as.numeric(as.factor(grdat_env_fst$species)) # intercept
grdat_env_fst$sp_id_s <- grdat_env_fst$sp_id + max(grdat_env_fst$sp_id) # slope

grdat_env_ne$sp_id <- as.numeric(as.factor(grdat_env_ne$species)) # intercept
grdat_env_ne$sp_id_s <- grdat_env_ne$sp_id + max(grdat_env_ne$sp_id) # slope

#add another vector of species identities
grdat_env$sp_id_2 = grdat_env$sp_id
grdat_env_fst$sp_id_2 = grdat_env_fst$sp_id
grdat_env_ne$sp_id_2 = grdat_env_ne$sp_id

n.species.ga <- max(grdat_env$sp_id)
n.species.fst <- max(grdat_env_fst$sp_id)
n.species.ne <- max(grdat_env_ne$sp_id)

## priors 
prior.prec <- list(prior = 'normal', param = c(0, 0.1)) # mean and precision
prior <- list(prec = prior.prec)

# Matrices to store results
DIC_result = matrix(NA, nrow = 9, ncol = 4)
WAIC_result = matrix(NA, nrow = 9, ncol = 4)
Moran_I_result = matrix(NA, nrow = 9, ncol = 4)
Moran_p_result = matrix(NA, nrow = 9, ncol = 4)
k_result = matrix(NA, nrow = 9, ncol = 4)
sum_fix_gd <- list()
sum_fix_ar <- list()
sum_fix_fst <- list()
sum_fix_ne <- list()

# Binary models ------
## non-spatial models ------
model_non_spat = function(response, data_set, n.species, xygd, xygd_ne, xygd_fst) {
  formulae =  as.formula(paste0(response, ' ~ PA_bin +
                        f(sp_id, model="iid2d", n=2*n.species, constr = T) +
                        f(sp_id_s, PA_bin, copy="sp_id")
                        '))
  model_pa = 0
  model_pa = try(inla(formula = formulae,
                       data = data_set,
                       control.predictor = list(compute = TRUE),
                       control.compute=list(dic = TRUE, waic = TRUE,
                                            cpo = TRUE, return.marginals.predictor=TRUE),
                       control.fixed = prior,
                       quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  
  if (length(model_pa)==1) {
    formulae =  as.formula(paste0(response, ' ~ PA_bin +
                        f(sp_id, model="iid2d", n=2*n.species, constr = T) +
                        f(sp_id_s, PA_bin, copy="sp_id")'))
    
    model_pa = try(inla(formula = formulae,
                         data = data_set,
                         control.predictor = list(compute = TRUE),
                         control.compute=list(dic = TRUE, waic = TRUE,
                                              cpo = TRUE, return.marginals.predictor=TRUE),
                         control.fixed = prior,
                         quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  }
  dic_temp = model_pa$dic$dic
  waic_temp = model_pa$waic$waic
  res = data_set[,colnames(data_set)==response]-model_pa$summary.fitted.values$mean
  
  moran_k = matrix(NA, 8, 2)
  
  for (i in 1:8) {
    if (response == 'scale_ne') {
      nb_matrix = chooseCN(xygd_ne, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    } else if(response == 'scale_fst'){
      nb_matrix = chooseCN(xygd_fst, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    }
      else {
      nb_matrix = chooseCN(xygd, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    }
    listw_temp = nb2listw(nb_matrix, style = 'W')
    moran_temp = moran.test(res, listw_temp)
    moran_k[i,1] = moran_temp[["estimate"]][["Moran I statistic"]]
    moran_k[i,2] = moran_temp[["p.value"]]
  }
  
  k_temp = which(moran_k[,2] < .05)
  k_temp = k_temp[which.max(moran_k[k_temp,1])]
  moran_I_temp = moran_k[k_temp,1]
  moran_p_temp = moran_k[k_temp,2]
  
  result = list(dic = dic_temp, waic = waic_temp, k = k_temp, moran_I = moran_I_temp, moran_p = moran_p_temp)
  
  return(result)
}

mod_phy1 = model_non_spat(response = 'scale_gd', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,1] = mod_phy1$dic
WAIC_result[1,1] = mod_phy1$waic
Moran_I_result[1,1] = mod_phy1$moran_I
Moran_p_result[1,1] = mod_phy1$moran_p
k_result[1,1] = mod_phy1$k

mod_phy1 = model_non_spat(response = 'scale_ar', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,2] = mod_phy1$dic
WAIC_result[1,2] = mod_phy1$waic
Moran_I_result[1,2] = mod_phy1$moran_I
Moran_p_result[1,2] = mod_phy1$moran_p
k_result[1,2] = mod_phy1$k

mod_phy1 = model_non_spat(response = 'scale_fst', data_set = grdat_env_fst, n.species = n.species.fst, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,3] = mod_phy1$dic
WAIC_result[1,3] = mod_phy1$waic
Moran_I_result[1,3] = mod_phy1$moran_I
Moran_p_result[1,3] = mod_phy1$moran_p
k_result[1,3] = mod_phy1$k

mod_phy1 = model_non_spat(response = 'scale_ne', data_set = grdat_env_ne, n.species = n.species.ne, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,4] = mod_phy1$dic
WAIC_result[1,4] = mod_phy1$waic
Moran_I_result[1,4] = mod_phy1$moran_I
Moran_p_result[1,4] = mod_phy1$moran_p
k_result[1,4] = mod_phy1$k



## Spatial models ---------
model_spat = function(response, data_set, n.species, xygd, xygd_ne, xygd_fst, i, Lattice.adj) {
  formulae =  as.formula(paste0(response, ' ~ PA_bin +
                                f(sp_id, model="iid2d", n=2*n.species, constr = T) +
                                f(sp_id_s, PA_bin, copy="sp_id") +
                                f(site_ID, model = "bym",
                                  graph = Lattice.adj)
                                '))
  model_pa = 0
  model_pa = try(inla(formula = formulae,
                       data = data_set,
                       control.predictor = list(compute = TRUE),
                       control.compute=list(dic = TRUE, waic = TRUE,
                                            cpo = TRUE, return.marginals.predictor=TRUE),
                       control.fixed = prior,
                       quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  
  if (length(model_pa) > 1) {
    
    dic_temp = model_pa$dic$dic
    waic_temp = model_pa$waic$waic
    model_summary_temp = summary_fixed(model_pa, resp_var = response, pred_var = 'PA_bin')
    res = data_set[,colnames(data_set)==response]-model_pa$summary.fitted.values$mean
    
    if (response == 'scale_ne') {
      nb_matrix = chooseCN(xygd_ne, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    } else if(response == 'scale_fst'){
      nb_matrix = chooseCN(xygd_fst, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    } else {
      nb_matrix = chooseCN(xygd, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    }
    
    listw_temp = nb2listw(nb_matrix, style = 'W')
    moran_temp = moran.test(res, listw_temp)
    moran_I_temp = moran_temp[["estimate"]][["Moran I statistic"]]
    moran_p_temp = moran_temp[["p.value"]]
    
    result = list(dic = dic_temp, waic = waic_temp, moran_I = moran_I_temp, moran_p = moran_p_temp, sf_temp = model_summary_temp)
  } else {
    result = list(dic = NA, waic = NA, moran_I = NA, moran_p = NA, sf_temp = NA)
  }
  
  return(result)
}

#loop over differen adjacency matrices
for (i in 1:8) {
  Lattice.adj = paste0("marine/networks/Lattice_knn_k", i, ".graph")
  LatticeFST.adj = paste0("marine/networks/Lattice_knn_k", i, "_FST.graph")
  LatticeNE.adj = paste0("marine/networks/Lattice_knn_k", i, "_NE.graph")
  
  mod_phy1 = model_spat(response = 'scale_gd', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = Lattice.adj)
  DIC_result[i+1,1] = mod_phy1$dic
  WAIC_result[i+1,1] = mod_phy1$waic
  Moran_I_result[i+1,1] = mod_phy1$moran_I
  Moran_p_result[i+1,1] = mod_phy1$moran_p
  sum_fix_gd[[i]] =  mod_phy1$sf_temp
  # k_result[i+1,1] = mod_phy1$k
  
  mod_phy1 = model_spat(response = 'scale_ar', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = Lattice.adj)
  DIC_result[i+1,2] = mod_phy1$dic
  WAIC_result[i+1,2] = mod_phy1$waic
  Moran_I_result[i+1,2] = mod_phy1$moran_I
  Moran_p_result[i+1,2] = mod_phy1$moran_p
  sum_fix_ar[[i]] =  mod_phy1$sf_temp
  # k_result[i+1,2] = mod_phy1$k
  
  mod_phy1 = model_spat(response = 'scale_fst', data_set = grdat_env_fst, n.species = n.species.fst, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = LatticeFST.adj)
  DIC_result[i+1,3] = mod_phy1$dic
  WAIC_result[i+1,3] = mod_phy1$waic
  Moran_I_result[i+1,3] = mod_phy1$moran_I
  Moran_p_result[i+1,3] = mod_phy1$moran_p
  sum_fix_fst[[i]] =  mod_phy1$sf_temp
  # k_result[i+1,3] = mod_phy1$k
  
  mod_phy1 = model_spat(response = 'scale_ne', data_set = grdat_env_ne, n.species = n.species.ne, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = LatticeNE.adj)
  DIC_result[i+1,4] = mod_phy1$dic
  WAIC_result[i+1,4] = mod_phy1$waic
  Moran_I_result[i+1,4] = mod_phy1$moran_I
  Moran_p_result[i+1,4] = mod_phy1$moran_p
  sum_fix_ne[[i]] =  mod_phy1$sf_temp
  # k_result[i+1,4] = mod_phy1$k
  
}

colnames(DIC_result) = c('scale_He', 'scale_AR', 'scale_FST', 'scale_NE_log')
rownames(DIC_result) = c('non_spat', paste0('spat_k_', 1:8))
#write.csv(DIC_result, file = 'marine/DIC_connectivity_matrix_PA.csv')

colnames(WAIC_result) = c('scale_He', 'scale_AR', 'scale_FST', 'scale_NE_log')
rownames(WAIC_result) = c('non_spat', paste0('spat_k_', 1:8))
#write.csv(WAIC_result, file = 'marine/WAIC_connectivity_matrix_PA.csv')

colnames(Moran_I_result) = c('scale_He', 'scale_AR', 'scale_FST', 'scale_NE_log')
rownames(Moran_I_result) = c('non_spat', paste0('spat_k_', 1:8))
#write.csv(Moran_I_result, file = 'marine/MoranI_conn_matrix_PA.csv')


sf_gdr <- do.call('rbind', sum_fix_gd)
rownames(sf_gdr) = c(paste0('GD_spat_k_', 1:8))
sf_arr <- do.call('rbind', sum_fix_ar)
rownames(sf_arr) = c(paste0('AR_spat_k_', 1:8))
sf_fstr <- do.call('rbind', sum_fix_fst)
rownames(sf_fstr) = c(paste0('FST_spat_k_', 1:8))
sf_ner <- do.call('rbind', sum_fix_ne)
rownames(sf_ner) = c(paste0('NE_spat_k_', 1:8))

sf_res <- rbind(sf_gdr, sf_arr, sf_fstr, sf_ner)
#write.csv(sf_res, 'marine/model_coef_compare_k.csv', quote = F)

# Distance models #####
DIC_result = matrix(NA, nrow = 9, ncol = 4)
WAIC_result = matrix(NA, nrow = 9, ncol = 4)
Moran_I_result = matrix(NA, nrow = 9, ncol = 4)
Moran_p_result = matrix(NA, nrow = 9, ncol = 4)
k_result = matrix(NA, nrow = 9, ncol = 4)
sum_fix_gd <- list()
sum_fix_ar <- list()
sum_fix_fst <- list()
sum_fix_ne <- list()

## non-spatial models ------
model_non_spat = function(response, data_set, n.species, xygd, xygd_ne, xygd_fst) {
  formulae =  as.formula(paste0(response, ' ~ scale_PA_dist +
                        f(sp_id, model="iid2d", n=2*n.species, constr = T) +
                        f(sp_id_s, PA_bin, copy="sp_id")
                        '))
  model_pa = 0
  model_pa = try(inla(formula = formulae,
                       data = data_set,
                       control.predictor = list(compute = TRUE),
                       control.compute=list(dic = TRUE, waic = TRUE,
                                            cpo = TRUE, return.marginals.predictor=TRUE),
                       control.fixed = prior,
                       quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  
  if (length(model_pa)==1) {
    formulae =  as.formula(paste0(response, ' ~ scale_PA_dist +
                        f(sp_id, model="iid2d", n=2*n.species, constr = T) +
                        f(sp_id_s, scale_PA_dist, copy="sp_id")'))
    
    model_pa = try(inla(formula = formulae,
                         data = data_set,
                         control.predictor = list(compute = TRUE),
                         control.compute=list(dic = TRUE, waic = TRUE,
                                              cpo = TRUE, return.marginals.predictor=TRUE),
                         control.fixed = prior,
                         quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  }
  dic_temp = model_pa$dic$dic
  waic_temp = model_pa$waic$waic
  model_summary_temp = summary_fixed(model_pa, resp_var = response, pred_var = 'scale_PA_dist')
  res = data_set[,colnames(data_set)==response]-model_pa$summary.fitted.values$mean
  
  moran_k = matrix(NA, 8, 2)
  
  for (i in 1:8) {
    if (response == 'scale_ne') {
      nb_matrix = chooseCN(xygd_ne, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    } else if(response == 'scale_fst'){
      nb_matrix = chooseCN(xygd_fst, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    }
      else {
      nb_matrix = chooseCN(xygd, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    }
    listw_temp = nb2listw(nb_matrix, style = 'W')
    moran_temp = moran.test(res, listw_temp)
    moran_k[i,1] = moran_temp[["estimate"]][["Moran I statistic"]]
    moran_k[i,2] = moran_temp[["p.value"]]
  }
  
  k_temp = which(moran_k[,2] < .05)
  k_temp = k_temp[which.max(moran_k[k_temp,1])]
  moran_I_temp = moran_k[k_temp,1]
  moran_p_temp = moran_k[k_temp,2]
  
  result = list(dic = dic_temp, waic = waic_temp, k = k_temp, moran_I = moran_I_temp, moran_p = moran_p_temp, sf_temp = model_summary_temp)
  
  return(result)
}

mod_phy1 = model_non_spat(response = 'scale_gd', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,1] = mod_phy1$dic
WAIC_result[1,1] = mod_phy1$waic
Moran_I_result[1,1] = mod_phy1$moran_I
Moran_p_result[1,1] = mod_phy1$moran_p
sum_fix_gd[[1]] =  mod_phy1$sf_temp
k_result[1,1] = mod_phy1$k

mod_phy1 = model_non_spat(response = 'scale_ar', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,2] = mod_phy1$dic
WAIC_result[1,2] = mod_phy1$waic
Moran_I_result[1,2] = mod_phy1$moran_I
Moran_p_result[1,2] = mod_phy1$moran_p
sum_fix_ar[[1]] =  mod_phy1$sf_temp
k_result[1,2] = mod_phy1$k

mod_phy1 = model_non_spat(response = 'scale_fst', data_set = grdat_env_fst, n.species = n.species.fst, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,3] = mod_phy1$dic
WAIC_result[1,3] = mod_phy1$waic
Moran_I_result[1,3] = mod_phy1$moran_I
Moran_p_result[1,3] = mod_phy1$moran_p
sum_fix_fst[[1]] =  mod_phy1$sf_temp
k_result[1,3] = mod_phy1$k

mod_phy1 = model_non_spat(response = 'scale_ne', data_set = grdat_env_ne, n.species = n.species.ne, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst)
DIC_result[1,4] = mod_phy1$dic
WAIC_result[1,4] = mod_phy1$waic
Moran_I_result[1,4] = mod_phy1$moran_I
Moran_p_result[1,4] = mod_phy1$moran_p
sum_fix_ne[[1]] =  mod_phy1$sf_temp
k_result[1,4] = mod_phy1$k



## Spatial models ---------
model_spat = function(response, data_set, n.species, xygd, xygd_ne, xygd_fst, i, Lattice.adj) {
  formulae =  as.formula(paste0(response, ' ~ scale_PA_dist +
                                f(sp_id, model="iid2d", n=2*n.species, constr = T) +
                                f(sp_id_s, scale_PA_dist, copy="sp_id") +
                                f(site_ID, model = "bym",
                                  graph = Lattice.adj)
                                '))
  model_pa = 0
  model_pa = try(inla(formula = formulae,
                       data = data_set,
                       control.predictor = list(compute = TRUE),
                       control.compute=list(dic = TRUE, waic = TRUE,
                                            cpo = TRUE, return.marginals.predictor=TRUE),
                       control.fixed = prior,
                       quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975)))
  
  if (length(model_pa) > 1) {
    
    dic_temp = model_pa$dic$dic
    waic_temp = model_pa$waic$waic
    model_summary_temp = summary_fixed(model_pa, resp_var = response, pred_var = 'scale_PA_dist')
    res = data_set[,colnames(data_set)==response]-model_pa$summary.fitted.values$mean
    
    if (response == 'scale_ne') {
      nb_matrix = chooseCN(xygd_ne, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    } else if(response == 'scale_fst'){
      nb_matrix = chooseCN(xygd_fst, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    } else {
      nb_matrix = chooseCN(xygd, type = 6, result.type = "nb", k = i, edit.nb = F, plot.nb = F)
    }
    
    listw_temp = nb2listw(nb_matrix, style = 'W')
    moran_temp = moran.test(res, listw_temp)
    moran_I_temp = moran_temp[["estimate"]][["Moran I statistic"]]
    moran_p_temp = moran_temp[["p.value"]]
    
    result = list(dic = dic_temp, waic = waic_temp, moran_I = moran_I_temp, moran_p = moran_p_temp, sf_temp = model_summary_temp)
  } else {
    result = list(dic = NA, waic = NA, moran_I = NA, moran_p = NA, sf_temp = NA)
  }
  
  return(result)
}

#loop over differen adjacency matrices
for (i in 1:8) {
  Lattice.adj = paste0("marine/networks/Lattice_knn_k", i, ".graph")
  LatticeFST.adj = paste0("marine/networks/Lattice_knn_k", i, "_FST.graph")
  LatticeNE.adj = paste0("marine/networks/Lattice_knn_k", i, "_NE.graph")
  
  mod_phy1 = model_spat(response = 'scale_gd', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = Lattice.adj)
  DIC_result[i+1,1] = mod_phy1$dic
  WAIC_result[i+1,1] = mod_phy1$waic
  Moran_I_result[i+1,1] = mod_phy1$moran_I
  Moran_p_result[i+1,1] = mod_phy1$moran_p
  sum_fix_gd[[i+1]] =  mod_phy1$sf_temp
  
  mod_phy1 = model_spat(response = 'scale_ar', data_set = grdat_env, n.species = n.species.ga, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = Lattice.adj)
  DIC_result[i+1,2] = mod_phy1$dic
  WAIC_result[i+1,2] = mod_phy1$waic
  Moran_I_result[i+1,2] = mod_phy1$moran_I
  Moran_p_result[i+1,2] = mod_phy1$moran_p
  sum_fix_ar[[i+1]] =  mod_phy1$sf_temp
  
  mod_phy1 = model_spat(response = 'scale_fst', data_set = grdat_env_fst, n.species = n.species.fst, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = LatticeFST.adj)
  DIC_result[i+1,3] = mod_phy1$dic
  WAIC_result[i+1,3] = mod_phy1$waic
  Moran_I_result[i+1,3] = mod_phy1$moran_I
  Moran_p_result[i+1,3] = mod_phy1$moran_p
  sum_fix_fst[[i+1]] =  mod_phy1$sf_temp
  
  mod_phy1 = model_spat(response = 'scale_ne', data_set = grdat_env_ne, n.species = n.species.ne, xygd = xygd, xygd_ne = xygd_ne, xygd_fst=xyfst, i = i, Lattice.adj = LatticeNE.adj)
  DIC_result[i+1,4] = mod_phy1$dic
  WAIC_result[i+1,4] = mod_phy1$waic
  Moran_I_result[i+1,4] = mod_phy1$moran_I
  Moran_p_result[i+1,4] = mod_phy1$moran_p
  sum_fix_ne[[i+1]] =  mod_phy1$sf_temp
  
}

colnames(DIC_result) = c('scale_He', 'scale_AR', 'scale_FST', 'scale_NE_log')
rownames(DIC_result) = c('non_spat', paste0('spat_k_', 1:8))
#write.csv(DIC_result, file = 'marine/DIC_connectivity_matrix_PA_DIST.csv')

colnames(WAIC_result) = c('scale_He', 'scale_AR', 'scale_FST', 'scale_NE_log')
rownames(WAIC_result) = c('non_spat', paste0('spat_k_', 1:8))
#write.csv(WAIC_result, file = 'marine/WAIC_connectivity_matrix_PA_DIST.csv')

colnames(Moran_I_result) = c('scale_He', 'scale_AR', 'scale_FST', 'scale_NE_log')
rownames(Moran_I_result) = c('non_spat', paste0('spat_k_', 1:8))
#write.csv(Moran_I_result, file = 'marine/MoranI_conn_matrix_PA_DIST.csv')




sf_gdr <- do.call('rbind', sum_fix_gd)
rownames(sf_gdr) = c(paste0('GD_spat_k_', 0:8))
sf_arr <- do.call('rbind', sum_fix_ar)
rownames(sf_arr) = c(paste0('AR_spat_k_', 0:8))
sf_fstr <- do.call('rbind', sum_fix_fst)
rownames(sf_fstr) = c(paste0('FST_spat_k_', 0:8))
sf_ner <- do.call('rbind', sum_fix_ne)
rownames(sf_ner) = c(paste0('NE_spat_k_', 0:8))

sf_res <- rbind(sf_gdr, sf_arr, sf_fstr, sf_ner)
#write.csv(sf_res, 'marine/model_coef_compare_k_DIST.csv', quote = F)


