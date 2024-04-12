# Summarise PA info --------
sites_to_PAs <- function(pop_id, intersect_list){
  
  # For each list item in intersect_list (right now just indexes of PAs overlapping genetic sites),
  # make a dataframe with columns: pop = ID of genetic site; PA = index of PA; PA_num = number of PAs overlapping for that pop (1 number per pop)
  PAlist <- Map(function(pop_id, pa_id){
    data.frame(pop = rep(pop_id, length(pa_id)),
               PA = pa_id,
               PA_num = rep(length(pa_id), length(pa_id)))
  }, pop_id = as.list(pop_id), pa_id = intersect_list)
  # Unlist the list to get all pops in 1 dataframe
  PAlist <- do.call('rbind', PAlist)
  
  # Sites that don't intersect were removed in previous steps, add them back
  popdf <- data.frame(pop = pop_id)
  fulldat <- merge(popdf, PAlist, by = 'pop', 
                   all = TRUE, incomparables = NA)
  fulldat$PA_num[is.na(fulldat$PA_num)] <- 0 # assign 0 to sites with no overlapping PAs
  
  # Create separate dataframe with 1 row per population & number of overlapping PAs (no PA IDs)
  sparsedat <- distinct(fulldat, pop, PA_num)
  
  # Return a list with (1) the number of PAs overlapping each site (sparse), and (2) this info + the index of PAs overlapping sites
  return(list(sparsedat, fulldat))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


coord_jitter <- function(dat){
  if('sf' %in% class(dat)){
    sfsites = dat
  } else(sfsites = st_as_sf(dat, coords = c("lon", "lat"), crs = 4326))
  
  # Jitter points at the same location:
  dupid = which(duplicated(st_geometry(sfsites)))
  dupes = sfsites[dupid,]
  dupesj = st_jitter(dupes, 0.1)
  jit = sfsites
  st_geometry(jit)[dupid] = st_geometry(dupesj)
  xygd = st_coordinates(jit)
  return(xygd)
}


run_INLA_model <- function(data, predictor, response, spatial=TRUE, Lattice.adj = NULL){
  
  data_temp <- data
  colnames(data_temp)[which(colnames(data_temp)==response)] = 'response'
  colnames(data_temp)[which(colnames(data_temp)==predictor)] = 'predictor'
  n.species <- length(unique(data_temp$species))
  
  if(spatial==FALSE){  
    formulae = response ~ 1 + predictor +                     # fixed effect
      f(sp_id, model="iid2d", n=2*n.species, constr = T) +    # random intercept (correlated with slope)
      f(sp_id_s, predictor, copy="sp_id") #+                  # random slope (correlated w/intercept)
    
    
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  control.fixed = prior,
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
  } else{
    formulae = response ~ 1 + predictor +                      # fixed effect
      f(sp_id, model="iid2d", n=2*n.species, constr = T) +     # random intercept (correlated with slope)
      f(sp_id_s, predictor, copy="sp_id") +                    # random slope (correlated w/intercept)
      f(site_ID, model = "bym",                                # spatial random effect
        graph = Lattice.adj) #+
   
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  control.fixed = prior,
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
    return(gd_mod)
  }
}

# marine PA IUCN category model (no random slope) 
mi_INLA_model <- function(data, predictor, response, spatial=TRUE, Lattice.adj = NULL){
  
  data_temp <- data
  colnames(data_temp)[which(colnames(data_temp)==response)] = 'response'
  colnames(data_temp)[which(colnames(data_temp)==predictor)] = 'predictor'
  n.species <- length(unique(data_temp$species))
  
  if(spatial==FALSE){  
    formulae = response ~ 1 + predictor +                     # fixed effect
      f(sp_id, model="iid", n=n.species, constr = T) 
    
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  control.fixed = prior,
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
  } else{
    formulae = response ~ 1 + predictor +                      # fixed effect
        f(site_ID, model = "besag",                                # spatial random effect
        graph = Lattice.adj) 

    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  control.fixed = prior,
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
    return(gd_mod)
  }
}


# Plotting ------

## Summarize fixed effects
summary_fixed <- function(model, resp_var, pred_var){
  dat <- model$summary.fixed
  dat$response_var <- resp_var
  dat$predictor <- pred_var
  dat <- dat[-1,]
  names(dat)[c(3:7)] <- c('lo_95CI', 'lo_90CI', '50CI', 'up_90CI', 'up_95CI')
  return(dat)
}

## Summarize species random slopes
summary_random <- function(model, resp_var, pred_var, summary_fixed, dat1){
  gdata <- dat1
  n.species <- length(unique(gdata$species))
  
  dat <- model$summary.random$sp_id_s[c((n.species+1):(2*n.species)),]
  dat$response_var <- resp_var
  dat$predictor <- pred_var
  names(dat)[c(1, 4:8)] <- c('species', 'lo_95CI', 'lo_90CI', '50CI', 'up_90CI', 'up_95CI')
  
  ## Add species slopes (offsets) to main effect:
  dat$overall_effect <- summary_fixed$mean
  dat$sp_fx <- dat$mean + dat$overall_effect
  names(dat)[1] <- 'sp_id_s'
  return(dat)
}