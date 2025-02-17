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

# Models ----------
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

nb_to_df = function(nb, coords){
  x = coords[, 1]
  y = coords[, 2]
  n = length(nb)
  cardnb = card(nb)
  i = rep(1:n, cardnb)
  j = unlist(nb)
  return(data.frame(x=x[i], xend=x[j],
                    y=y[i], yend=y[j]))
}

run_INLA_model <- function(data, predictor, response, spatial=TRUE, Lattice.adj = NULL){
  
  data_temp <- data
  colnames(data_temp)[which(colnames(data_temp)==response)] = 'response'
  colnames(data_temp)[which(colnames(data_temp)==predictor)] = 'predictor'
  n.species <- length(unique(data_temp$species))
  
  if(spatial==FALSE){  
    formulae = response ~ 1 + predictor +                     # fixed effect
      f(sp_id, model="iid2d", n=2*n.species, constr = T) +    # random intercept (correlated with slope)
      f(sp_id_s, predictor, copy="sp_id")                     # random slope (correlated w/intercept)
    
    
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
  } else{
    formulae = response ~ 1 + predictor +                      # fixed effect
      f(sp_id, model="iid2d", n=2*n.species, constr = T) +     # random intercept (correlated with slope)
      f(sp_id_s, predictor, copy="sp_id") +                    # random slope (correlated w/intercept)
       f(site_ID, model = "besag",
        graph = Lattice.adj, constr = T,
        scale.model = T,
        hyper = hyper_besag)
    
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
    return(gd_mod)
  }
}

run_INLA_model_r1 <- function(data, predictor, response, spatial=TRUE, Lattice.adj = NULL){
  
  data_temp <- data
  colnames(data_temp)[which(colnames(data_temp)==response)] = 'response'
  colnames(data_temp)[which(colnames(data_temp)==predictor)] = 'predictor'
  n.species <- length(unique(data_temp$species))
  
  if(spatial==FALSE){
    if (response=='scale_arsp'){
      formulae =  response ~ predictor +
        f(sp_id, model="iid", constr = T) + 
        f(sp_id_s, predictor, model="iid", hyper = pcprior)
      
    } else{
      formulae = response ~ 1 + predictor +                     # fixed effect
        f(sp_id, model="iid2d", n=2*n.species, constr = T) +    # random intercept (correlated with slope)
        f(sp_id_s, predictor, copy="sp_id")
    }
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
  } else{
    if (response=='scale_arsp'){
      formulae =  response ~ predictor +
        f(sp_id, model="iid", constr = T) + 
        f(sp_id_s, predictor, model="iid", hyper = pcprior) +
        f(site_ID, model = "besag",
          graph = Lattice.adj, constr = T,
          scale.model = T,
          hyper = hyper_besag)
      
    } else{
      formulae =  response ~ predictor +
        f(sp_id, model="iid2d", n=2*n.species, constr = T) +
        f(sp_id_s, predictor, copy="sp_id") +
        f(site_ID, model = "besag",
          graph = Lattice.adj, constr = T,
          scale.model = T,
          hyper = hyper_besag)
      
    }
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
    return(gd_mod)
  }
}

# PA count buffer model
run_INLA_model_buf_r1 <- function(data, predictor, response, spatial=TRUE, Lattice.adj = NULL){
  
  data_temp <- data
  colnames(data_temp)[which(colnames(data_temp)==response)] = 'response'
  colnames(data_temp)[which(colnames(data_temp)==predictor)] = 'predictor'
  n.species <- length(unique(data_temp$species))
  
  if(spatial==FALSE){
    
    formulae = response ~ 1 + PA_bin + predictor +            # fixed effect
      f(sp_id, model="iid2d", n=2*n.species, constr = T) +    # random intercept (correlated with slope)
      f(sp_id_s2, predictor, copy="sp_id")                    # random slope (correlated w/intercept)
    
    
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
                  quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
  } 
  else{
    formulae = response ~ 1 + PA_bin + predictor +
      f(sp_id, model="iid2d", n=2*n.species, constr = T) +    # random intercept (correlated with slope)
      f(sp_id_s2, predictor, copy="sp_id") +                  # random slope (correlated w/intercept)
      f(site_ID, model = "besag",
        graph = Lattice.adj, constr = T,
        scale.model = T,
        hyper = hyper_besag)
    
    
    gd_mod = inla(formulae,
                  data = as.data.frame(data_temp),
                  control.predictor = list(compute = TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE,
                                         cpo = TRUE, return.marginals.predictor=TRUE),
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

summary_random_buf <- function(model, resp_var, pred_var, summary_fixed, dat1){
  gdata <- dat1
  n.species <- length(unique(gdata$species))
  summary_fixed_fil <- summary_fixed[2,]
  
  dat <- model$summary.random$sp_id_s2[c((n.species+1):(2*n.species)),]
  dat$response_var <- resp_var
  dat$predictor <- pred_var
  names(dat)[c(1, 4:8)] <- c('species', 'lo_95CI', 'lo_90CI', '50CI', 'up_90CI', 'up_95CI')
  
  ## Add species slopes (offsets) to main effect:
  dat$overall_effect <- summary_fixed_fil$mean
  dat$sp_fx <- dat$mean + dat$overall_effect
  names(dat)[1] <- 'sp_id_s'
  return(dat)
}

## Plot function ------
orchard_plot <- function(random_summary, fixed_summary, trunkcol, fruitcol){
  p <- ggplot() +
    # y axis lines that appear btwn groups:
    geom_vline(xintercept=seq(1.5, length(unique(fixed_summary$response_var))-0.5, 1),
               lwd=1, colour="gray90") +
    
    # Zero line
    geom_hline(yintercept = 0, colour = "black", lty = 2) +
    
    ## Fruit layer:
    geom_jitter(data = random_summary, aes(x = response_var, y = sp_fx, size = num_sites, color = sp_fx), 
                width = 0.3, alpha = 0.6) +
    
    # Make the legends the same so they are merged:
    scale_colour_gradientn(
      limits  = range(random_summary$sp_fx),
      colours = fruitcol[c(1, seq_along(fruitcol), length(fruitcol))],
      guide="none") +
    ## 90% CIs on overall effect size:
    geom_linerange(data = fixed_summary, aes(x = response_var, ymin = lo_90CI, ymax = up_90CI), 
                   lwd = 2.5, position = position_dodge(width = 1), color = trunkcol) +
    
    ## 95% CIs + overall coefficient:
    geom_pointrange(data = fixed_summary, aes(x = response_var, y = mean, ymin = lo_95CI,
                                              ymax = up_95CI),
                    lwd = 1, position = position_dodge(width = 1),
                    shape = 21, fill = "white", stroke = 3, color = trunkcol) +
    
    # Change axis text
    scale_x_discrete(labels = c(expression(population-specific~F[ST]),
                                'allelic richness', 'gene diversity', 'effective population size')) +
    
    # Everything else
    coord_flip() + 
    theme_minimal(base_size = 14) +
    theme(axis.ticks.y = element_blank()) +
    labs(x= "", y = "model coefficients", title = "", size = "sample sites") +
    theme(text=element_text(family="Roboto Medium"))
  
  return(p)
}