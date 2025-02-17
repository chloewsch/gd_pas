# 1b. Additional data prep ------
## Allelic richness rarefied to a different minimum per species
## Proportion of all alleles that are within PAs

library(adegenet)
library(hierfstat)
library(tidyr)
library(dplyr)
library(stringr)


# Terrestrial ------
# Re-estimate AR with species minimum
OGwd <- getwd()
# Directory with microsat data:
str <- 'str'
setwd(str)

file_list <- list.files(full.names = FALSE)
species_list <- list()

# Read in all files as genind objects
for (i in file_list){
  if(grepl(".str", i)){
    junk <- read.table(i)
    s <- read.structure(i, n.ind = (nrow(junk)), n.loc = ((ncol(junk)-1)/2),
                        onerowperind = TRUE, col.lab = 0, col.pop = 1, row.marknames = 0, ask = FALSE)}
  else {
    k <- try(read.genepop(i, ncode = 2))
    s <- if(inherits(k, "try-error")){read.genepop(i, ncode = 3)}
    else {
      s <- (read.genepop(i, ncode = 2))}
  }
  
  species_list[[i]] <- s
  
}

setwd(OGwd)

## Load metadata ------
# Population names, coordinates, and other metadata for merging
popn <- read.csv('genetdata_PA_terr.csv', header = T)

## List of names for all populations -----
allpoplist <- data.frame(unlist(lapply(species_list, popNames)))

## Allelic richness ------
# Find minimum number of alleles per species
# Species name from file
fi <- sub('\\..[^\\.]*$', '', file_list)
species_names_fi <- lapply(fi, function(f){
  s1 <- strsplit(f, split="_")[[1]][c(2,3)]
  sp <- paste(s1[1], s1[2], sep='_')
})
species_names_fi <- unlist(species_names_fi)
species_names_file_order <- data.frame(index = c(1:length(species_names_fi)),
                                       filename = list.files(str),
                                       species = species_names_fi)

# Minimum sample size per species:
min.n <- popn %>% 
  group_by(species) %>% 
  slice_min(num_individuals) %>% 
  summarise(min_n = mean(num_individuals))

min.n.fi <- merge(species_names_file_order, min.n, by='species', all.x = TRUE, all.y = FALSE)
min.n.fi <- min.n.fi[order(min.n.fi$index),]

# Re-rarefy with species minimum

arich0 <- list()

for(i in 1:length(species_list)){ 
  in_slice = min.n.fi %>% filter(index == i)
  minn = in_slice$min_n
  ar <- colMeans(allelic.richness(species_list[[i]], min.n = minn, diploid = TRUE)$Ar, 
                 na.rm = TRUE)
  arich0[[i]] <- data.frame(allelic_richness = ar,
                            filename = in_slice$filename)
}
arich1 <- do.call('rbind', arich0)

arich1 <- na.omit(arich1) # allelic.richness function creates dummy populations for data with only 1 population; delete
arich <- cbind.data.frame(allpoplist, arich1)
names(arich) <- c("pop","allelic_richness_rbysp", 'filename')
rownames(arich) <- NULL

## Count alleles -----
# Of the total number of alleles in ALL populations of a species, how many are in PAs?

pops_in_pas <- popn %>% 
  filter(PA_bin==1) %>% 
  select(filename, species, pop)

filelist <- list.files(str)
pa_allele_df <- list()

for(i in 1:length(species_list)){
  num_pops <- length(unique(species_list[[i]]$pop))
  print(i)
  # Count total alleles:
  pop_in_data <- popn %>% 
    filter(filename == filelist[i])
  
  if(num_pops==1 & nrow(pop_in_data)==1){
    print(paste('1/', i))
    total_alleles <- sum(summary(species_list[[i]])$loc.n.all)
    print(paste("1/ Total alleles = ", total_alleles, i))
  } else if(num_pops>1 & nrow(pop_in_data)==1){
    print(paste('2/', i))
    allpops_separated <- seppop(species_list[[i]], drop=TRUE)
    tot_alleles <- sum(summary(allpops_separated[[pop_in_data$pop]])$loc.n.all)
    total_alleles <- tot_alleles
    print(paste("2/ Total alleles = ", total_alleles, i))
  } else if(num_pops>1 & nrow(pop_in_data)>1){
    print(paste('3/', i))
    allpops_separated <- seppop(species_list[[i]], drop=FALSE)
    all_pops <- repool(allpops_separated[pop_in_data$pop])
    tot_alleles <- sum(summary(all_pops)$loc.n.all)
    total_alleles <- tot_alleles
    print(paste("3/ Total alleles = ", total_alleles, i))
  }
  
  
  # Count alleles in PAs:
  pipa <- pops_in_pas %>% 
    filter(filename == filelist[i])
  
  ## If no populations in PAs, number of alleles in PAs is 0:
  if(nrow(pipa)==0){ 
    pa_allele_count <- 0
    
    ## If only one population in PAs, count number of alleles in this population:  
  } else if(num_pops==1 & nrow(pipa)==1){
    pa_allele_count <- sum(summary(species_list[[i]])$loc.n.all)
    
    ## If more than one population in PAs, count number of alleles across all populations in PAs:  
  } else if(num_pops>1 & nrow(pipa)==1){
    pops_separated <- seppop(species_list[[i]], drop=TRUE)
    pa_alleles <- sum(summary(pops_separated[[pipa$pop]])$loc.n.all)
    pa_allele_count <- pa_alleles
    
  } else if(num_pops>1 & nrow(pipa)>1){
    pops_separated <- seppop(species_list[[i]], drop=FALSE)
    pa_pops <- repool(pops_separated[pipa$pop])
    pa_alleles <- sum(summary(pa_pops)$loc.n.all)
    pa_allele_count <- pa_alleles
  }
  
  pa_allele_df[[i]] <- data.frame(filename = filelist[i],
                                  total_alleles = total_alleles,
                                  pa_allele_count = pa_allele_count)
}

pa_allele_df <- do.call('rbind', pa_allele_df)

# Check that PA allele count is never greater than total allele count:
which(pa_allele_df$pa_allele_count>pa_allele_df$total_alleles)

# Attach to data:
gdat <- merge(popn, arich %>% select(pop, allelic_richness_rbysp), all.x = TRUE, all.y = FALSE, by = 'pop')
gdat <- merge(gdat, pa_allele_df, by = 'filename', all.x = TRUE, all.y = FALSE) 

# Check for NAs in allele count:
gdat[which(is.na(gdat$pa_allele_count)),]
gdat[which(is.na(gdat$total_alleles)),]

# Proportion of alleles in PAs:
gdat$prop_alleles_inPA <- gdat$pa_allele_count/gdat$total_alleles

#write.csv(gdat, 'genetdata_PA_terr.csv', quote = F, row.names = F)

# Marine -----
rm(list = ls())
# Re-estimate AR with species minimum
OGwd <- getwd()
# Directory with microsat data:
str <- 'marine_str'
setwd(str)

file_list <- list.files(full.names = FALSE)
species_list <- list()

# Read in all files as genind objects
for (i in file_list){
  if(grepl(".str", i)){
    junk <- read.table(i)
    s <- read.structure(i, n.ind = (nrow(junk)), n.loc = ((ncol(junk)-1)/2),
                        onerowperind = TRUE, col.lab = 0, col.pop = 1, row.marknames = 0, ask = FALSE)}
  else {
    k <- try(read.genepop(i, ncode = 2))
    s <- if(inherits(k, "try-error")){read.genepop(i, ncode = 3)}
    else {
      s <- (read.genepop(i, ncode = 2))}
  }
  
  species_list[[i]] <- s
  
}

setwd(OGwd)

# filenames with populations
popdf <- data.frame(filename = rep(file_list, unlist(lapply(species_list, function(x) length(levels(pop(x)))))),
                    pop = unlist(lapply(species_list, function(x) levels(pop(x)))))

## Load metadata ------
# Population names, coordinates, and other metadata for merging
grdat_env <- read.csv('genetdata_PA_mar.csv',h=T)

popn <- merge(grdat_env, popdf, by='pop', all.x=TRUE, all.y= FALSE)

# Attach filenames to populations
popdf <- data.frame(filename = rep(file_list, unlist(lapply(species_list, function(x) length(levels(pop(x)))))),
                    pop = unlist(lapply(species_list, function(x) levels(pop(x)))))
# Add filename to popn:
popn <- merge(grdat_env, popdf, by='pop', all.x=TRUE, all.y= FALSE)

# Check:
setdiff(popn$filename, list.files(str))
setdiff(list.files(str), popn$filename)

## List of names for all populations -----
allpoplist <- data.frame(unlist(lapply(species_list, popNames)))

## Allelic richness ------

# Link file names to species names:
species_names_fi <- popn %>% 
  select(filename, species) %>%
  distinct()
setdiff(species_names_fi$filename, list.files(str))

# Put df rows in the same order as file list
index <- match(file_list, species_names_fi$filename)

species_names_file_order <- species_names_fi[index,]
all.equal(species_names_file_order$filename, file_list)

species_names_file_order <- species_names_file_order %>% 
  mutate(index = index) %>% 
  select(index, filename, species)

# Minimum sample size per species:
min.n <- popn %>% 
  group_by(species) %>% 
  slice_min(num_individuals) %>% 
  summarise(min_n = mean(num_individuals))

min.n.fi <- merge(species_names_file_order, min.n, by='species', all.x = TRUE, all.y = FALSE)
min.n.fi <- min.n.fi[order(min.n.fi$index),]

# Re-rarefy with species minimum

arich0 <- list()

for(i in 1:length(species_list)){ 
  in_slice = min.n.fi %>% filter(index == i)
  minn = in_slice$min_n
  ar <- colMeans(allelic.richness(species_list[[i]], min.n = minn, diploid = TRUE)$Ar, 
                 na.rm = TRUE)
  arich0[[i]] <- data.frame(allelic_richness = ar,
                            filename = in_slice$filename)
}
arich1 <- do.call('rbind', arich0)

arich1 <- na.omit(arich1) # allelic.richness function creates dummy populations for data with only 1 population; delete
arich <- cbind.data.frame(allpoplist, arich1)
names(arich) <- c("pop","allelic_richness_rbysp", 'filename')
rownames(arich) <- NULL

## Count alleles -----
# Of the total number of alleles in ALL populations of a species, how many are in PAs?

pops_in_pas <- popn %>% 
  filter(PA_bin==1) %>% 
  select(filename, species, pop)

filelist <- list.files(str)
pa_allele_df <- list()

for(i in 1:length(species_list)){
  num_pops <- length(unique(species_list[[i]]$pop))
  print(i)
  # Count total alleles:
  pop_in_data <- popn %>% 
    filter(filename == filelist[i])
  
  if(num_pops==1 & nrow(pop_in_data)==1){
    print(paste('1/', i))
    total_alleles <- sum(summary(species_list[[i]])$loc.n.all)
    print(paste("1/ Total alleles = ", total_alleles, i))
  } else if(num_pops>1 & nrow(pop_in_data)==1){
    print(paste('2/', i))
    allpops_separated <- seppop(species_list[[i]], drop=TRUE)
    tot_alleles <- sum(summary(allpops_separated[[pop_in_data$pop]])$loc.n.all)
    total_alleles <- tot_alleles
    print(paste("2/ Total alleles = ", total_alleles, i))
  } else if(num_pops>1 & nrow(pop_in_data)>1){
    print(paste('3/', i))
    allpops_separated <- seppop(species_list[[i]], drop=FALSE)
    all_pops <- repool(allpops_separated[pop_in_data$pop])
    tot_alleles <- sum(summary(all_pops)$loc.n.all)
    total_alleles <- tot_alleles
    print(paste("3/ Total alleles = ", total_alleles, i))
  }
  
  
  # Count alleles in PAs:
  pipa <- pops_in_pas %>% 
    filter(filename == filelist[i])
  
  ## If no populations in PAs, number of alleles in PAs is 0:
  if(nrow(pipa)==0){ 
    pa_allele_count <- 0
    
    ## If only one population in PAs, count number of alleles in this population:  
  } else if(num_pops==1 & nrow(pipa)==1){
    pa_allele_count <- sum(summary(species_list[[i]])$loc.n.all)
    
    ## If more than one population in PAs, count number of alleles across all populations in PAs:  
  } else if(num_pops>1 & nrow(pipa)==1){
    pops_separated <- seppop(species_list[[i]], drop=TRUE)
    pa_alleles <- sum(summary(pops_separated[[pipa$pop]])$loc.n.all)
    pa_allele_count <- pa_alleles
    
  } else if(num_pops>1 & nrow(pipa)>1){
    pops_separated <- seppop(species_list[[i]], drop=FALSE)
    pa_pops <- repool(pops_separated[pipa$pop])
    pa_alleles <- sum(summary(pa_pops)$loc.n.all)
    pa_allele_count <- pa_alleles
  }
  
  pa_allele_df[[i]] <- data.frame(filename = filelist[i],
                                  total_alleles = total_alleles,
                                  pa_allele_count = pa_allele_count)
}

pa_allele_df <- do.call('rbind', pa_allele_df)

# Check that PA allele count is never greater than total allele count:
which(pa_allele_df$pa_allele_count>pa_allele_df$total_alleles)

# Attach to data:
gdat <- merge(popn, arich %>% select(pop, allelic_richness_rbysp), all.x = TRUE, all.y = FALSE, by = 'pop')
gdat <- merge(gdat, pa_allele_df, by = 'filename', all.x = TRUE, all.y = FALSE) 

# Check for NAs in allele count:
gdat[which(is.na(gdat$pa_allele_count)),]
gdat[which(is.na(gdat$total_alleles)),]

# Proportion of alleles in PAs:
gdat$prop_alleles_inPA <- gdat$pa_allele_count/gdat$total_alleles

#write.csv(gdat, 'genetdata_PA_mar.csv', quote = F, row.names = F)

