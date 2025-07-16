library(rstan)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)


# keep = 300
nchains = 1
iter = 500

##############################################
## 1. Load data
##############################################

# load built dataset 
dat = readRDS('data/tree_data_HARVARD.RDS')

# dat$sig_d_obs = sig_d_obs
sig_d_obs = readRDS('data/sig_d_obs_RWC.RDS')
dat$sig_d_obs = sig_d_obs

##############################################
## 2. compile and run model: RW 
##############################################

model_compiled_RW = stan_model('models/growth_model_RW.stan')

# fit and extract values
model_fit_RW = sampling(model_compiled_RW, 
                        data = dat, 
                        iter = iter, 
                        chains = nchains,
                        verbose=TRUE)


# get organized values for MCMC chains so we can make trace plots and also save
model_post_RW = rstan::extract(model_fit_RW, permuted = FALSE)

saveRDS(model_post_RW, file = file.path('output/model_post_RW_HARVARD.RDS'))

##############################################
## 2. compile and run model: CENSUS
##############################################

model_compiled_CENSUS = stan_model('models/growth_model_CENSUS.stan')

# fit and extract values
model_fit_CENSUS = sampling(model_compiled_RW_CENSUS, 
                            data = dat, 
                            iter = iter, 
                            chains = nchains,
                            verbose=TRUE)


# get organized values for MCMC chains so we can make trace plots and also save
model_post_CENSUS = rstan::extract(model_fit_CENSUS, permuted = FALSE)


saveRDS(model_post_CENSUS, file = file.path('output/model_post_CENSUS_HARVARD.RDS'))


##############################################
## 2. compile and run model: RW + CENSUS
##############################################

model_compiled_RW_CENSUS = stan_model('models/growth_model_RW_CENSUS.stan')

# fit and extract values
model_fit_RW_CENSUS = sampling(model_compiled_RW_CENSUS, 
                               data = dat, 
                               iter = iter, 
                               chains = nchains,
                               verbose=TRUE)


# get organized values for MCMC chains so we can make trace plots and also save
model_post_RW_CENSUS = rstan::extract(model_fit_RW_CENSUS, permuted = FALSE)


saveRDS(model_post_RW_CENSUS, file = file.path('output/model_post_RW_CENSUS_HARVARD.RDS'))

##############################################
## 3. organize model output
##############################################

fnames = paste0(c('model_post_RW',
                  'model_post_CENSUS',
                  'model_post_RW_CENSUS'), '_', 'HARVARD.RDS')
models = c('RW', 'Census', 'RW + Census')

nmodels = length(fnames)

# we only need "D" for diameter and we only need the iterations that we are planning to keep
post_d = list()
post_rw = list()
post_bt = list()
post_sx = list()
post_sxobs = list()
post_sdobs = list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  out = readRDS(paste0('output/', fname_model))
  
  # get all array slices for diameters
  variables = names(out[1,1,])
  allDs = grep('D\\[',variables)
  allRWs = grep('X\\[',variables)
  allBTs = grep('beta_t\\[',variables)
  allSX = grep('sig_x\\[',variables)
  allSXOBS = grep('sig_x_obs',variables)
  
  # we need to put into matrix for use in processing, some compile info from all chains
  out.temp.d = out[,,allDs]
  out.temp.rw = out[,,allRWs]
  out.temp.bt = out[,,allBTs]
  out.temp.sx = out[,,allSX]
  out.temp.sxobs = out[,,allSXOBS]
  # out.temp.sdobs = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allSDOBS]
  
  # if(nchains>1){
  #   out = out.temp[,1,]
  #   out.rw = out.temp.rw[,1,]
  #   out.bt = out.temp.bt[,1,]
  #   out.sx = out.temp.sx[,1,]
  #   
  #   for (j in 2:ncol(out.temp)){
  #     out = rbind(out, out.temp[,j,])
  #     out.rw = rbind(out.rw, out.temp.rw[,j,])
  #     out.sx = rbind(out.rw, out.temp.sx[,j,])
  #   }
  #   post[[i]] = out
  #   post_rw[[i]] = out.rw
  # } else {
  post_d[[i]] = out.temp.d
  post_rw[[i]] = out.temp.rw
  post_bt[[i]] = out.temp.bt
  post_sx[[i]] = out.temp.sx
  post_sx[[i]] = out.temp.sx
  post_sxobs[[i]] = out.temp.sxobs
  # }
  
  # if (models[i] == 'Census'){
  #   allSDOBS = grep('sig_d_obs',variables)
  #   out.temp.sdobs = out[,,allSDOBS]
  #   # if(nchains>1){
  #   #   out.sdobs = out.temp.sdobs[,1,]
  #   #   for (j in 2:ncol(out.temp.sdobs)){
  #   #     out.sdobs = rbind(out.sdobs, out.temp.sdobs[,j,])
  #   #   }
  #   # } else {
  #     out.sdobs = out.temp.sdobs
  #   # }
  #   post_sdobs_C = out.sdobs
  #   sigma_d_obs_C = median(post_sdobs_C)
  # }
  # 
  # if (models[i] == 'RW + Census'){
  #   allSDOBS = grep('sig_d_obs',variables)
  #   out.temp.sdobs = out[,,allSDOBS]
  #   if(nchains>1){
  #     out.sdobs = out.temp.sdobs[,1,]
  #     for (j in 2:ncol(out.temp.sdobs)){
  #       out.sdobs = rbind(out.sdobs, out.temp.sdobs[,j,])
  #     }
  #   } else {
  #     out.sdobs = out.temp.sdobs
  #   }
  #   post_sdobs_RWC = out.sdobs
  #   sigma_d_obs_RWC = median(post_sdobs_RWC)
  # }
  
}

# 
# # load built data for site 
# dat = readRDS(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input',paste0('tree_data_', site ,'_STAN_',mvers,'_', dvers, '.RDS')))
N_years = dat$N_years
N_Tr = dat$N_Tr
N_taxa = dat$N_taxa
X2Tr = dat$X2Tr
X2year = dat$X2year
Tr = dat$Tr %>% arrange(stat_id)
taxon = Tr$taxon
plot = Tr$plot
years = dat$years
year_lo = min(years)
year_hi = max(years)
# sig_d_obs = dat$sig_d_obs
pdbh = exp(dat$logTr)
# 
list2env(dat, envir = globalenv())

N_C = dat$N_C
X2C = dat$X2C
X2year_C = dat$X2year_C

allTrees = dat$allTrees %>% arrange(stat_id)
taxon_C = allTrees$taxon
plot_C = allTrees$plot
distance = allTrees$distance
# }
# 
# match species acronyms to level3a available species/pfts
choj = read.csv('data/acronym_to_chojnacky_v0.1.csv', stringsAsFactors = FALSE)

# use HAVI (average of all hardwoods) for those species not found in chojnacky equations
gen = choj[which(choj$acronym == 'HAVI'),]
choj = choj %>% filter(acronym %in% unique(taxon_C))

# #####################################################
# ################ 1a. Plot model and data ############
# #####################################################
# 
# 
# variables = names(model_post_RW[1,1,])
# allDs = grep('D\\[',variables)
# 
# post_D = list()
# post_D[[1]] = model_post_RW[,1,allDs]

#####################################################
################ 1a. Plot model and data ############
#####################################################

pdf('figures/tree_growth_model_ALL.pdf', width=10, height=6)
for (tree in 1:N_C){
  
  print(tree)
  in.RW = tree %in% X2Tr
  
  if (in.RW){
    inds = which(X2Tr == tree)
    yrinds = X2year[inds]
    
    dbh_iter = t(post_d[[1]][,inds])
    dbh_iter = data.frame(dbh_iter)
    dbh_iter = data.frame(year=years[yrinds], dbh_iter)
    
    dbh_mean = apply(dbh_iter[,2:ncol(dbh_iter)], 1, mean, na.rm=TRUE)
    dbh_quant = t(apply(dbh_iter[,2:ncol(dbh_iter)], 1, 
                        function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)))
    
    dbh_tree = data.frame(d_mean = dbh_mean, 
                          d_median = dbh_quant[,2], 
                          d_lo = dbh_quant[,1], 
                          d_hi = dbh_quant[,3], 
                          year = years[yrinds],
                          model = 'RW')
    
  }
  
  # determine which estimates correspond to this tree
  inds_C = which(X2Tr_C == tree)
  yrinds_C = X2year_C[inds_C]

  # C
  dbh_iter_C = t(post_d[[2]][,inds_C])
  dbh_iter_C = data.frame(dbh_iter_C)
  dbh_iter_C = data.frame(year=years[yrinds_C], dbh_iter_C)
  
  dbh_mean_C = apply(dbh_iter_C[,2:ncol(dbh_iter_C)], 1, mean, na.rm=TRUE)
  dbh_quant_C = t(apply(dbh_iter_C[,2:ncol(dbh_iter_C)], 1, 
                        function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)))
  
  dbh_tree_C = data.frame(d_mean = dbh_mean_C, 
                          d_median = dbh_quant_C[,2], 
                          d_lo = dbh_quant_C[,1], 
                          d_hi = dbh_quant_C[,3], 
                          year = years[yrinds_C],
                          model = 'Census')
  
  # RWC
  dbh_iter_RWC = t(post_d[[3]][,inds_C])
  dbh_iter_RWC = data.frame(dbh_iter_RWC)
  dbh_iter_RWC = data.frame(year=years[yrinds_C], dbh_iter_RWC)
  
  dbh_mean_RWC = apply(dbh_iter_RWC[,2:ncol(dbh_iter_RWC)], 1, mean, na.rm=TRUE)
  dbh_quant_RWC = t(apply(dbh_iter_RWC[,2:ncol(dbh_iter_RWC)], 1, 
                          function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)))
  
  dbh_tree_RWC = data.frame(d_mean = dbh_mean_RWC, 
                            d_median = dbh_quant_RWC[,2], 
                            d_lo = dbh_quant_RWC[,1], 
                            d_hi = dbh_quant_RWC[,3], 
                            year = years[yrinds_C],
                            model = 'RW + Census')
  
  if (in.RW){
    dbh_tree = rbind(dbh_tree,
                     dbh_tree_C,
                     dbh_tree_RWC)
  } else {
    dbh_tree = dbh_tree_C
  }
  
  idx_d_obs_C = which(Dobs$stat_id == tree)
  
  dbh_obs_C = data.frame(d_obs = Dobs$dbh[idx_d_obs_C],
                         year = years[Dobs$year[idx_d_obs_C]])
  
  stem_id = Dobs$ID[idx_d_obs_C][1]
  
  if (in.RW){
    idx_d_obs = which(Tr$stat_id == tree)
    
    dbh_obs = data.frame(d_obs = Tr$dbh[idx_d_obs],
                         year = years[Tr$year[idx_d_obs]])
    
  } else {
    dbh_obs = data.frame(d_obs = numeric(0),
                         year = numeric(0))
  }
  
  # Create a text
  grob = grobTree(textGrob(paste0('Tree ', tree, '; Stem ID ', stem_id, '; Species ', taxon[tree] ), x=0.05,  y=0.9, hjust=0,
                           gp=gpar(col="black", fontsize=22)))
  
  p1 = ggplot() +
    geom_ribbon(data=dbh_tree, aes(x=year, ymin=d_lo, ymax=d_hi, fill=model), alpha=0.5) +
    geom_line(data=dbh_tree, aes(x=year, y=d_median, colour=model)) +
    geom_point(data=dbh_obs, aes(x=year, y=d_obs), size=2) +
    geom_point(data=dbh_obs_C, aes(x=year, y=d_obs), size=2, shape=1) +
    xlab('year') +
    ylab('dbh (cm)') +
    xlim(c(year_lo, year_hi)) +
    theme_bw(16)  +
    annotation_custom(grob)
  
  print(p1)
  
  if (in.RW){
    inds = which(X2Tr == tree)
    yrinds = X2year[inds]
    
    rw_iter = t(post_rw[[1]][,inds])
    rw_iter = data.frame(rw_iter)
    rw_iter = data.frame(year=years[yrinds], rw_iter)
    
    rw_mean = apply(rw_iter[,2:ncol(rw_iter)], 1, mean, na.rm=TRUE)
    rw_quant = t(apply(rw_iter[,2:ncol(rw_iter)], 1, 
                       function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)))
    
    rw_tree = data.frame(x_mean = rw_mean, 
                         x_median = rw_quant[,2], 
                         x_lo = rw_quant[,1], 
                         x_hi = rw_quant[,3], 
                         year = years[yrinds],
                         model = 'RW')
  }
  
  
  inds_C = which(X2Tr_C == tree)
  yrinds_C = X2year_C[inds_C]
  
  rw_iter_C = t(post_rw[[2]][,inds_C])
  rw_iter_C = data.frame(rw_iter_C)
  rw_iter_C = data.frame(year=years[yrinds_C], rw_iter_C)
  
  rw_mean_C = apply(rw_iter_C[,2:ncol(rw_iter_C)], 1, mean, na.rm=TRUE)
  rw_quant_C = t(apply(rw_iter_C[,2:ncol(rw_iter_C)], 1, 
                       function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)))
  
  rw_tree_C = data.frame(x_mean = rw_mean_C, 
                         x_median = rw_quant_C[,2], 
                         x_lo = rw_quant_C[,1], 
                         x_hi = rw_quant_C[,3], 
                         year = years[yrinds_C],
                         model = 'Census')
  
  rw_iter_RWC = t(post_rw[[3]][,inds_C])
  rw_iter_RWC = data.frame(rw_iter_RWC)
  rw_iter_RWC = data.frame(year=years[yrinds_C], rw_iter_RWC)
  
  rw_mean_RWC = apply(rw_iter_RWC[,2:ncol(rw_iter_RWC)], 1, mean, na.rm=TRUE)
  rw_quant_RWC = t(apply(rw_iter_RWC[,2:ncol(rw_iter_RWC)], 1, 
                         function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm=TRUE)))
  
  rw_tree_RWC = data.frame(x_mean = rw_mean_RWC, 
                           x_median = rw_quant_RWC[,2], 
                           x_lo = rw_quant_RWC[,1], 
                           x_hi = rw_quant_RWC[,3], 
                           year = years[yrinds_C],
                           model = 'RW + Census')
  
  
  if (in.RW){
    rw_tree = rbind(rw_tree,
                    rw_tree_C,
                    rw_tree_RWC)
  } else {
    rw_tree = rw_tree_C
  }
  
  
  if (in.RW){
    idx_rw_obs = which(Xobs$stat_id == tree)
    
    rw_obs = data.frame(x_obs = Xobs$incr[idx_rw_obs],
                        year = years[Xobs$year[idx_rw_obs]])
  } else {
    rw_obs = data.frame(x_obs = numeric(0),
                        year = numeric(0))
  }
  
  # Create a text
  grob <- grobTree(textGrob(paste0('Tree ', tree, '; Stem ID ', stem_id, '; Species ', taxon[tree] ), x=0.05,  y=0.9, hjust=0,
                            gp=gpar(col="black", fontsize=22)))
  
  p2 = ggplot() +
    # geom_line(data=dbh_tree, aes(x=year, y=d_mean)) +
    # geom_ribbon(data=rw_tree, aes(x=year, ymin=x_lo, ymax=x_hi), fill='lightgrey') +
    geom_ribbon(data=rw_tree, aes(x=year, ymin=x_lo, ymax=x_hi, fill=model), alpha=0.5) +
    geom_line(data=rw_tree, aes(x=year, y=x_median, colour=model)) +
    geom_point(data=rw_obs, aes(x=year, y=x_obs), size=2, alpha=0.4) +
    # geom_dog(data=rw_obs, aes(x=year, y=x_obs, dog='glasses'), size=2) +
    # ylim(c(0,500)) +
    xlab('year') +
    ylab('rw (mm)') +
    xlim(c(year_lo, year_hi)) +
    theme_bw(16)  #+
  # ggtitle(paste0('Tree ', i)) +
  # annotation_custom(grob)
  
  # print(p2)
  
  grid.arrange(p1, p2, nrow = 2)
  
  
}
dev.off()

#####################################################
################ 2. Estimate biomass ################
#####################################################

# first for RW only model (no census)
dbh_array_RW = array(NA, dim = c(N_Tr, N_years, keep))
agb_array_RW = array(NA, dim = c(N_Tr, N_years, keep))

for (t in 1:N_Tr){
  
  # determine which estimates correspond to this tree
  inds = which(X2Tr == t)
  yrinds = X2year[inds]
  
  # extract diameter data
  dbh_array_RW[t,yrinds,] = t(post_d[[1]][,inds])
  
  # get equation coefficients based on taxon
  beta0 = choj$beta0[which(choj$acronym == taxon[t])]
  beta1 = choj$beta1[which(choj$acronym == taxon[t])]
  
  if (length(beta0)==0){
    beta0 = gen$beta0
    beta1 = gen$beta1
  }  
  
  # use biomass equation to estimate biomass from diameter
  agb_array_RW[t,,] = exp(beta0 + beta1 * log(dbh_array_RW[t,,]))
}

# next for RW + Census Model (if applicable)
dbh_array_C = array(NA, dim = c(N_C, N_years, keep))
agb_array_C = array(NA, dim = c(N_C, N_years, keep))

dbh_array_RWC = array(NA, dim = c(N_C, N_years, keep))
agb_array_RWC = array(NA, dim = c(N_C, N_years, keep))

for (t in 1:N_C){
  
  # determine which estimates correspond to this tree
  inds = which(X2C == t)
  yrinds = X2year_C[inds]
  
  # extract diameter data
  dbh_array_C[t,yrinds,] = t(post_d[[2]][,inds])
  dbh_array_RWC[t,yrinds,] = t(post_d[[3]][,inds])
  
  # get equation coefficients based on taxon
  beta0 = choj$beta0[which(choj$acronym == taxon_C[t])]
  beta1 = choj$beta1[which(choj$acronym == taxon_C[t])]
  
  # use biomass equation to estimate biomass from diameter
  agb_array_C[t,,] = exp(beta0 + beta1 * log(dbh_array_C[t,,]))
  agb_array_RWC[t,,] = exp(beta0 + beta1 * log(dbh_array_RWC[t,,]))
}
####################################################################
################ 3. Remove biomass from small trees ################
####################################################################

# first for RW only model
for (tree in 1:N_Tr){
  for (year in 1:N_years){
    
    # determine mean DBH for this year and tree
    dbh_mean = mean(dbh_array_RW[tree, year, ], na.rm=TRUE)
    
    # if smaller than 5 cm., eliminate the data 
    if (is.na(dbh_mean) | dbh_mean >= 5) next
    dbh_array_RW[tree, year, ] = rep(NA, keep)
    agb_array_RW[tree, year, ] = rep(NA, keep)
  }
}

# second for RW + Census model
for (tree in 1:N_C){
  for (year in 1:N_years){
    
    # determine mean DBH for this year and tree
    dbh_mean_C = mean(dbh_array_C[tree, year, ], na.rm=TRUE)
    
    # if smaller than 5 cm., eliminate the data 
    if (is.na(dbh_mean_C) | dbh_mean_C >= 5) next
    dbh_array_C[tree, year, ] = rep(NA, keep)
    agb_array_C[tree, year, ] = rep(NA, keep)
  }
}



for (tree in 1:N_C){
  for (year in 1:N_years){
    
    # determine mean DBH for this year and tree
    dbh_mean_RWC = mean(dbh_array_RWC[tree, year, ], na.rm=TRUE)
    
    # if smaller than 5 cm., eliminate the data 
    if (is.na(dbh_mean_RWC) | dbh_mean_RWC >= 5) next
    dbh_array_RWC[tree, year, ] = rep(NA, keep)
    agb_array_RWC[tree, year, ] = rep(NA, keep)
  }
}

###############################################################################
################ 4. Apply census smoothing (RW + CENSUS MODEL) ################
###############################################################################

# If the census is the final record of a tree and the tree was not alive for all of the censuses, 
# we need to determine the year in which the tree died stochastically since it could have been anytime 
# between the censuses.
# determine which trees we need to consider for smoothing (all those marked as "dead" in finalCond)
smoothIDs = unique(dat$Dobs$stat_id[which(dat$Dobs$finalCond == 'dead')])

# loop through all of the trees 
for (id in smoothIDs){
  
  # get last census year
  cYr = max(dat$Dobs$year[which(dat$Dobs$stat_id == id)]) 
  if (cYr == -Inf) next
  
  # get last data year 
  dYr = X2year_C[dat$idx_C$lastidx[which(dat$idx_C$stat_id == id)]]
  
  timeRange = seq(cYr, dYr)
  
  # loop through all iterations and stochastically choose the last year the tree lived  
  for (k in 1:keep){
    lYr = sample(timeRange, 1)
    if (lYr != dYr){
      agb_array_C[id,((lYr+1):dYr),k] = rep(NA, length(c((lYr+1):dYr)))
      agb_array_RWC[id,((lYr+1):dYr),k] = rep(NA, length(c((lYr+1):dYr)))
    } 
  }
}

################################################################
################ 5. Calculate biomass increment ################
################################################################

# determine biomass increment

# first for RW only model
abi = apply(agb_array_RW, c(1,3), function(x) diff(x))
abi = aperm(abi, c(2, 1, 3))
abi_melt = melt(abi)
colnames(abi_melt) = c('tree', 'year', 'iter', 'value')
abi_melt = abi_melt %>% filter(!is.na(value))
abi_melt$year = years[abi_melt$year]
abi_melt$plot = plot[abi_melt$tree]
abi_melt$taxon = taxon[abi_melt$tree]
abi_melt$model = rep("RW", nrow(abi_melt))
abi_melt$type = rep('abi',nrow(abi_melt))
# rm(abi)

# then for RW + census model 
abi_C = apply(agb_array_C, c(1,3), function(x) diff(x))
abi_C = aperm(abi_C, c(2, 1, 3))
abi_melt_C = melt(abi_C)
colnames(abi_melt_C) = c('tree', 'year', 'iter', 'value')
abi_melt_C = abi_melt_C %>% filter(!is.na(value))
abi_melt_C$year = years[abi_melt_C$year]
abi_melt_C$plot = plot_C[abi_melt_C$tree]
abi_melt_C$taxon = taxon_C[abi_melt_C$tree]
abi_melt_C$model = rep("Census", nrow(abi_melt_C))
abi_melt_C$type = rep('abi',nrow(abi_melt_C))
abi_melt = rbind(abi_melt, abi_melt_C)

abi_RWC = apply(agb_array_RWC, c(1,3), function(x) diff(x))
abi_RWC = aperm(abi_RWC, c(2, 1, 3))
abi_melt_RWC = melt(abi_RWC)
colnames(abi_melt_RWC) = c('tree', 'year', 'iter', 'value')
abi_melt_RWC = abi_melt_RWC %>% filter(!is.na(value))
abi_melt_RWC$year = years[abi_melt_RWC$year]
abi_melt_RWC$plot = plot_C[abi_melt_RWC$tree]
abi_melt_RWC$taxon = taxon_C[abi_melt_RWC$tree]
abi_melt_RWC$model = rep("RW + Census", nrow(abi_melt_RWC))
abi_melt_RWC$type = rep('abi',nrow(abi_melt_RWC))
abi_melt = rbind(abi_melt, abi_melt_RWC)

###################################################################
################ 6. Organize data into data frames ################
###################################################################

# melt down dbh_array to data frame
dbh_melt = melt(dbh_array_RW)
colnames(dbh_melt) = c('tree','year','iter','value')
dbh_melt = dbh_melt %>% filter(!is.na(value))
dbh_melt$year = years[dbh_melt$year]
dbh_melt$plot = plot[dbh_melt$tree]
dbh_melt$taxon = taxon[dbh_melt$tree]
dbh_melt$model = rep("RW", nrow(dbh_melt))
dbh_melt$type = rep('dbh',nrow(dbh_melt))
# rm(dbh_array)

# melt down agb_array to data frame
agb_melt = melt(agb_array_RW)
colnames(agb_melt) = c('tree','year','iter','value')
agb_melt = agb_melt %>% filter(!is.na(value))
agb_melt$year = years[agb_melt$year]
agb_melt$plot = plot[agb_melt$tree]
agb_melt$taxon = taxon[agb_melt$tree]
agb_melt$model = rep("RW", nrow(agb_melt))
agb_melt$type = rep('ab',nrow(agb_melt))

if (census_site){
  
  # melt down dbh_array_C to data frame
  dbh_melt_C = melt(dbh_array_C)
  colnames(dbh_melt_C) = c('tree','year','iter','value')
  dbh_melt_C = dbh_melt_C %>% filter(!is.na(value))
  dbh_melt_C$year = years[dbh_melt_C$year]
  dbh_melt_C$plot = plot_C[dbh_melt_C$tree]
  dbh_melt_C$taxon = taxon_C[dbh_melt_C$tree]
  dbh_melt_C$model = rep("Census", nrow(dbh_melt_C))
  dbh_melt_C$type = rep('dbh',nrow(dbh_melt_C))
  dbh_melt = rbind(dbh_melt, dbh_melt_C)
  # rm(dbh_array_C,dbh_melt_C)
  
  # melt down agb_array_C to data frame
  agb_melt_C = melt(agb_array_C)
  colnames(agb_melt_C) = c('tree','year','iter','value')
  agb_melt_C = agb_melt_C %>% filter(!is.na(value))
  agb_melt_C$year = years[agb_melt_C$year]
  agb_melt_C$plot = plot_C[agb_melt_C$tree]
  agb_melt_C$taxon = taxon_C[agb_melt_C$tree]
  agb_melt_C$model = rep("Census", nrow(agb_melt_C))
  agb_melt_C$type = rep('ab',nrow(agb_melt_C))
  agb_melt = rbind(agb_melt, agb_melt_C)
  
  # melt down dbh_array_C to data frame
  dbh_melt_RWC = melt(dbh_array_RWC)
  colnames(dbh_melt_RWC) = c('tree','year','iter','value')
  dbh_melt_RWC = dbh_melt_RWC %>% filter(!is.na(value))
  dbh_melt_RWC$year = years[dbh_melt_RWC$year]
  dbh_melt_RWC$plot = plot_C[dbh_melt_RWC$tree]
  dbh_melt_RWC$taxon = taxon_C[dbh_melt_RWC$tree]
  dbh_melt_RWC$model = rep("RW + Census", nrow(dbh_melt_RWC))
  dbh_melt_RWC$type = rep('dbh',nrow(dbh_melt_RWC))
  dbh_melt = rbind(dbh_melt, dbh_melt_RWC)
  
  # melt down agb_array_C to data frame
  agb_melt_RWC = melt(agb_array_RWC)
  colnames(agb_melt_RWC) = c('tree','year','iter','value')
  agb_melt_RWC = agb_melt_RWC %>% filter(!is.na(value))
  agb_melt_RWC$year = years[agb_melt_RWC$year]
  agb_melt_RWC$plot = plot_C[agb_melt_RWC$tree]
  agb_melt_RWC$taxon = taxon_C[agb_melt_RWC$tree]
  agb_melt_RWC$model = rep("RW + Census", nrow(agb_melt_RWC))
  agb_melt_RWC$type = rep('ab',nrow(agb_melt_RWC))
  agb_melt = rbind(agb_melt, agb_melt_RWC)
}

# remove incomplete rw/census years if applicable 
agb_melt = agb_melt %>% filter(year <= finalyr)
abi_melt = abi_melt %>% filter(year <= finalyr)
dbh_melt = dbh_melt %>% filter(year <= finalyr)


if (census_site){
  pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','tree_agb_model.pdf'), width=10, height=6)
  # for (tree in 1:N_Tr){
  for (tree in 1:N_C){
    
    print(tree)
    
    agb_tree = agb_melt[which(agb_melt$tree == tree),]
    agb_tree_quants = agb_tree %>% 
      dplyr::group_by(tree, year, model) %>%
      dplyr::summarize(agb_mean = mean(value, na.rm=TRUE),
                       agb_median = quantile(value, c(0.5)),
                       agb_lo = quantile(value, c(0.025)),
                       agb_hi = quantile(value, c(0.975)), .groups = 'keep') 
    
    species_id = agb_tree$taxon[1]  
    stem_id = Tr$id[which(Tr$stat_id == tree)]  
    
    grob <- grobTree(textGrob(paste0('Tree ', tree, '; Stem ID ', stem_id, '; Species ', species_id ), x=0.05,  y=0.9, hjust=0,
                              gp=gpar(col="black", fontsize=22)))
    
    p1 = ggplot() +
      geom_ribbon(data = agb_tree_quants, aes(x = year, ymin = agb_lo, ymax = agb_hi, fill = model), alpha=0.5) +
      geom_line(data=agb_tree_quants, aes(x=year, y=agb_median, colour = model)) +
      xlab('year') +
      ylab('agb (kg)') +
      xlim(c(year_lo, year_hi)) +
      theme_bw(16)  +
      # ggtitle(paste0('Tree ', i)) +
      annotation_custom(grob) 
    
    
    abi_tree = abi_melt[which(abi_melt$tree == tree),]
    abi_tree_quants = abi_tree %>% 
      dplyr::group_by(tree, year, model) %>%
      dplyr::summarize(abi_mean = mean(value, na.rm=TRUE),
                       abi_median = quantile(value, c(0.5)),
                       abi_lo = quantile(value, c(0.025)),
                       abi_hi = quantile(value, c(0.975)), .groups = 'keep') 
    
    p2 = ggplot() +
      geom_ribbon(data = abi_tree_quants, aes(x = year, ymin = abi_lo, ymax = abi_hi, fill = model), alpha=0.5) +
      geom_line(data=abi_tree_quants, aes(x=year, y=abi_median, colour = model)) +
      xlab('year') +
      ylab('abi (kg / year)') +
      xlim(c(year_lo, year_hi)) +
      theme_bw(16)  #+
    # ggtitle(paste0('Tree ', i)) +
    # annotation_custom(grob) 
    # print(p1)
    
    grid.arrange(p1, p2, nrow = 2)
    
  }
  dev.off()
}


####################################################################
################ 7. Save individual-level RDS files ################
####################################################################

# save AGB RDS and CSV files
filename = file.path(output_dir,paste0('AGB_STAN_',site,'_',mvers,'_',dvers))
saveRDS(agb_melt, paste0(filename,'.RDS'))

# save AGBI RDS and CSV files
filename2 = file.path(output_dir,paste0('AGBI_STAN_',site,'_',mvers,'_',dvers))
saveRDS(abi_melt, paste0(filename2,'.RDS'))

# save DBH RDS and CSV files
filename3 = file.path(output_dir,paste0('DBH_STAN_',site,'_',mvers,'_',dvers))
saveRDS(dbh_melt, paste0(filename3,'.RDS'))

################################################################################
################ 8. Perform sampling correction (RW ONLY MODEL) ################
################################################################################

# sampling correction adjusts biomass values to account for the PalEON sampling method 

# determine the measured diameter of trees at time of coring
pdbh = exp(dat$logTr)

# option to not apply the sampling correction
if (nest == 'nofix'){
  agb_array_RW = agb_array_RW * (1/(pi*plot_radius^2)) * (1/0.0001) * (1/1000)
}

if (nest == 'single'){
  agb_array_RW = agb_array_RW * (1/(pi*plot_radius^2)) * (1/0.0001) * (1/1000)
}

# double-nested plots 
if (nest == 'double'){
  idx_small  = which(pdbh<20)
  idx_large = which(pdbh>=20)
  
  # we need to adjust the biomass units from kg/plot to Mg/ha
  inner_factor = (1 / (pi*13^2)) * (1/0.0001) * (1/1000)
  outer_factor = (1 / (pi*20^2)) * (1/0.0001) * (1/1000)
  agb_array_RW[idx_small,,] = agb_array_RW[idx_small,,] * inner_factor
  agb_array_RW[idx_large,,] = agb_array_RW[idx_large,,] * outer_factor
}

# triple-nested plots
if (nest == 'triple'){
  idx_small = which(pdbh<20)
  idx_med = which((pdbh>=20) & (pdbh<30))
  idx_large = which(pdbh>=30)
  
  # we need to adjust the biomass units from kg/plot to Mg/ha
  inner_factor = (1 / (pi*13^2)) * (1/0.0001) * (1/1000)
  mid_factor = (1 / (pi*20^2)) * (1/0.0001) * (1/1000)
  outer_factor = (1 / (pi*30^2)) * (1/0.0001) * (1/1000)
  agb_array_RW[idx_small,,] = agb_array_RW[idx_small,,] * inner_factor
  agb_array_RW[idx_med,,] = agb_array_RW[idx_med,,] * mid_factor
  agb_array_RW[idx_large,,] = agb_array_RW[idx_large,,] * outer_factor
}

if (census_site){
  if (nest == 'double') {
    
    # kg/plot to Mg/ha
    
    # kg / plot * 1 Mg / 1000 Kg * 1 plot / 2.88 ha
    
    # agb_array_C = agb_array_C * (1 / 2.88) * (1/1000)
    agb_array_C = agb_array_C * (1 / (pi*20^2)) * (1/0.0001) * (1/1000)
    agb_array_RWC = agb_array_RWC * (1 / (pi*20^2)) * (1/0.0001) * (1/1000)
    
    # idx_small  = which(pdbh<20)
    # idx_large = which(pdbh>=20)
    # 
    # # we need to adjust the biomass units from kg/plot to Mg/ha
    # inner_factor = (1 / (pi*13^2)) * (1/0.0001) * (1/1000)
    # outer_factor = (1 / (pi*20^2)) * (1/0.0001) * (1/1000)
    # agb_array[idx_small,,] = agb_array[idx_small,,] * inner_factor
    # agb_array[idx_large,,] = agb_array[idx_large,,] * outer_factor
  }
  if (nest == 'triple') {
    agb_array_C = agb_array_C * (1 / (pi*30^2)) * (1/0.0001) * (1/1000)
    agb_array_RWC = agb_array_RWC * (1 / (pi*30^2)) * (1/0.0001) * (1/1000)
    
  }
  if (nest == 'nofix') {
    agb_array_C = agb_array_C * (1 / (pi*plot_radius^2)) * (1/0.0001) * (1/1000)
    agb_array_RWC = agb_array_RWC * (1 / (pi*plot_radius^2)) * (1/0.0001) * (1/1000)
    
  }
  if (nest == 'single') {
    agb_array_C = agb_array_C * (1 / (pi*plot_radius^2)) * (1/0.0001) * (1/1000)
    agb_array_RWC = agb_array_RWC * (1 / (pi*plot_radius^2)) * (1/0.0001) * (1/1000)
    
  }
}

# recreate agb data frame using corrected array 
agb_melt = melt(agb_array_RW)
colnames(agb_melt) = c('tree','year','iter','value')
agb_melt = agb_melt %>% filter(!is.na(value))
agb_melt$year = years[agb_melt$year]
agb_melt$plot = plot[agb_melt$tree]
agb_melt$taxon = taxon[agb_melt$tree]
agb_melt$model = rep("RW", nrow(agb_melt))
agb_melt$type = rep('ab',nrow(agb_melt))
#rm(agb_array)

if (census_site){
  # melt down agb_array_C to data frame
  agb_melt_C = melt(agb_array_C)
  colnames(agb_melt_C) = c('tree','year','iter','value')
  agb_melt_C = agb_melt_C %>% filter(!is.na(value))
  agb_melt_C$year = years[agb_melt_C$year]
  agb_melt_C$plot = plot_C[agb_melt_C$tree]
  agb_melt_C$taxon = taxon_C[agb_melt_C$tree]
  agb_melt_C$model = rep("Census", nrow(agb_melt_C))
  agb_melt_C$type = rep('ab',nrow(agb_melt_C))
  agb_melt = rbind(agb_melt, agb_melt_C)
  
  
  # melt down agb_array_C to data frame
  agb_melt_RWC = melt(agb_array_RWC)
  colnames(agb_melt_RWC) = c('tree','year','iter','value')
  agb_melt_RWC = agb_melt_RWC %>% filter(!is.na(value))
  agb_melt_RWC$year = years[agb_melt_RWC$year]
  agb_melt_RWC$plot = plot_C[agb_melt_RWC$tree]
  agb_melt_RWC$taxon = taxon_C[agb_melt_RWC$tree]
  agb_melt_RWC$model = rep("RW + Census", nrow(agb_melt_RWC))
  agb_melt_RWC$type = rep('ab',nrow(agb_melt_RWC))
  agb_melt = rbind(agb_melt, agb_melt_RWC)
  
}


# first for RW only model
abi = apply(agb_array_RW, c(1,3), function(x) diff(x))
abi = aperm(abi, c(2, 1, 3))
abi_melt = melt(abi)
colnames(abi_melt) = c('tree', 'year', 'iter', 'value')
abi_melt = abi_melt %>% filter(!is.na(value))
abi_melt$year = years[abi_melt$year]
abi_melt$plot = plot[abi_melt$tree]
abi_melt$taxon = taxon[abi_melt$tree]
abi_melt$model = rep("RW", nrow(abi_melt))
abi_melt$type = rep('abi',nrow(abi_melt))
# rm(abi)

# then for RW + census model
if (census_site){
  abi_C = apply(agb_array_C, c(1,3), function(x) diff(x))
  abi_C = aperm(abi_C, c(2, 1, 3))
  abi_melt_C = melt(abi_C)
  colnames(abi_melt_C) = c('tree', 'year', 'iter', 'value')
  abi_melt_C = abi_melt_C %>% filter(!is.na(value))
  abi_melt_C$year = years[abi_melt_C$year]
  abi_melt_C$plot = plot_C[abi_melt_C$tree]
  abi_melt_C$taxon = taxon_C[abi_melt_C$tree]
  abi_melt_C$model = rep("Census", nrow(abi_melt_C))
  abi_melt_C$type = rep('abi',nrow(abi_melt_C))
  abi_melt = rbind(abi_melt, abi_melt_C)
  
  
  abi_RWC = apply(agb_array_RWC, c(1,3), function(x) diff(x))
  abi_RWC = aperm(abi_RWC, c(2, 1, 3))
  abi_melt_RWC = melt(abi_RWC)
  colnames(abi_melt_RWC) = c('tree', 'year', 'iter', 'value')
  abi_melt_RWC = abi_melt_RWC %>% filter(!is.na(value))
  abi_melt_RWC$year = years[abi_melt_RWC$year]
  abi_melt_RWC$plot = plot_C[abi_melt_RWC$tree]
  abi_melt_RWC$taxon = taxon_C[abi_melt_RWC$tree]
  abi_melt_RWC$model = rep("RW + Census", nrow(abi_melt_RWC))
  abi_melt_RWC$type = rep('abi',nrow(abi_melt_RWC))
  abi_melt = rbind(abi_melt, abi_melt_RWC)
}

if (census_site){
  pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','tree_agb_model_fix.pdf'), width=10, height=6)
  for (tree in 1:N_C){
    
    print(tree)
    
    agb_tree = agb_melt[which(agb_melt$tree == tree),]
    agb_tree_quants = agb_tree %>% 
      dplyr::group_by(tree, year, model) %>%
      dplyr::summarize(agb_mean = mean(value, na.rm=TRUE),
                       agb_median = quantile(value, c(0.5)),
                       agb_lo = quantile(value, c(0.025)),
                       agb_hi = quantile(value, c(0.975)), .groups = 'keep') 
    
    species_id = agb_tree$taxon[1]  
    stem_id = Tr$id[which(Tr$stat_id == tree)]  
    
    grob <- grobTree(textGrob(paste0('Tree ', tree, '; Stem ID ', stem_id, '; Species ', species_id ), x=0.05,  y=0.9, hjust=0,
                              gp=gpar(col="black", fontsize=22)))
    
    p1 = ggplot() +
      geom_ribbon(data = agb_tree_quants, aes(x = year, ymin = agb_lo, ymax = agb_hi, colour=model, fill=model), alpha=0.5) +
      geom_line(data=agb_tree_quants, aes(x=year, y=agb_median, colour=model)) +
      xlab('year') +
      ylab('agb (Mg/ha)') +
      xlim(c(year_lo, year_hi)) +
      theme_bw(16)  +
      # ggtitle(paste0('Tree ', i)) +
      annotation_custom(grob)
    
    # print(p1)
    
    abi_tree = abi_melt[which(abi_melt$tree == tree),]
    abi_tree_quants = abi_tree %>% 
      dplyr::group_by(tree, year, model) %>%
      dplyr::summarize(abi_mean = mean(value, na.rm=TRUE),
                       abi_median = quantile(value, c(0.5)),
                       abi_lo = quantile(value, c(0.025)),
                       abi_hi = quantile(value, c(0.975)), .groups = 'keep') 
    
    p2 = ggplot() +
      geom_ribbon(data = abi_tree_quants, aes(x = year, ymin = abi_lo, ymax = abi_hi, fill = model), alpha=0.5) +
      geom_line(data=abi_tree_quants, aes(x=year, y=abi_median, colour = model)) +
      xlab('year') +
      ylab('abi (kg / year)') +
      xlim(c(year_lo, year_hi)) +
      theme_bw(16)  #+
    # ggtitle(paste0('Tree ', i)) +
    # annotation_custom(grob) 
    # print(p2)
    
    grid.arrange(p1, p2, nrow = 2)
    
    
  }
  dev.off()
}

##############################################################################
################ 9. Save taxon-level total aboveground biomass ################
##############################################################################

# sum annual biomass across taxa for each year, iteration, and plot 
# also remove incomplete final years if applicable
agb_taxa = agb_melt %>%
  filter(year <= finalyr) %>%
  dplyr::group_by(taxon, year, plot, iter, model) %>% 
  dplyr::summarize(ab = sum(value), .groups='keep')

# sum annual biomass across taxa for each year, iteration, and plot 
# also remove incomplete final years if applicable
abi_taxa = abi_melt %>%
  filter(year <= finalyr) %>%
  dplyr::group_by(taxon, year, plot, iter, model) %>% 
  dplyr::summarize(abi = sum(value), .groups='keep')

# save file
filename4 = file.path(output_dir,paste0('AGB_TAXA_STAN_',site,'_',mvers,'_',dvers))
saveRDS(agb_taxa, paste0(filename4,'.RDS'))

# save file
filename4b = file.path(output_dir,paste0('AGBI_TAXA_STAN_',site,'_',mvers,'_',dvers))
saveRDS(abi_taxa, paste0(filename4b,'.RDS'))
