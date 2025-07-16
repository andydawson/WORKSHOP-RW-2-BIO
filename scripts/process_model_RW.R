## Step 3: Process model output into biomass estimates 

## This script processes the output of the stat model into a format that is generally useful in PalEON and for the PEcAn workflow
## It will automatically process the data with and without the sampling correction depending on the site configuration file.

# Need to add plot_radius if no sampling correction 

process_rw_model <- function(census_site, mvers, dvers, site, nest,
                             finalyr = NULL, plot_radius = NULL, 
                             keep = 250, pool = 250, nchains = 3){
  
  ###############################################################
  ################ 1. Prepare workspace and data ################
  ###############################################################

  # decide how many models you need
  if (census_site){
    # fnames = paste0(c('ring_model_t_pdbh_sigd_species_sigxk_STAN', 'ring_model_t_pdbh_species_sigxk_STAN'), '_', 
    #                 site, '_',mvers,'_',dvers,'.RDS')
    # models = c('Model RW', 'Model RW + Census')
    fnames = paste0(c('ring_model_t_pdbh_sigd_species_sigxk_RW_STAN', 
                      'ring_model_t_pdbh_species_sigxk_C_STAN',
                      'ring_model_t_pdbh_species_sigxk_RWC_STAN'), '_', 
                    site, '_',mvers,'_',dvers,'.RDS')
    models = c('Model RW', 'Model Census', 'Model RW + Census')
  }else{
    # fnames = paste0(c('ring_model_t_pdbh_sigd_STAN'), '_', site, '_',mvers,'_',dvers,'.RDS')
    # fnames = paste0(c('ring_model_t_pdbh_sigd_species_STAN'), '_', site, '_',mvers,'_',dvers,'.RDS')
    fnames = paste0(c('ring_model_t_pdbh_sigd_species_sigxk_STAN'), '_', site, '_',mvers,'_',dvers,'.RDS')
    models = c('Model RW')
  }
  
  # extract stat model output 
  output_dir = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output')
  nmodels = length(fnames)
  
  # we only need "D" for diameter and we only need the iterations that we are planning to keep
  post = list()
  post_rw = list()
  post_bt = list()
  post_sx = list()
  post_sxobs = list()
  post_sdobs = list()
  for (i in 1:length(fnames)) {
    fname_model = fnames[i]
    out = readRDS(paste0(output_dir,'/', fname_model))
    
    # get all array slices for diameters 
    variables = names(out[1,1,])
    allDs = grep('D\\[',variables)
    allRWs = grep('X\\[',variables)
    allBTs = grep('beta_t\\[',variables)
    allSX = grep('sig_x\\[',variables)
    allSXOBS = grep('sig_x_obs',variables)
    
    # we need to put into matrix for use in processing, some compile info from all chains
    out.temp = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allDs]
    out.temp.rw = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allRWs]
    out.temp.bt = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allBTs]
    out.temp.sx = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allSX]
    out.temp.sxobs = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allSXOBS]
    # out.temp.sdobs = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allSDOBS]
    
    if(nchains>1){
      out = out.temp[,1,]
      out.rw = out.temp.rw[,1,]
      out.bt = out.temp.bt[,1,]
      out.sx = out.temp.sx[,1,]
      
      for (j in 2:ncol(out.temp)){
        out = rbind(out, out.temp[,j,])
        out.rw = rbind(out.rw, out.temp.rw[,j,])
        out.sx = rbind(out.rw, out.temp.sx[,j,])
      }
      post[[i]] = out
      post_rw[[i]] = out.rw
    } else {
      post[[i]] = out.temp
      post_rw[[i]] = out.temp.rw
      post_bt[[i]] = out.temp.bt
      post_sx[[i]] = out.temp.sx
      post_sx[[i]] = out.temp.sx
      post_sxobs[[i]] = out.temp.sxobs
    }
    
    if (models[i] == 'Model Census'){
      allSDOBS = grep('sig_d_obs',variables)
      out.temp.sdobs = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allSDOBS]
      if(nchains>1){
        out.sdobs = out.temp.sdobs[,1,]
        for (j in 2:ncol(out.temp.sdobs)){
          out.sdobs = rbind(out.sdobs, out.temp.sdobs[,j,])
        }
      } else {
        out.sdobs = out.temp.sdobs
      }
      post_sdobs_C = out.sdobs
      sigma_d_obs_C = median(post_sdobs_C)
    }
    
    if (models[i] == 'Model RW + Census'){
      allSDOBS = grep('sig_d_obs',variables)
      out.temp.sdobs = out[seq(dim(out)[1]-pool+1, dim(out)[1], pool/(keep/nchains)),,allSDOBS]
      if(nchains>1){
        out.sdobs = out.temp.sdobs[,1,]
        for (j in 2:ncol(out.temp.sdobs)){
          out.sdobs = rbind(out.sdobs, out.temp.sdobs[,j,])
        }
      } else {
        out.sdobs = out.temp.sdobs
      }
      post_sdobs_RWC = out.sdobs
      sigma_d_obs_RWC = median(post_sdobs_RWC)
    }
    
  }  
  rm(out, out.temp, allDs, variables)
  
  # load built data for site 
  dat = readRDS(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input',paste0('tree_data_', site ,'_STAN_',mvers,'_', dvers, '.RDS')))
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
  
  list2env(dat, envir = globalenv())
  
  if (is.null(finalyr)) finalyr = max(years)
  
  if (census_site){
    N_C = dat$N_C
    X2C = dat$X2C
    X2year_C = dat$X2year_C
    
    allTrees = dat$allTrees %>% arrange(stat_id)
    taxon_C = allTrees$taxon
    plot_C = allTrees$plot
    distance = allTrees$distance
  }
  
  # match species acronyms to level3a available species/pfts 
  choj = read.csv('data/acronym_to_chojnacky_v0.1.csv', stringsAsFactors = FALSE)
  
  # use HAVI (average of all hardwoods) for those species not found in chojnacky equations
  gen = choj[which(choj$acronym == 'HAVI'),]
  
  if (!census_site){
    choj = choj %>% filter(acronym %in% unique(taxon))
  }else{
    choj = choj %>% filter(acronym %in% unique(taxon_C))
  }
  
  #####################################################
  ################ 1a. Plot model and data ############
  #####################################################
  
  if (census_site){
    pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','tree_growth_model_data.pdf'), width=10, height=6)
    # pdf(paste0('figures/dbh_vs_year_estimated_', model, '.pdf'), width=10, height=6)
    
    # first for RW only model (no census)
    # dbh_array = array(NA, dim = c(N_Tr, N_years, keep))
    
    for (tree in 1:N_C){
      
      print(tree)
      
      # determine which estimates correspond to this tree
      
      in.RW = tree %in% X2Tr
      
      if (in.RW){
        inds = which(X2Tr == tree)
        yrinds = X2year[inds]
        
        # extract diameter data
        # dbh_array[t,yrinds,] = t(post[[1]][,inds])
        
        dbh_iter = t(post[[1]][,inds])
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
      
      # extract diameter data
      # dbh_array[t,yrinds,] = t(post[[1]][,inds])
      
      # C
      dbh_iter_C = t(post[[2]][,inds_C])
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
      dbh_iter_RWC = t(post[[3]][,inds_C])
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
      
      # idx_d_obs_C = which(X2Tr_C[Dobs2X] == tree)
      idx_d_obs_C = which(Dobs$stat_id == tree)
      
      dbh_obs_C = data.frame(d_obs = Dobs$dbh[idx_d_obs_C],
                             year = years[Dobs$year[idx_d_obs_C]])
      
      stem_id = Dobs$ID[idx_d_obs_C][1]
      
      if (in.RW){
        idx_d_obs = which(Tr$stat_id == tree)
        
        dbh_obs = data.frame(d_obs = Tr$dbh[idx_d_obs],
                             year = years[Tr$year[idx_d_obs]])
        
        # stem_id = Tr$id[idx_d_obs[1]]
      } else {
        dbh_obs = data.frame(d_obs = numeric(0),
                             year = numeric(0))
      }
      
      # Create a text
      grob = grobTree(textGrob(paste0('Tree ', tree, '; Stem ID ', stem_id, '; Species ', taxon[tree] ), x=0.05,  y=0.9, hjust=0,
                               gp=gpar(col="black", fontsize=22)))
      
      p1 = ggplot() +
        geom_ribbon(data=dbh_tree, aes(x=year, ymin=d_lo, ymax=d_hi, fill=model), alpha=0.5) +
        # geom_ribbon(data=dbh_tree, aes(x=year, ymin=d_lo, ymax=d_hi), fill='lightgrey') +
        geom_line(data=dbh_tree, aes(x=year, y=d_median, colour=model)) +
        geom_point(data=dbh_obs, aes(x=year, y=d_obs), size=2) +
        geom_point(data=dbh_obs_C, aes(x=year, y=d_obs), size=2, shape=1) +
        xlab('year') +
        ylab('dbh (cm)') +
        xlim(c(year_lo, year_hi)) +
        theme_bw(16)  +
        annotation_custom(grob)
      
      # print(p1)
      
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
  } 
  
 
  # # next for RW + Census Model (if applicable)
  # if (census_site){
  #   
  #   dbh_array_C = array(NA, dim = c(N_C, N_years, keep))
  #   agb_array_C = array(NA, dim = c(N_C, N_years, keep))
  #   
  #   for (tree in 1:N_C){
  #     
  #     # determine which estimates correspond to this tree
  #     inds = which(X2C == tree)
  #     yrinds = X2year_C[inds]
  #     
  #     # extract diameter data
  #     dbh_array_C[tree,yrinds,] = t(post[[2]][,inds])
  #     
  #     
  #     # extract diameter data
  #     # dbh_array[t,yrinds,] = t(post[[1]][,inds])
  #     
  #     # dbh_iter = t(post[[2]][,inds])
  #     # dbh_iter = data.frame(dbh_iter)
  #     # dbh_iter = data.frame(year=years[yrinds], dbh_iter)
  #     
  #     # get equation coefficients based on taxon
  #     beta0 = choj$beta0[which(choj$acronym == taxon_C[tree])]
  #     beta1 = choj$beta1[which(choj$acronym == taxon_C[tree])]
  #     
  #     # use biomass equation to estimate biomass from diameter
  #     agb_array_C[tree,,] = exp(beta0 + beta1 * log(dbh_array_C[tree,,]))
  #   }
  # }
  # 
  
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
    dbh_array_RW[t,yrinds,] = t(post[[1]][,inds])
    
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
  if (census_site){
    
    dbh_array_C = array(NA, dim = c(N_C, N_years, keep))
    agb_array_C = array(NA, dim = c(N_C, N_years, keep))
    
    dbh_array_RWC = array(NA, dim = c(N_C, N_years, keep))
    agb_array_RWC = array(NA, dim = c(N_C, N_years, keep))
    
    for (t in 1:N_C){
      
      # determine which estimates correspond to this tree
      inds = which(X2C == t)
      yrinds = X2year_C[inds]
      
      # extract diameter data
      dbh_array_C[t,yrinds,] = t(post[[2]][,inds])
      dbh_array_RWC[t,yrinds,] = t(post[[3]][,inds])
      
      # get equation coefficients based on taxon
      beta0 = choj$beta0[which(choj$acronym == taxon_C[t])]
      beta1 = choj$beta1[which(choj$acronym == taxon_C[t])]
      
      # use biomass equation to estimate biomass from diameter
      agb_array_C[t,,] = exp(beta0 + beta1 * log(dbh_array_C[t,,]))
      agb_array_RWC[t,,] = exp(beta0 + beta1 * log(dbh_array_RWC[t,,]))
    }
  }
  
  # rm(post)
  
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
  if (census_site){
    
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
  }
  
  ###############################################################################
  ################ 4. Apply census smoothing (RW + CENSUS MODEL) ################
  ###############################################################################
  
  # If the census is the final record of a tree and the tree was not alive for all of the censuses, 
  # we need to determine the year in which the tree died stochastically since it could have been anytime 
  # between the censuses.
  
  if (census_site){
    
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
  abi_melt$model = rep("Model RW", nrow(abi_melt))
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
    abi_melt_C$model = rep("Model Census", nrow(abi_melt_C))
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
    abi_melt_RWC$model = rep("Model RW + Census", nrow(abi_melt_RWC))
    abi_melt_RWC$type = rep('abi',nrow(abi_melt_RWC))
    abi_melt = rbind(abi_melt, abi_melt_RWC)
    # rm(abi_C,abi_melt_C)
    
    
    
  }
  
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
  dbh_melt$model = rep("Model RW", nrow(dbh_melt))
  dbh_melt$type = rep('dbh',nrow(dbh_melt))
  # rm(dbh_array)
  
  # melt down agb_array to data frame
  agb_melt = melt(agb_array_RW)
  colnames(agb_melt) = c('tree','year','iter','value')
  agb_melt = agb_melt %>% filter(!is.na(value))
  agb_melt$year = years[agb_melt$year]
  agb_melt$plot = plot[agb_melt$tree]
  agb_melt$taxon = taxon[agb_melt$tree]
  agb_melt$model = rep("Model RW", nrow(agb_melt))
  agb_melt$type = rep('ab',nrow(agb_melt))
  
  if (census_site){
    
    # melt down dbh_array_C to data frame
    dbh_melt_C = melt(dbh_array_C)
    colnames(dbh_melt_C) = c('tree','year','iter','value')
    dbh_melt_C = dbh_melt_C %>% filter(!is.na(value))
    dbh_melt_C$year = years[dbh_melt_C$year]
    dbh_melt_C$plot = plot_C[dbh_melt_C$tree]
    dbh_melt_C$taxon = taxon_C[dbh_melt_C$tree]
    dbh_melt_C$model = rep("Model Census", nrow(dbh_melt_C))
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
    agb_melt_C$model = rep("Model Census", nrow(agb_melt_C))
    agb_melt_C$type = rep('ab',nrow(agb_melt_C))
    agb_melt = rbind(agb_melt, agb_melt_C)
    
    # melt down dbh_array_C to data frame
    dbh_melt_RWC = melt(dbh_array_RWC)
    colnames(dbh_melt_RWC) = c('tree','year','iter','value')
    dbh_melt_RWC = dbh_melt_RWC %>% filter(!is.na(value))
    dbh_melt_RWC$year = years[dbh_melt_RWC$year]
    dbh_melt_RWC$plot = plot_C[dbh_melt_RWC$tree]
    dbh_melt_RWC$taxon = taxon_C[dbh_melt_RWC$tree]
    dbh_melt_RWC$model = rep("Model RW + Census", nrow(dbh_melt_RWC))
    dbh_melt_RWC$type = rep('dbh',nrow(dbh_melt_RWC))
    dbh_melt = rbind(dbh_melt, dbh_melt_RWC)
    
    # melt down agb_array_C to data frame
    agb_melt_RWC = melt(agb_array_RWC)
    colnames(agb_melt_RWC) = c('tree','year','iter','value')
    agb_melt_RWC = agb_melt_RWC %>% filter(!is.na(value))
    agb_melt_RWC$year = years[agb_melt_RWC$year]
    agb_melt_RWC$plot = plot_C[agb_melt_RWC$tree]
    agb_melt_RWC$taxon = taxon_C[agb_melt_RWC$tree]
    agb_melt_RWC$model = rep("Model RW + Census", nrow(agb_melt_RWC))
    agb_melt_RWC$type = rep('ab',nrow(agb_melt_RWC))
    agb_melt = rbind(agb_melt, agb_melt_RWC)
  }
  
  # remove incomplete rw/census years if applicable 
  agb_melt = agb_melt %>% filter(year <= finalyr)
  abi_melt = abi_melt %>% filter(year <= finalyr)
  dbh_melt = dbh_melt %>% filter(year <= finalyr)
  
  # if ()
  #   pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','tree_agb_model.pdf'), width=10, height=6)
  # for (tree in 1:N_Tr){
  #   
  #   print(tree)
  #   
  #   agb_tree = agb_melt[which(agb_melt$tree == tree),]
  #   agb_tree_quants = agb_tree %>% 
  #     group_by(tree, year, model) %>%
  #     summarize(agb_mean = mean(value, na.rm=TRUE),
  #               agb_median = quantile(value, c(0.5)),
  #               agb_lo = quantile(value, c(0.025)),
  #               agb_hi = quantile(value, c(0.975)), .groups = 'keep') 
  #   
  #   species_id = agb_tree$taxon[1]  
  #   stem_id = Tr$id[which(Tr$stat_id == tree)]  
  #   
  #   grob <- grobTree(textGrob(paste0('Tree ', tree, '; Stem ID ', stem_id, '; Species ', species_id ), x=0.05,  y=0.9, hjust=0,
  #                             gp=gpar(col="black", fontsize=22)))
  #   
  #   p1 = ggplot() +
  #     geom_ribbon(data = agb_tree_quants, aes(x = year, ymin = agb_lo, ymax = agb_hi, fill = model), alpha=0.5) +
  #     geom_line(data=agb_tree_quants, aes(x=year, y=agb_median, colour = model)) +
  #     xlab('year') +
  #     ylab('agb (kg)') +
  #     xlim(c(year_lo, year_hi)) +
  #     theme_bw(16)  +
  #     # ggtitle(paste0('Tree ', i)) +
  #     annotation_custom(grob) 
  #   
  #   print(p1)
  #   
  # }
  # dev.off()
  
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
  agb_melt$model = rep("Model RW", nrow(agb_melt))
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
    agb_melt_C$model = rep("Model Census", nrow(agb_melt_C))
    agb_melt_C$type = rep('ab',nrow(agb_melt_C))
    agb_melt = rbind(agb_melt, agb_melt_C)
    
    
    # melt down agb_array_C to data frame
    agb_melt_RWC = melt(agb_array_RWC)
    colnames(agb_melt_RWC) = c('tree','year','iter','value')
    agb_melt_RWC = agb_melt_RWC %>% filter(!is.na(value))
    agb_melt_RWC$year = years[agb_melt_RWC$year]
    agb_melt_RWC$plot = plot_C[agb_melt_RWC$tree]
    agb_melt_RWC$taxon = taxon_C[agb_melt_RWC$tree]
    agb_melt_RWC$model = rep("Model RW + Census", nrow(agb_melt_RWC))
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
  abi_melt$model = rep("Model RW", nrow(abi_melt))
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
    abi_melt_C$model = rep("Model Census", nrow(abi_melt_C))
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
    abi_melt_RWC$model = rep("Model RW + Census", nrow(abi_melt_RWC))
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
  
  #####################################################################################################
  ################ 10. Calculate approximate diameters from measured DBH and increments ################
  #####################################################################################################
  
  # this section calculates aboveground biomass based on the RW and the measured DBH values 
  dbh_data = matrix(NA, N_Tr, N_years)
  agb_data = matrix(NA, N_Tr, N_years)
  
  for (i in 1:N_Tr){
    print(i)
    
    # get data for this tree
    data_now = dat$Xobs %>% filter(stat_id == i)
    
    if (nrow(data_now) == 0){
      print(paste0('No RW for tree with stat id: ', i))
      next
    }
    
    dbh_now = pdbh[i]
    taxon_now = taxon[i]
    yrs = sort(unique(data_now$year), decreasing = TRUE)
    
    # set diameter at time of coring
    dbh_data[i,yrs[1]] = dbh_now
    dbh_last = dbh_now
    
    # loop through years with available data
    for (j in 2:length(yrs)){
      yr = yrs[j]
      # incr = min((data_now %>% filter(year == yr))$incr, na.rm = TRUE)
      incr = mean((data_now %>% filter(year == yr))$incr, na.rm = TRUE)
      dbh_temp = dbh_last - (2 * incr/10)
      
      # make sure diameter is not less than 5 cm
      if (dbh_temp < 5){
        dbh_last = 0 
      }else{
        dbh_data[i,yr] = dbh_temp 
        dbh_last = dbh_temp
      }
    }
    
    # then determine biomass based on chojnacky equations for this species
    beta0 = choj$beta0[which(choj$acronym == taxon_now)]
    beta1 = choj$beta1[which(choj$acronym == taxon_now)]
    
    # use general if value doesn't exist for the taxon
    if (length(beta0)==0){
      beta0 = gen$beta0
      beta1 = gen$beta1
    }  
    
    agb_data[i,] = exp(beta0 + beta1 * log(dbh_data[i,]))
  }
  
  # rescale for sampling correction
  
  # if we do not want to apply sampling correction
  if (nest == 'nofix'){
    agb_data = agb_data * (1/(pi*plot_radius^2)) * (1/0.0001) * (1/1000)
  }
  
  # double-nested plots 
  if (nest == 'double'){
    agb_data[idx_small,] = agb_data[idx_small,] * inner_factor
    agb_data[idx_large,] = agb_data[idx_large,] * outer_factor
  }
  
  # triple-nested plots (most other plots)
  if (nest == 'triple'){
    agb_data[idx_small,] = agb_data[idx_small,] * inner_factor
    agb_data[idx_med,] = agb_data[idx_med,] * mid_factor
    agb_data[idx_large,] = agb_data[idx_large,] * outer_factor
  }
  
  if (nest == 'single') {
    agb_data = agb_data * (1 / (pi*plot_radius^2)) * (1/0.0001) * (1/1000)
  }
  
  # melt to data frames
  data_melt = melt(agb_data)
  colnames(data_melt) = c("tree", "year", "value")
  data_melt = data_melt %>% filter(!is.na(value))
  data_melt$year = years[data_melt$year]
  data_melt$plot = plot[data_melt$tree]
  data_melt$taxon = taxon[data_melt$tree]
  data_melt$type = rep('ab',nrow(data_melt))
  # rm(agb_data, dbh_data, data_now)
  
  # remove incomplete years 
  data_melt = data_melt %>% filter(year <= finalyr)
  
  #############################################
  ################ 11. Figures ################
  #############################################
  
  N_plots = length(unique(plot))
  
  # determine quantiles for graphing 
  agb_plot = agb_taxa %>%
    dplyr::group_by(model, plot, taxon, year) %>%
    dplyr::summarize(ab025 = quantile(ab, 0.025),
                     ab50 = quantile(ab, 0.5),
                     ab975 = quantile(ab, 0.975), 
                     .groups = 'keep') %>% 
    ungroup()
  
  agb_plot$plot =  as.numeric(agb_plot$plot)
  
  sum_plot = agb_taxa %>%
    dplyr::group_by(model, plot, year, iter) %>%
    dplyr::summarize(ab = sum(ab), .groups = 'keep') %>%
    ungroup() %>% 
    group_by(model, plot, year) %>%
    dplyr::summarize(ab025 = quantile(ab, 0.025),
                     ab50 = quantile(ab, 0.5),
                     ab975 = quantile(ab, 0.975), 
                     .groups = 'keep')
  
  data_pft_plot = data_melt %>% 
    group_by(year, plot, taxon) %>% 
    dplyr::summarize(ab = sum(value), .groups = 'keep') %>%
    ungroup()
  data_pft_plot$plot = as.numeric(data_pft_plot$plot) 
  
  data_plot = data_melt %>%
    group_by(year, plot) %>%
    dplyr::summarize(ab = sum(value), .groups = 'keep') %>% 
    ungroup()
  
  agb_pft_plot_data = left_join(agb_plot, data_pft_plot, by = c('plot','taxon','year'))
  
  # figure to compare biomass by PFT for each plot 
  pl1 = ggplot(data = left_join(agb_plot, data_pft_plot, by = c('plot','taxon','year'))) + 
    facet_grid(~plot~as.factor(model)) + 
    geom_line(aes(x = year, y = ab50, group = taxon, color = taxon)) + 
    geom_ribbon(aes(x = year, ymin = ab025, ymax = ab975, 
                    fill = taxon, group = taxon, color = taxon), alpha = 0.4) +
    geom_line(aes(x = year, y = ab, group = taxon)) +
    theme_bw() + 
    labs(x = 'Year', y = 'Biomass (Mg/ha)', color = 'Species', fill = 'Species',
         title = 'Aboveground Biomass by PFT')
  print(pl1)
  ggsave(pl1, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_AGB.pdf'), width=10, height=8)
  # ggsave(pl1, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_AGB.jpg'))
  # ggsave(pl1, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_AGB.jpg'))
  
  # determine quantiles for graphing 
  abi_plot = abi_taxa %>%
    dplyr::group_by(model, plot, taxon, year) %>%
    dplyr::summarize(abi025 = quantile(abi, 0.025),
                     abi50 = quantile(abi, 0.5),
                     abi975 = quantile(abi, 0.975), 
                     .groups = 'keep') %>% 
    ungroup()
  
  abi_sum_plot = abi_taxa %>%
    dplyr::group_by(model, plot, year, iter) %>%
    dplyr::summarize(abi = sum(abi), .groups = 'keep') %>%
    ungroup() %>% 
    group_by(model, plot, year) %>%
    dplyr::summarize(abi025 = quantile(abi, 0.025),
                     abi50 = quantile(abi, 0.5),
                     abi975 = quantile(abi, 0.975), 
                     .groups = 'keep')
  
  # data_pft_plot = data_melt %>% 
  #   group_by(year, plot, taxon) %>% 
  #   summarize(ab = sum(value), .groups = 'keep') %>%
  #   ungroup()
  # 
  # data_plot = data_melt %>%
  #   group_by(year, plot) %>%
  #   summarize(ab = sum(value), .groups = 'keep') %>% 
  #   ungroup()
  
  # figure to compare biomass increment by PFT for each plot 
  # pl1 = ggplot(data = left_join(abi_plot, data_pft_plot, by = c('plot','taxon','year'))) + 
  pl1 = ggplot(data = abi_plot) + 
    facet_grid(~plot~as.factor(model)) + 
    geom_line(aes(x = year, y = abi50, group = taxon, color = taxon)) + 
    geom_ribbon(aes(x = year, ymin = abi025, ymax = abi975, 
                    fill = taxon, group = taxon, color = taxon), alpha = 0.4) +
    # geom_line(aes(x = year, y = abi, group = taxon)) +
    theme_bw() + 
    labs(x = 'Year', y = 'Biomass Increment (Mg/ha)', color = 'Species', fill = 'Species',
         title = 'Aboveground Biomass Increment by PFT')
  print(pl1)
  ggsave(pl1, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_ABI.pdf'), width=10, height=8)
  # ggsave(pl1, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_ABI.jpg'))
  
  
  # figure to compare total biomass for each plot 
  pl2 = ggplot() + 
    geom_line(data = abi_sum_plot, aes(x = year, y = abi50, 
                                       group = model, color = model)) + 
    geom_ribbon(data = abi_sum_plot, aes(x = year, ymin = abi025, ymax = abi975, 
                                         group = model, color = model, fill = model), alpha = 0.4) +
    # geom_line(data = data_plot, aes(x = year, y = abi)) + 
    facet_wrap(~plot) + 
    theme_bw() +
    labs(x = 'Year', y = 'Biomass Increment (Mg/ha)', color = 'Model', fill = 'Model', 
         title = "Total Aboveground Biomass Increment by Plot") 
  print(pl2)
  ggsave(pl2, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_total_plot_ABI.pdf'), width=10, height=8)
  # ggsave(pl2, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_total_plot_ABI.jpg'))
  
  # figure to compare total site biomass 
  sum_site = agb_taxa %>%
    group_by(model, plot, year, iter) %>%
    dplyr::summarize(ab = sum(ab), .groups = 'keep') %>% 
    ungroup() %>% 
    group_by(model, iter,  year) %>%
    dplyr::summarize(ab = mean(ab), .groups = 'keep') %>% 
    ungroup() %>% 
    group_by(model, year) %>% 
    dplyr::summarize(ab025 = quantile(ab, 0.025),
                     ab50 = quantile(ab, 0.5),
                     ab975 = quantile(ab, 0.975), .groups = 'keep')
  
  data_site = data_plot %>% 
    group_by(year) %>%
    dplyr::summarize(ab = mean(ab), .groups = 'keep')
  
  pl3 = ggplot(sum_site) +
    geom_line(data = sum_site, aes(x = year, y = ab50, 
                                   group = model, color = model)) + 
    geom_ribbon(data = sum_site, aes(x = year, ymin = ab025, ymax = ab975, 
                                     group = model, color = model, fill = model), alpha = 0.4) +
    # geom_line(data = data_site, aes(x = year, y = ab)) + 
    theme_bw() +
    labs(x = 'Year', y = 'Biomass (Mg/ha)', color = 'Model', fill = 'Model', 
         title = "Total Aboveground Biomass") 
  print(pl3)
  ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_total_site_AGB.pdf'), width=10, height=8)
  # ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_total_site_AGB.jpg'))
  
  # figure to compare total site biomass 
  abi_sum_site = abi_taxa %>%
    group_by(model, plot, year, iter) %>%
    dplyr::summarize(abi = sum(abi), 
                     .groups = 'keep') %>% 
    ungroup() %>% 
    group_by(model, iter,  year) %>%
    dplyr::summarize(abi = mean(abi), 
                     .groups = 'keep') %>% 
    ungroup() %>% 
    group_by(model, year) %>% 
    dplyr::summarize(abi025 = quantile(abi, 0.025),
                     abi50 = quantile(abi, 0.5),
                     abi975 = quantile(abi, 0.975), 
                     .groups = 'keep')
  
  # data_site = data_plot %>% 
  #   group_by(year) %>%
  #   summarize(ab = mean(ab))
  
  pl3 = ggplot(abi_sum_site) +
    geom_line(data = abi_sum_site, aes(x = year, y = abi50, 
                                       group = model, color = model)) + 
    geom_ribbon(data = abi_sum_site, aes(x = year, ymin = abi025, ymax = abi975, 
                                         group = model, color = model, fill = model), alpha = 0.4) +
    # geom_line(data = data_site, aes(x = year, y = ab)) + 
    theme_bw() +
    labs(x = 'Year', y = 'Biomass (Mg/ha)', color = 'Model', fill = 'Model', 
         title = "Total Aboveground Biomass Increment") 
  print(pl3)
  ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_total_site_ABI.pdf'), width=10, height=8)
  # ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_total_site_ABI.jpg'))
  
  
  
  # if (N_plots > 1){
  
  # figure to compare pft site biomass
  sum_site_pft = agb_taxa %>%
    group_by(model, taxon, year, iter) %>%
    dplyr::summarize(ab = mean(ab), .groups = 'keep') %>%
    ungroup() %>%
    group_by(model, taxon, year) %>%
    dplyr::summarize(ab_mean = mean(ab),
                     ab025 = quantile(ab, 0.025),
                     ab50 = quantile(ab, 0.5),
                     ab975 = quantile(ab, 0.975), .groups = 'keep')
  
  # data_site = data_plot %>%
  #   group_by(year) %>%
  #   summarize(ab = mean(ab), .groups = 'keep')
  
  pl3 = ggplot(sum_site_pft) +
    geom_line(data = sum_site_pft, aes(x = year, y = ab50,
                                       group = taxon, color = taxon)) +
    geom_ribbon(data = sum_site_pft, aes(x = year, ymin = ab025, ymax = ab975,
                                         group = taxon, color = taxon, fill = taxon), alpha = 0.4) +
    # geom_line(data = data_site, aes(x = year, y = ab)) +
    theme_bw() +
    labs(x = 'Year', y = 'Biomass (Mg/ha)', color = 'Model', fill = 'Model',
         title = "Total Aboveground Biomass") +
    facet_grid(.~model)
  print(pl3)
  ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_site_AGB.pdf'), width=10, height=8)
  # ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_site_AGB.jpg'))
  
  # figure to compare total site biomass
  abi_sum_site_pft = abi_taxa %>%
    group_by(model, taxon, year, iter) %>%
    dplyr::summarize(abi = sum(abi),
                     .groups = 'keep') %>%
    ungroup() %>%
    group_by(model, taxon, year) %>%
    dplyr::summarize(abi_mean = mean(abi),
                     abi025 = quantile(abi, 0.025),
                     abi50 = quantile(abi, 0.5),
                     abi975 = quantile(abi, 0.975),
                     .groups = 'keep')
  # 
  # # data_site = data_plot %>% 
  # #   group_by(year) %>%
  # #   summarize(ab = mean(ab))
  # 
  pl3 = ggplot(abi_sum_site_pft) +
    geom_line(data = abi_sum_site_pft, aes(x = year, y = abi50,
                                           group = taxon, color = taxon)) +
    geom_ribbon(data = abi_sum_site_pft, aes(x = year, ymin = abi025, ymax = abi975,
                                             group = taxon, color = taxon, fill = taxon), alpha = 0.4) +
    # geom_line(data = data_site, aes(x = year, y = ab)) +
    theme_bw() +
    labs(x = 'Year', y = 'Biomass Increment (Mg/ha)', color = 'Model', fill = 'Model',
         title = "Total Aboveground Biomass Increment") +
    facet_grid(.~model)
  print(pl3)
  ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_site_ABI.pdf'), width=10, height=8)
  # ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','processed_pft_site_ABI.jpg'))
  
  # }
  
  # figure to show cumulative biomass contribution across species and demonstrate which species are most important to site 
  prior_mat <- agb_taxa %>% 
    group_by(year, model, plot, taxon) %>%
    dplyr::summarize(ab = mean(ab), .groups = 'keep') %>%
    ungroup() %>% 
    group_by(year, model, taxon) %>% 
    dplyr::summarize(ab = mean(ab), .groups = 'keep') %>%
    ungroup() %>%
    group_by(model, taxon) %>%
    dplyr::summarize(contr = sum(ab), .groups = 'keep') %>%
    arrange(model, desc(contr))
  
  # first,  let's look at RW model 
  prior_mat1 = prior_mat %>% filter(model == 'Model RW')
  prior_mat1$perc = prior_mat1$contr/sum(prior_mat1$contr,na.rm=T)
  prior_mat1$cumsum = cumsum(prior_mat1$perc)
  prior_mat1$taxon = factor(x = prior_mat1$taxon, levels = prior_mat1$taxon)
  pl4 = ggplot(prior_mat1) +
    geom_point(aes(x=taxon, y=cumsum)) + 
    geom_hline(yintercept = 0.98, col = 'red') + 
    labs(title = 'overall cumulative proportion of biomass by species - Model RW', 
         y = 'cumulative proportion of biomass', 
         x = 'species')
  print(pl4)
  ggsave(pl4, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','98_biomass_modelRW.pdf'), width=10, height=8)
  # ggsave(pl4, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','98_biomass_modelRW.jpg'))
  
  # then RW + Census model if applicable 
  if (census_site){
    prior_mat2 = prior_mat %>% filter(model == 'Model RW + Census')
    prior_mat2$perc = prior_mat2$contr/sum(prior_mat2$contr,na.rm=T)
    prior_mat2$cumsum = cumsum(prior_mat2$perc)
    prior_mat2$taxon = factor(x = prior_mat2$taxon, levels = prior_mat2$taxon)
    pl5 = ggplot(prior_mat2) +
      geom_point(aes(x=taxon, y=cumsum)) + 
      geom_hline(yintercept = 0.98, col = 'red') + 
      labs(title = 'overall cumulative proportion of biomass by species - Model RW + Census', 
           y = 'cumulative proportion of biomass', 
           x = 'species')
    ggsave(pl5, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','98_biomass_modelRW_census.pdf'), width=10, height=8)
    # ggsave(pl5, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','98_biomass_modelRW_census.jpg'))
  }
}
