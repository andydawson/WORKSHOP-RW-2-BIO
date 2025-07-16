## Step 2: Run model 

## This script takes the built data from the last step and fits the STAN model, giving annual DBH estimates for each individual. 

# run_rw_model <- function(census_site, site, mvers, 
#                          dvers, keep = 300, nchains = 3, pool = 2500, iter = iter){
  run_rw_model <- function(census_site, site, mvers, 
                           dvers, keep = 300, nchains = 3, iter = iter){
  
  ##############################################
  ################ 1. Load data ################
  ##############################################
  
  #library(rstan)
  #library(gridExtra)
  #library(ggplotify)
  
  # Create directories if needed
  site_dir <- file.path('sites',site)
  
  # load built dataset 
  dat = readRDS(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input',paste0('tree_data_', site ,'_STAN_',mvers,'_', dvers, '.RDS')))
  
  ########################################################################
  ################ 2A. Run RW + census model if applicable ################
  ########################################################################
  
  if (census_site){
    
    # compile STAN model
    #compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_STAN.stan'))
    # compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_STAN.stan'))
    
    compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_species_sigxk_STAN.stan'))
    
    # fit and extract values
    fit <- sampling(compiled, 
                    data = dat, 
                    iter = iter, 
                    chains = nchains,
                    verbose=TRUE)
    # rm(compiled)
    
  ########################################################################
  ################ 3. RW + Census Model MCMC diagnostics  ################
  ########################################################################
    
    # diagnostic 1 : rhat
    mcmc_diags = list()
    trk_ind = 1
    rhat_val = rstan::stan_rhat(fit)
    mcmc_diags[[trk_ind]] = rhat_val
    trk_ind = trk_ind + 1
    rm(rhat_val)
    
    # diagnostic 2 : ess
    summary_val = summary(fit)
    ess = summary_val$summary[,9]
    
    # add bar plot for singular parameters
    inds = which(names(ess) %in% c('beta0','beta_sd','beta_t_sd','sig_x','sig_x_obs','sig_d_obs'))
    ess_singular = data.frame(names = names(ess)[inds], ess = ess[inds])
    pl1 = ggplot(ess_singular) + 
      geom_col(aes(x = names, y = ess)) + 
      geom_hline(yintercept = 100, color = 'red') + 
      labs(x = 'parameters', y = 'effective sample size', title = 'Effective Sample Size for Model Parameters')
    mcmc_diags[[trk_ind]] = pl1
    trk_ind = trk_ind + 1
    rm(pl1, ess_singular)
    
    # add histogram of ess for all beta_trees 
    inds = grep('beta\\[', names(ess))
    ess_betas = data.frame(ess = ess[inds])
    pl2 = ggplot(ess_betas) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for beta trees')
    mcmc_diags[[trk_ind]] = pl2
    trk_ind = trk_ind + 1
    rm(ess_betas, pl2)
    
    # add histogram of ess for all beta_years
    inds = grep('beta_t\\[', names(ess))
    ess_beta_ts = data.frame(ess = ess[inds])
    pl3 = ggplot(ess_beta_ts) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for beta years')
    mcmc_diags[[trk_ind]] = pl3
    trk_ind = trk_ind + 1
    rm(ess_beta_ts, pl3)
    
    # add histogram for increment estimates 
    inds = grep('X\\[', names(ess))
    ess_Xs = data.frame(ess = ess[inds])
    pl4 = ggplot(ess_Xs) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for increments')
    mcmc_diags[[trk_ind]] = pl4
    trk_ind = trk_ind + 1
    rm(ess_Xs, pl4)
    
    # add histogram for D0 values 
    inds = grep('D0\\[', names(ess))
    ess_D0s = data.frame(ess = ess[inds])
    pl5 = ggplot(ess_D0s) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for D0 values')
    mcmc_diags[[trk_ind]] = pl5
    trk_ind = trk_ind + 1
    rm(ess_D0s, pl5)
    
    # get organized values for MCMC chains so we can make trace plots and also save
    post=rstan::extract(fit, permuted = FALSE)
    variables = names(post[1,1,])
    # rm(fit)
    
    # diagnostic 3 : trace plots 
    # first, for the singulars 
    nchains = dim(post)[2]
    niter = dim(post)[1]
    for (par in c('beta0','beta_sd','beta_t_sd','sig_x_obs', 'sig_d_obs')){
      temp.ind = which(variables == par)
      temp.data = reshape2::melt(post[,,temp.ind])
      temp.data$iterations = seq(1, nrow(temp.data))
      if(nchains == 1){
        temp.data$chains = rep(1, nrow(temp.data))
      }
      pl.temp = ggplot(temp.data) + 
        geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
        labs(x = 'iteration', y = 'value', title = paste0('trace plot for ',par), 
             color = 'chain')
      print(pl.temp)
      mcmc_diags[[trk_ind]] = pl.temp
      trk_ind = trk_ind + 1
    }
    
    nchains = dim(post)[2]
    niter = dim(post)[1]
    for (k in 1:dat$N_taxa){
      temp.ind = which(variables == paste0('sig_x[',k,']'))
      temp.data = reshape2::melt(post[,,temp.ind])
      temp.data$iterations = seq(1, nrow(temp.data))
      if(nchains == 1){
        temp.data$chains = rep(1, nrow(temp.data))
      }
      pl.temp = ggplot(temp.data) + 
        geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
        labs(x = 'iteration', y = 'value', title = paste0('trace plot for sig_x ', k), 
             color = 'chain')
      print(pl.temp)
      mcmc_diags[[trk_ind]] = pl.temp
      trk_ind = trk_ind + 1
    }
    
    # then, for all beta trees 
    for (btr in 1:dat$N_Tr){
      temp.ind = which(variables == paste0('beta[',btr,']'))
      temp.data = reshape2::melt(post[,,temp.ind])
      temp.data$iterations = seq(1, nrow(temp.data))
      if (nchains==1){
        temp.data$chains = rep(1)
      }
      pl.temp = ggplot(temp.data) + 
        geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
        labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta ',btr), 
             color = 'chain')
      mcmc_diags[[trk_ind]] = pl.temp
      trk_ind = trk_ind + 1
    }
    
    # # lastly, for all beta years 
    # for (byr in 1:dat$N_years){
    #   for (tree in 1:dat$N_Tr){
    #   print(byr)
    #   temp.ind = which(substr(variables, 1, 7) == 'beta_t[')#which(variables == paste0('beta_t[', tree, ',' ,byr, ']'))
    #   temp.data = reshape2::melt(post[,,temp.ind])
    #   # temp.data$iterations = seq(1, nrow(temp.data))
    #   var.split = strsplit(as.vector(temp.data$parameters), "\\[|,|\\]")
    #   temp.data$year = as.numeric(lapply(var.split, function(x) x[[2]]))
    #   temp.data$taxon_id = as.numeric(lapply(var.split, function(x) x[[3]]))
    #   if (nchains==1){
    #     temp.data$chains = rep(1)
    #   }
    #   # if (nchains==1){
    #   #   pl.temp = ggplot(temp.data) + 
    #   #     geom_line(aes(x = iterations, y = value)) + 
    #   #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #   #          color = 'chain')
    #   # } else {
    #     pl.temp = ggplot(temp.data) + 
    #       geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
    #       labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #            color = 'chain')
    #   # }
    #   mcmc_diags[[trk_ind]] = pl.temp
    #   trk_ind = trk_ind + 1
    # }
    
    # lastly, for all beta years 
    # for (byr in 1:dat$N_years){
    #   print(byr)
    #   temp.ind = which(variables == paste0('beta_t[',byr,']'))
    #   temp.data = reshape2::melt(post[,,temp.ind])
    #   temp.data$iterations = seq(1, nrow(temp.data))
    #   if (nchains==1){
    #     temp.data$chains = rep(1)
    #   }
    #   # if (nchains==1){
    #   #   pl.temp = ggplot(temp.data) + 
    #   #     geom_line(aes(x = iterations, y = value)) + 
    #   #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #   #          color = 'chain')
    #   # } else {
    #   pl.temp = ggplot(temp.data) + 
    #     geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
    #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #          color = 'chain')
    #   # }
    #   mcmc_diags[[trk_ind]] = pl.temp
    #   trk_ind = trk_ind + 1
    # }
    
    pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','MCMC_diagnostics_RWC.pdf'), onefile = TRUE)
    for (i in seq(length(mcmc_diags))) {
      print(i)
      grid.arrange(mcmc_diags[[i]])
    }
    dev.off()
    rm(mcmc_diags)
    
    #################################################################
    ################ 4. Saving RW + Census RDS file  ################
    #################################################################
    
    # this is a big file, so let's just save the iterations we need 
    #postTemp = post[seq(dim(post)[1]-pool+1, dim(post)[1], pool/(keep/nchains)),,]
    #post = postTemp
    #rm(postTemp,ess)
    rm(ess)
    
    # save as RDS file 
    # saveRDS(post, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output', paste0('ring_model_t_pdbh_STAN_', site, '_', mvers, '_', dvers, '.RDS')))
    saveRDS(post, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output', paste0('ring_model_t_pdbh_species_sigxk_RWC_STAN_', site, '_', mvers, '_', dvers, '.RDS')))
    
    # generate a quick figure to roughly check model output  
    X2taxon_C = sapply(dat$X2C, function(id){dat$allTrees$taxon[which(dat$allTrees$stat_id == id)]})
    variables = names(post[1,1,])
    allDs = grep('D\\[',variables)
    D = matrix(NA,1,length(allDs))
    sig_d_obs = c()
    for (i in 1:nchains){
      D = rbind(D,post[,i,allDs])
      sig_d_obs = c(sig_d_obs, post[,i,which(variables == 'sig_d_obs')])
    }
    D = D[-1,]
    output = data.frame(D = apply(D, 2, mean), year = dat$X2year_C, tree = dat$X2C, taxon = X2taxon_C)
    
    
    # median is better measure of center here due to distribution skewness:
    sig_d_obs = median(sig_d_obs)
    # rm(post)
    
    saveRDS(sig_d_obs, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input','sig_d_obs_RWC.RDS'))
    
    
    # plot estimated diameter for all individuals over time 
    pl = ggplot(output) + 
      geom_line(aes(x = year, y = D, group = tree, color = tree)) +
      facet_wrap(~taxon) + 
      labs(x = 'Year', y = 'Diameter (cm)', title = 'Estimated Diameter over Time') + 
      theme(legend.position = 'none')
    ggsave(pl, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','estimated_species_growth_RWC.jpg'))
    
    rm(output, D, pl, variables, X2taxon_C, allDs)
  }
  
  ########################################################################
  ################ 2B. Run census model if applicable ####################
  ########################################################################
  
  if (census_site){
    
    compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_species_sigxk_CENSUS_STAN.stan'))
    
    # fit and extract values
    fit <- sampling(compiled, 
                    data = dat, 
                    iter = iter, 
                    chains = nchains,
                    verbose=TRUE)
    # rm(compiled)
    
    ########################################################################
    ################ 3. RW + Census Model MCMC diagnostics  ################
    ########################################################################
    
    # diagnostic 1 : rhat
    mcmc_diags = list()
    trk_ind = 1
    rhat_val = rstan::stan_rhat(fit)
    mcmc_diags[[trk_ind]] = rhat_val
    trk_ind = trk_ind + 1
    rm(rhat_val)
    
    # diagnostic 2 : ess
    summary_val = summary(fit)
    ess = summary_val$summary[,9]
    
    # add bar plot for singular parameters
    inds = which(names(ess) %in% c('beta0','beta_sd','beta_t_sd','sig_x','sig_x_obs','sig_d_obs'))
    ess_singular = data.frame(names = names(ess)[inds], ess = ess[inds])
    pl1 = ggplot(ess_singular) + 
      geom_col(aes(x = names, y = ess)) + 
      geom_hline(yintercept = 100, color = 'red') + 
      labs(x = 'parameters', y = 'effective sample size', title = 'Effective Sample Size for Model Parameters')
    mcmc_diags[[trk_ind]] = pl1
    trk_ind = trk_ind + 1
    rm(pl1, ess_singular)
    
    # add histogram of ess for all beta_trees 
    inds = grep('beta\\[', names(ess))
    ess_betas = data.frame(ess = ess[inds])
    pl2 = ggplot(ess_betas) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for beta trees')
    mcmc_diags[[trk_ind]] = pl2
    trk_ind = trk_ind + 1
    rm(ess_betas, pl2)
    
    # add histogram of ess for all beta_years
    inds = grep('beta_t\\[', names(ess))
    ess_beta_ts = data.frame(ess = ess[inds])
    pl3 = ggplot(ess_beta_ts) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for beta years')
    mcmc_diags[[trk_ind]] = pl3
    trk_ind = trk_ind + 1
    rm(ess_beta_ts, pl3)
    
    # add histogram for increment estimates 
    inds = grep('X\\[', names(ess))
    ess_Xs = data.frame(ess = ess[inds])
    pl4 = ggplot(ess_Xs) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for increments')
    mcmc_diags[[trk_ind]] = pl4
    trk_ind = trk_ind + 1
    rm(ess_Xs, pl4)
    
    # add histogram for D0 values 
    inds = grep('D0\\[', names(ess))
    ess_D0s = data.frame(ess = ess[inds])
    pl5 = ggplot(ess_D0s) + 
      geom_histogram(aes(x = ess)) + 
      geom_vline(xintercept = 100, color = 'red') + 
      labs(x = 'effective sample size', title = 'effective sample size for D0 values')
    mcmc_diags[[trk_ind]] = pl5
    trk_ind = trk_ind + 1
    rm(ess_D0s, pl5)
    
    # get organized values for MCMC chains so we can make trace plots and also save
    post=rstan::extract(fit, permuted = FALSE)
    variables = names(post[1,1,])
    # rm(fit)
    
    # diagnostic 3 : trace plots 
    # first, for the singulars 
    nchains = dim(post)[2]
    niter = dim(post)[1]
    for (par in c('beta0','beta_sd','beta_t_sd', 'sig_d_obs')){
    # for (par in c('beta0','beta_sd','beta_t_sd','sig_x_obs', 'sig_d_obs')){
      print(par)
      temp.ind = which(variables == par)
      temp.data = reshape2::melt(post[,,temp.ind])
      temp.data$iterations = seq(1, nrow(temp.data))
      if(nchains == 1){
        temp.data$chains = rep(1, nrow(temp.data))
      }
      pl.temp = ggplot(temp.data) + 
        geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
        labs(x = 'iteration', y = 'value', title = paste0('trace plot for ',par), 
             color = 'chain')
      print(pl.temp)
      mcmc_diags[[trk_ind]] = pl.temp
      trk_ind = trk_ind + 1
    }
    
    nchains = dim(post)[2]
    niter = dim(post)[1]
    for (k in 1:dat$N_taxa){
      temp.ind = which(variables == paste0('sig_x[',k,']'))
      temp.data = reshape2::melt(post[,,temp.ind])
      temp.data$iterations = seq(1, nrow(temp.data))
      if(nchains == 1){
        temp.data$chains = rep(1, nrow(temp.data))
      }
      pl.temp = ggplot(temp.data) + 
        geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
        labs(x = 'iteration', y = 'value', title = paste0('trace plot for sig_x ', k), 
             color = 'chain')
      print(pl.temp)
      mcmc_diags[[trk_ind]] = pl.temp
      trk_ind = trk_ind + 1
    }
    
    # then, for all beta trees 
    for (btr in 1:dat$N_Tr){
      temp.ind = which(variables == paste0('beta[',btr,']'))
      temp.data = reshape2::melt(post[,,temp.ind])
      temp.data$iterations = seq(1, nrow(temp.data))
      if (nchains==1){
        temp.data$chains = rep(1)
      }
      pl.temp = ggplot(temp.data) + 
        geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
        labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta ',btr), 
             color = 'chain')
      mcmc_diags[[trk_ind]] = pl.temp
      trk_ind = trk_ind + 1
    }
    
    # # lastly, for all beta years 
    # for (byr in 1:dat$N_years){
    #   for (tree in 1:dat$N_Tr){
    #   print(byr)
    #   temp.ind = which(substr(variables, 1, 7) == 'beta_t[')#which(variables == paste0('beta_t[', tree, ',' ,byr, ']'))
    #   temp.data = reshape2::melt(post[,,temp.ind])
    #   # temp.data$iterations = seq(1, nrow(temp.data))
    #   var.split = strsplit(as.vector(temp.data$parameters), "\\[|,|\\]")
    #   temp.data$year = as.numeric(lapply(var.split, function(x) x[[2]]))
    #   temp.data$taxon_id = as.numeric(lapply(var.split, function(x) x[[3]]))
    #   if (nchains==1){
    #     temp.data$chains = rep(1)
    #   }
    #   # if (nchains==1){
    #   #   pl.temp = ggplot(temp.data) + 
    #   #     geom_line(aes(x = iterations, y = value)) + 
    #   #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #   #          color = 'chain')
    #   # } else {
    #     pl.temp = ggplot(temp.data) + 
    #       geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
    #       labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #            color = 'chain')
    #   # }
    #   mcmc_diags[[trk_ind]] = pl.temp
    #   trk_ind = trk_ind + 1
    # }
    
    # lastly, for all beta years 
    # for (byr in 1:dat$N_years){
    #   print(byr)
    #   temp.ind = which(variables == paste0('beta_t[',byr,']'))
    #   temp.data = reshape2::melt(post[,,temp.ind])
    #   temp.data$iterations = seq(1, nrow(temp.data))
    #   if (nchains==1){
    #     temp.data$chains = rep(1)
    #   }
    #   # if (nchains==1){
    #   #   pl.temp = ggplot(temp.data) + 
    #   #     geom_line(aes(x = iterations, y = value)) + 
    #   #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #   #          color = 'chain')
    #   # } else {
    #   pl.temp = ggplot(temp.data) + 
    #     geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
    #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
    #          color = 'chain')
    #   # }
    #   mcmc_diags[[trk_ind]] = pl.temp
    #   trk_ind = trk_ind + 1
    # }
    
    pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','MCMC_diagnostics_C.pdf'), onefile = TRUE)
    for (i in seq(length(mcmc_diags))) {
      print(i)
      grid.arrange(mcmc_diags[[i]])
    }
    dev.off()
    rm(mcmc_diags)
    
    #################################################################
    ################ 4. Saving RW + Census RDS file  ################
    #################################################################
    
    # this is a big file, so let's just save the iterations we need 
    #postTemp = post[seq(dim(post)[1]-pool+1, dim(post)[1], pool/(keep/nchains)),,]
    #post = postTemp
    #rm(postTemp,ess)
    rm(ess)
    
    # save as RDS file 
    # saveRDS(post, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output', paste0('ring_model_t_pdbh_STAN_', site, '_', mvers, '_', dvers, '.RDS')))
    saveRDS(post, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output', paste0('ring_model_t_pdbh_species_sigxk_C_STAN_', site, '_', mvers, '_', dvers, '.RDS')))
    
    # generate a quick figure to roughly check model output  
    X2taxon_C = sapply(dat$X2C, function(id){dat$allTrees$taxon[which(dat$allTrees$stat_id == id)]})
    variables = names(post[1,1,])
    allDs = grep('D\\[',variables)
    D = matrix(NA,1,length(allDs))
    sig_d_obs = c()
    for (i in 1:nchains){
      D = rbind(D,post[,i,allDs])
      sig_d_obs = c(sig_d_obs, post[,i,which(variables == 'sig_d_obs')])
    }
    D = D[-1,]
    output = data.frame(D = apply(D, 2, mean), year = dat$X2year_C, tree = dat$X2C, taxon = X2taxon_C)
    
    
    # median is better measure of center here due to distribution skewness:
    sig_d_obs = median(sig_d_obs)
    # rm(post)
    
    saveRDS(sig_d_obs, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input','sig_d_obs_C.RDS'))
    
    
    # plot estimated diameter for all individuals over time 
    pl = ggplot(output) + 
      geom_line(aes(x = year, y = D, group = tree, color = tree)) +
      facet_wrap(~taxon) + 
      labs(x = 'Year', y = 'Diameter (cm)', title = 'Estimated Diameter over Time') + 
      theme(legend.position = 'none')
    ggsave(pl, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','estimated_species_growth_C.jpg'))
    
    rm(output, D, pl, variables, X2taxon_C, allDs)
  }
  

  ######################################################
  ################ 5. Run RW only model ################
  ######################################################
  
  # if we do not have census data for this site, we have to obtain mean measurement error for DBH 
  # from a site that has both of these datasets (i.e. Harvard Forest)
  # TO DO: this needs to be automated I think so that we can adjust which site we want to use
  if (!census_site){
    # out = readRDS('sites/HARVARD/runs/v2.0_102020/output/ring_model_t_pdbh_STAN_HARVARD_v2.0_102020.RDS')
    # 
    # col_names = sapply(strsplit(colnames(out[,1,]), '\\['), function(x) x[[1]])
    # hist(out[,1,which(col_names=="sig_d_obs")])
    # sig_d_obs_chains = apply(out, 2, function(x) x[, which(col_names=="sig_d_obs")])
    # sig_d_obs = mean(sig_d_obs_chains)
    
    out = readRDS('sites/HARVARD/runs/v3.1_102020/output/ring_model_t_pdbh_species_sigxk_RWC_STAN_HARVARD_v3.1_102020.RDS')
    
    col_names = sapply(strsplit(colnames(out[,1,]), '\\['), function(x) x[[1]])
    hist(out[,1,which(col_names=="sig_d_obs")])
    sig_d_obs_chains = apply(out, 2, function(x) x[, which(col_names=="sig_d_obs")])
    sig_d_obs = mean(sig_d_obs_chains)
    
    # col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
    # hist(out[,which(col_names=="sig_d_obs")])
    # sig_d_obs = mean(out[,which(col_names=="sig_d_obs")])
    

    
    # # extract measurement error from dataset and find median (better measure of center due to skewness of distribution)
    # needed = which(names(out[1,1,]) == 'sig_d_obs')
    # all = c(out[,1,needed], out[,2,needed], out[,3,needed])
    # dat$sig_d_obs = median(all)
    dat$sig_d_obs = sig_d_obs
    # print(needed)
    # print(all) 
    rm(out)
    
  # otherwise we use the value found in the model above 
  }else{
    # dat$sig_d_obs = sig_d_obs
    sig_d_obs = readRDS(file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input','sig_d_obs_RWC.RDS'))
    dat$sig_d_obs = sig_d_obs
  }

  # compile STAN model
  # compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_sigd_STAN.stan'))
  # compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_sigd_species_STAN.stan'))
  compiled <- stan_model(file = paste0('models/ring_model_t_pdbh_sigd_species_sigxk_STAN.stan'))
  

  # fit and extract values
  fit <- sampling(compiled, 
                  data = dat, 
                  iter = iter, 
                  chains = nchains,
                  verbose = TRUE)
  # rm(compiled)
  
  ##############################################################
  ################ 6. RW only MCMC diagnostics  ################
  ##############################################################
  
  # diagnostic 1 : rhat
  mcmc_diags = list()
  trk_ind = 1
  rhat_val = rstan::stan_rhat(fit)
  mcmc_diags[[trk_ind]] = rhat_val
  trk_ind = trk_ind + 1
  rm(rhat_val)
  
  # diagnostic 2 : ess
  summary_val = summary(fit)
  ess = summary_val$summary[,9]
  
  # add bar plot for singular parameters
  inds = which(names(ess) %in% c('beta0','beta_sd','beta_t_sd','sig_x','sig_x_obs'))
  ess_singular = data.frame(names = names(ess)[inds], ess = ess[inds])
  pl1 = ggplot(ess_singular) + 
    geom_col(aes(x = names, y = ess)) + 
    geom_hline(yintercept = 100, color = 'red') + 
    labs(x = 'parameters', y = 'effective sample size', title = 'Effective Sample Size for Model Parameters')
  mcmc_diags[[trk_ind]] = pl1
  trk_ind = trk_ind + 1
  rm(pl1, ess_singular)
  
  # add histogram of ess for all beta_trees 
  inds = grep('beta\\[', names(ess))
  ess_betas = data.frame(ess = ess[inds])
  pl2 = ggplot(ess_betas) + 
    geom_histogram(aes(x = ess)) + 
    geom_vline(xintercept = 100, color = 'red') + 
    labs(x = 'effective sample size', title = 'effective sample size for beta trees')
  mcmc_diags[[trk_ind]] = pl2
  trk_ind = trk_ind + 1
  rm(ess_betas, pl2)
  
  # add histogram of ess for all beta_years
  inds = grep('beta_t\\[', names(ess))
  ess_beta_ts = data.frame(ess = ess[inds])
  pl3 = ggplot(ess_beta_ts) + 
    geom_histogram(aes(x = ess)) + 
    geom_vline(xintercept = 100, color = 'red') + 
    labs(x = 'effective sample size', title = 'effective sample size for beta years')
  mcmc_diags[[trk_ind]] = pl3
  trk_ind = trk_ind + 1
  rm(ess_beta_ts, pl3)
  
  # add histogram for increment estimates 
  inds = grep('X\\[', names(ess))
  ess_Xs = data.frame(ess = ess[inds])
  pl4 = ggplot(ess_Xs) + 
    geom_histogram(aes(x = ess)) + 
    geom_vline(xintercept = 100, color = 'red') + 
    labs(x = 'effective sample size', title = 'effective sample size for increments')
  mcmc_diags[[trk_ind]] = pl4
  trk_ind = trk_ind + 1
  rm(ess_Xs, pl4)
  
  # add histogram for D0 values 
  inds = grep('D0\\[', names(ess))
  ess_D0s = data.frame(ess = ess[inds])
  pl5 = ggplot(ess_D0s) + 
    geom_histogram(aes(x = ess)) + 
    geom_vline(xintercept = 100, color = 'red') + 
    labs(x = 'effective sample size', title = 'effective sample size for D0 values')
  mcmc_diags[[trk_ind]] = pl5
  trk_ind = trk_ind + 1
  rm(ess_D0s, pl5)

  # get organized values for MCMC chains so we can make trace plots and also save
  post=rstan::extract(fit, permuted = FALSE)
  variables = names(post[1,1,])
  # rm(fit)
  
  # diagnostic 3 : trace plots 
  # first, for the singulars 
  # nchains = dim(post)[2]
  # niter = dim(post)[1]
  # for (par in c('beta0','beta_sd','beta_t_sd','sig_x','sig_x_obs')){
  #   temp.ind = which(variables == par)
  #   temp.data = reshape2::melt(post[,,temp.ind])
  #   temp.data$iterations = seq(1, nrow(temp.data))
  #   if(nchains == 1){
  #     temp.data$chains = rep(1, nrow(temp.data))
  #   }
  #   pl.temp = ggplot(temp.data) + 
  #     geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
  #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for ',par), 
  #          color = 'chain')
  #   print(pl.temp)
  #   mcmc_diags[[trk_ind]] = pl.temp
  #   trk_ind = trk_ind + 1
  # }
  
  nchains = dim(post)[2]
  niter = dim(post)[1]
  for (par in c('beta0','beta_sd','beta_t_sd','sig_x_obs')){
    temp.ind = which(variables == par)
    temp.data = reshape2::melt(post[,,temp.ind])
    temp.data$iterations = seq(1, nrow(temp.data))
    if(nchains == 1){
      temp.data$chains = rep(1, nrow(temp.data))
    }
    pl.temp = ggplot(temp.data) + 
      geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
      labs(x = 'iteration', y = 'value', title = paste0('trace plot for ',par), 
           color = 'chain')
    print(pl.temp)
    mcmc_diags[[trk_ind]] = pl.temp
    trk_ind = trk_ind + 1
  }
  
  nchains = dim(post)[2]
  niter = dim(post)[1]
  for (k in 1:dat$N_taxa){
    temp.ind = which(variables == paste0('sig_x[',k,']'))
    temp.data = reshape2::melt(post[,,temp.ind])
    temp.data$iterations = seq(1, nrow(temp.data))
    if(nchains == 1){
      temp.data$chains = rep(1, nrow(temp.data))
    }
    pl.temp = ggplot(temp.data) + 
      geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
      labs(x = 'iteration', y = 'value', title = paste0('trace plot for sig_x ', k), 
           color = 'chain')
    print(pl.temp)
    mcmc_diags[[trk_ind]] = pl.temp
    trk_ind = trk_ind + 1
  }
  
  # then, for all beta trees 
  for (btr in 1:dat$N_Tr){
    print(btr)
    temp.ind = which(variables == paste0('beta[',btr,']'))
    temp.data = reshape2::melt(post[,,temp.ind])
    temp.data$iterations = seq(1, nrow(temp.data))
    if(nchains == 1){
      temp.data$chains = rep(1, nrow(temp.data))
    }
    pl.temp = ggplot(temp.data) + 
      geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
      labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta ',btr), 
           color = 'chain')
    print(pl.temp)
    mcmc_diags[[trk_ind]] = pl.temp
    trk_ind = trk_ind + 1
  }
  
  
  # # lastly, for all beta years 
  # for (byr in 1:dat$N_years){
  #   for (tree in 1:dat$N_Tr){
  #   print(byr)
  #   temp.ind = which(substr(variables, 1, 7) == 'beta_t[')#which(variables == paste0('beta_t[', tree, ',' ,byr, ']'))
  #   temp.data = reshape2::melt(post[,,temp.ind])
  #   # temp.data$iterations = seq(1, nrow(temp.data))
  #   var.split = strsplit(as.vector(temp.data$parameters), "\\[|,|\\]")
  #   temp.data$year = as.numeric(lapply(var.split, function(x) x[[2]]))
  #   temp.data$taxon_id = as.numeric(lapply(var.split, function(x) x[[3]]))
  #   if (nchains==1){
  #     temp.data$chains = rep(1)
  #   }
  #   # if (nchains==1){
  #   #   pl.temp = ggplot(temp.data) + 
  #   #     geom_line(aes(x = iterations, y = value)) + 
  #   #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
  #   #          color = 'chain')
  #   # } else {
  #     pl.temp = ggplot(temp.data) + 
  #       geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
  #       labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
  #            color = 'chain')
  #   # }
  #   mcmc_diags[[trk_ind]] = pl.temp
  #   trk_ind = trk_ind + 1
  # }
  
  # # lastly, for all beta years 
  # for (byr in 1:dat$N_years){
  #   print(byr)
  #   temp.ind = which(variables == paste0('beta_t[',byr,']'))
  #   temp.data = reshape2::melt(post[,,temp.ind])
  #   temp.data$iterations = seq(1, nrow(temp.data))
  #   if(nchains == 1){
  #     temp.data$chains = rep(1, nrow(temp.data))
  #   }
  #   pl.temp = ggplot(temp.data) + 
  #     geom_line(aes(x = iterations, y = value, group = as.factor(chains), color = as.factor(chains))) + 
  #     labs(x = 'iteration', y = 'value', title = paste0('trace plot for beta year ',byr), 
  #          color = 'chain')
  #   print(pl.temp)
  #   mcmc_diags[[trk_ind]] = pl.temp
  #   trk_ind = trk_ind + 1
  # }
  # 
  pdf(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','MCMC_diagnostics_RW.pdf'), onefile = TRUE)
  for (i in seq(length(mcmc_diags))) {
    print(i)
    grid.arrange(mcmc_diags[[i]])
  }
  dev.off()
  rm(mcmc_diags)
  
  #############################################################
  ################ 7. Saving RW only RDS file  ################
  #############################################################
  
  # this is a big file, so let's just save the iterations we need 
  #postTemp = post[seq(dim(post)[1]-pool+1, dim(post)[1], pool/(keep/nchains)),,]
  #post = postTemp
  #rm(postTemp)
  
  # save as RDS file 
  # saveRDS(post, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output',paste0('ring_model_t_pdbh_sigd_species_STAN_', site, '_', mvers,'_', dvers, '.RDS')))
  saveRDS(post, file = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output',paste0('ring_model_t_pdbh_sigd_species_sigxk_RW_STAN_', site, '_', mvers,'_', dvers, '.RDS')))
  
  # generate a quick figure to roughly check model output 
  allDs = grep('D\\[',variables)
  D = matrix(NA,1,length(allDs))
  for (i in 1:nchains){
    D = rbind(D,post[,i,allDs])
  }
  D = D[-1,]
  output = data.frame(D = apply(D, 2, mean), year = dat$X2year, tree = dat$X2Tr, taxon = dat$Tr$taxon[dat$X2Tr])
  # rm(post) 
  
  # plot estimated diameter for all individuals over time
  pl = ggplot(output) + 
    geom_line(aes(x = year, y = D, group = tree, color = tree)) +
    facet_wrap(~taxon) + 
    labs(x = 'Year', y = 'Diameter (cm)', title = 'Estimated Diameter over Time') + 
    theme(legend.position = 'none')
  ggsave(pl, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','estimated_species_growth_RW.jpg'))
}
