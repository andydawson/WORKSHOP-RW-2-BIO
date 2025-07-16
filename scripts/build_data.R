## Step 1: Build and format ring width data 

## This script takes the ring width files (rwl format) and the stand-level inputs from the year the trees were cored and formats the data
## to run in the STAN ring width model. In getting to the correct format, this script also matches cores with tree species, 
## DBH, and location within the plot. The CSV file should be in the same format as that of the sample CSV file available at (INSERT REPO URL). 

build_data <- function(site, dvers, mvers, prefix, 
                       census_site, cutoff = 1900){
   
  # Prepare workspace 
  # Had to comment these out for submitting on ND campus cluster (stupid issues with loading packages)
  #library(plotrix)
  #library(dplR)
  #library(fields)
  #library(reshape2)
  #library(dplyr)
  #library(plyr)
  #library(ggplot2)
  
  # Create save folders for data 
  site_dir <- file.path('sites',site)
  if (!file.exists(file.path(site_dir,'runs')))   dir.create(file.path(site_dir,'runs'))
  if (!file.exists(file.path(site_dir,'runs',paste0(mvers,'_',dvers))))   dir.create(file.path(site_dir,'runs',paste0(mvers,'_',dvers)))
  if (!file.exists(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output')))   dir.create(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'output'))
  if (!file.exists(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input')))   dir.create(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input'))
  if (!file.exists(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures')))   dir.create(file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures'))
  
  rwl_dir = file.path('sites',site,'data','raw','rwl')
  meta_loc = file.path('sites',site,'data','raw',paste0(site,'_treeMeta_',dvers,'.csv'))
  RDS_loc = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'input',paste0('tree_data_', site ,'_STAN_',mvers,'_', dvers, '.RDS'))
  
  if (census_site){
    census_loc = file.path('sites',site,'data','raw',paste0(site,'_census_',dvers,'.csv'))
  }
  
  #################################################
  ################ 1. Extract data ################
  #################################################
  
  # Load CSV with tree data
  treeMeta = read.csv(meta_loc, stringsAsFactors = FALSE)
  treeMeta$species[which(treeMeta$species == 'ASCA')] = 'ACSA'
  if (site == 'HMC'){
    treeMeta$site = treeMeta$plot
    # treeMeta$ID = treeMeta$tree_number
    # treeMeta$ID[which(nchar(treeMeta$ID) == 1)] = paste0('00',  treeMeta$ID[which(nchar(treeMeta$ID) == 1)])
    # treeMeta$ID[which(nchar(treeMeta$ID) == 2)] = paste0('0',  treeMeta$ID[which(nchar(treeMeta$ID) == 2)])
  }
  
  len = nchar(treeMeta$ID[1])
  
  # List of files with tree ring increments, extract only RWL files 
  rwFiles <- list.files(rwl_dir)
  rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
  rwData <- list()
  for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    # Insert the contents of each file into the rwData list
    rwData[[id]] <- t((read.tucson(file.path(rwl_dir, fn))))  # rows are tree, cols are times
  }
  
  # Bind into one dataframe and insert core IDs into the matrix 
  # Dimensions: core x year
  incr = ldply(rwData, rbind)
  incr = incr[,c(".id", sort(colnames(incr)[2:ncol(incr)]))]
  rownames(incr) = as.vector(unlist(lapply(rwData, rownames)))
  incr[,1] = rownames(incr)
  
  incr_missing = incr
  incr_missing_data = melt(incr_missing)
  colnames(incr_missing_data) = c('id', 'year', 'incr')
  incr_missing_data$year = as.vector(incr_missing_data$year)
  # Assign core ID to each core
  incr_missing_data$id = substr(incr_missing_data$id, 1, len)
  
  # foo = incr_missing_data %>%
  #   group_by(id, year) %>%
  #   mutate(any_missing = any(incr==0, na.rm=TRUE),
  #          not_missing = any(incr>0, na.rm=TRUE))
  
  foo = incr_missing_data %>%
    group_by(id, year) %>%
    dplyr::summarize(any_missing = any(incr==0, na.rm=TRUE),
           not_missing = any(incr>0, na.rm=TRUE))
  
  head(foo)
  
  sum(foo$not_missing)
  sum(foo$any_missing)
  
  bar = foo[which(foo$any_missing),]
  bar2 = foo[which((foo$any_missing)),]
  
  incr = replace(incr, incr==0, NA)
  
  # Create incr_data data frame 
  incr_data = melt(incr)
  colnames(incr_data) = c('id', 'year', 'incr')
  incr_data$year = as.vector(incr_data$year)
  
  # Remove NA values right away to quicken processing (only keep the series for which trees were alive/existed)
  # Also remove all increment values for before the determined cut off year, which has a default value of 1900
  # Also remove all data that comes after the first year of full data (some sites have multiple years of coring)
  incr_data = incr_data %>% 
    mutate(year = as.numeric(year)) %>% 
    filter(!is.na(incr), year >= cutoff)
  
  # Assign core ID to each core
  incr_data$id = substr(incr_data$id, 1, len)
  
  # Assign plot to each core in data frame (if only one plot, all 1)
  if (!('site' %in% colnames(treeMeta))){
    incr_data$plot   = rep(1, nrow(incr_data))
  }else{
    incr_data$plot = sapply(1:length(incr_data$id), 
                            function(ind){
                              choose = treeMeta$site[which(treeMeta$ID == incr_data$id[ind])]
                              return(as.numeric(substr(choose, nchar(prefix)+1, nchar(choose))))
                            })
  }
  
  # Assign stat IDs
  #for HMC, there are some sampled trees without RWs
  treeMeta[which(!(treeMeta$ID %in% unique(incr_data$id))),]
  
  
  treeMeta = treeMeta %>% filter(ID %in% unique(incr_data$id))
  treeMeta$stat_id = seq(1,nrow(treeMeta))
  incr_data$stat_id = as.numeric(plyr::mapvalues(incr_data$id, 
                                                   from = treeMeta$ID,
                                                   to = treeMeta$stat_id,
                                                   warn_missing = FALSE))
  
  #########################################################
  ################ 2. Organize census data ################
  #########################################################
  
  # Load and organize census data if census available
  if (census_site){
    census = read.csv(census_loc, stringsAsFactors = FALSE) %>% 
      filter(species %in% unique(treeMeta$species))
    
    if (site == 'HMC'){
      census$site = as.numeric(substr(census$plot, 4, 4))
      census$finalCond = NA
      census$ID = census$tree_number
    }
    # if (site == 'HARVARD'){
    if (site == 'HMC'){
      census_long = census[, c('site', 'ID', 'species', 'distance','finalCond', 'year', 'dbh')]
      census_long = census_long %>% 
        filter(!is.na(dbh))
    } else {
      census_long = melt(census, 
                         id=c('site', 'ID', 'species', 'distance', 'finalCond'))
      colnames(census_long) = c('site', 'ID', 'species', 'distance','finalCond', 'year', 'dbh')
      # remove NA years and reformat census year column
      census_long = census_long %>% 
        filter(!is.na(dbh)) %>%
        mutate(year = substr(year, 2, 3)) %>%
        mutate(year = ifelse(as.numeric(year) < 30, paste0('20',year), paste0('19',year)))
    }
 # }
    

    
    # get census years
    census_years = as.numeric(sort(unique(census_long$year)))
    
    # create a list of all trees that are present in available data
    allTrees = full_join(census, treeMeta, by = c('ID')) %>% 
      mutate(taxon = ifelse(is.na(species.x), species.y, species.x),
             distance = ifelse(is.na(distance.y), distance.x, distance.y), 
             plot = as.numeric(substr(ID,nchar(prefix)+1,nchar(prefix)+1))) %>% 
      select(ID, stat_id, taxon, distance, plot)
    allTrees$stat_id[which(is.na(allTrees$stat_id))] = seq(nrow(treeMeta)+1, nrow(treeMeta)+length(which(is.na(allTrees$stat_id))))
    
    census_long$stat_id = as.numeric(plyr::mapvalues(census_long$ID, from = allTrees$ID, to = allTrees$stat_id, warn_missing = FALSE))
  }
  
  ############################################################
  ################ 3. Organize increment data ################
  ############################################################
  
  # Assign species to each observation
  incr_data$taxon = sapply(c(1:length(incr_data$id)), 
                           function(ind){
                             treeMeta$species[which(treeMeta$ID == incr_data$id[ind])]
                           })
  
  # Find year range for data
  year_end = max(as.numeric(incr_data$year), na.rm=TRUE)
  year_start = min(as.numeric(incr_data$year), na.rm=TRUE)
  years = seq(year_start, year_end)
  
  # Make units of time in "years since first recorded tree ring"
  incr_data$year = as.numeric(incr_data$year) - year_start + 1
  
  # Order by tree and year 
  incr_data = incr_data %>% arrange(stat_id, year)
  
  # Data frame of the first and last year of observation for each tree
  # First year for all years is the first year of data (which will tend to be the cutoff year)
  year_idx = data.frame(stat_id = as.numeric(aggregate(year~stat_id, data=incr_data, FUN=min, na.rm=TRUE)[,1]),
                        year_end=as.numeric(aggregate(year~stat_id, incr_data, max)[,2]))
  year_idx$year_start = rep(1, length(year_idx$stat_id))
  
  
  taxa = unique(incr_data$taxon)
  N_taxa = length(taxa)
  # X2taxon = Tr2taxon[X2Tr]
  
  #########################################################
  ################ 4. Organize RW DBH data ################
  #########################################################
  
  # create dataframe that gives year index corresponding to RW DBH measurement 
  pdbh = aggregate(year~stat_id+plot+id, incr_data, max, na.rm=TRUE) %>% arrange(stat_id)
  pdbh$dbh = as.numeric(plyr::mapvalues(pdbh$id, from = treeMeta$ID, to = treeMeta$dbh,
                                        warn_missing = FALSE))
  pdbh$distance = as.numeric(plyr::mapvalues(pdbh$id, from = treeMeta$ID, to = treeMeta$distance,
                                             warn_missing = FALSE))
  pdbh$taxon = plyr::mapvalues(pdbh$id, from = treeMeta$ID, to = treeMeta$species, 
                               warn_missing = FALSE)
  
  Tr2taxon = match(pdbh$taxon, taxa)
  #####################################################################
  ################ 5. Organize RW only Model estimates ################
  #####################################################################
  
  # create dataframe that would hold all RW values we need to estimate based only on RW data
  # AKA loop through each tree and create a row of info for each RW value we need to estimate
  X_data = data.frame(meas=numeric(0), stat_id=numeric(0), year=numeric(0))
  n = 1
  for (tree in 1:length(year_idx$stat_id)){
    stat = year_idx$stat_id[tree]
    year = seq(year_idx$year_start[tree], year_idx$year_end[tree])
    meas = seq(n, n+length(year)-1)
    n = n + length(year)
    X_data = rbind(X_data, data.frame(meas=meas, stat_id=rep(stat, length(year)), year=year))
  }
  
  #########################################################################
  ################ 6. Organize RW + Census Model estimates ################
  #########################################################################
  
  if (census_site){
    
    # Need to adjust dates of census data 
    census_long$year = as.numeric(census_long$year) - year_start + 1
    
    # create dataframe that would hold all RW values we need to estimate based on both RW + Census data 
    # AKA loop through all trees and create a row of info for each RW value we need to estimate
    X_data_C = data.frame(meas=numeric(0), 
                          stat_id=numeric(0), 
                          species_id = numeric(0), 
                          year=numeric(0))
    n = 1
    for (tree in 1:length(allTrees$stat_id)){
      print(tree)
      stat = allTrees$stat_id[tree]
      
      if (stat==23){
        print(tree)
      }
      
      # check to see if available in RW data
      in.RW = ifelse(stat %in% year_idx$stat_id, TRUE, FALSE)
      
      # check to see if in census 
      in.C = ifelse(stat %in% census_long$stat_id, TRUE, FALSE)
      
      # if in both, determine latest date with available data 
      # if RW ends before the last census measurement, we need to enable smoothing 
      # last year will be either:
      # (1) year before first census without that tree if it is not in all censuses or
      # (2) last year with recorded data if after last census 
      if (in.RW & in.C){
        
        # first year is always first year of available data (tends to be cutoff)
        firstyr = 1
        
        # what was the last census this tree was in?
        lastCensus = max(which(which(years %in% census_years) %in% 
                                 as.numeric(census_long$year[which(census_long$stat_id == stat)])))
        lastData = year_idx$year_end[year_idx$stat_id == stat]
        
        # if RW ends after census, then we just use last RW year
        if ((which(years %in% census_years)[lastCensus]) < lastData){
          lastyr = lastData
          
        # we will stochastically pick death year for trees in processing, always run up until first year of coring 
        # which in this case will just be the last year of data 
        }else{
          # if in last census, use last year of data (deal with in processing)
          if(lastCensus == length(census_years)){
            lastyr = length(years)
          # if not in last census, use year before first census where the tree is missing and mark for smoothing 
          }else{
            lastyr = which(years == (census_years[lastCensus+1] - 1))
            census_long$finalCond[which(census_long$stat_id == stat)] = 'dead'
          }
        }
        
        species = census_long$species[which(census_long$stat_id == stat)][1]
        species_id = match(species, taxa)
      }
      
      # if only in RW, use year_idx 
      if (in.RW & !in.C){
        firstyr = 1
        lastyr = year_idx$year_end[year_idx$stat_id == stat]
        # if (stat == 23){
        #   species_id = match('BEAL', taxa)
        # } else if (stat == 31) {
        #   
        # } else{
          species = pdbh$taxon[match(stat, pdbh$stat_id)]
          species_id = match(species, taxa)
        # }
      }
      
      # if only in census, last year is either:
        # (1) year before first census date without that tree if it is not in all censuses or 
        # (2) final year if in last census (will use firstCond later to determine if smoothing needed) 
      # this allows us to perform mortality smoothing in processing the model output
      if (!in.RW & in.C){
        
        # first year is going to be the first year with any data (tends to be cutoff)
        firstyr = 1
        
        # determine the last year with census data 
        lastCensus = max(which(which(years %in% census_years) %in% unique(census_long$year[which(census_long$stat_id == stat)])))
        
        # if tree was in most recent census
        if (lastCensus == length(census_years)){
          lastyr = length(years)
        # otherwise, run it until the year before the next census and mark for smoothing
        }else{
          lastyr = which(years == (census_years[lastCensus+1] - 1))
          census_long$finalCond[which(census_long$stat_id == stat)] = 'dead'
        }
        
        species = census_long$species[which(census_long$stat_id == stat)][1]
        species_id = match(species, taxa)
      }
      
      if (!in.RW & !in.C){
        print(paste(tree, 'is missing from both datasets'))
        next
      }
      
      # add those estimates needed for this stat ID 
      year = seq(firstyr, lastyr)
      meas = seq(n, n+length(year)-1)
      n = n + length(year)
      X_data_C = rbind(X_data_C, 
                       data.frame(meas=meas, stat_id=rep(stat, length(year)), 
                                  species_id =rep(species_id, length(year)), year=year))
    }
  }
  
  ######################################################################################
  ################ 7. Save variables as needed for model and processing ################
  ######################################################################################
  
  # if there are stat IDs without data, we just need to make sure they are removed
  pdbh = pdbh %>% filter(stat_id %in% unique(X_data$stat_id))
  
  # obtain required values for rw model and processing 
  # Number of measurements
  N_Xobs = nrow(incr_data)
  # Number of years where we have RW data
  N_years = length(years)
  # Number of trees with RW data 
  N_Tr = nrow(pdbh)
  # Number of values to estimate with RW only model
  N_X = nrow(X_data)

  if (census_site){
    
    # Number of trees that will need to be estimated
    N_C = length(allTrees$stat_id)
    
    # if there are stat IDs without data, we just need to make sure they are removed
    allTrees = allTrees %>% filter(stat_id %in% unique(X_data_C$stat_id))
    
    # Number of census DBH measurements
    N_Dobs = nrow(census_long)

    # Number of values to estimate with the RW + Census model 
    N_X_C = length(X_data_C$meas)
    
    # Dataframe that lists first and last RW estimate index for each tree for RW only model 
    idx_C = data.frame(stat_id = allTrees$stat_id)
    idx_C$firstidx = sapply(c(1:length(idx_C$stat_id)), 
                             function(ind){
                               val = min(which((X_data_C$stat_id == idx_C$stat_id[ind])))
                               if (val == -Inf){
                                 return(NA)
                               }else{
                                 return(val)
                               }
                             })
    idx_C$lastidx = sapply(c(1:length(idx_C$stat_id)), 
                            function(ind){
                              val = max(which((X_data_C$stat_id == idx_C$stat_id[ind])))
                              if (val == -Inf){
                                return(NA)
                              }else{
                                return(val)
                              }
                            })
    
    # Maps diameter values to their respective RW estimate index 
    Dobs2X = sapply(c(1:length(census_long$site)),
                    function(ind){
                      inds = which(((X_data_C$stat_id == census_long$stat_id[ind]) & (X_data_C$year == census_long$year[ind])))
                      return(inds)
                    })
      
    # Determine log of all census diameter measurements to fit with RW STAN model
    logDobs = log(census_long$dbh)  
    
    # Maps RW diameter measurements to respective RW estimate index for RW + Census model 
    Tr2X_C = sapply((1:length(pdbh$stat_id)),
                    function(ind){
                      which((X_data_C$stat_id == pdbh$stat_id[ind]) & (X_data_C$year == pdbh$year[ind]))
                    })
    # Tr2X_C = unlist(Tr2X_C)
    
    # Maps RW estimates to year index and stat_id for RW + Census model 
    X2year_C = X_data_C$year
    X2C = X_data_C$stat_id
    X2Tr_C = X_data_C$stat_id 
    
    # Tr2taxon_C = as.numeric(X_data_C$species_id[match(seq(1, N_C), X_data_C$stat_id)])
    Tr2taxon_C = match(allTrees$taxon, taxa)

    # Maps observed values of RW to respective RW estimate index for RW + Census model 
    # Xobs2X_C = sapply((1:length(incr_data$id)),
    #                 function(ind){
    #                   which((X_data_C$stat_id == incr_data$stat_id[ind]) & (X_data_C$year == incr_data$year[ind]))
    #                 })
    # Xobs2X_C = unlist(Xobs2X_C)
    
    
    # Xobs2X_C = lapply((1:length(incr_data$id)),
    #                   function(ind){
    #                     which((X_data_C$stat_id == incr_data$stat_id[ind]) & (X_data_C$year == incr_data$year[ind]))
    #                   })
    
    Xobs2X_C = rep(NA, length=N_Xobs)
    for (i in 1:N_Xobs){
      idx = which((X_data_C$stat_id == incr_data$stat_id[i]) & (X_data_C$year == incr_data$year[i]))
      
      if (length(idx)>1){
        print(paste0("Ring width obs ", i, "; match ", idx))
        stop()
      }

      Xobs2X_C[i] =  idx
    }
    
    # Data frame that contains all census measurements/data
    Dobs = census_long
  }
  
  # Maps RW diameter measurements to respective RW estimate index for RW only model 
  Tr2X = sapply((1:length(pdbh$stat_id)),
                function(ind){
                  which((X_data$stat_id == pdbh$stat_id[ind]) & (X_data$year == pdbh$year[ind]))
                })
  
  # Maps RW estimates to year index and stat_ids for RW only model 
  X2year = X_data$year
  X2Tr = X_data$stat_id
  
  # Data frame that lists first and last RW estimate index for each tree for RW only model 
  idx_Tr = year_idx
  idx_Tr$firstidx = sapply(c(1:length(idx_Tr$stat_id)), 
                           function(ind){
                             min(which((X_data$stat_id == idx_Tr$stat_id[ind])))
                           })
  idx_Tr$lastidx = sapply(c(1:length(idx_Tr$stat_id)), 
                           function(ind){
                             max(which((X_data$stat_id == idx_Tr$stat_id[ind])))
                           })
  idx_Tr = idx_Tr %>% dplyr::select(stat_id, firstidx, lastidx)
  
  # Maps observed values of RW to respective RW estimate index for RW only model 
  Xobs2X = sapply((1:length(incr_data$id)),
                  function(ind){
                    which((X_data$stat_id == incr_data$stat_id[ind]) & (X_data$year == incr_data$year[ind]))
                  })

  # Determine log of all RW diameter measurements to fit with RW STAN model
  logTr = log(pdbh$dbh)
  
  

  
  # Determine log of all RW increment measurements to fit with RW STAN model
  Xobs    = incr_data$incr
  Xobs[Xobs==0] = 0.0001
  logXobs = log(Xobs)
  
  # X2taxon = match(incr_data$taxon, taxa)
  X2taxon = Tr2taxon[X2Tr]
  X2taxon_C = Tr2taxon_C[X2Tr_C]
  
  # Larger save matrices
  # All RW measurement data
  Xobs = incr_data
  # All RW diameter measurement data
  Tr = pdbh

  #####################################################
  ################ 4. Save as RDS file ################
  #####################################################
   
  if (census_site){
    saveRDS(list(N_Xobs = N_Xobs, 
                 N_years = N_years, 
                 N_X = N_X, 
                 N_X_C = N_X_C, 
                 N_Tr = N_Tr,
                 N_Dobs = N_Dobs,
                 N_C = N_C, 
                 N_taxa = N_taxa,
                 Tr2X = Tr2X,
                 Tr2X_C = Tr2X_C,
                 Tr2taxon = Tr2taxon,
                 Tr2taxon_C = Tr2taxon_C,
                 X2C = X2C, 
                 X2Tr = X2Tr,
                 X2Tr_C = X2Tr_C,
                 X2year = X2year, 
                 X2year_C = X2year_C,
                 X2taxon = X2taxon,
                 X2taxon_C = X2taxon_C,
                 idx_C = idx_C, 
                 idx_Tr = idx_Tr, 
                 Dobs2X = Dobs2X, 
                 Xobs2X = Xobs2X, 
                 Xobs2X_C = Xobs2X_C,
                 logDobs = logDobs,
                 logTr = logTr, 
                 logXobs = logXobs,
                 Dobs = Dobs, 
                 Tr = Tr,
                 Xobs = Xobs,
                 allTrees = allTrees,
                 years = years,
                 taxa = taxa
    ),
    file=RDS_loc)
  }else{
    saveRDS(list(N_Xobs = N_Xobs, 
                 N_years = N_years, 
                 N_X = N_X, 
                 N_Tr = N_Tr,
                 N_taxa = N_taxa,
                 Tr2X = Tr2X,
                 Tr2taxon = Tr2taxon,
                 X2Tr = X2Tr,
                 X2year = X2year, 
                 X2taxon = X2taxon,
                 idx_Tr = idx_Tr, 
                 Xobs2X = Xobs2X, 
                 logTr = logTr, 
                 logXobs = logXobs,
                 Xobs = Xobs, 
                 Tr = Tr,
                 years = years,
                 taxa = taxa
    ),
    file=RDS_loc) 
  }
  
  ###############################################
  ################ 5. Check data ################
  ###############################################
  
  # Organize diameter estimates
  trees = unique(incr_data$stat_id)
  D0s = rep(NA,length(trees))
  D = rep(NA, length(incr_data$id))
  for (t in seq_along(trees)){
    
    stat = trees[t] 
    yrs = rev(seq(1,max((incr_data %>% filter(stat_id == stat))$year)))
    Dlast = Tr$dbh[which(Tr$stat_id == stat)]
    
    D[which((incr_data$stat_id == stat) & (incr_data$year == yrs[1]))] = Dlast
    
    Dnow = Dlast
    incr_now = min(incr_data$incr[which((incr_data$stat_id == stat) & (incr_data$year == yrs[1]))])
    
    # loop through years and calculate approximate diameter value 
    for (y in 2:length(yrs)){
      
      yrnow = yrs[y]
      Dnow = Dnow - (2 * incr_now / 10)
      
      # if year has increment data we need to save it 
      if (length(which((incr_data$stat_id == stat) & (incr_data$year == yrnow))) > 0){
        D[which((incr_data$stat_id == stat) & (incr_data$year == yrnow))] = Dnow
        incr_now = min(incr_data$incr[which((incr_data$stat_id == stat) & (incr_data$year == yrnow))])
      }else{
        incr_now = median(incr_data$incr)
      }
    }
  
    # get approximate D0 value for tree
    D0s[t] = Dnow - (2 * incr_now / 10)
  }
  
  incr_data$D = D
  D0s = data.frame(D0 = D0s, stat_id = trees)
  
  # check range of predicted D0s 
  if (min(D0s$D0) < -30 | max(D0s$D0) > 80) print('warning: estimated D0 value outside of model range')
  
  # PLOT 1: Look at values of D0 to see if diameter prior range is reasonable
  pl1 = ggplot(D0s) + 
    geom_histogram(aes(x = D0), binwidth = 1) + 
    labs(x = 'Estimated D0 value', 
         title = paste0('Histogram of D0 Values (', 
                        min(D0s$D0), ', ', max(D0s$D0), ')')) + 
    xlim(-40,100) + 
    geom_vline(xintercept = 80, color = 'red') + 
    geom_vline(xintercept = -30, color = 'red')
  ggsave(pl1, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','D0_histogram.jpg'))
  
  # PLOT 2: Check to make sure growth makes sense 
  pl2 = ggplot(incr_data) + 
    geom_line(aes(x = year, y = D, group = stat_id, color = stat_id)) + 
    facet_wrap(~as.factor(taxon)) + 
    labs(y = 'diameter (cm)', title = 'Species Growth over Time') + 
    theme(legend.position = 'none')
  ggsave(pl2, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','species_growth_check.jpg'))
  
  # Now, we need to check the D0 assumption for the census data
  if (census_site){
    
    # Then, we are going to determine what the D0 value would be in 1960 for each of the trees
    # based on their most recent census measurement
    medInc = median(incr_data$incr)
    
    # Let's just consider those trees without RW data because we already checked the other trees
    census_only = allTrees$stat_id[which(!(allTrees$stat_id %in% unique(incr_data$stat_id)))]
    D0s_census = rep(NA, length(census_only))
    
    for (t in seq_along(census_only)){
      
      stat = census_only[t] 
      datanow = census_long %>% filter(stat_id == stat)
      firstyr = min((X_data_C %>% filter(stat_id == stat))$year)
      
      if (length(datanow$site) < 1) next
      
      # Find earliest census measurement
      C1 = datanow$dbh[which.min(datanow$year)]
      C1year = min(as.numeric(datanow$year))
      
      # Approximate a D0 value in first year based on average increment 
      nyears = C1year - firstyr
      D0s_census[t] = C1 - (2 * nyears * medInc / 10)
    }
    
    D0s_census = data.frame(D0 = D0s_census, stat_id = census_only)
    
    # check range of predicted D0s 
    if (min(D0s_census$D0) < -30 | max(D0s_census$D0) > 80) print('warning: estimated D0 value outside of model range for census')
    
    pl3 = ggplot(D0s_census) + 
      geom_histogram(aes(x = D0), binwidth = 1) + 
      labs(x = 'Estimated D0 value', 
           title = paste0('Histogram of D0 Values - Census (', 
                          min(D0s_census$D0), ', ', max(D0s_census$D0), ')')) + 
      xlim(-40,100) + 
      geom_vline(xintercept = 80, color = 'red') + 
      geom_vline(xintercept = -30, color = 'red')
    ggsave(pl3, filename = file.path(site_dir,'runs',paste0(mvers,'_',dvers),'figures','D0_histogram_census.jpg'))
  }
}

