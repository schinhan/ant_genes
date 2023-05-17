#!/usr/bin/env Rscript
# ====================================================================================
# === Calculate correlations of randomly rotated ANT maps and gene expression maps ===
# ====================================================================================

# set working directory
setwd('/.../02_data')

# attach required packages
for (pkg in c('doParallel','data.table','bigmemory')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# read data
ant_data = read.csv(file = '/03_fMRI/ant_lausanne.csv') # use /home/schinhan instead of data
genes_data = data.frame(fread('/02_parcellation/02_lausanne_parc_atlas.csv', sep = ',')) # use /home/schinhan instead of data
genes_data = genes_data[,-1] # delete column 'label'
perm_id_data = data.frame(fread('gzip -dc /02_parcellation/Lausanne_perm_ids.txt.gz'))

# calculate observed correlations and save them
for (cond in c('alert', 'orient', 'control')) {
  tmp = data.frame(t(cor(ant_data[,cond], genes_data, use = 'pairwise.complete.obs')))
  tmp$genes = row.names(tmp)
  row.names(tmp) = NULL
  names(tmp) = c('cor.estimate','genes')
  tmp = tmp[,c('genes','cor.estimate')]
  assign(paste0(cond,'.obs.corr'),tmp)
  write.table(tmp, sprintf('results/%s.obs.corr.txt',cond), sep = '\t', row.names = F, quote = F)
}

# set up parallel pool
n.cores = 50 # parallel::detectCores() - 1
my.cluster = parallel::makeCluster(n.cores, type = "FORK")
doParallel::registerDoParallel(cl = my.cluster) # foreach::getDoParRegistered()

# create big matrices that can be accessed by multiple workers (without copying them to each worker)
perm_big_matrix = as.big.matrix(x = as.matrix(perm_id_data), type = "integer", 
                                separated = FALSE, 
                                backingfile = "perms.bin", 
                                descriptorfile = "perms.desc")
genes_big_matrix = as.big.matrix(x = as.matrix(genes_data), type = "double", 
                                 separated = FALSE, 
                                 backingfile = "genes.bin", 
                                 descriptorfile = "genes.desc")

# calculate correlations
n.permutations = ncol(perm_id_data)
iterationsPerWorker = n.permutations/n.cores # must be an integer, check by typing: iterationsPerWorker == as.integer(iterationsPerWorker)
geneLabels = names(genes_data)

start_time = Sys.time()
for (cond in c('alert', 'orient', 'control')) {
  message(paste0('Starting with ', cond, ' | n.permutations: ', n.permutations, ' | n.cores: ', n.cores, ' | iterations per worker: ', iterationsPerWorker)) 
  system(paste0('rm -f /03_fMRI/',cond, '/*; mkdir -p ','results/',cond))
  
  # start parallel loop
  x = foreach(i = 1:n.cores, .packages = c('bigmemory', 'data.table')) %dopar% { # .combine = 'rbind', 
    
    # get big matrices from shared memory
    perms = attach.big.matrix("perms2.desc")
    genes = attach.big.matrix("genes2.desc")
    
    # define block for worker
    j.start = i*iterationsPerWorker-iterationsPerWorker+1
    j.end = i*iterationsPerWorker
    
    # calculate correlations
    tmp = data.frame(matrix(NA, nrow = 1, ncol = ncol(genes[]))) # mp = data.frame(matrix(NA, nrow = iterationsPerWorker, ncol = ncol(genes[])))
    names(tmp) = geneLabels
      
    # k = 0
    for (j in j.start:j.end) {
      # k = k+1
      tmp[1,] = cor(ant_data[perms[,j],cond], genes[], use = 'pairwise.complete.obs') # tmp[k,] = cor(ant_data[perms[,j],cond], genes[], use = 'pairwise.complete.obs')
      if (j == j.start) { appendLogical = F } else { appendLogical = T }
      data.table::fwrite(tmp, file = paste0('/03_fMRI/',cond,'/exp.corr.',sprintf("%04d", i),'.txt'), sep = '\t', quote = F, row.names = F, append = appendLogical)
    }   
  }
}
Sys.time() - start_time

# delete temporary big matrix files
system('rm -f genes.bin; rm -f genes.desc; rm -f perms.bin; rm -f perms.desc')

# shut down cluster
stopCluster(my.cluster)

