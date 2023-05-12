#!/usr/bin/env Rscript
# ==============================================
# === randomly rotate Lausanne parcellations ===
# ==============================================

# set working directory
setwd('/slow/projects/coco_genes')

# attach required packages
for (pkg in c('data.table','doRNG','foreach','matrixStats')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# source function rotate.parcellation
source('code/functions/rotate.parcellation.R')

# load data
coord.l = as.matrix(read.table('results/Lausanne_Centroids_left.csv', sep = ',', header = F))
coord.r = as.matrix(read.table('results/Lausanne_Centroids_right.csv', sep = ',', header = F))

# set up parallel pool
n.cores = 50 # parallel::detectCores() - 1
my.cluster = parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster) # foreach::getDoParRegistered()

# create random rotations
nrotations = 1000000
start_time = Sys.time()
set.seed(7990)
x = foreach(i = 1:n.cores, .combine = 'cbind', .packages='matrixStats') %dorng% {
  rotate.parcellation(coord.l,coord.r,nrot=nrotations/n.cores)
}
end_time = Sys.time()
end_time - start_time

# save rotations
write.table(x, 'results/Lausanne_perm_ids.txt', quote = F, col.names = F, row.names = F)
system('gzip -f results/Lausanne_perm_ids.txt')
system('chmod 770 results/Lausanne_perm_ids.txt.gz')

