#!/usr/bin/env Rscript
# ====================================================================
# === combine observed correlations and permutation-based p-values ===
# ====================================================================

# set working directory
setwd('/slow/projects/coco_genes')

# attach required packages to current R session
for (pkg in c('doParallel','data.table','bigmemory')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load observed correlations
alert.cor = read.delim('results/alert.obs.corr.txt', sep = '\t', header = T)
orient.cor = read.delim('results/orient.obs.corr.txt', sep = '\t', header = T)
control.cor = read.delim('results/control.obs.corr.txt', sep = '\t', header = T)

# load permutation-based p-valyes of observed correlations
alert.p = as.vector(unlist(read.table('results/alert.pval.twotailed.txt')))
orient.p = as.vector(unlist(read.table('results/orient.pval.twotailed.txt'))) 
control.p = as.vector(unlist(read.table('results/control.pval.twotailed.txt')))
pvals = data.frame(gene = alert.corr$genes, alert.p = alert.p, control.p = control.p, orient.p = orient.p) 

# merge with pvals
df = data.frame(gene = alert.cor$genes,
  alert.p = alert.p, control.p = control.p, orient.p = orient.p,
  alert.cor = alert.cor$cor.estimate, control.cor = control.cor$cor.estimate, orient.cor = orient.cor$cor.estimate) 

# calculate t-values
df$alert.t = qnorm(df$alert.p)*-sign(df$alert.cor)
df$control.t = qnorm(df$control.p)*-sign(df$control.cor)
df$orient.t = qnorm(df$orient.p)*-sign(df$orient.cor)

# ======================
# === get some stats ===
# ======================

# how many genes are nominally significant?
sapply(data.frame(df[,c('alert.p','control.p','orient.p')]<0.05),sum)

# how many genes are FDR-significant?
sapply(data.frame(sapply(data.frame(df[,c('alert.p','control.p','orient.p')]),p.adjust,method='BH')<0.05),sum)

# how many genes are Bonferroni-significant? caveat: number of permutations must be larger than 1/(0.05/nrow(results))
sapply(data.frame(df[,c('alert.p','control.p','orient.p')]<(0.05/nrow(df))),sum)

# what is the lowest p-value
sapply(data.frame(df[,c('alert.p','control.p','orient.p')]),min)

# how many unique p-values
sapply(sapply(data.frame(df[,c('alert.p','control.p','orient.p')]),unique),length)

# how many unique t-values
sapply(sapply(data.frame(df[,c('alert.t','control.t','orient.t')]),unique),length)

# save results
save(df, file = 'results/perm.results.RData')
write.table(df, file = 'results/perm.results.txt', sep = '\t', quote = F, col.names = T, row.names = F)

# save pvals
pvals = df[,c('gene','alert.p','control.p','orient.p')]
save(pvals, file = 'results/perm.pvals.RData')
write.table(pvals, file = 'results/perm.pvals.txt', sep = '\t', quote = F, col.names = T, row.names = F)
