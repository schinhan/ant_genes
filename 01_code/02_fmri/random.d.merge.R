#!/usr/bin/env Rscript
# ====================================================================
# === combine observed correlations and permutation-based p-values ===
# ====================================================================

# set working directory
setwd('/.../02_data/03_fMRI')

# attach required packages to current R session
for (pkg in c('dplyr','data.table')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load observed correlations
alert.cor = read.delim('alert.obs.corr.txt', sep = '\t', header = T)
orient.cor = read.delim('orient.obs.corr.txt', sep = '\t', header = T)
control.cor = read.delim('control.obs.corr.txt', sep = '\t', header = T)

# load permutation-based p-valyes of observed correlations
alert.p = as.vector(unlist(read.table('alert.pval.twotailed.txt')))
orient.p = as.vector(unlist(read.table('orient.pval.twotailed.txt'))) 
control.p = as.vector(unlist(read.table('control.pval.twotailed.txt')))
pvals = data.frame(gene = alert.cor$genes, alert.p = alert.p, control.p = control.p, orient.p = orient.p) 

# merge with pvals
df = data.frame(gene = alert.cor$genes,
  alert.p = alert.p, control.p = control.p, orient.p = orient.p,
  alert.cor = alert.cor$cor.estimate, control.cor = control.cor$cor.estimate, orient.cor = orient.cor$cor.estimate) 

# calculate z-values
df$alert.z = qnorm(df$alert.p)*-sign(df$alert.cor)
df$control.z = qnorm(df$control.p)*-sign(df$control.cor)
df$orient.z = qnorm(df$orient.p)*-sign(df$orient.cor)

# ======================
# === get some stats ===
# ======================

# how many genes are nominally significant?
sapply(data.frame(df[,c('alert.p','control.p','orient.p')]<0.05),sum)

# how many genes are FDR-significant?
df$alert.fdr = p.adjust(df$alert.p,method='BH')
df$control.fdr = p.adjust(df$control.p,method='BH')
df$orient.fdr = p.adjust(df$orient.p,method='BH')
sapply(data.frame(df[,c('alert.fdr','control.fdr','orient.fdr')]<0.05),sum)

# how many genes are Bonferroni-significant? caveat: number of permutations must be larger than 1/(0.05/nrow(results))
sapply(data.frame(df[,c('alert.p','control.p','orient.p')]<(0.05/nrow(df))),sum)

# what is the lowest p-value
sapply(data.frame(df[,c('alert.p','control.p','orient.p')]),min)

# how many unique p-values
sapply(sapply(data.frame(df[,c('alert.p','control.p','orient.p')]),unique),length)

# how many unique z-values
sapply(sapply(data.frame(df[,c('alert.z','control.z','orient.z')]),unique),length)

# save results
save(df, file = 'perm.results.RData')
write.table(df, file = 'perm.results.txt', sep = '\t', quote = F, col.names = T, row.names = F)

# save pvals
pvals = df[,c('gene','alert.p','control.p','orient.p')]
save(pvals, file = 'perm.pvals.RData')
write.table(pvals, file = 'results/perm.pvals.txt', sep = '\t', quote = F, col.names = T, row.names = F)

# ========================================
# === save gene results for supplement ===
# ========================================

# get gene description from prepared RefSeq file)
refseq = data.frame(fread(cmd=paste0("gzip -dc .../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/GRCh37_latest_genomic.edit.gff.gz"), sep='\t', header=F, stringsAsFactors=FALSE))

# get info
message('Querying gene infos.')
info = data.frame(matrix(NA, nrow = nrow(df), ncol = 6))
names(info) = c('gene','chr','start','stop','symbol','description')
pb = txtProgressBar(min = 1, max = nrow(info), style = 3)
for (i in 1:nrow(df)) {
  setTxtProgressBar(pb, i)
  info$gene[i] = df$gene[i]
  tmp = refseq[grep(paste0("\\b",df$gene[i],"\\b"),refseq$V10),]
  if (nrow(tmp) == 0) { tmp = refseq[grep(paste0("\\b",df$gene[i],"\\b"),refseq$V13),] }
  if (nrow(tmp) == 0) { tmp = refseq[grep(paste0("\\b",df$gene[i],"\\b"),refseq$V14),] }
  if (nrow(tmp) > 0) {
    info$chr[i] = tmp$V9[1]
    info$start[i] = tmp$V4[1]
    info$stop[i] = tmp$V5[1]
    info$symbol[i] = tmp$V10[1]
    info$description[i] = tmp$V12[1]
  }
}

# merge refseq and df
df = left_join(df, info, by = 'gene')

# create output
df$top.corAbs =  df[,grep('.cor',names(df))] %>% abs() %>% apply(1, FUN = max)
df$top.p = df[,grep('[.]p',names(df))] %>% apply(1, FUN = min)
df$top.fdr = df[,grep('[.]fdr',names(df))] %>% apply(1, FUN = min)
df$`NA` = ''
output = df[,c('gene','chr','start','stop','description','top.corAbs','top.p','top.fdr',
    'NA','alert.cor','alert.z','alert.p','alert.fdr',
    'NA','control.cor','control.z','control.p','control.fdr',
    'NA','orient.cor','orient.z','orient.p','orient.fdr')]
output = output[order(output$top.fdr,output$top.p,-output$top.corAbs),]
write.table(output, file = 'perm.results.suppl.txt', sep = '\t', quote = F, col.names = T, row.names = F)
