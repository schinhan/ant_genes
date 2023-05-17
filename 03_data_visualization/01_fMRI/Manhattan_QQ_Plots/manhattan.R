#!/usr/bin/env Rscript

# ======================================
# === draw gene-based manhattan plot ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
pFile = args[1] # pFile='results/perm.pvals.txt'
refseqFile = args[2] # refseqFile = "/.../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/GRCh37_latest_genomic.edit.gff.gz"
targetDir = args[3] # targetDir='/.../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/' # annotationThresh = args[4] # annotationThresh = 'bonferroni' # 'fdr'
yend = as.numeric(args[4]) # yend = 10
ysteps = as.numeric(args[5]) # ysteps = 2
width = as.numeric(args[6]) # width = 11
height = as.numeric(args[7]) # height = 4

message(paste0('\n--- Manhatten plot settings ---',
               '\npFile: ', pFile,
               '\nrefseqFile: ', refseqFile,
               '\ntargetDir: ', targetDir,
             # '\nannotationThresh: ', annotationThresh,
               '\nyend: ', yend,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# load packages
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'ggrepel', 'stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# create targetDir
system(paste0('mkdir -p ', targetDir))

# import data
message(paste0('Importing data.'))
pvals = read.delim(pFile, sep='\t', header=T, stringsAsFactors=FALSE)
refseq = data.frame(fread(cmd=paste0("gzip -dc ", refseqFile), sep='\t', header=F, stringsAsFactors=FALSE))

# get info
message('Querying gene infos.')
info = data.frame(matrix(NA, nrow = nrow(pvals), ncol = 6))
names(info) = c('gene','chr','start','stop','symbol','description')
pb = txtProgressBar(min = 1, max = nrow(info), style = 3)
for (i in 1:nrow(pvals)) {
  setTxtProgressBar(pb, i)
  info$gene[i] = pvals$gene[i]
  tmp = refseq[grep(paste0("\\b",pvals$gene[i],"\\b"),refseq$V10),]
  if (nrow(tmp) == 0) { tmp = refseq[grep(paste0("\\b",pvals$gene[i],"\\b"),refseq$V13),] }
  if (nrow(tmp) == 0) { tmp = refseq[grep(paste0("\\b",pvals$gene[i],"\\b"),refseq$V14),] }
  if (nrow(tmp) > 0) {
    info$chr[i] = tmp$V9[1]
    info$start[i] = tmp$V4[1]
    info$stop[i] = tmp$V5[1]
    info$symbol[i] = tmp$V10[1]
    info$description[i] = tmp$V12[1]
  }
}

# merge refseq and df
df = left_join(pvals, info, by = 'gene')

# 74 without match, set manually
idx = is.na(df$symbol)
df$symbol[idx] = df$gene[idx]
df$chr[idx] = 99
df$start[idx] = df$stop[idx] = 1:sum(idx)
df$description[idx] = 'no description'

# create function for plot
manhplot = function(GWAS = NULL, pCol = NULL) {

  # prepare for manhattan plot
  names(GWAS)[names(GWAS)=='chr'] = 'CHR'
  names(GWAS)[names(GWAS)=='start'] = 'BP'
  names(GWAS)[names(GWAS)==pCol] = 'P'
  GWAS$CHR[GWAS$CHR=='X'] = 23
  GWAS$CHR[GWAS$CHR=='Y'] = 24
  #GWAS$CHR[GWAS$CHR==25] = 23 # GWAS$CHR[GWAS$CHR=='XY'] = 23
  #GWAS$CHR[GWAS$CHR==26] = 25 # GWAS$CHR[GWAS$CHR=='MT'] = 25
  GWAS$CHR = suppressWarnings(as.numeric(GWAS$CHR))
  GWAS$CHR[is.na(GWAS$CHR)] = 99
  
  # sort by chromosome
  idx = sort(GWAS$CHR, index.return = T)$ix
  GWAS = GWAS[idx,]
  
  # get cumulative base pair positions and center positions
  # credits to Daniël Roelfs (http://www.danielroelfs.com/coding/manhattan_plots/)
  message(paste0('Geting cumulative base pair positions.'))
  nCHR = length(unique(GWAS$CHR))
  GWAS$BPcum = NA
  s = 0
  nbp = c()
  for (i in unique(GWAS$CHR)){
    nbp[i] = max(GWAS[GWAS$CHR == i,]$BP)
    GWAS[GWAS$CHR == i,'BPcum'] = GWAS[GWAS$CHR == i,'BP'] + s
    s = s + nbp[i]
  }
  axis.set = GWAS %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  
  # edit GWAS findings exceeding plot
  GWAS$P_exceeding = GWAS$P
  GWAS$BPcum_exceeding = GWAS$BPcum
  GWAS$Gene_exceeding = GWAS$gene
  
  if (sum(!is.na(GWAS$Gene) & GWAS$P < 10^-yend) > 0) {
    GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-yend,]$P_exceeding = 10^-yend
    GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-yend,]$BPcum_exceeding = GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-yend,]$BPcum + 60000000
    GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-yend,]$Gene_exceeding = paste0(GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-yend,]$Gene, '\n(p = ', sprintf('%0.e',GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-yend,]$P), ')')
  }
  
  # set yend significance threshold
  GWAS$FDR = p.adjust(GWAS$P, method = 'BH')
  sigFDR = max(GWAS$P[GWAS$FDR < 0.05])
  sigBonf = 0.05/nrow(GWAS)
  
  # # define Genes that shall be annotated
  GWAS$GENE_COUNT = 1 
  
  # create manhattan plots
  message('Creating manhattan plot.')
  set.seed(8245)
  manhplot = ggplot(data = subset(GWAS, P >= 10^-yend)) +
    geom_point(aes(x=BPcum, y=-log10(P), color=as.factor(CHR)), alpha = 1, size = 1.2, stroke = 0, shape = 16) +
    #geom_point(data=subset(GWAS, GENE_COUNT == 1 & P <= sigFDR & P > sigBonf), aes(x=BPcum, y=-log10(P)), color='black', shape=1, size=1.5) +
    geom_point(data=subset(GWAS, GENE_COUNT == 1 & P <= sigBonf & P >= 10^-yend), aes(x=BPcum, y=-log10(P)), color='black', shape=5, size=2) +
    scale_color_manual(values = rep(c('#282873', '#6e90ca'), nCHR)) +
    scale_x_continuous(expand = expansion(mult = c(0.03,0.03), add = c(0,0)), label = c(1:18,'', 20, '', 22, 'X', 'Y U', ''), breaks = axis.set$center) + # label = axis.set$CHR % label = c(1:22, 'X', 'Y MT', '') label = c(1:18,'', 20, '', 22, 'X', 'Y MT', '')
    scale_y_continuous(expand = expansion(mult = c(0,0), add = c(0,0.5)), limits = c(0,yend), breaks = seq(0,yend,ysteps)) +
    geom_hline(yintercept = -log10(sigBonf), color = 'black', linetype = 'dashed', size = 0.25) +
    geom_hline(yintercept = -log10(sigFDR), color = 'black', linetype = 'solid', size = 0.25) + 
    labs(x = 'Chromosome', y = expression('-log'[10]*'('*italic(p)*')')) + 
    geom_segment(aes(x=min(axis.set$center),xend=max(axis.set$center),y=-Inf,yend=-Inf), colour = 'black', size = 0.25)+
    geom_segment(aes(y=0,yend=yend,x=-Inf,xend=-Inf), colour = 'black', size = 0.25) +
    theme_bw() +
    theme( 
      legend.position = 'none',
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = 'black', size = 0.25),
      axis.ticks.length=unit(.15, 'cm'),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
      axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 5, b = 0, l = 5)),
      axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      plot.margin=unit(c(0.25,0.75,0,0),'cm')
    )
  manhplot
}

# save plots
message('Saving manhattan plots.')
pCols = names(df)[grep('\\.p',names(df))]
traits = str_replace(pCols, '.p', '')
for (i in 1:length(traits)) {
  tmp = manhplot(df, pCols[i])
  ggsave(paste0(targetDir, '/', traits[i], '.manhattan.png'), tmp, width = width, height = height, units = 'in', dpi = 300)
  system(paste0('chmod -R 770 ', targetDir, '/', traits[i], '.manhattan.png'))
}
message(paste0('--- Manhatten plots finished ---\n'))
