#!/usr/bin/env Rscript

# ===========================
# === draw manhattan plot ===
# ===========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=10) {
  stop(paste0('expected 10 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
abagenFile = args[1] # abagenFile='results/perm.pvals.txt'
targetDir = args[2] # targetDir="results/violin/"
prune = args[3] # prune = TRUE
drawCI = args[4] # drawCI = FALSE
xend = as.numeric(args[5]) # xend = 6
xsteps = as.numeric(args[6]) # xsteps = 2
yend = as.numeric(args[7]) # yend = 8
ysteps = as.numeric(args[8]) # ysteps = 2
width = as.numeric(args[9]) # width = 3
height = as.numeric(args[10]) # height = 4

message(paste0('\n--- qq-plot settings ---',
               '\nabagenFile: ', abagenFile,
               '\ntargetDir: ', targetDir,
               '\nprune: ', prune,
               '\ndrawCI: ', drawCI,
               '\nxend: ', xend,
               '\nxsteps: ', xsteps,
               '\nyend: ', yend,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('data.table','dplyr','ggplot2', 'stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# create folder
system(paste0('mkdir -p ', targetDir))

# import data
message(paste0('[1/2] Importing data.'))
df = read.delim(abagenFile, sep='\t', header=T, stringsAsFactors=FALSE)

# -----------------------------
# -- 'makeviolin' function ---
# -----------------------------

# prepare matrix for ggplot
df.m = reshape2::melt(df, id="gene", variable.name="ant", value.name="pvalue", na.rm = F)
df.m$ant = factor(df.m$ant, levels=c('p.alert','p.orient','p.control'))

# draw plot
violin = ggplot(df.m, aes(x = ant, y = -log10(pvalue), fill = ant)) + # aes_string(x = 'sex', fill = 'sex', y = paste0(col.name) ))
  geom_violin(lwd = 0.25) +
  #geom_point(alpha = 0.04) +
  geom_boxplot(width=0.025, fill = 'white', outlier.shape = NA, alpha = 1) + # 
  # ggtitle(header) +
  # ylab(y.title) +
  # scale_y_continuous(limits = c(ylower,yupper), breaks = seq(ylower,yupper,ysteps)) +
  theme_bw(base_size=10) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = -0.5),
        axis.text.x = element_text(size = 10 ),  # rotate x-axis labels so they don't overlap,
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        legend.position = 'none',
        axis.ticks = element_line(size = 0.25),
        plot.margin = unit(c(3, 3, 3, 3), "mm")) +
  # show.y +
  geom_segment(aes(x=min(as.numeric(ant)),xend=max(as.numeric(ant)),y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
  geom_segment(aes(y=0,yend=6,x=-Inf,xend=-Inf), colour = "black", size = 0.25)

hist_large = ggplot(df.m, aes(x = pvalue, fill = ant)) + # aes_string(x = 'sex', fill = 'sex', y = paste0(col.name) ))
  geom_histogram(lwd = 0.25) +
  facet_wrap(~ant) +
  theme_bw(base_size=10) +
  guides(fill="none")

hist_small = 
ggplot(df.m, aes(x = pvalue, fill = ant)) + # aes_string(x = 'sex', fill = 'sex', y = paste0(col.name) ))
  geom_histogram(lwd = 0.25) +
  scale_x_continuous(lim = c(0,0.0001), labels = function(x) sprintf("%g", x)) +
  facet_wrap(~ant) +
  theme_bw(base_size=10) +
  guides(fill="none")

library(patchwork)
violin / hist_large / hist_small

+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = -0.5),
          axis.text.x = element_text(size = 10 ),  # rotate x-axis labels so they don't overlap,
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          #axis.title.y = element_blank(),
          legend.position = 'none',
          axis.ticks = element_line(size = 0.25),
          plot.margin = unit(c(3, 3, 3, 3), "mm")) +
    # show.y +
    geom_segment(aes(x=min(as.numeric(ant)),xend=max(as.numeric(ant)),y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
    geom_segment(aes(y=0,yend=8,x=-Inf,xend=-Inf), colour = "black", size = 0.25)

  #geom_point(alpha = 0.04) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  # ggtitle(header) +
  # ylab(y.title) +
  # scale_y_continuous(limits = c(ylower,yupper), breaks = seq(ylower,yupper,ysteps)) +
  


hist(df$p.orient)



