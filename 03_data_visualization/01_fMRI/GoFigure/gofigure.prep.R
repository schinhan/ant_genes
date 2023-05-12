# convert output files from gene set enrichment analysis 
# with Panther to correct input file format
# to make summary visualizations with GO-Figure

# required package: stringr

# set wd to data folder
setwd("../../03_data_visualization/01_fMRI/GoFigure")

library(stringr)

#read in signifcant molecular functions from gene set enrichment analysis with corresponding p-values
alert_pval <- read.table("../../02_data/03_fMRI/GSEA.alert.index.05.GO.txt", sep = '\t', header = TRUE)
alert <- alert_pval[,c(1,5)]              #delete unnecessary cols
names(alert)[1] <- "%  GO term"           #rename cols
names(alert)[2] <- "enrichment_P-value"
alert[,1] <- str_sub(alert[,1],-11,-1)    #set characters as input
alert[,1] <- substr(alert[,1], 1, 10)  
write.table(alert, file = 'gofigure.input.fdr.alert.txt', sep = '\t', quote = F, col.names = T, row.names = F)

control_pval <- read.table("../../02_data/03_fMRI/GSEA.control.index.05.GO.txt", sep = '\t', header = TRUE)
control <- control_pval[,c(1,5)]              #delete unnecessary cols
names(control)[1] <- "%  GO term"           #rename cols
names(control)[2] <- "enrichment_P-value"
control[,1] <- str_sub(control[,1],-11,-1)    #set characters as input
control[,1] <- substr(control[,1], 1, 10)  
write.table(control, file = 'gofigure.input.fdr.control.txt', sep = '\t', quote = F, col.names = T, row.names = F)

orient_pval <- read.table("../../02_data/03_fMRI/GSEA.orient.index.05.GO.txt", sep = '\t', header = TRUE)
orient <- orient_pval[,c(1,5)]              #delete unnecessary cols
names(orient)[1] <- "%  GO term"           #rename cols
names(orient)[2] <- "enrichment_P-value"
orient[,1] <- str_sub(orient[,1],-11,-1)    #set characters as input
orient[,1] <- substr(orient[,1], 1, 10)  
write.table(orient, file = 'gofigure.input.fdr.orient.txt', sep = '\t', quote = F, col.names = T, row.names = F)
