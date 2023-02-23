###############################################
############ GO-Figure PREPARATION ############
###############################################

#set wd to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#packages
library(stringr)

#read in signifcant molecular functions from gene set enrichment analysis with corresponding p-values
alert_pval <- read.table("GSEA.alert.index.05.GO.txt", sep = '\t', header = TRUE)
alert <- alert_pval[,c(1,4)]              #delete unnecessary cols
names(alert)[1] <- "%  GO term"           #rename cols
names(alert)[2] <- "enrichment_P-value"
alert[,1] <- str_sub(alert[,1],-11,-1)    #set characters as input
alert[,1] <- substr(alert[,1], 1, 10)  
write.table(alert, file = 'alert.gofigure.input.txt', sep = '\t', quote = F, col.names = T, row.names = F)

control_pval <- read.table("GSEA.control.index.05.GO.txt", sep = '\t', header = TRUE)
control <- control_pval[,c(1,4)]              #delete unnecessary cols
names(control)[1] <- "%  GO term"           #rename cols
names(control)[2] <- "enrichment_P-value"
control[,1] <- str_sub(control[,1],-11,-1)    #set characters as input
control[,1] <- substr(control[,1], 1, 10)  
write.table(control, file = 'control.gofigure.input.txt', sep = '\t', quote = F, col.names = T, row.names = F)

orient_pval <- read.table("GSEA.orient.index.05.GO.txt", sep = '\t', header = TRUE)
orient <- orient_pval[,c(1,4)]              #delete unnecessary cols
names(orient)[1] <- "%  GO term"           #rename cols
names(orient)[2] <- "enrichment_P-value"
orient[,1] <- str_sub(orient[,1],-11,-1)    #set characters as input
orient[,1] <- substr(orient[,1], 1, 10)  
write.table(orient, file = 'orient.gofigure.input.txt', sep = '\t', quote = F, col.names = T, row.names = F)
