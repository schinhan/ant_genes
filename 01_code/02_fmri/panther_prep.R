#################################
### prepare files for panther ###
#################################

# purpose: arrange genes according to their z-values obtained from permutation testing
# and prepare file for panther 

setwd("~/.../02_data/03_fMRI/")

#read and prepare data 
pvals <- read.table("04_perm.pvals.txt", sep = '\t',header = TRUE)
genes_data <- read.csv(file = "lausanne_parc_atlas.csv", sep = ',')
genes_dat <- genes_data[,-c(1)]
ant_data <- read.csv(file = "ant_lausanne.csv", sep = ',')
dat <- cbind(ant_data,genes_dat) #merge dfs
dat <- na.omit(dat)              #delete NAs
ant_dat <- dat[,2:4]             #seperate for correlations
genes_dat <- dat[,5:15636]

#calculate descriptive statistics
fdr.alert <- p.adjust(pvals$p.alert, method = "fdr", n = length(pvals$p.alert))
fdr.control <- p.adjust(pvals$p.control, method = "fdr", n = length(pvals$p.control))
fdr.orient <- p.adjust(pvals$p.orient, method = "fdr", n = length(pvals$p.orient))

#calculate observed correlations
alert.corr.obs <- cor (ant_dat[,1], genes_dat, method = "pearson") #calculate correlations 
alert.corr.obs <- as.data.frame(t(alert.corr.obs))  #format df 
colnames(alert.corr.obs) <- "alert.corr.obs" #rename to prevent confusion
orient.corr.obs <- cor (ant_dat[,2], genes_dat, method = "pearson")
orient.corr.obs <- as.data.frame(t(orient.corr.obs))
colnames(orient.corr.obs) <- "orient.corr.obs"
control.corr.obs <- cor (ant_dat[,3], genes_dat, method = "pearson")
control.corr.obs <- as.data.frame(t(control.corr.obs))
colnames(control.corr.obs) <- "control.corr.obs"

#calculate z-values
alert.zval = qnorm(pvals$p.alert/2) * -sign(ant.genes.corr.obs$alert.corr.obs)
orient.zval = qnorm(pvals$p.orient/2) * -sign(ant.genes.corr.obs$orient.corr.obs)
control.zval = qnorm(pvals$p.control/2) * -sign(ant.genes.corr.obs$control.corr.obs)

#### arrange genes for panther ###
library(dplyr)
library(tibble)
ant.genes.corr.obs <- tibble::rownames_to_column(ant.genes.corr.obs, "genes")

#merge data frames
genes <- merge(ant.genes.corr.obs, perm.zvals, by='genes')

#arrange by zval and obs.corr
alert <- genes[,c("genes", "alert.corr.obs", "alert.zval")]
alert <- alert %>% arrange(alert.zval, alert.corr.obs)
alert$range.index <- 1:15632
alert <- alert[,c(1,4)]

orient <- genes[,c("genes", "orient.corr.obs", "orient.zval")]
orient <- orient %>% arrange(orient.zval, orient.corr.obs)
orient$range.index <- 1:15632
orient <- orient[,c(1,4)]

control <- genes[,c("genes", "control.corr.obs", "control.zval")]
control <- control %>% arrange(control.zval, control.corr.obs)
control$range.index <- 1:15632
control <- control[,c(1,4)]

#write data table 
write.table(alert, file = 'alert.range.index.txt', sep = '\t', quote = F, col.names = F, row.names = F)
write.table(orient, file = 'orient.range.index.txt', sep = '\t', quote = F, col.names = F, row.names = F)
write.table(control, file = 'control.range.index.txt', sep = '\t', quote = F, col.names = F, row.names = F)


