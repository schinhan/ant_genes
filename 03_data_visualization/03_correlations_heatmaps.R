###########################################
######## Correlation Visualisation ########
###########################################

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ClusterR)
library(cluster)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets the wd to the folder where this script is in

#####################################
#######  3 ANT x 15.000 Genes #######
#####################################

#read data 
ant.genes.corr.obs <- read.table("ant.genes.corr.obs.txt", sep = '\t',header = TRUE)

#change variable names
colnames(ant.genes.corr.obs)[1] <- "Alerting    "
colnames(ant.genes.corr.obs)[2] <- "Orienting    "
colnames(ant.genes.corr.obs)[3] <- "Control    "

#switch col order to match with other figures
ant.genes.corr.obs <- ant.genes.corr.obs[,c(3,2,1)]

#melt data frame
ant.obs.corr.melt <- melt(as.matrix(ant.genes.corr.obs))

#create heat map 
ant.obs.corr.plt <-
  ggplot(data = ant.obs.corr.melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  xlab("Genes") +
  theme(text=element_text(size=16,family="serif"),
        plot.margin = unit(c(7.6,4,7.6,4), "pt"),
        legend.title=element_blank(),                 #set legend
        legend.key.height= unit(1, 'cm'),
        legend.box.margin=margin(0,0,0,13),
        axis.text.x=element_blank(),
        axis.title.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color = "black", size=16)) +
  scale_fill_gradient2(limit = c(-0.55,0.55), low = "dodgerblue3", high =  "red4", mid = "white", midpoint = 0)

#save plot
ggsave('ant.obs.corr.png', ant.obs.corr.plt, width = 7.11, height = 2.89, units = 'in', dpi = 300)

#######################################
######## 15.000 Genes x 219 ROIs ######
#######################################

ant_data <- read.csv(file = "lausanne_parc_atlas.csv", sep = ',')
ant_data <- na.omit(ant_data)
ant_data <- ant_data[,-1]
rownames(ant_data) <- 1:nrow(ant_data)
df <- data.frame(t(ant_data)) #rotate df for clustering

#determine optimal number of clusters 
library(factoextra)
# Elbow method
fviz_nbclust(ant_data, kmeans, method = "wss", k.max = 20) +
  labs(subtitle = "Elbow method")

#result: 2 clusters

#cluster ROIs
set.seed(240) # Setting seed
cluster.roi <- kmeans(x=ant_data, centers=2) #build clusters
cluster.index.roi <- cluster.roi$cluster  #extract cluster indices
df.roi <- cbind(cluster.index.roi, ant_data) #add indices to df
df.roi <- df.roi[order(cluster.index.roi),] #order df by cluster index
df.roi <- df.roi[,-1] #delete index col


#cluster genes 
df <- data.frame(t(df.roi)) #rotate df for clustering
set.seed(240) # Setting seed
cluster <- kmeans(x=df, centers=2) #build clusters
cluster.index <- cluster$cluster  #extract cluster indices
df <- cbind(cluster.index, df) #add indices to df
df <- df[order(cluster.index),] #order df by cluster index
df <- df[,-1] #delete index col
df.new = df[seq(1, nrow(df), 1), ] #only keep every 8th gene
df.melt <- melt(as.matrix(df.new)) #melt df for plotting

#create heat map 
gene.expr.plt <-
ggplot(data = df.melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  xlab("AHBA genes") +
  ylab("Lausanne ROIs") +
  theme(text=element_text(size=16,family="serif"),
        plot.margin = unit(c(7.6,4,7.6,4), "pt"),
        legend.title=element_blank(),                 #set legend
        legend.key.height= unit(1, 'cm'),
        legend.box.margin=margin(0,0,0,13),
        axis.text.x=element_blank(),
        axis.title.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_text(size=16, margin = margin(r = 0))) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "red4")

#save plot
ggsave('gene.expr.plt.png', gene.expr.plt, width = 7.11, height = 2.89, units = 'in', dpi = 300)



#########################################
###### 219 ROIs x 3 ANT Conditions ######
#########################################

ant.activation.obs <- read.csv(file = "ant_lausanne.csv", sep = ',')
ant.activation.obs <- ant.activation.obs[,-1]
ant.activation.obs <- ant.activation.obs[, c(2,3,1)] #reorder coloumns
rownames(ant.activation.obs) <- 1:219 #add ascending index

#rename cols
colnames(ant.activation.obs)[3] <- "Alerting"
colnames(ant.activation.obs)[2] <- "Control"
colnames(ant.activation.obs)[1] <- "Orienting"

#switch col order to match with other figures
ant.activation.obs <- ant.activation.obs[,c(2,1,3)]

#melt df
ant.activation.melt <- melt(as.matrix(ant.activation.obs))

#create heat map 
ant.activation <- 
ggplot(data = ant.activation.melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  xlab("Lausanne ROIs") +
  theme(text=element_text(size=14,family="serif"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title=element_blank(),                 #set legend
        legend.key.height= unit(1, 'cm'),
        legend.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(color = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color = "black", size=16)) +
  scale_fill_gradient2(limit = c(-6.2,6.2), low = "dodgerblue3", high =  "red4", mid = "white", midpoint = 0)
 
#save plot
ggsave('ant.activation.plt.png', ant.activation, width = 7.11, height = 2.89, units = 'in', dpi = 300)
