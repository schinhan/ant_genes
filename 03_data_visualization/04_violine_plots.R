#############################
####### Violine Plots #######
#############################


#set the wd to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#packages
library(ggplot2)
library(tidyr)
library(data.table)
library(yarrr)


#################################
####### Violine Plots ANT #######
#################################

##prepare data
ant.exp.corr <- read.table("ant.exp.corr.txt", sep = '\t', header = TRUE)
exp.dat <- ant.exp.corr

#rename cols for correct order
colnames(exp.dat) <- c("control.alert.exp",
                       "alert.orient.exp",
                       "control.orient.exp")

#transform df to long format
exp.dat <- gather(exp.dat, condition, corr.exp, control.alert.exp:control.orient.exp)

#observed correlations
ant_data = read.csv(file = 'ant_lausanne.csv')

alert.control.obs <- cor(ant_data$alert, ant_data$control)
alert.orient.obs <- cor(ant_data$alert, ant_data$orient)
control.orient.obs <- cor(ant_data$orient, ant_data$control)

#calculate p-cutoff
#sort in ascending order and select value at .05 cut-off
p.alert.control   <- sort(abs(ant.exp.corr$alert.control.exp), decreasing = TRUE)[round(0.05*nrow(ant.exp.corr))] 
p.alert.orient  <- sort(abs(ant.exp.corr$alert.orient.exp), decreasing = TRUE)[round(0.05*nrow(ant.exp.corr))] 
p.control.orient <- sort(abs(ant.exp.corr$control.orient.exp), decreasing = TRUE)[round(0.05*nrow(ant.exp.corr))] 

# Basic violin
violine.ant <-
ggplot(exp.dat, aes(x=condition, y=corr.exp, fill=condition)) +   #violine
  scale_x_discrete(labels = c("Alert vs. Orient", "Alert vs. Control", "Control vs. Orient")) +
  scale_y_continuous(limits = c(-0.45, 1), breaks = seq(-0.4, 1, by = 0.2)) +
  geom_point(aes(x="control.alert.exp", y=alert.control.obs),     #observed corr
             col = yarrr::transparent("lightgreen", trans.val = .98), 
             shape = 20, 
             size = 5) +
  geom_point(aes(x="alert.orient.exp", y=alert.orient.obs),
             col = yarrr::transparent("lightpink", trans.val = .98), 
             shape = 20, 
             size = 5) +
  geom_point(aes(x="control.orient.exp", y=control.orient.obs),
             col = yarrr::transparent("lightskyblue1", trans.val = .98), 
             shape = 20, 
             size = 5) +
  geom_violin()+
  geom_segment(aes(x = 1.7, y = p.alert.control, xend = 2.3, yend = p.alert.control),
               colour = "#333333") +   #p-cut-offl lines
  geom_segment(aes(x = 0.7, y = p.alert.orient, xend = 1.3, yend = p.alert.orient), 
               colour="#333333") + 
  geom_segment(aes(x = 2.7, y = p.control.orient, xend = 3.3, yend = p.control.orient),
               colour="#333333") + 
  labs(y = "Pearson's r") +
  theme(text=element_text(family="serif"),
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

#save plot
ggsave('violine.ant.png', violine.ant, width = 7, height = 4, units = 'in', dpi = 300)


###################################
####### Violine Plots Genes #######
###################################

#read data
ant.genes.corr.obs <- read.table("ant.genes.corr.obs.txt", sep = '\t',header = TRUE)

ant.exp.corr <- read.table(compare.exp.corr.txt", sep = '\t', header = TRUE)
ant.genes.corr.exp <- ant.exp.corr
colnames(ant.genes.corr.exp) <- c("exp.cor.genes.alert.orient",
                                  "exp.cor.genes.control.alert",
                                  "exp.cor.genes.control.orient")

#transform df to long format
exp.dat <- ant.genes.corr.exp 
exp.dat <- gather(exp.dat, condition, corr.exp, exp.cor.genes.alert.orient:exp.cor.genes.control.orient)

#observed correlations ANT x ANT 
alert.control.corr.obs <- cor(ant.genes.corr.obs$alert.corr.obs, ant.genes.corr.obs$control.corr.obs)
alert.orient.corr.obs <- cor(ant.genes.corr.obs$alert.corr.obs, ant.genes.corr.obs$orient.corr.obs)
orient.control.corr.obs <- cor(ant.genes.corr.obs$orient.corr.obs, ant.genes.corr.obs$control.corr.obs)

#calculate p-cutoff
#sort in ascending order and select value at .05 cut-off
p.cutoff.alert.control   <- sort(abs(ant.genes.corr.exp$exp.cor.genes.control.alert),
                          decreasing = TRUE)[round(0.05*nrow(ant.genes.corr.exp))] 

p.cutoff.alert.orient  <- sort(abs(ant.genes.corr.exp$exp.cor.genes.alert.orient),
                        decreasing = TRUE)[round(0.05*nrow(ant.genes.corr.exp))] 

p.cutoff.control.orient <- sort(abs(ant.genes.corr.exp$exp.cor.genes.control.orient),
                         decreasing = TRUE)[round(0.05*nrow(ant.genes.corr.exp))] 

# Basic violin
violine <- ggplot(exp.dat, aes(x=condition, y=corr.exp, fill=condition)) +   #violine
  scale_x_discrete(labels = c("Alert vs. Orient", "Alert vs. Control", "Control vs. Orient")) +
  scale_y_continuous(limits = c(-0.45, 1), breaks = seq(-0.4, 1, by = 0.2)) +
  geom_violin()+
  geom_point(aes(x="exp.cor.genes.alert.orient", y=alert.orient.corr.obs),
             col = "lightpink1", 
             shape = 20, 
             size = 5) +
  geom_point(aes(x="exp.cor.genes.control.alert", y=alert.control.corr.obs),     #observed corr
             col = "lightgreen", 
             shape = 20, 
             size = 5) +
  geom_point(aes(x="exp.cor.genes.control.orient", y=orient.control.corr.obs),
             col = "lightskyblue1", 
             shape = 20, 
             size = 5) +
  geom_segment(aes(x = 0.7, y = p.cutoff.alert.orient, xend = 1.3, yend = p.cutoff.alert.orient),
               colour = "#333333") +   #p-cut-offl lines
  geom_segment(aes(x = 1.7, y = p.cutoff.alert.control, xend = 2.3, yend = p.cutoff.alert.control), 
               colour="#333333") + 
  geom_segment(aes(x = 2.7, y = p.cutoff.control.orient, xend = 3.3, yend = p.cutoff.control.orient),
               colour="#333333") + 
  labs(y = "Pearson's r") +
  theme(text=element_text(size=18,family="serif"),
        plot.title = element_text(hjust = 0.5, face = "bold"), 
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 18, color = "black"),
        legend.position="none",
        legend.title=element_text(size=18),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

#save plot
ggsave("violine.genes.png", violine, width = 7, height = 4, units = 'in', dpi = 300)
 