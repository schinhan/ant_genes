###############################################
############ GO-Figures self-made #############
###############################################

library(ggplot2)
library(RColorBrewer)
library(autoimage)

#set wd to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#read tsv
alert.data <- read.table(file = 'molecular_function_alerting_results.tsv', sep = '\t', header = TRUE)
alert.df <- alert.data
control.data <- read.table(file = 'molecular_function_control_results.tsv', sep = '\t', header = TRUE)
control.df <- control.data
orient.data <- read.table(file = 'molecular_function_orienting_results.tsv', sep = '\t', header = TRUE)
orient.df <- orient.data
df <- rbind.data.frame(alert.df, control.df, orient.df)

#define breaks for color scale
pvals <- c(alert.df$pval, orient.df$pval, control.df$pval)
ggplot2:::breaks(c(min(pvals),max(pvals)),"width",n = 6) #calculate breaks

#set categories for continous color scale 
df$p <- cut(df$pval,
                     breaks=c(-Inf, 1e-10, 1e-08, 1e-06, 1e-04, Inf),
                     labels=c('<1e-10', '1e-10 - 1e-08', '1e-08 - 1e-06', '1e-06 - 1e-04', '>1e-04'))

#make plot 
ggplot(df, aes(x, y)) +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) + #set axis limits and breaks
  scale_y_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) +
  geom_point(aes(colour = p),
             size = 8*df$size)+
  scale_color_manual(values=c("#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7"),
                     drop=FALSE) +
  labs(x = "Semantic Space X",
       y = "Semantic Space Y") +
  ggtitle("Alerting") +
  geom_text(aes(label=dotCount), 
            size=3.5,
            family = "serif")+
  theme(text=element_text(family="serif"),
        plot.title = element_text(size=16, hjust = 0.5, face = "bold"), 
        axis.title.x=element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(color = guide_legend(override.aes = list(size=6)),
         fill=guide_legend(title="New Legend Title"))

#Alerting
alert.df$p <- cut(alert.df$pval,
                  breaks=c(-Inf, 1e-10, 1e-08, 1e-06, 1e-04, Inf),
                  labels=c('<1e-10', '1e-10 - 1e-08', '1e-08 - 1e-06', '1e-06 - 1e-04', '>1e-04'))

alert.gofigure <-
ggplot(alert.df, aes(x, y)) +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) + #set axis limits and breaks
  scale_y_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) +
  geom_point(aes(colour = p),
             size = 8*alert.df$size)+
  scale_color_manual(values=c("#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7"),
                     drop=FALSE) +
  labs(x = "Semantic Space X",
       y = "Semantic Space Y") +
  ggtitle("Alerting, si = 0.5") +
  geom_text(aes(label=dotCount), 
            size=3.5,
            family = "serif")+
  theme(text=element_text(family="serif"),
        plot.title = element_text(size=16, hjust = 0.5, face = "bold"), 
        axis.title.x=element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(color = guide_legend(override.aes = list(size=6)))

#Orienting
orient.df$p <- cut(orient.df$pval,
                  breaks=c(-Inf, 1e-10, 1e-08, 1e-06, 1e-04, Inf),
                  labels=c('<1e-10', '1e-10 - 1e-08', '1e-08 - 1e-06', '1e-06 - 1e-04', '>1e-04'))


orient.gofigure <-
  ggplot(orient.df, aes(x, y)) +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) + #set axis limits and breaks
  scale_y_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) +
  geom_point(aes(colour = p),
             size = 8*orient.df$size)+
  scale_color_manual(values=c("#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7"),
                     drop=FALSE) +
  labs(x = "Semantic Space X",
       y = "Semantic Space Y") +
  ggtitle("Orienting, si = 0.5") +
  geom_text(aes(label=dotCount), 
            size=3.5,
            family = "serif")+
  theme(text=element_text(family="serif"),
        plot.title = element_text(size=16, hjust = 0.5, face = "bold"), 
        axis.title.x=element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(color = guide_legend(override.aes = list(size=6)))

#Control
control.df$p <- cut(control.df$pval,
                   breaks=c(-Inf, 1e-10, 1e-08, 1e-06, 1e-04, Inf),
                   labels=c('<1e-10', '1e-10 - 1e-08', '1e-08 - 1e-06', '1e-06 - 1e-04', '>1e-04'))


control.gofigure <- 
ggplot(control.df, aes(x, y)) +
  scale_x_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) + #set axis limits and breaks
  scale_y_continuous(limits = c(-1.2, 1.2), breaks = seq(-1, 1, by = 0.5)) +
  geom_point(aes(colour = p),
             size = 8*control.df$size)+
  scale_color_manual(values=c("#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7"),
                     drop=FALSE) +
  labs(x = "Semantic Space X",
       y = "Semantic Space Y") +
  ggtitle("Control, si = 0.5") +
  geom_text(aes(label=dotCount), 
            size=3.5,
            family = "serif")+
  theme(text=element_text(family="serif"),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        axis.title.x=element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), 
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  guides(color = guide_legend(override.aes = list(size=6)))

#save plots
ggsave('alert.gofigure.png', alert.gofigure, width = 1525, height = 1125, units = 'px')
ggsave('orient.gofigure.png', orient.gofigure, width = 1525, height = 1125, units = 'px')
ggsave('control.gofigure.png', control.gofigure, width = 1525, height = 1125, units = 'px')

