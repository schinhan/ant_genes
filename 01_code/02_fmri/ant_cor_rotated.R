# Calculate correlations between observed and randomly rotated ANT maps 

# required packages: data.table

# set wd to data folder
setwd("../../02_data/03_fMRI/")

# read data
library(data.table)
ant_data = read.csv(file = 'ant_lausanne.csv')
perm_id_data = data.frame(fread('/../../02_data/02_parcellation/Lausanne_perm_ids.txt.gz'))

# calculate observed correlations
alert.control.corr.obs <- cor(ant_data$alert, ant_data$control)
alert.orient.corr.obs <- cor(ant_data$alert, ant_data$orient)
orient.control.corr.obs <- cor(ant_data$orient, ant_data$control)


###### rotate ANT activation #####
perm.alert <- as.data.frame(1:219) #create empty data frame

for(i in 1:10) {       # loop over cols to permutate ant act
  rot = ant_data$alert[order(perm_id_data[,i])]
  perm.alert = cbind(perm.alert, rot)
}
perm.alert <- perm.alert[,-1] #delete index col


#do this for control
perm.control <- as.data.frame(1:219) #create empty data frame
for(i in 1:10000) {       # loop over cols to permutate ant act
  rot = ant_data$control[order(perm_id_data[,i])]
  perm.control = cbind(perm.control, rot)
}
perm.control <- perm.control[,-1] #delete index col

#do this for orient
perm.orient <- as.data.frame(1:219) #create empty data frame
for(i in 1:10000) {       # loop over cols to permutate ant act
  rot = ant_data$orient[order(perm_id_data[,i])]
  perm.orient = cbind(perm.orient, rot)
}
perm.orient <- perm.orient[,-1] #delete index col

## correlate observed ant_dat from alerting with perm_dat from control
#alert.obs vs. control.perm
alert.obs.control.perm <- vector("numeric", ncol(perm.control))    # create empty vector 
for(i in 1:ncol(perm.control)) {       # loop over cols to permutate ant act
  alert.obs.control.perm[i] = cor(ant_data$alert, perm.control[,i])
}

#control.obs vs. alert.perm
control.obs.alert.perm <- vector("numeric", ncol(perm.alert))    # create empty vector 
for(i in 1:ncol(perm.alert)) {       # loop over cols to permutate ant act
  control.obs.alert.perm[i] = cor(ant_data$control, perm.alert[,i])
}

#alert.obs vs. orient.perm
alert.obs.orient.perm <- vector("numeric", ncol(perm.orient))    # create empty vector 
for(i in 1:ncol(perm.orient)) {       # loop over cols to permutate ant act
  alert.obs.orient.perm[i] = cor(ant_data$alert, perm.orient[,i])
}

#orient vs. alert.perm
orient.obs.alert.perm <- vector("numeric", ncol(perm.orient))    # create empty vector 
for(i in 1:ncol(perm.orient)) {       # loop over cols to permutate ant act
  orient.obs.alert.perm[i] = cor(ant_data$orient, perm.alert[,i])
}

#control.obs vs. orient.perm
control.obs.orient.perm <- vector("numeric", ncol(perm.orient))    # create empty vector 
for(i in 1:ncol(perm.orient)) {       # loop over cols to permutate ant act
  control.obs.orient.perm[i] = cor(ant_data$control, perm.orient[,i])
}

#orient vs. control.perm
orient.obs.control.perm <- vector("numeric", ncol(perm.control))    # create empty vector 
for(i in 1:ncol(perm.control)) {       # loop over cols to permutate ant act
  orient.obs.control.perm[i] = cor(ant_data$orient, perm.control[,i])
}


#################
##### pvals #####
#################

## alerting ###
p.alert.obs.control.perm <- sum(abs(alert.obs.control.perm) >= abs(alert.control.corr.obs))/length(alert.obs.control.perm)
p.alert.obs.orient.perm <- sum(abs(alert.obs.orient.perm) >= abs(alert.orient.corr.obs))/length(alert.obs.orient.perm)

## control ##
p.control.obs.alert.perm <- sum(abs(control.obs.alert.perm) >= abs(alert.control.corr.obs))/length(control.obs.alert.perm)
p.control.obs.orient.perm <- sum(abs(control.obs.orient.perm) >= abs(orient.control.corr.obs))/length(control.obs.orient.perm)

## orient ##
p.orient.obs.alert.perm <- sum(abs(orient.obs.alert.perm) >= abs(alert.orient.corr.obs))/length(orient.obs.alert.perm)
p.orient.obs.control.perm <- sum(abs(orient.obs.control.perm) >= abs(orient.control.corr.obs))/length(orient.obs.control.perm)




