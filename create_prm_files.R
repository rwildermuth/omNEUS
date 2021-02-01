library(stringr)
source("change_param.R")

outfile <- "newtable.out"
origfile <- "harvest_sandbox.prm"

param.name <- "mFC_"
name2 <- NULL
param.type=2
newval = 2  #multiplier on original F
change.prm3(path,origfile,param.name,outfile,mapfile,param.type,name2, newval, type = 1)

#fishing scenarios
param.name <- ""
name2 <- "_mFC_changes"
param.type=2
newval = 1
change.prm3(path,origfile,param.name,outfile,mapfile,param.type,name2, newval,type = 2)

param.name <- "mFCchange_start_"
name2 <- ""
param.type=2
newval = 18251
change.prm3(path,origfile,param.name,outfile,mapfile,param.type,name2, newval,type = 2)

param.name <- "mFCchange_mult_"
name2 <- ""
param.type=2
newval = 0.5
change.prm3(path,origfile,param.name,outfile,mapfile,param.type,name2, newval,type = 2)
newval = 2.0
change.prm3(path,origfile,param.name,outfile,mapfile,param.type,name2, newval,type = 2)

#wind
param.name <- "mFC_"
name2 <- NULL
param.type=2
newval = 2  #multiplier on original F
fishmort <- change.prm3(path,origfile,param.name,outfile,mapfile,param.type,name2, newval, type = 3)

wind <- data.frame(groups = c("BD", "BC", "BFS", "BFF", "BMS", "BML", "BG"),
                   add_f = c(0.000000005,
                             0.000004,
                             0.000000015,
                             0.000000025,
                             0.00000015,
                             0.0000001,
                             0.00000005))  #0.05 * linear mortality
xx <- fishmort$output2[match(wind$groups,fishmort$Child),]

yy <- rep(0,nrow(xx))
for (i in 1:nrow(xx)) {
  yy[i] <- (xx[i,xx[i,]>0]+wind$add_f[i])/xx[i,xx[i,]>0]
}
# new F multipliers for the wind
#yy
#[1] 1.051043       NA 1.003655 1.000115 1.014543 1.002589 1.000151

#BC
(0.0000000000001+0.000004)/0.0000000000001
# 4e+07

#mod _mFC_changes by hand
#also "mFCchange_start_" 18251
#use the new multipliers for "mFCchange_mult_"

