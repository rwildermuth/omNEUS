# Data prep for NEUS MSE of BDN
# Created: 1/18/2021, Robert Wildermuth

library(tidyverse)

# Load Atlantis NEUS pseudo-data
load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusHistoricData.RData")

neusHistoricData


# Fit Fox-Gompertz model and get MSY -------------------------------------------

#Use lab 8 code instead
source('NEUS_MSE/lab8_functions.r') # loads write.datfile() and get.pars()
#load(file="C:/Users/rwildermuth/Dropbox/PhD_UMass/GB_BBN/MSYcalcs/thetaMLEs.RData")

# scan over range of initial biomass proportions of K and record estimates of MSY, Bmsy, and Biomass (and objective fxn val)
intThetas <- seq(0.05, 1, by=0.05)

# Need a conversion/sampling factor for GF and FF to help estimate catchability?

# Shellfish
SFobjFxn <- data.frame(theta=intThetas, objFxnVal=NA, convg=NA)

SFthetaMLEs <- list()

# assume a single survey
BObs <- cbind(neusHistoricData[, c("time", "inverts")], survey = 1) 
BObs <- BObs[, c("time", "survey", "inverts")]
for(i in 1:length(intThetas)){
  write.datfile(Nsp=1,BioObs=BObs,
                CatObs=as.data.frame(neusHistoricData[, c("time", "commSF")]), 
                zVal = 1.00005, #Syear=11, yrBreak1=26,yrBreak2=50,
                thetaPhase = -3, theta = intThetas[i], K=800000, r=0.05, outfile=paste0("NEUS_MSE/SF_",i,".dat"))
  command <- paste("C: & cd ",getwd(),"/NEUS_MSE & prodmodel -ind SF_",i,".dat", sep="")
  shell(command,wait=TRUE,invisible=TRUE)
  
  SFobjFxn[i,"objFxnVal"] <- as.numeric(scan(file="NEUS_MSE/prodmodel.par", what=character())[11])
  SFobjFxn[i,"convg"] <- as.numeric(scan(file="NEUS_MSE/prodmodel.par", what=character())[16])
  
  if(SFobjFxn[i,"convg"] < 1){
    res <- read.delim("NEUS_MSE/prodmodel.cor", skip = 1, sep="")
    res <- res[,1:4]
    SFthetaMLEs[[i]] <- res
  } else {
    SFthetaMLEs[[i]] <- "no convergence"
  }
  
}

SFobjFxn
plot(SFobjFxn[SFobjFxn$convg < 1, "theta"], 
     SFobjFxn[SFobjFxn$convg < 1, "objFxnVal"], type = "p", main = "Shellfish theta profile")

# edit theta to be estimated in .dat file with minimum 'objFxnVal' and re-estimate MLEs
command <- paste("C: & cd ",getwd(),"/NEUS_MSE & prodmodel -ind SF_5.dat",sep="") 
shell(command,wait=TRUE,invisible=TRUE)

# Evaluate model fit
par(mfrow=c(2,3))
resSF <- read.delim("NEUS_MSE/prodmodel_SF20210118_NoObsErr.cor", skip = 1, sep="")
resSF <- resSF[,1:4]
resSF$lower <- resSF$value * exp(-1.96*resSF$std.dev/resSF$value)
resSF$upper <- resSF$value * exp(1.96*resSF$std.dev/resSF$value)

plot(1:51, resSF[resSF$name %in% "Bio", "value"], type="l", main = "Estimated Biomass",
     ylim=range(resSF[resSF$name %in% "Bio", c("value", "lower", "upper")],
                neusHistoricData$inverts/resSF$value[7], na.rm = TRUE, 
                resSF[resSF$name %in% "Bmsy", "value"]), ylab = "Biomass (mt)")
lines(1:51, resSF[resSF$name %in% "Bio", "lower"], type="l", lty=2)
lines(1:51, resSF[resSF$name %in% "Bio", "upper"], type="l", lty=2)
points(1:50, neusHistoricData$inverts/resSF$value[7])
abline(h=resSF[resSF$name %in% "Bmsy", "value"], col=2)
abline(h=resSF[resSF$name %in% "Bmsy", "value"]*0.5, lty = 2, col=2)
abline(h=resSF[resSF$name %in% "K", "value"], col="darkblue")

# try plotting in log space
plot(1:51, log(resSF[resSF$name %in% "Bio", "value"]), type="l", main = "Estimated log-Biomass",
     ylim=range(log(resSF[resSF$name %in% "Bio", c("value", "lower", "upper")]),
                log(neusHistoricData$inverts) - log(resSF$value[7]), na.rm = TRUE), ylab = "Biomass (ln-mt)")
lines(1:51, log(resSF[resSF$name %in% "Bio", "lower"]), type="l", lty=2)
lines(1:51, log(resSF[resSF$name %in% "Bio", "upper"]), type="l", lty=2)
points(1:50, log(neusHistoricData$inverts) - log(resSF$value[7]))
abline(h=log(resSF[resSF$name %in% "Bmsy", "value"]), col=2)
abline(h=log(resSF[resSF$name %in% "Bmsy", "value"]*0.5), lty=2, col=2)
abline(h=log(resSF[resSF$name %in% "K", "value"]), col="darkblue")

plot(1:51, resSF[resSF$name %in% "ratioBmsy", "value"], type="l", main = "B/Bmsy",
     ylim=range(resSF[resSF$name %in% "ratioBmsy", c("value", "lower", "upper")], 1), ylab = "Ratio (B/Bmsy)")
lines(1:51, resSF[resSF$name %in% "ratioBmsy", "lower"], type="l", lty=2)
lines(1:51, resSF[resSF$name %in% "ratioBmsy", "upper"], type="l", lty=2)
abline(h = 1, col=2)
abline(h = 0.5, lty=2, col=2)
abline(h=resSF[resSF$name %in% "K", "value"]/resSF[resSF$name %in% "Bmsy", "value"], col="darkblue")

plot(1:50, neusHistoricData$commSF, type='l', main = "Catch", ylab = "Landings (mt)")
abline(h=resSF[resSF$name %in% "MSY", "value"], col=2)

residSF <- read.delim("NEUS_MSE/prodmodel_SF20210118_NoObsErr.rep", skip = 58, sep="", nrows = 50, header = FALSE)
names(residSF) <- c("Yr", "obsI", "estI")
plot(residSF$obsI - residSF$estI, type="p", main = "Residuals")
abline(h=0)
plot(residSF$estI, log(neusHistoricData$inverts), type="p", main = "QQ Plot")
abline(a=0,b=1)

# calculate the index as the probability of being below the minimum stock size threshold B_t/B_MSY = 0.5 
# as estimated from the MLE fit of the Fox-Gompertz model to the catch and survey data. 
fi <- read.delim("NEUS_MSE/prodmodel_SF20210118_NoObsErr.cor", skip = 1, sep="")
fi <- fi[,1:4]
fiMSY <- fi[fi$name == "MSY", "value"]
fi <- fi[fi$name == "ratioBmsy",]
fi$CV <- fi$std.dev/fi$value
for(i in 1:nrow(fi)){
  fi$pLess1[i] <- plnorm(0.5, mean = log(fi$value[i]), sd = fi$CV[i]) 
}


# Use Sean Lucey's catchability (q) values to create survey sample without error for groundfish and forage fish indices
load("C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/Species_codes.RData")

# Bring in lookup of Atlantis codes and names for BDN node indices
indBDNLookup <- read.csv("C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/NeusGroups_modforusewNEUS1.0.csv", 
                         stringsAsFactors = FALSE)
indBDNLookup <- indBDNLookup %>% dplyr::select(Code, Name, Long.Name, BNindexGroup) %>%
  filter(BNindexGroup != "")

# Get groundfish spp qs
indBDNLookup %>% filter(BNindexGroup == "groundfish")
# match Atlantis NEUS groups to qs using Link et al. 2011 tech memo (Table 2, pg 55)
spp %>% filter(grepl("HERRING", spp$COMNAME, fixed = TRUE))
spp %>% filter(grepl(toupper("Scophthalmus"), spp$SCINAME, fixed = TRUE))

qGF <- spp %>% filter(COMNAME %in%
                c("OFFSHORE HAKE", "SILVER HAKE", "WHITE HAKE", "RED HAKE", "SPOTTED HAKE", "ARMORED SEAROBIN", "NORTHERN SEAROBIN", 
                  "STRIPED SEAROBIN", "CONGER EEL", "AMERICAN EEL", "FOURBEARD ROCKLING", "THREEBEARD ROCKLING", "ATLANTIC COD", "ATLANTIC MOONFISH", 
                  "ATLANTIC CROAKER", "ATLANTIC NEEDLEFISH", "ATLANTIC SALMON", "ATLANTIC HALIBUT", "SUMMER FLOUNDER", "FOURSPOT FLOUNDER", 
                  "YELLOWTAIL FLOUNDER", "WINTER FLOUNDER", "WITCH FLOUNDER", "GULF STREAM FLOUNDER", "WINDOWPANE", "BLACK SEA BASS", "SEA RAVEN", 
                  "AMERICAN PLAICE", "ATLANTIC ANGEL SHARK", "NURSE SHARK", "PORBEAGLE SHARK", "BULL SHARK", "LEMON SHARK", "ATLANTIC SHARPNOSE SHARK",
                  "SCALLOPED HAMMERHEAD SHARK", "THRESHER SHARK", "SANDBAR SHARK", "SAND TIGER", "MOUSTACHE SCULPIN", "SHORTHORN SCULPIN", 
                  "LONGHORN SCULPIN", "BARNDOOR SKATE", "WINTER SKATE", "CLEARNOSE SKATE", "ROSETTE SKATE", "LITTLE SKATE", "SMOOTH SKATE", 
                  "THORNY SKATE", "ROUGHTAIL STINGRAY", "COWNOSE RAY", "SPOTTED EAGLE RAY", "BULLNOSE RAY", "BLUNTNOSE STINGRAY", "YELLOW STINGRAY",
                  "ATLANTIC STINGRAY", "OYSTER TOADFISH", "LARGESCALE LIZARDFISH", "SHORTJAW LIZARDFISH", "INSHORE LIZARDFISH", "OFFSHORE LIZARDFISH", 
                  "RED LIZARDFISH", "ACADIAN REDFISH", "GOOSEFISH", "TILEFISH", "LONGSPINE SNIPEFISH", "NORTHERN PIPEFISH", "BLACKBELLY ROSEFISH",
                  "SMOOTH DOGFISH", "CHAIN DOGFISH", "SPINY DOGFISH", "NORTHERN KINGFISH", "ATLANTIC HAGFISH", "HOGFISH", "WEAKFISH", 
                  "LONGSPINE SCORPIONFISH", "SCORPIONFISH AND ROCKFISH UNCL", "TRUMPETFISH", "LUMPFISH", "HADDOCK", "POLLOCK", "ATLANTIC WOLFFISH", 
                  "HOGCHOKER", "CUSK", "GRENADIER UNCL", "JOHN DORY", "THREESPINE STICKLEBACK", "LOOKDOWN", "SCUP", 'SILVER PERCH', "SAND PERCH", 
                  "STRIPED BASS", 'BLACK DRUM', "SPOT", "TAUTOG", "CUNNER", 'NORTHERN PUFFER', "OCEAN POUT", "WRYMOUTH", "SHEEPSHEAD", 
                  "ATLANTIC STURGEON", "ATLANTIC TOMCOD", 'ATLANTIC TORPEDO')) %>%
          dplyr::select(COMNAME, SCINAME, Fall.q, Spring.q)
any(is.na(qGF))
qGF <- mean(c(qGF$Fall.q, qGF$Spring.q))

# Groundfish
GFobjFxn <- data.frame(theta=intThetas, objFxnVal=NA, convg=NA)

GFthetaMLEs <- list()

# assume a single survey
BObs <- cbind(neusHistoricData[, c("time", "ground")], survey = 1) 
BObs <- BObs[, c("time", "survey", "ground")]
# apply catchability and 75% survey coverage to biomass timeseries
BObs <- BObs %>% mutate(ground = ground * qGF * 0.75)

for(i in 1:length(intThetas)){
  write.datfile(Nsp=1,BioObs=BObs,
                CatObs=as.data.frame(neusHistoricData[, c("time", "commGF")]), 
                zVal = 1.00005, #Syear=11, yrBreak1=26,yrBreak2=50,
                thetaPhase = -3, theta = intThetas[i], K=800000, r=0.2, outfile=paste0("NEUS_MSE/GF_",i,".dat"))
  command <- paste("C: & cd ",getwd(),"/NEUS_MSE & prodmodel -ind GF_",i,".dat", sep="")
  shell(command,wait=TRUE,invisible=TRUE)
  
  GFobjFxn[i,"objFxnVal"] <- as.numeric(scan(file="NEUS_MSE/prodmodel.par", what=character())[11])
  GFobjFxn[i,"convg"] <- as.numeric(scan(file="NEUS_MSE/prodmodel.par", what=character())[16])
  
  if(GFobjFxn[i,"convg"] < 1){
    res <- read.delim("NEUS_MSE/prodmodel.cor", skip = 1, sep="")
    res <- res[,1:4]
    GFthetaMLEs[[i]] <- res
  } else {
    GFthetaMLEs[[i]] <- "no convergence"
  }
  
}

GFobjFxn
plot(GFobjFxn[GFobjFxn$convg < 1, "theta"], 
     GFobjFxn[GFobjFxn$convg < 1, "objFxnVal"], type = "p", main = "Groundfish theta profile")

# edit theta to be estimated in .dat file with minimum 'objFxnVal' and re-estimate MLEs
command <- paste("C: & cd ",getwd(),"/NEUS_MSE & prodmodel -ind GF_8.dat",sep="") # better convirgence than 3 param where init theta = 0.45
shell(command,wait=TRUE,invisible=TRUE)

# Evaluate model fit
par(mfrow=c(2,3))
resGF <- read.delim("NEUS_MSE/prodmodel_GF20210120_NoObsErr.cor", skip = 1, sep="")
resGF <- resGF[,1:4]
resGF$lower <- resGF$value * exp(-1.96*resGF$std.dev/resGF$value)
resGF$upper <- resGF$value * exp(1.96*resGF$std.dev/resGF$value)

plot(1:51, resGF[resGF$name %in% "Bio", "value"], type="l", main = "Estimated Biomass",
     ylim=range(resGF[resGF$name %in% "Bio", c("value", "lower", "upper")],
                BObs$ground/resGF$value[7], na.rm = TRUE, 
                resGF[resGF$name %in% "Bmsy", "value"]), ylab = "Biomass (mt)")
lines(1:51, resGF[resGF$name %in% "Bio", "lower"], type="l", lty=2)
lines(1:51, resGF[resGF$name %in% "Bio", "upper"], type="l", lty=2)
points(1:50, BObs$ground/resGF$value[7])
abline(h=resGF[resGF$name %in% "Bmsy", "value"], col=2)
abline(h=resGF[resGF$name %in% "Bmsy", "value"]*0.5, lty = 2, col=2)
abline(h=resGF[resGF$name %in% "K", "value"], col="darkblue")

# try plotting in log space
plot(1:51, log(resGF[resGF$name %in% "Bio", "value"]), type="l", main = "Estimated log-Biomass",
     ylim=range(log(resGF[resGF$name %in% "Bio", c("value", "lower", "upper")]),
                log(BObs$ground) - log(resGF$value[7]), na.rm = TRUE), ylab = "Biomass (ln-mt)")
lines(1:51, log(resGF[resGF$name %in% "Bio", "lower"]), type="l", lty=2)
lines(1:51, log(resGF[resGF$name %in% "Bio", "upper"]), type="l", lty=2)
points(1:50, log(BObs$ground) - log(resGF$value[7]))
abline(h=log(resGF[resGF$name %in% "Bmsy", "value"]), col=2)
abline(h=log(resGF[resGF$name %in% "Bmsy", "value"]*0.5), lty=2, col=2)
abline(h=log(resGF[resGF$name %in% "K", "value"]), col="darkblue")

plot(1:51, resGF[resGF$name %in% "ratioBmsy", "value"], type="l", main = "B/Bmsy",
     ylim=range(resGF[resGF$name %in% "ratioBmsy", c("value", "lower", "upper")], 1), ylab = "Ratio (B/Bmsy)")
lines(1:51, resGF[resGF$name %in% "ratioBmsy", "lower"], type="l", lty=2)
lines(1:51, resGF[resGF$name %in% "ratioBmsy", "upper"], type="l", lty=2)
abline(h = 1, col=2)
abline(h = 0.5, lty=2, col=2)
abline(h=resGF[resGF$name %in% "K", "value"]/resGF[resGF$name %in% "Bmsy", "value"], col="darkblue")

plot(1:50, neusHistoricData$commGF, type='l', main = "Catch", ylab = "Landings (mt)")
abline(h=resGF[resGF$name %in% "MSY", "value"], col=2)

residGF <- read.delim("NEUS_MSE/prodmodel_GF20210120_NoObsErr.rep", skip = 58, sep="", nrows = 50, header = FALSE)
names(residGF) <- c("Yr", "obsI", "estI")
plot(residGF$obsI - residGF$estI, type="p", main = "Residuals")
abline(h=0)
plot(residGF$estI, log(BObs$ground), type="p", main = "QQ Plot")
abline(a=0,b=1)

# calculate the index as the probability of being below the minimum stock size threshold B_t/B_MSY = 0.5 
# as estimated from the MLE fit of the Fox-Gompertz model to the catch and survey data. 
gf <- read.delim("NEUS_MSE/prodmodel_GF20210120_NoObsErr.cor", skip = 1, sep="")
gf <- gf[,1:4]
gfMSY <- gf[gf$name == "MSY", "value"]
gf <- gf[gf$name == "ratioBmsy",]
gf$CV <- gf$std.dev/gf$value
for(i in 1:nrow(gf)){
  gf$pLess1[i] <- plnorm(0.5, mean = log(gf$value[i]), sd = gf$CV[i]) 
}


# Forage fish
indBDNLookup %>% filter(BNindexGroup == "forage")

qFF <- spp %>% filter(COMNAME %in%
                c("ATLANTIC MACKEREL", "ATLANTIC HERRING", "ATLANTIC SILVERSIDE", "AMERICAN SHAD", "BLUEBACK HERRING", "ROUND HERRING",
                  "BLUEBACK HERRING", "HICKORY SHAD", "GIZZARD SHAD", "ATLANTIC ARGENTINE", "BAY ANCHOVY", "STRIPED ANCHOVY",
                  "BUTTERFISH", "NORTHERN SAND LANCE", "HARVESTFISH", "CHUB MACKEREL", "SILVERSTRIPE HALFBEAK", "FLYING HALFBEAK",
                  "ALEWIFE", "MENHADEN", "RAINBOW SMELT", "NORTHERN SHORTFIN SQUID", "LONGFIN SQUID")) %>%
                dplyr::select(COMNAME, SCINAME, Fall.q, Spring.q)
any(is.na(qFF))
qFF <- mean(c(qFF$Fall.q, qFF$Spring.q))


FFobjFxn <- data.frame(theta=intThetas, objFxnVal=NA, convg=NA)

FFthetaMLEs <- list()

# assume a single survey
BObs <- cbind(neusHistoricData[, c("time", "forage")], survey = 1) 
BObs <- BObs[, c("time", "survey", "forage")]
# apply catchability and 75% survey coverage to biomass timeseries
BObs <- BObs %>% mutate(forage = forage * qFF * 0.75)

for(i in 1:length(intThetas)){
  write.datfile(Nsp=1,BioObs=BObs,
                CatObs=as.data.frame(neusHistoricData[, c("time", "commPel")]), 
                zVal = 1.00005, #Syear=11, yrBreak1=26,yrBreak2=50,
                thetaPhase = -3, theta = intThetas[i], K=500000, r=0.35, outfile=paste0("NEUS_MSE/FF_",i,".dat"))
  command <- paste("C: & cd ",getwd(),"/NEUS_MSE & prodmodel -ind FF_",i,".dat", sep="")
  shell(command,wait=TRUE,invisible=TRUE)
  
  FFobjFxn[i,"objFxnVal"] <- as.numeric(scan(file="NEUS_MSE/prodmodel.par", what=character())[11])
  FFobjFxn[i,"convg"] <- as.numeric(scan(file="NEUS_MSE/prodmodel.par", what=character())[16])
  
  if(FFobjFxn[i,"convg"] < 1){
    res <- read.delim("NEUS_MSE/prodmodel.cor", skip = 1, sep="")
    res <- res[,1:4]
    FFthetaMLEs[[i]] <- res
  } else {
    FFthetaMLEs[[i]] <- "no convergence"
  }
  
}

FFobjFxn
plot(FFobjFxn[FFobjFxn$convg < 1, "theta"], 
     FFobjFxn[FFobjFxn$convg < 1, "objFxnVal"], type = "p", main = "Forage fish theta profile")

# edit theta to be estimated in .dat file with minimum 'objFxnVal' and re-estimate MLEs
command <- paste("C: & cd ",getwd(),"/NEUS_MSE & prodmodel -ind FF_3.dat",sep="") 
shell(command,wait=TRUE,invisible=TRUE)

# Evaluate model fit
par(mfrow=c(2,3))
resFF <- read.delim("NEUS_MSE/prodmodel_FF20210123_NoObsErr.cor", skip = 1, sep="")
resFF <- resFF[,1:4]
resFF$lower <- resFF$value * exp(-1.96*resFF$std.dev/resFF$value)
resFF$upper <- resFF$value * exp(1.96*resFF$std.dev/resFF$value)

plot(1:51, resFF[resFF$name %in% "Bio", "value"], type="l", main = "Estimated Biomass",
     ylim=range(resFF[resFF$name %in% "Bio", c("value", "lower", "upper")],
                BObs$forage/resFF$value[7], na.rm = TRUE, 
                resFF[resFF$name %in% "Bmsy", "value"]), ylab = "Biomass (mt)")
lines(1:51, resFF[resFF$name %in% "Bio", "lower"], type="l", lty=2)
lines(1:51, resFF[resFF$name %in% "Bio", "upper"], type="l", lty=2)
points(1:50, BObs$forage/resFF$value[7])
abline(h=resFF[resFF$name %in% "Bmsy", "value"], col=2)
abline(h=resFF[resFF$name %in% "Bmsy", "value"]*0.5, lty = 2, col=2)
abline(h=resFF[resFF$name %in% "K", "value"], col="darkblue")

# try plotting in log space
plot(1:51, log(resFF[resFF$name %in% "Bio", "value"]), type="l", main = "Estimated log-Biomass",
     ylim=range(log(resFF[resFF$name %in% "Bio", c("value", "lower", "upper")]),
                log(BObs$forage) - log(resFF$value[7]), na.rm = TRUE), ylab = "Biomass (ln-mt)")
lines(1:51, log(resFF[resFF$name %in% "Bio", "lower"]), type="l", lty=2)
lines(1:51, log(resFF[resFF$name %in% "Bio", "upper"]), type="l", lty=2)
points(1:50, log(BObs$forage) - log(resFF$value[7]))
abline(h=log(resFF[resFF$name %in% "Bmsy", "value"]), col=2)
abline(h=log(resFF[resFF$name %in% "Bmsy", "value"]*0.5), lty=2, col=2)
abline(h=log(resFF[resFF$name %in% "K", "value"]), col="darkblue")

plot(1:51, resFF[resFF$name %in% "ratioBmsy", "value"], type="l", main = "B/Bmsy",
     ylim=range(resFF[resFF$name %in% "ratioBmsy", c("value", "lower", "upper")], 1), ylab = "Ratio (B/Bmsy)")
lines(1:51, resFF[resFF$name %in% "ratioBmsy", "lower"], type="l", lty=2)
lines(1:51, resFF[resFF$name %in% "ratioBmsy", "upper"], type="l", lty=2)
abline(h = 1, col=2)
abline(h = 0.5, lty=2, col=2)
abline(h=resFF[resFF$name %in% "K", "value"]/resFF[resFF$name %in% "Bmsy", "value"], col="darkblue")

plot(1:50, neusHistoricData$commPel, type='l', main = "Catch", ylab = "Landings (mt)")
abline(h=resFF[resFF$name %in% "MSY", "value"], col=2)

residFF <- read.delim("NEUS_MSE/prodmodel_FF20210123_NoObsErr.rep", skip = 58, sep="", nrows = 50, header = FALSE)
names(residFF) <- c("Yr", "obsI", "estI")
plot(residFF$obsI - residFF$estI, type="p", main = "Residuals")
abline(h=0)
plot(residFF$estI, log(BObs$forage), type="p", main = "QQ Plot")
abline(a=0,b=1)

# calculate the index as the probability of being below the minimum stock size threshold B_t/B_MSY = 0.5 
# as estimated from the MLE fit of the Fox-Gompertz model to the catch and survey data. 
ff <- read.delim("NEUS_MSE/prodmodel_FF20210123_NoObsErr.cor", skip = 1, sep="")
ff <- ff[,1:4]
ffMSY <- ff[ff$name == "MSY", "value"]
ff <- ff[ff$name == "ratioBmsy",]
ff$CV <- ff$std.dev/ff$value
for(i in 1:nrow(ff)){
  ff$pLess1[i] <- plnorm(0.5, mean = log(ff$value[i]), sd = ff$CV[i]) 
}


# Colate and scale dataset for Netica -------------------------------------
# use probability of being less than 1/2 B/Bmsy as index for managed groups w/out last year estimate
bdnMSEData <- neusHistoricData %>% mutate(inverts = fi$pLess1[-51],
                                          forage = ff$pLess1[-51],
                                          ground = gf$pLess1[-51])

# Log-transform biomass  and count indices
bdnMSEData$gfMAB <- log(bdnMSEData$gfMAB)
bdnMSEData$GZ <- log(bdnMSEData$GZ)
bdnMSEData$copepod <- log(bdnMSEData$copepod)
bdnMSEData$benthos <- log(bdnMSEData$benthos)
bdnMSEData$commPel <- log(bdnMSEData$commPel)
bdnMSEData$commGF <- log(bdnMSEData$commGF)
bdnMSEData$commSF <- log(bdnMSEData$commSF)
bdnMSEData$recParticip <- log(bdnMSEData$recParticip)

#write.csv(bdnMSEData, "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/GB_BDN_MSE_data_20210125.csv", row.names = FALSE)
bdnMSEData <- read.csv("C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/GB_BDN_MSE_data_20210125.csv",
                      stringsAsFactors = FALSE)

# Scale and center raw data indices (not anomalies, percentages, or categorical variables)
# Need to create anomalies for Temp and Sal data and copepods
bdnMSEData$BT <- scale(bdnMSEData$BT, center = TRUE, scale = TRUE)
bdnMSEData$BS <- scale(bdnMSEData$BS, center = TRUE, scale = TRUE)
bdnMSEData$SST <- scale(bdnMSEData$SST, center = TRUE, scale = TRUE)
bdnMSEData$SSS <- scale(bdnMSEData$SSS, center = TRUE, scale = TRUE)
bdnMSEData$copepod <- scale(bdnMSEData$copepod, center = TRUE, scale = TRUE)
bdnMSEData$gfMAB <- scale(bdnMSEData$gfMAB, center = TRUE, scale = TRUE)
bdnMSEData$GZ <- scale(bdnMSEData$GZ, center = TRUE, scale = TRUE)
bdnMSEData$strat <- scale(bdnMSEData$strat, center = TRUE, scale = TRUE)
bdnMSEData$habPel <- scale(bdnMSEData$habPel, center = TRUE, scale = TRUE)
bdnMSEData$habDem <- scale(bdnMSEData$habDem, center = TRUE, scale = TRUE)
bdnMSEData$PP <- scale(bdnMSEData$PP, center = TRUE, scale = TRUE)
bdnMSEData$detBac <- scale(bdnMSEData$detBac, center = TRUE, scale = TRUE)
bdnMSEData$benthos <- scale(bdnMSEData$benthos, center = TRUE, scale = TRUE)
bdnMSEData$commPel <- scale(bdnMSEData$commPel, center = TRUE, scale = TRUE)
bdnMSEData$commGF <- scale(bdnMSEData$commGF, center = TRUE, scale = TRUE)
bdnMSEData$commSF <- scale(bdnMSEData$commSF, center = TRUE, scale = TRUE)
bdnMSEData$recParticip <- scale(bdnMSEData$recParticip, center = TRUE, scale = TRUE)

# Need to turn habNear into percentage change data
bdnMSEData$habNearPct <- NA
for(i in 2:nrow(bdnMSEData)){
  bdnMSEData$habNearPct[i] <- (((bdnMSEData$habNear[i]/bdnMSEData$habNear[i-1])^(1/2))-1)*100
}
bdnMSEData %>% dplyr::select(time, habNear, habNearPct)

bdnMSEData$habNear <- bdnMSEData$habNearPct

nodes <- names(bdnMSEData[,-1])
names(bdnMSEData) <- c("Year",paste0(nodes, "2")) # For t=2 period

# Include columns with 1-year lag offsets
lag1 <- rbind(NA, bdnMSEData[,-1])
names(lag1) <- paste0(nodes, "1")
bdnMSEData <- rbind(bdnMSEData, NA)
bdnMSEData <- cbind(bdnMSEData,lag1)


# write.table(bdnMSEData, "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/GB_BDN_MSE_data_20210125.txt", sep="\t", row.names = FALSE)


# Scaled thresholds:
# MSY estimates
sfScale <- scale(log(neusHistoricData$commSF), center = TRUE, scale = TRUE)
(log(fiMSY)-attr(sfScale, "scaled:center"))/attr(sfScale, "scaled:scale")

ffScale <- scale(log(neusHistoricData$commPel), center = TRUE, scale = TRUE)
(log(ffMSY)-attr(ffScale, "scaled:center"))/attr(ffScale, "scaled:scale")

gfScale <- scale(log(neusHistoricData$commGF), center = TRUE, scale = TRUE)
(log(gfMSY)-attr(gfScale, "scaled:center"))/attr(gfScale, "scaled:scale")
