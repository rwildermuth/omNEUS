
library(LoopAnalyst)
library(tidyverse)
library(igraph)

#---------------------------------------------------------
# Author: Robert Wildermuth, rwildermuth@umassd.edu
# Created: 2/24/2016
# Last Modified: 1/12/2021

# Description:
# Imports the adjacency matrix for the Georges Bank conceptual model 2.0
# and adapts the data to a community matrix for loop analysis. 
# This code also performs correlation analysis on Atlantis pseudo-data
# to modify community links based on observed correlations.
# After the community matrix ('commMat') is created, loops and 
# paths are defined and enumerated, feedback and stability are evaluated, 
# and press perturbation is performed using the make.adjoint() function.
# The total feedback and weighted community matrices are also calculated 
# to determine effect determinacy (essentially number of + or - effects 
# devided by the total number of paths possible). All matrices are then 
# saved to file.

#---------------------------------------------------------



# Load the GB adjacency matrix
mmMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/QNMproject/GB/GB_QNM_adjacencyMat.csv', 
                  header = TRUE)

rNames <- names(mmMat)[-1]
rNames <- gsub(".", "", x = rNames, fixed = TRUE)
mmMat <- as.matrix(mmMat[, -1])
# 31 x 31 matrix
rownames(mmMat) <- rNames

# Translate into a (-1,0,1) community matrix
# intermediate function to assess values for positive or negative influence
influenceFxn <- function(x){
  if(is.na(x)){
    x <- 0
  } else if(x < 0){
    x <- -1
  } else if(x > 0){
    x <- 1
  } else if(x == 0){
    x <- 0
  } else {
    x <- NA # double check that all values are being captured
  }
  return(x)
}

test2 <- sapply(mmMat, influenceFxn)
class(test2)

# Must transpose to keep actors in columns and receivers in rows
commMat <- matrix(test2, nrow = nrow(mmMat), ncol = ncol(mmMat), byrow=TRUE)
rownames(commMat) <- rNames
colnames(commMat) <- rNames
# 31 x 31 matrix


# Add self-limitation to each component
for(i in 1:nrow(commMat)){
  commMat[i,i] <- -1
}
any(is.na(commMat))
# 31 x 31 matrix


# Atlantis pseudo-data correlation analysis -------------------------------

load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusHistoricData.RData")

# Spearman correlations without categorical PS index
atlGBCorr <- neusHistoricData %>% select( -c("time", "PS")) %>% cor(method = "spearman")

graphGB <- graph.adjacency(commMat)
edgesGB <- as.data.frame(get.edgelist(graphGB), stringsAsFactors = FALSE)
colnames(edgesGB) <- c("commMat_To", "commMat_From")

# match community matrix names with Atlantis/BBN names
colnames(atlGBCorr)
namesLookupGB <- data.frame(commMatNames = rNames, 
                            bnNames = c(NA, NA, "SST",
                                        "SSS", "BT", "BS",
                                        "strat", "habPel", "habNear",
                                        "habDem", NA, NA,
                                        "forage", "ground", "inverts",          
                                        "copepod", "PP", "detBac",             
                                        "commGF", "PS", "benthos",
                                        "GZ", "commPel", "commSF",   
                                        NA, "recParticip", NA,         
                                        "gfMAB", NA, NA, NA),
                            stringsAsFactors = FALSE)

edgesGB$corrMat_To <- edgesGB$corrMat_From <- NA

for(i in 1:nrow(edgesGB)){
  edgesGB[i, "corrMat_From"] <-  namesLookupGB[which(namesLookupGB$commMatNames == edgesGB[i, "commMat_From"]), 
                                               "bnNames"]
  edgesGB[i, "corrMat_To"] <-  namesLookupGB[which(namesLookupGB$commMatNames == edgesGB[i, "commMat_To"]), 
                                               "bnNames"]
}

edgesGB <- na.omit(edgesGB)
edgesGB <- edgesGB[edgesGB$corrMat_From != "PS" &
                     edgesGB$corrMat_To != "PS", ]

edgesGB$corrVal <- NA
for(i in 1:nrow(edgesGB)){
  edgesGB$corrVal[i] <- atlGBCorr[edgesGB$corrMat_From[i], edgesGB$corrMat_To[i]]
}


# remove links with abs(correlation coefficient) < 0.21, except for physical drivers

for(j in seq(0, 0.6, by = 0.01)){
  remLinks <- edgesGB %>% filter(abs(corrVal) < j,
                                 !corrMat_From %in% c("SST", "SSS", "BT", "BS", "strat"))

  if(nrow(remLinks) == 0){ next}
  atlAssessMat <- commMat

  for(i in 1:nrow(remLinks)){
    #print(commMat[remLinks$commMat_To[i], remLinks$commMat_From[i]])
    atlAssessMat[remLinks$commMat_To[i], remLinks$commMat_From[i]] <- 0
  }
  
  # Check stability
  print(j)
  print(det(-atlAssessMat))
}

remLinks <- edgesGB %>% filter(abs(corrVal) < 0.21,
                               !corrMat_From %in% c("SST", "SSS", "BT", "BS", "strat"))

atlAssessMat <- commMat

for(i in 1:nrow(remLinks)){
  #print(commMat[remLinks$commMat_To[i], remLinks$commMat_From[i]])
  atlAssessMat[remLinks$commMat_To[i], remLinks$commMat_From[i]] <- 0
}

# Perform Loop Analysis ---------------------------------------------------

# gives number of loops in the matrix (number of elements in the list)
# each list item contains the loop pathway described by the community 
# members (row or column number in 'atlAssessMat')
enumerate.loops(atlAssessMat)

# gives number of pathways between community members 'i' and 'j'
#enumerate.paths(atlAssessMat, i = 6, j = 8) # from pelagic habitat to commercial fishery 

# determine feedback and stability
feedback(atlAssessMat)
det(-atlAssessMat)
#(-1^(nrow(atlAssessMat)+1))*det(atlAssessMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(atlAssessMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)


# negative inverse matrix
#solve(-atlAssessMat)

# Get adjoint matrix
adjMat <- make.adjoint(atlAssessMat, status = TRUE)

# write.csv(adjMat, 
#           'C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/adjMat2.0.csv')

# Get absolute feedback matrix
totMat <- make.T(atlAssessMat)

# write.csv(totMat, 
#           'C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/totalFeedback2.0.csv')

# get the weighted feedback matrix
weightedFB <- make.wfm(atlAssessMat, status = TRUE)

# write.csv(weightedFB, 
#           'C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/weightedFeedback2.0.csv')


# Management Scenarios ----------------------------------------------------
adjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/adjMat2.0.csv')
totMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/totalFeedback2.0.csv')
weightMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/weightedFeedback2.0.csv')

# calculate simultaneous presses as sum of individual columns
qnmScenarios <- data.frame(node = adjMat$X,
                           baselineSSTAdj = adjMat$SurfaceTemperature,
                           baselineSSTTot = totMat$SurfaceTemperature,
                           baselineSSTWeight = weightMat$SurfaceTemperature,
                           fishingNode1Adj = adjMat[, "CommercialGroundfishFishery"],
                           fishingNode2Adj = adjMat[, "CommercialShellfishFishery"],
                           fishingNode3Adj = adjMat[, "CommercialPelagicFishery"],
                           fishingNode1Tot = totMat[, "CommercialGroundfishFishery"],
                           fishingNode2Tot = totMat[, "CommercialShellfishFishery"],
                           fishingNode3Tot = totMat[, "CommercialPelagicFishery"],
                           energyNode1Adj = adjMat[, "HabitatNearshore"],
                           energyNode2Adj = adjMat[, "HabitatSeafloorDemersal"],
                           energyNode1Tot = totMat[, "HabitatNearshore"],
                           energyNode2Tot = totMat[, "HabitatSeafloorDemersal"])

qnmScenarios$fishingSumAdj <- rowSums(qnmScenarios[, c("baselineSSTAdj",
                                                       "fishingNode1Adj",
                                                       "fishingNode2Adj",
                                                       "fishingNode3Adj")])
qnmScenarios$fishingSumTot <- rowSums(qnmScenarios[, c("baselineSSTTot",
                                                       "fishingNode1Tot",
                                                       "fishingNode2Tot",
                                                       "fishingNode3Tot")])
qnmScenarios$fishingWeight <- abs(qnmScenarios$fishingSumAdj)/qnmScenarios$fishingSumTot

# negate adjacency values for decreases
qnmScenarios$decFishSumAdj <- rowSums(-qnmScenarios[, c("fishingNode1Adj",
                                                       "fishingNode2Adj",
                                                       "fishingNode3Adj")])
qnmScenarios$decFishSumAdj <- qnmScenarios$decFishSumAdj + qnmScenarios$baselineSSTAdj
qnmScenarios$decFishWeight <- abs(qnmScenarios$decFishSumAdj)/qnmScenarios$fishingSumTot

qnmScenarios$energySumAdj <- rowSums(-qnmScenarios[, c("energyNode1Adj",
                                                       "energyNode2Adj")])
qnmScenarios$energySumTot <- rowSums(qnmScenarios[, c("baselineSSTTot",
                                                      "energyNode1Tot",
                                                       "energyNode2Tot")])
qnmScenarios$energySumAdj <- qnmScenarios$energySumAdj + qnmScenarios$baselineSSTAdj
qnmScenarios$energyWeight <- abs(qnmScenarios$energySumAdj)/qnmScenarios$energySumTot

# write.csv(qnmScenarios, file = 'C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/NEUS_MSE/qnmScenarios_Truth_20210310.csv')
