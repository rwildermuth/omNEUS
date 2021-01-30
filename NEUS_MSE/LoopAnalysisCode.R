
library(LoopAnalyst)
library(DiagrammeR)

#---------------------------------------------------------
# Author: Robert Wildermuth, rwildermuth@umassd.edu
# Created: 2/24/2016
# Last Modified: 9/11/2017

# Description:
# Imports the interaction matrix for the Georges Bank conceptual model 2.0
# and adapts the data to a (-1,0,1)-specified community matrix for loop 
# analysis. After the community matrix ('commMat') is created, loops and 
# paths are defined and enumerated, feedback and stability are evaluated, 
# and press perturbation is performed using the make.adjoint() function.
# The total feedback and weighted community matrices are also calculated 
# to determine effect determinacy (essentially number of + or - effects 
# devided by the total number of paths possible). All matrices are then 
# saved to file.

#---------------------------------------------------------



# Load the Mental Modeler matrix
mmMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/GBmerged_RemoveProfit_Culture_Feedbacks_2HabPPlinks.csv', 
                  header = TRUE)
# mmMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/GB_Merged_WGNARS_baseline1.csv', 
#                    header = TRUE)
# 32 x 33 data.frame

# Remove 'Climate' because it has no confirmed affects or influences
mmMat <- mmMat[!(mmMat$X %in% 'Climate'), !(colnames(mmMat) %in% c('Climate',"X.1"))]
# 31x 32 data.frame

# 2/26/2016 try removing weather drivers (start with surface and bottom temp/salinity, etc.)
# 3/1/2016  removing all results in det(-A)=-15 and feedback=NA
#           removing Tidal Forcing only has same det(-A) and feedback
# -> no real effect on stability or feedback
#mmMat <- mmMat[!(mmMat$X %in% c('Winds', 
#                                'Source Water Proportions', 
#                                'Tidal Forcing', 
#                                'Air Temperature', 
#                                'Precipitation'
#                                )), 
#               !(colnames(mmMat) %in% c('Winds', 
#                                        'Source.Water.Proportions', 
#                                        'Tidal.Forcing', 
#                                        'Air.Temperature', 
#                                        'Precipitation'
#                                        ))]

# and scup and striped bass
#mmMat <- mmMat[!(mmMat$X %in% 'Scup & Striped Bass'), !(colnames(mmMat) %in% 'Scup...Striped.Bass')]

# Try removing physical forcing
# feedback = NA, det(-A) = -15
#mmMat <- mmMat[!(mmMat$X %in% c('Surface Temperature', 'Surface Salinity', 'Bottom Temperature', 'Bottom Salinity', 'Stratification')), 
#               !(colnames(mmMat) %in% c('Surface.Temperature', 'Surface.Salinity', 'Bottom.Temperature', 'Bottom.Salinity', 'Stratification'))]

# Try just representing the trophic web without fisheries
#mmMat <- mmMat[!(mmMat$X %in% c('Habitat: Pelagic', 'Habitat: Nearshore', 'Habitat: Seafloor & Demersal')), 
#               !(colnames(mmMat) %in% c('Habitat..Pelagic', 'Habitat..Nearshore', 'Habitat..Seafloor...Demersal'))]
#mmMat <- mmMat[!(mmMat$X %in% c('Recreational Fishery', 'Commercial Fishery', 'Cultural Practices & Attachments')), 
#               !(colnames(mmMat) %in% c('Recreational.Fishery', 'Commercial.Fishery', 'Cultural.Practices...Attachments'))]


rNames <- mmMat[, 'X']
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
rownames(commMat) <- colnames(mmMat)#rNames
colnames(commMat) <- colnames(mmMat)
# 31 x 31 matrix

# Try adding self limitation to the trophic web only matrix
#rind <- c('Forage Fish', 'Protected Species', 'Primary Production', 'Groundfish', 
#          'Fished Invertebrates', 'Copepods & Micronekton', 'Benthos', 
#          'Detritus & Bacteria', 'Gelatinous Zooplankton', 'Scup & Striped Bass')
#cind <- c('Forage.Fish', 'Protected.Species','Primary.Production', 'Groundfish', 
#          'Fished.Invertebrates', 'Copepods...Micronekton', 'Benthos', 
#          'Detritus...Bacteria', 'Gelatinous.Zooplankton', 'Scup...Striped.Bass')
#for(i in 1:length(rind)){
#  commMat[rind[i],cind[i]] <- -1
#}

# Add self-limitation to each component
for(i in 1:nrow(commMat)){
  commMat[i,i] <- -1
}
any(is.na(commMat))
# 31 x 31 matrix

# gives number of loops in the matrix (number of elements in the list)
# each list item contains the loop pathway described by the community 
# members (row or column number in 'commMat')
enumerate.loops(commMat)

# gives number of pathways between community members 'i' and 'j'
#enumerate.paths(commMat, i = 6, j = 8) # from pelagic habitat to commercial fishery 

# determine feedback and stability
feedback(commMat)
det(-commMat)
#(-1^(nrow(commMat)+1))*det(commMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(commMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)

# Stability depends on order of column-row pairs for v2.0
indx <- sample(colnames(commMat), length(colnames(commMat)))
indx
#commMat[indx,indx]
det(-commMat)
det(-commMat[indx,indx])
eigen(commMat[indx,indx], symmetric=FALSE, only.values = TRUE)$values

# Save the stable order
commMat <- commMat[indx,indx]
#write.csv(commMat, 
#          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/commMat2.0.csv')
commMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/commMat2.0.csv', 
                      header = TRUE)
# simpleMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/simpCommMat.csv', 
#                       header = TRUE)
rNames <- commMat[, 'X']
commMat <- as.matrix(commMat[, -1])
commMat <- apply(commMat, 2, as.numeric)
rownames(commMat) <- rNames

# negative inverse matrix
#solve(-commMat)

# Get adjoint matrix
adjMat <- make.adjoint(commMat, status = TRUE)

write.csv(adjMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/adjMat2.0.csv')

# Get absolute feedback matrix
totMat <- make.T(commMat)

write.csv(totMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/totalFeedback2.0.csv')

# get the weighted feedback matrix
weightedFB <- make.wfm(commMat, status = TRUE)

write.csv(weightedFB, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/weightedFeedback2.0.csv')

# make the community effects matrix
#commEffects <- make.cem(commMat)
#any(commEffects != 0)

#write.csv(commEffects, 
#          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/commEffects2.0.csv')

#----------------------------------------------------------
# Calculations for simplified model comparison 

# Bring in simplified comunity matrix with single Commercial Fishery and no other Human dimensions nodes
simpleMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpCommMat2.0.csv', 
                      header = TRUE)
# simpleMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/simpCommMat.csv', 
#                       header = TRUE)
rNames <- simpleMat[, 'X']
simpleMat <- as.matrix(simpleMat[, -1])
simpleMat <- apply(simpleMat, 2, as.numeric)
rownames(simpleMat) <- rNames

# determine feedback and stability
feedback(simpleMat)
det(-simpleMat)
#(-1^(nrow(simpleMat)+1))*det(simpleMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(simpleMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)

# negative inverse matrix
#solve(-simpleMat)

# Get adjoint matrix
simpAdjMat <- make.adjoint(simpleMat, status = TRUE)
# testAdj <- function(A) {det(A)*solve(A)}
test5 <- testAdj(-simpleMat)

write.csv(simpAdjMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpAdjMat2.0.csv')

# Get absolute feedback matrix
simpTotMat <- make.T(simpleMat)

write.csv(simpTotMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpTotMat2.0.csv')

# get the weighted feedback matrix
simpWeight <- make.wfm(simpleMat, status = TRUE)

write.csv(simpWeight, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpWeight2.0.csv')

#----------------------------------------------------------
# Calculations for simplified model with lower trophic connections comparison 

# Bring in simplified comunity matrix with single Commercial Fishery and no other Human dimensions nodes
trophSimpMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/simpCommMat_detail2.0.csv', 
                      header = TRUE)
# trophSimpMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/simpCommMat_trophic.csv', 
#                          header = TRUE)
rNames <- trophSimpMat[, 'X']
trophSimpMat <- as.matrix(trophSimpMat[, -1])
trophSimpMat <- apply(trophSimpMat, 2, as.numeric)
rownames(trophSimpMat) <- rNames

# determine feedback and stability
feedback(trophSimpMat)
det(-trophSimpMat)
#(-1^(nrow(trophSimpMat)+1))*det(trophSimpMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(trophSimpMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)

# negative inverse matrix
#solve(-trophSimpMat)

# Get adjoint matrix
trsmpAdjMat <- make.adjoint(trophSimpMat, status = TRUE)
# testAdj <- function(A) {det(A)*solve(A)}
#test5 <- testAdj(-trophSimpMat)

write.csv(trsmpAdjMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/trsmpAdjMat2.0.csv')

# Get absolute feedback matrix
trsmpTotMat <- make.T(trophSimpMat)

write.csv(trsmpTotMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/trsmpTotMat2.0.csv')

# get the weighted feedback matrix
trsmpWeight <- make.wfm(trophSimpMat, status = TRUE)

write.csv(trsmpWeight, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/trsmpWeight2.0.csv')

#----------------------------------------------------------
# Calculations for complex model with lower trophic connections comparison 

trophCommMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/commMat_detail2.0.csv', 
                         header = TRUE)
# trophCommMat <- read.csv('/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/commMat_trophic.csv', 
#                          header = TRUE)
rNames <- trophCommMat[, 'X']
trophCommMat <- as.matrix(trophCommMat[, -1])
trophCommMat <- apply(trophCommMat, 2, as.numeric)
rownames(trophCommMat) <- rNames

# determine feedback and stability
feedback(trophCommMat)
det(-trophCommMat)
#(-1^(nrow(trophCommMat)+1))*det(trophCommMat)
# assess Lyapunov stability (Dambacher et al. 2003 Am Nat)
eigVals <- eigen(trophCommMat, symmetric=FALSE, only.values = TRUE)$values
all(Re(eigVals)<=0)

# negative inverse matrix
#solve(-trophCommMat)

# Get adjoint matrix
trcmplxAdjMat <- make.adjoint(trophCommMat, status = TRUE)
# testAdj <- function(A) {det(A)*solve(A)}
#test5 <- testAdj(-trophCommMat)

write.csv(trcmplxAdjMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/trcmplxAdjMat2.0.csv')

# Get absolute feedback matrix
trcmplxTotMat <- make.T(trophCommMat)

write.csv(trcmplxTotMat, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/trcmplxTotMat2.0.csv')

# get the weighted feedback matrix
trcmplxWeight <- make.wfm(trophCommMat, status = TRUE)

write.csv(trcmplxWeight, 
          'C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_2.0/trcmplxWeight2.0.csv')

