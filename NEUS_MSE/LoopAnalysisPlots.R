library(shape)
library(RColorBrewer)
library(circlize)
#---------------------------------------------------------
# Author: Robert Wildermuth (rwildermuth@umassd.edu), Gavin Fay (gfay@umassd.edu)
# Created: 6/6/2016/2016
# Last Modified: 12/28/2016

# Description:
# Code from Gavin Fay is modified to create heat maps of the positive 
# and negative influences and weights for each loop analysis. Plots 
# are also constructed to summarize positive, negative, and neutral 
# effects on management objectives.

# Code is also provided to save the .dot language code to create a network
# figure using graphVis().
#---------------------------------------------------------

# image.scale() from http://menugget.blogspot.de/2013/12/new-version-of-imagescale-function.html
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}

#---------------------------------------------------------
# Generalized code from Gavin to plot weights and trends together
# RPW: still need to figure out appropriate legend code
PlotLoopResults <- function(adjMat, wMat, rNames){
  #image showing sign of adjoint and weights
  colfunc <- colorRampPalette(c("white", "steelblue"))
  layout(matrix(c(1,1,2,3), nrow=2, ncol=2), widths=c(4,lcm(3)), heights=c(1, 4))
  layout.show(3)
  par(mar=c(1.1,20,16,1.1))
  image(1:nrow(adjMat),1:nrow(adjMat),t(wMat[nrow(adjMat):1,]),col = colfunc(10), ylab='', xlab='', xaxt='n', yaxt='n')
  #colorlegend(col = colfunc(10), zlim = c(0,1), zlevels = 11, dz = 0.1,
  #            zval = seq(0,1,by=0.1), log = FALSE, posx = c(0.9, 0.93), 
  #            posy = c(0.05, 0.9), main = NULL, main.cex = 1.0, 
  #            main.col = "black", lab.col = "black", 
  #            digit = 1, left = FALSE)
  # add axes with component lables
  axis(3, at=1:nrow(adjMat), labels=FALSE, tck=0)
  text(x=1:nrow(adjMat), y=par("usr")[2]+0.5,
       labels=rNames, srt=-45, pos = 2, offset = -0.25, xpd=NA, cex = 1.5)
  #axis(side = 1,
  #   at = 1:nrow(adjMat),
  #   labels = rNames,
  #   tck=0, las=2, cex.axis=1.5)
  axis(side = 2,
     at = nrow(adjMat):1,
     labels = rNames,
     tck=0, las=2, cex.axis=1.5)
  #text(1,0,"+",col="black",cex=2)
  # index of entries in 'adjMat' that are positive
  pick <- which(adjMat>0)
  # add the positive effects
  x <- ceiling(pick/nrow(adjMat))
  y <- (nrow(adjMat)+1)-pick%%nrow(adjMat)#31
  y[y==(nrow(adjMat)+1)] <- 1
  text(x,y,"+",col="black",cex=1.5)
  # index of entries in 'adjMat' that are negative
  pick <- which(adjMat<0)
  # add the negative effects
  x <- ceiling(pick/nrow(adjMat))
  y <- (nrow(adjMat)+1)-pick%%nrow(adjMat)#31
  y[y==(nrow(adjMat)+1)] <- 1
  text(x,y,"-",col="black",cex=1.5)
  
  # index of entries in 'adjMat' that are neutral
  #pick <- which(adjMat==0)
  # add the negative effects
  #x <- ceiling(pick/nrow(adjMat))
  #y <- (nrow(adjMat)+1)-pick%%nrow(adjMat)#31
  #y[y==(nrow(adjMat)+1)] <- 1
  #text(x,y,".",col="black",cex=1.5)
  par(mar=c(0.5,0.5,0.5,0.5))
  plot.new()
  #text(x=0.5, y=0.75, labels="+  Positive\n   Response", cex=1.5)
  #text(x=0.5, y=0.5, labels="-  Negative\n    Response", cex = 1.5)
  
  par(mar=c(1.1,0.5,3,6))
  image.scale(z=t(wMat[nrow(adjMat):1,]), col = colfunc(10), 
              breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), axis.pos=4, add.axis=TRUE)
  mtext("Reliability Weight", side=4, line=3, cex=1.25)
 
}

#---------------------------------------------------------

# create plot with the complex model
compAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/adjMat.csv',
                   header=TRUE)
rNames <- compAdjMat[,1]
compAdjMat <- as.matrix(compAdjMat[,2:32])
row.names(compAdjMat) <- rNames
compWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/weightedFeedback.csv',
                 header=TRUE)
compWMat <- as.matrix(compWMat[,2:32])
row.names(compWMat) <- rNames

PlotLoopResults(adjMat = compAdjMat, wMat = compWMat)

# create plot with the simple model
#simpAdjMat <- read.csv('~/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/simpAdjMat.csv',
#                   header=TRUE)
simpAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/IMCC poster/simpAdjMat.csv',
                       header=TRUE)
rNames <- simpAdjMat[,1]
simpAdjMat <- as.matrix(simpAdjMat[,2:27])
row.names(simpAdjMat) <- rNames
#simpWMat <- read.csv('~/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/simpWeight.csv',
#                 header=TRUE)
simpWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/IMCC poster/simpWeight.csv',
                     header=TRUE)
simpWMat <- as.matrix(simpWMat[,2:27])
row.names(simpWMat) <- rNames

#png("~/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/SimpleAdjPlot.png",
png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/SimpleAdjPlot.png",
    width = 12.5, height = 11, units = 'in', res = 800)
PlotLoopResults(adjMat = simpAdjMat, wMat = simpWMat, rNames = rNames)
dev.off()
# for poster: 11.5" wide x 9.75" high

# create plot with the complex trophic model
compAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/trcmplxAdjMat.csv',
                       header=TRUE)
rNames <- compAdjMat[,1]
compAdjMat <- as.matrix(compAdjMat[,2:32])
row.names(compAdjMat) <- rNames
compWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/trcmplxWeight.csv',
                     header=TRUE)
compWMat <- as.matrix(compWMat[,2:32])
row.names(compWMat) <- rNames

PlotLoopResults(adjMat = compAdjMat, wMat = compWMat)

# create plot with the simple trophic model
compAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/trsmpAdjMat.csv',
                       header=TRUE)
rNames <- compAdjMat[,1]
compAdjMat <- as.matrix(compAdjMat[,2:27])
row.names(compAdjMat) <- rNames
compWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GB_GOM_1.5/trsmpWeight.csv',
                     header=TRUE)
compWMat <- as.matrix(compWMat[,2:27])
row.names(compWMat) <- rNames

PlotLoopResults(adjMat = compAdjMat, wMat = compWMat)
#---------------------------------------------------------
# Plot like Reum et al. 2015 Fig 2, but showing whole adjoint matrix
simpAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/IMCC poster/simpAdjMat.csv',
                       header=TRUE)
rNames <- simpAdjMat[,1]
simpAdjMat <- as.matrix(simpAdjMat[,2:27])
row.names(simpAdjMat) <- rNames

simpWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/IMCC poster/simpWeight.csv',
                     header=TRUE)
simpWMat <- as.matrix(simpWMat[,2:27])
row.names(simpWMat) <- rNames

# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- c("#8da0cb", "#fc8d62", "#66c2a5")
xMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)
colMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)
pchMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)

# create matrix indicating character types and colors
for(i in 1:ncol(simpWMat)){
  for(j in 1:nrow(simpWMat)){
    if(simpWMat[j,i]<0.5){
      pchMat[j,i] <- 21
    } else {
      pchMat[j,i] <- 16
    }
    
    if(simpAdjMat[j,i] < 0){
      colMat[j,i] <- colKey[1]
    } else if(simpAdjMat[j,i] > 0){
      colMat[j,i] <- colKey[3]
    } else {
      colMat[j,i] <- colKey[2]
    }

  }
  xMat[, i] <- i
}
png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/SimpleAdjPlot.png",
    width = 14, height = 11, units = 'in', res = 800)
layout(matrix(c(1,1,2,3), nrow=2, ncol=2), widths=c(4,1), heights=c(1, 4))
layout.show(3)
par(mar=c(1,20.5,16,0.1))
plot(as.vector(t(xMat)), rev(rep(c(1:nrow(simpWMat)), each = nrow(simpWMat))), pch = as.vector(t(pchMat)), 
     col = as.vector(t(colMat)), bg = "white", cex = as.vector(t(simpWMat))+1,
     ylab='', xlab='', xaxt='n', yaxt='n')
axis(3, at=1:nrow(simpAdjMat), labels=FALSE, tck=0)
text(x=1:nrow(simpAdjMat), y=par("usr")[2]+0.5,
     labels=rNames, srt=-45, pos = 2, offset = -0.25, xpd=NA, cex = 1.5)
axis(side = 2,
     at = nrow(simpAdjMat):1,
     labels = c("Tidal Forcing", "Winds", "Air Temperature", "Source Water Proportions", "Precipitation", 
                "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity",
                "Stratification", "* Habitat: Pelagic", "* Habitat: Nearshore", "* Habitat: Seafloor & Demersal", 
                "* Protected Species", "* Forage Fish", "* Groundfish", "* Fished Invertebrates", "Copepods & Micronekton",
                "Gelatinous Zooplankton", "Benthos", "Mid Atlantic Groundfish", "Primary Production", 
                "Detritus & Bacteria", "* Recreational Groundfish Fishery", "* Commercial Fishery", 
                "* Cultural Practices & Attachments"),
     tck=0, las=2, cex.axis=1.5)
abline(v = c(11.5, 13.5, 24.5, 25.5))
par(mar=c(0.5,0.5,0.5,0.5))
plot.new()
par(mar=c(20,0.1,20,6))
plot(x=rep(c(1,2,3),3), y=rep(c(1.5,2.5,3.5), each=3), pch=rep(c(16,16,21), 3), 
     col=rep(colKey, each=3), bg="white", cex=rep(c(2,1.5,1), 3), ylab='', xlab='', xaxt='n', yaxt='n', bty="n", 
     xlim = c(0,4), ylim=c(1,5))
axis(4, at=c(1.5, 2.5, 3.5), labels=c("Decrease", "No change", "Increase"), tck=0, las=1, cex.axis=1.25)
axis(1, at=1:3, labels=c(1, 0.5, 0), tck=-0.05, las=1, cex.axis=1.25)
mtext("Reliability\nWeight", side=1, line=4, cex=1.25)
dev.off()
#---------------------------------------------------------
# Revised Fig 2 plot
# Modified 9/14/2017
simpAdjMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/simpAdjMat2.0.csv',
                       header=TRUE)
rNames <- simpAdjMat[,1]
simpAdjMat <- as.matrix(simpAdjMat[,2:27])
row.names(simpAdjMat) <- rNames

simpWMat <- read.csv('C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/simpWeight2.0.csv',
                     header=TRUE)
simpWMat <- as.matrix(simpWMat[,2:27])
row.names(simpWMat) <- rNames

# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- c("#beaed4", "#fdc086", "#7fc97f")
xMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)
colMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)
pchMat <- matrix(NA, nrow = nrow(simpWMat), ncol = nrow(simpWMat), byrow = FALSE)

# create matrix indicating character types and colors
for(i in 1:ncol(simpWMat)){
  for(j in 1:nrow(simpWMat)){
    # if(simpWMat[j,i]<0.5){
    #   pchMat[j,i] <- 21
    # } else {
    #   pchMat[j,i] <- 16
    # }
    
    if(simpAdjMat[j,i] < 0){
      colMat[j,i] <- colKey[1]
      pchMat[j,i] <- 25
    } else if(simpAdjMat[j,i] > 0){
      colMat[j,i] <- colKey[3]
      pchMat[j,i] <- 24
    } else {
      colMat[j,i] <- colKey[2]
      pchMat[j,i] <- 21
    }
    
  }
  xMat[, i] <- i
}
png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/SimpleAdjPlot.png",
    width = 15, height = 11, units = 'in', res = 800)
layout(matrix(c(1,1,2,3), nrow=2, ncol=2), widths=c(5,1), heights=c(1, 4))
layout.show(3)
par(mar=c(1,20.5,16,0.1))
plot(as.vector(t(xMat)), rev(rep(c(1:nrow(simpWMat)), each = nrow(simpWMat))), pch = as.vector(t(pchMat)), 
     col = as.vector(t(colMat)), bg = as.vector(t(colMat)), cex = as.vector(t(simpWMat))+1,
     ylab='', xlab='', xaxt='n', yaxt='n')
axis(3, at=1:nrow(simpAdjMat), labels=FALSE, tck=0)
text(x=1:nrow(simpAdjMat), y=par("usr")[2]+0.5,
     labels=c("Tidal Forcing", "Winds", "Air Temperature", "Source Water Proportions", "Precipitation", 
              "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity",
              "Stratification", "Habitat: Pelagic", "Habitat: Nearshore", "Habitat: Seafloor & Demersal", 
              "Protected Species", "* Forage Fish", "Groundfish", "Fished Invertebrates", "Copepods & Micronekton",
              "Gelatinous Zooplankton", "Benthos", "Mid Atlantic Groundfish", "Primary Production", 
              "Detritus & Bacteria", "Recreational Groundfish Fishery", "Commercial Fishery", 
              "Cultural Practices & Attachments"), 
     srt=-45, pos = 2, offset = -0.25, xpd=NA, cex = 1.5)
axis(side = 2,
     at = nrow(simpAdjMat):1,
     labels = c("Tidal Forcing", "Winds", "Air Temperature", "Source Water Proportions", "Precipitation", 
                "Surface Temperature", "Surface Salinity", "Bottom Temperature", "Bottom Salinity",
                "Stratification", "* Habitat: Pelagic", "* Habitat: Nearshore", "* Habitat: Seafloor & Demersal", 
                "* Protected Species", "* Forage Fish", "* Groundfish", "* Fished Invertebrates", "Copepods & Micronekton",
                "Gelatinous Zooplankton", "Benthos", "Mid Atlantic Groundfish", "Primary Production", 
                "Detritus & Bacteria", "* Recreational Groundfish Fishery", "* Commercial Fishery", 
                "* Cultural Practices & Attachments"),
     tck=0, las=2, cex.axis=1.5)
polygon(x = c(11.5, 11.5, 13.5, 13.5), y=c(0.5,26.5,26.5,0.5), lty = 2)
polygon(x = c(24.5, 24.5, 25.5, 25.5), y=c(0.5,26.5,26.5,0.5), lty = 2)
#brackets(11.5, -2, 13.5, -2, lwd=2, type=1)
par(mar=c(0.5,0.5,0.5,0.5))
plot.new()
par(mar=c(20,0.1,20,6.5))
plot(x=rep(c(1,2,3),3), y=rep(c(1.5,2.5,3.5), each=3), pch=rep(c(25,21,24), each=3), 
     col=rep(colKey, each=3), bg=rep(colKey, each=3), cex=rep(c(2,1.5,1), 3), ylab='', xlab='', xaxt='n', yaxt='n', bty="n", 
     xlim = c(0,4), ylim=c(1,5))
axis(4, at=c(1.5, 2.5, 3.5), labels=c("Decrease", "No change", "Increase"), tck=0, las=1, cex.axis=1.25)
axis(1, at=1:3, labels=c(1, 0.5, 0), tck=-0.05, las=1, cex.axis=1.25)
mtext("Reliability\nWeight", side=1, line=4, cex=1.25)
dev.off()
#---------------------------------------------------------
# plots of cumulative sums of positive, negative, and neutral outcomes for different strategies
# Strategies are increased fishing and reduced nearshore and demersal habitat from energy extraction
loopResults <- data.frame(Model = c(rep("Simple", 2), rep("Complex", 2)), 
                          Strategy = c("Fishing", "Energy Extraction", "Fishing", "Energy Extraction"),
                          Positives = c(4, 6, 9, 9),
                          Negatives = c(4, 3, 4, 5),
                          Neutrals  = c(2, 1, 2, 1))

par(mfrow=c(1,2))
barplot(as.matrix(loopResults[1:2,3:5]), beside = TRUE, 
        ylab = "Outcomes of objectives", main = "Simple model", col = c("grey40", "grey80"))
barplot(as.matrix(loopResults[3:4,3:5]), beside = TRUE,
        main = "Complex model", col = c("grey40", "grey80"))
legend("topright", legend = c("Fishing", "Energy\nExtraction"), 
       fill = c("grey40", "grey80"), box.col = FALSE)

#---------------------------------------------------------
# heatplot of objectives vs strategies and models

## Function for desaturating colors by specified proportion from http://stackoverflow.com/questions/26314701/r-reducing-colour-saturation-of-a-colour-palette
desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- matrix(c("#7bccc4", "#feb24c", "#d95f02", # high certainty cols
                   "#bae4bc", "#ffeda0", "#fc8d62"), # low certainty cols
                 nrow = 2, ncol = 3, byrow = TRUE)

# matrix of values indicating good - 1, neutral - 2, and bad - 3 objective outcomes
# Columns in order: Simple-FocalSpp-Fishing, Complex-FocalSpp-Fishing, Simple-Trophic-Fishing, Complex-Trophic-Fishing, Simple-FocalSpp-Energy, Complex-FocalSpp-Energy, Simple-Trophic-Energy, Complex-Trophic-Energy
objMat <- matrix(data = c(3, 3, 3, 1, 3, 1, 3, 3, # seafloor and demersal habitat
                          2, 2, 2, 2, 1, 1, 1, 1, # nearshore habitat
                          2, 2, 2, 2, 2, 2, 2, 2, # pelagic habitat
                          3, 1, 3, 1, 3, 3, 3, 3, # Fished inverts
                          1, 1, 1, 1, 3, 3, 3, 3, # Forage fish
                          3, 3, 3, 1, 1, 1, 3, 1, # Groundfish
                          3, 1, 3, 1, 1, 1, 2, 1, # TEP species
                          1, 1, 1, 1, 1, 3, 1, 1, # overall wellbeing (Cultural Practices & Attachments)
                          1, 1, 1, 1, 1, 1, 1, 1, # profits
                          1, 1, 1, 1, 1, 1, 1, 1, # Employment
                          1, 1, 1, 1, 1, 1, 1, 1, # food production
                          1, 1, 1, 1, 1, 1, 1, 1), # recreation
                 nrow = 12, ncol = 8, byrow = TRUE)

# matrix of predictability weights, rounded to one decimal place
weightMat <- matrix(data = c(0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # seafloor and demersal habitat
                             1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, # nearshore habitat
                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, # pelagic habitat
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, # Fished inverts
                             0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, # Forage fish
                             0.2, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, # Groundfish
                             0.2, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, # TEP species
                             1.0, 0.1, 0.1, 0.0, 0.1, 0.0, 0.0, 0.0, # overall wellbeing (Cultural Practices & Attachments)
                             0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # profits
                             0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # Employment
                             0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # food production
                             1.0, 0.1, 0.1, 0.0, 0.1, 0.0, 0.0, 0.0), # recreation
                 nrow = 12, ncol = 8, byrow = TRUE)
# shift weight up by 0.05 so weights of 0 are colored in the plot
#weightMat <- weightMat + 0.25
#weightMat <- matrix(sapply(weightMat, function(x){if(x == 1.05) x <- 1 else x <- x}), 
#                    nrow = 11, ncol = 8, byrow = FALSE)

#colMat <- matrix(sapply(objMat, function(x){x <- rev(brewer.pal(n = 3, name = "Spectral"))[x]}), 
#                 nrow = 12, ncol = 8, byrow = FALSE)
colMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)

# reduce saturation of cell based on weighted feedback value
for(i in 1:ncol(weightMat)){
  for(j in 1:nrow(weightMat)){
    if(weightMat[j,i]<0.5){
      colMat[j,i] <- colKey[2, objMat[j,i]]
      # colMat[j,i] <- desat(cols=colMat[j,i], sat=0.35)
    } else {
      colMat[j,i] <- colKey[1, objMat[j,i]]
    }
#    colMat[j,i] <- desat(cols=colMat[j,i], sat=weightMat[j,i])
  }
}

index <- as.factor(colMat)
indexMat <- matrix(as.numeric(index), 
                   nrow = 12, ncol = 8, byrow = FALSE)

#dev.new(width=10.5, height=6.75)
png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ObjectivesPlot2.png",
    width = 12, height = 10, units = 'in', res = 800)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(7,2))
layout.show(2)
par(mar=c(3,17.5,5,0.5)) 

image(1:ncol(objMat),1:nrow(objMat),t(indexMat),
      col = as.character(sort(unique(index))), 
      ylab='', xlab='', xaxt='n', yaxt='n')

# add axes with objective and model/strategy lables
#axis(side = 1,
#     at = 1:ncol(objMat),
#     labels = c('Simple\nFocal', 'Complex\nFocal', # Fishing
#                'Simple\nTrophic', 'Complex\nTrophic', # Fishing
#                'Simple\nFocal', 'Complex\nFocal', # energy
#                'Simple\nTrophic', 'Complex\nTrophic'), # energy
#     tck=0, las=1, cex.axis=1.5)
mtext('Simple', SOUTH<-1, line=1, at=1, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=2, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=2, at=1.5, cex=1.1, col=1)
mtext('Simple', SOUTH<-1, line=1, at=3, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=4, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=2, at=3.5, cex=1.1, col=1)
mtext('Simple', SOUTH<-1, line=1, at=5, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=6, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=2, at=5.5, cex=1.1, col=1)
mtext('Simple', SOUTH<-1, line=1, at=7, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=8, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=2, at=7.5, cex=1.1, col=1)
mtext("Socioeconomics:", SOUTH<-1, line=1, at=-0.56, cex=1.1, col=1, adj = 1)
mtext("Trophic interactions:", SOUTH<-1, line=2, at=-0.56, cex=1.1, col=1, adj = 1)

axis(side = 2,
     at = nrow(objMat):1,
     labels = c("Improve recreation",
                "Optimize seafood *",
               "Increase employment *",
               "Increase profits *",
               "Improve human wellbeing",
               "Increase protected species",
               "Increase groundfish stocks",
               "Increase forage fish stocks",
               "Increase fished invert. stocks",
               "Maintain pelagic habitat",
               "Maintain nearshore habitat",
               "Maintain demersal habitat"),
     tck=0, las=2, cex.axis=1.5)
#abline(v=2.5)
abline(v=4.5)
#abline(v=6.5)
mtext("Fishing", NORTH<-3, line=1, at=2.5, cex=1.5,col=1)
mtext("Energy Production", NORTH<-3, line=1, at=6.5, cex=1.5,col=1)
mtext("Objectives", NORTH<-3, line=3, at=-0.56, cex=1.5, col=1, adj = 1)
mtext("Management Strategy", NORTH<-3, line=3, at = 4.5, cex=1.5, col=1)
for(i in 1:nrow(weightMat)){
  for(j in 1:ncol(weightMat)){
    text(j, i, t(weightMat[i,j]), col = "black", cex=1.5)
  }
}

# index of entries in 'objMat' that are derived from commercial fishery node
#pick <- c(1, 2, 3, 23, 24, 25, 34, 35, 36, 45, 46, 47)
#pick <- 2
# add the negative effects
#x <- ceiling(pick/nrow(objMat))
#x <- c(1,1,1,3,3,3,5,5,5,7,7,7)
#y <- c(9,10,11)
#text(x,y,"*",col="black",cex=1.5)
# for poster: 10.5" wide x 6.75" high
#dev.off()

# make a key for outcomes plot
#scaledWeights <- seq(0,1, by = 0.1)#+0.25

#colKey <- matrix(rep(rev(brewer.pal(n = 3, name = "Spectral")), each = 11), 
#                 nrow = 11, ncol = 3, byrow = FALSE)
#colKey <- matrix(rep(rev(brewer.pal(n = 3, name = "Spectral")), each = 2), 
#                 nrow = 2, ncol = 3, byrow = FALSE)

#for(i in 1:3){
#  for(j in 1:11){
#    colKey[j,i] <- desat(cols=colKey[j,i], sat=scaledWeights[j])
#  }
#}
#for(i in 1:3){
#  colKey[2,i] <- desat(cols=colKey[2,i], sat=0.35)
#}


index <- as.factor(colKey)
indexMat <- matrix(as.numeric(index), 
                   nrow = 2, ncol = 3, byrow = FALSE)
par(mar=c(15,1.1,15,6))
image(1:nrow(colKey),1:ncol(colKey),indexMat,
      col = as.character(sort(unique(index))), 
      ylab='', xlab='', xaxt='n', yaxt='n')
axis(side = 1,
     at = 1:2,
     labels = c("> 0.5", "< 0.5"),
     tck=0, las=1, cex.axis=1.1)
mtext("Weight", SOUTH<-1, line=3, at=1.5, cex=1.5, col=1)
axis(side = 4,
     at = 1:3,
     labels = c("Positive", "Neutral", "Negative"),
     tck=0, las=2, cex.axis=1.5)
dev.off()
#---------------------------------------------------------
# Plot similar to Reum et al. 2015 Fig 2
# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- matrix(c("#7bccc4", "#feb24c", "#d95f02", # high certainty cols
                   "#bae4bc", "#ffeda0", "#fc8d62"), # low certainty cols
                 nrow = 2, ncol = 3, byrow = TRUE)

# matrix of values indicating good - 1, neutral - 2, and bad - 3 objective outcomes
# Columns in order: Simple-FocalSpp-Fishing, Complex-FocalSpp-Fishing, Simple-Trophic-Fishing, Complex-Trophic-Fishing, Simple-FocalSpp-Energy, Complex-FocalSpp-Energy, Simple-Trophic-Energy, Complex-Trophic-Energy
objMat <- matrix(data = c(3, 3, 3, 1, 3, 1, 3, 3, # seafloor and demersal habitat
                          2, 2, 2, 2, 1, 1, 1, 1, # nearshore habitat
                          2, 2, 2, 2, 2, 2, 2, 2, # pelagic habitat
                          3, 1, 3, 1, 3, 3, 3, 3, # Fished inverts
                          1, 1, 1, 1, 3, 3, 3, 3, # Forage fish
                          3, 3, 3, 1, 1, 1, 3, 1, # Groundfish
                          3, 1, 3, 1, 1, 1, 2, 1, # TEP species
                          1, 1, 1, 1, 1, 3, 1, 1, # overall wellbeing (Cultural Practices & Attachments)
                          1, 1, 1, 1, 1, 1, 1, 1, # profits
                          1, 1, 1, 1, 1, 1, 1, 1, # Employment
                          1, 1, 1, 1, 1, 1, 1, 1, # food production
                          1, 1, 1, 1, 1, 1, 1, 1), # recreation
                 nrow = 12, ncol = 8, byrow = TRUE)

# matrix of predictability weights, rounded to one decimal place
weightMat <- matrix(data = c(0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # seafloor and demersal habitat
                             1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, # nearshore habitat
                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, # pelagic habitat
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, # Fished inverts
                             0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, # Forage fish
                             0.2, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, # Groundfish
                             0.2, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, # TEP species
                             1.0, 0.1, 0.1, 0.0, 0.1, 0.0, 0.0, 0.0, # overall wellbeing (Cultural Practices & Attachments)
                             0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # profits
                             0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # Employment
                             0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, # food production
                             1.0, 0.1, 0.1, 0.0, 0.1, 0.0, 0.0, 0.0), # recreation
                    nrow = 12, ncol = 8, byrow = TRUE)

xMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)
colMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)
pchMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)
#bgMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)

# create matrix indicating character types
for(i in 1:ncol(weightMat)){
  for(j in 1:nrow(weightMat)){
    if(weightMat[j,i]<0.5){
      pchMat[j,i] <- 21
    } else {
      pchMat[j,i] <- 19
    }
    
    colMat[j,i] <- colKey[1, objMat[j,i]]
    #bgMat[j,i] <- colKey[2, objMat[j,i]]
  }
  xMat[, i] <- i
}
png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/ObjPlot_v2.png",
    width = 12, height = 10, units = 'in', res = 800)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(7,2))
layout.show(2)
par(mar=c(3,17.5,5,0.5)) 

plot(as.vector(t(xMat)), rep(c(1:12), each = 8), pch = as.vector(t(pchMat)), 
     col = as.vector(t(colMat)), bg = "white", cex = as.vector(t(weightMat))+1,
     ylab='', xlab='', xaxt='n', yaxt='n')
mtext('Simple', SOUTH<-1, line=1, at=1, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=2, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=2, at=1.5, cex=1.1, col=1)
mtext('Simple', SOUTH<-1, line=1, at=3, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=4, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=2, at=3.5, cex=1.1, col=1)
mtext('Simple', SOUTH<-1, line=1, at=5, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=6, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=2, at=5.5, cex=1.1, col=1)
mtext('Simple', SOUTH<-1, line=1, at=7, cex=1.1, col=1)
mtext('Complex', SOUTH<-1, line=1, at=8, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=2, at=7.5, cex=1.1, col=1)
mtext("Socioeconomics:", SOUTH<-1, line=1, at=-0.56, cex=1.1, col=1, adj = 1)
mtext("Trophic interactions:", SOUTH<-1, line=2, at=-0.56, cex=1.1, col=1, adj = 1)

axis(side = 2,
     at = nrow(objMat):1,
     labels = c("Improve recreation",
                "Optimize seafood *",
                "Increase employment *",
                "Increase profits *",
                "Improve human wellbeing",
                "Increase protected species",
                "Increase groundfish stocks",
                "Increase forage fish stocks",
                "Increase fished invert. stocks",
                "Maintain pelagic habitat",
                "Maintain nearshore habitat",
                "Maintain demersal habitat"),
     tck=0, las=2, cex.axis=1.5)
#abline(v=2.5)
abline(v=4.5)
#abline(v=6.5)
mtext("Fishing", NORTH<-3, line=1, at=2.5, cex=1.5,col=1)
mtext("Energy Production", NORTH<-3, line=1, at=6.5, cex=1.5,col=1)
mtext("Objectives", NORTH<-3, line=3, at=-0.56, cex=1.5, col=1, adj = 1)
mtext("Management Strategy", NORTH<-3, line=3, at = 4.5, cex=1.5, col=1)

par(mar=c(20,0.1,20,6))
plot(x=rep(c(1,2,3),3), y=rep(c(1.5,2.5,3.5), each=3), pch=rep(c(16,16,21), 3), 
     col=rep(colKey[1,], each=3), bg="white", cex=rep(c(2,1.5,1), 3), ylab='', xlab='', xaxt='n', yaxt='n', bty="n", 
     xlim = c(0,4), ylim=c(1,5))
axis(4, at=c(1.5, 2.5, 3.5), labels=c("Positive", "Neutral", "Negative"), tck=0, las=1, cex.axis=1.25)
axis(1, at=1:3, labels=c(1, 0.5, 0), tck=-0.05, las=1, cex.axis=1.25)
mtext("Reliability\nWeight", side=1, line=4, cex=1.25)
dev.off()
#---------------------------------------------------------
# Try a radial barchart to show objectives
circos.clear()
circos.par("track.height" = 0.7, gap.after=0)
objScores <- objMat[,1]
objScores[objScores == 3] <- -1
objScores[objScores == 2] <- 0.25

objWt <- weightMat[,1]
objWt[objWt == 0] <- 0.05

# Initialize chart
factors = letters[1:length(objScores)]
circos.initialize(factors = factors, xlim = c(-1, 1))

# Function for deciding color
colFxn <- function(y){
  if(y==1){ 
    return("#7bccc4")
  } else if(y==-1){
    return("#d95f02")
  } else { return("#feb24c")}
}

circos.track(ylim = c(-1,1), y=objScores, x=objWt, panel.fun = function(x, y) {
  #circos.axis(labels.cex=0.5, labels.font=1, lwd=0.8, h="bottom", direction="inside")
  circos.rect(xleft=-x, ybottom=ifelse(y==0.25,-0.25,0), xright=x, ytop=y, col=colFxn(y))
  circos.lines(x=c(-1,1), y=c(0,0), lwd = 2)
}, bg.border=NA)
text(x=c(0.3, 0.75, 1.1, 1.1, 0.75, 0.3),
     y=c(1, 0.8, 0.3, -0.3, -0.8, -1), 
     labels = c("Increase employment *",
                "Optimize seafood *",
                "Improve recreation",
                "Maintain demersal habitat",
                "Maintain nearshore habitat",
                "Maintain pelagic habitat"), adj = 0)
text(x=-c(0.3, 0.75, 1.1, 1.1, 0.75, 0.3),
     y=c(1, 0.8, 0.3, -0.3, -0.8, -1), 
     labels = c("Increase profits *",
                "Improve human wellbeing",
                "Increase protected species",
                "Increase groundfish stocks",
                "Increase forage fish stocks",
                "Increase fished invert. stocks"), adj = 1)
text(0,0,"Simple\nReduced")

#---------------------------------------------------------
# Revised Figure 3
# Modified 9/14/2017

# matrix of values indicating good - 1, neutral - 2, and bad - 3 objective outcomes
# Columns in order: Simple-FocalSpp-Fishing, Complex-FocalSpp-Fishing, Simple-Trophic-Fishing, Complex-Trophic-Fishing, Simple-FocalSpp-Energy, Complex-FocalSpp-Energy, Simple-Trophic-Energy, Complex-Trophic-Energy
objMat <- matrix(data = c(3, 3, 3, 3, 3, 3, 3, 3, # seafloor and demersal habitat
                          2, 2, 2, 2, 3, 2, 3, 2, # nearshore habitat
                          2, 2, 2, 2, 2, 2, 2, 2, # pelagic habitat
                          3, 1, 3, 1, 3, 1, 3, 1, # Fished inverts
                          1, 1, 1, 1, 1, 1, 1, 1, # Forage fish
                          3, 3, 3, 3, 3, 3, 3, 3, # Groundfish
                          3, 3, 3, 3, 3, 3, 3, 3, # TEP species
                          1, 1, 1, 1, 1, 1, 1, 1, # overall wellbeing (Cultural Practices & Attachments)
                          1, 1, 1, 1, 1, 1, 1, 1, # profits
                          1, 1, 1, 1, 1, 1, 1, 1, # Employment
                          1, 1, 1, 1, 1, 1, 1, 1, # food production
                          2, 2, 1, 1, 2, 2, 1, 1), # recreation
                 nrow = 12, ncol = 8, byrow = TRUE)

# matrix of predictability weights, rounded to one decimal place
weightMat <- matrix(data = c(1.0, 0.3, 1.0, 0.3, 0.1, 0.1, 0.1, 0.1, # seafloor and demersal habitat
                             1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, # nearshore habitat
                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, # pelagic habitat
                             0.2, 0.1, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, # Fished inverts
                             0.4, 0.3, 0.4, 0.3, 0.1, 0.1, 0.1, 0.1, # Forage fish
                             0.8, 0.3, 0.8, 0.3, 0.1, 0.1, 0.1, 0.1, # Groundfish
                             0.5, 0.5, 0.5, 0.5, 0.1, 0.2, 0.1, 0.2, # TEP species
                             1.0, 0.4, 0.4, 0.3, 0.1, 0.1, 0.0, 0.1, # overall wellbeing (Cultural Practices & Attachments)
                             1.0, 0.2, 1.0, 0.2, 0.1, 0.1, 0.1, 0.1, # profits
                             1.0, 0.2, 1.0, 0.2, 0.1, 0.1, 0.1, 0.1, # Employment
                             1.0, 0.1, 1.0, 0.1, 0.1, 0.0, 0.1, 0.1, # food production
                             1.0, 1.0, 0.1, 0.2, 1.0, 1.0, 0.0, 0.1), # recreation
                    nrow = 12, ncol = 8, byrow = TRUE)

# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- c("#fc8d62", "#ffeda0", "#bae4bc")
xMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)
colMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)
pchMat <- matrix(NA, nrow = 12, ncol = 8, byrow = FALSE)

# create matrix indicating character types and colors
for(i in 1:ncol(objMat)){
  for(j in 1:nrow(objMat)){
    
    if(objMat[j,i] == 3){
      colMat[j,i] <- colKey[1]
      pchMat[j,i] <- 7
    } else if(objMat[j,i] == 1){
      colMat[j,i] <- colKey[3]
      pchMat[j,i] <- 15
    } else {
      colMat[j,i] <- colKey[2]
      pchMat[j,i] <- 16
    }
    
  }
  xMat[, i] <- i
}

png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/ObjPlot_v3.png",
    width = 12, height = 8, units = 'in', res = 800)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,2))
layout.show(2)
par(mar=c(5,20,5,1)) 

plot(as.vector(t(xMat)), rep(c(1:12), each = 8), pch = as.vector(t(pchMat)), 
     col = as.vector(t(colMat)), bg = "white", cex = (as.vector(t(weightMat))*3)+1.5,
     ylab='', xlab='', xaxt='n', yaxt='n', xlim=c(0.75,8.25))
text(x=seq(1,8)-0.4, y=par("usr")[3]-0.25, adj=0, labels=rep(c('Simple','Complex'), 4), cex=1.1, srt=-45, xpd = TRUE)
# mtext('Complex', SOUTH<-1, line=1, at=2, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=4, at=1.5, cex=1.1, col=1)
# mtext('Simple', SOUTH<-1, line=1, at=3, cex=1.1, col=1)
# mtext('Complex', SOUTH<-1, line=1, at=4, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=4, at=3.5, cex=1.1, col=1)
# mtext('Simple', SOUTH<-1, line=1, at=5, cex=1.1, col=1)
# mtext('Complex', SOUTH<-1, line=1, at=6, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=4, at=5.5, cex=1.1, col=1)
# mtext('Simple', SOUTH<-1, line=1, at=7, cex=1.1, col=1)
# mtext('Complex', SOUTH<-1, line=1, at=8, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=4, at=7.5, cex=1.1, col=1)
mtext("Socioeconomics:", SOUTH<-1, line=1, at=-0.56, cex=1.1, col=1, adj = 1)
mtext("Trophic interactions:", SOUTH<-1, line=4, at=-0.56, cex=1.1, col=1, adj = 1)

axis(side = 2,
     at = nrow(objMat):1,
     labels = c("Improve recreation",
                "Optimize seafood *",
                "Increase employment *",
                "Increase profits *",
                "Improve human wellbeing",
                "Increase protected species",
                "Increase groundfish stocks",
                "Increase forage fish stocks",
                "Increase fished invert. stocks",
                "Maintain pelagic habitat",
                "Maintain nearshore habitat",
                "Maintain demersal habitat"),
     tck=0, las=2, cex.axis=1.5)
#abline(v=2.5)
abline(v=4.5)
#abline(v=6.5)
mtext("Fishing", NORTH<-3, line=1, at=2.5, cex=1.5,col=1)
mtext("Energy Production", NORTH<-3, line=1, at=6.5, cex=1.5,col=1)
mtext("Objectives", NORTH<-3, line=3, at=-0.56, cex=1.5, col=1, adj = 1)
mtext("Management Strategy", NORTH<-3, line=3, at = 4.5, cex=1.5, col=1)

par(mar=c(15,0.1,15,6.5))
plot(x=rep(c(1,2,3),3), y=rep(c(1.5,2.5,3.5), each=3), pch=rep(c(7,16,15), each=3), 
     col=rep(colKey, each=3), bg=rep(colKey, each=3), cex=(rep(c(1,0.5,0), 3)*3)+1.5, ylab='', xlab='', xaxt='n', yaxt='n', bty="n", 
     xlim = c(0.5,3.5), ylim=c(1,5))
axis(4, at=c(1.5, 2.5, 3.5), labels=c("Negative", "Neutral", "Positive"), tck=0, las=1, cex.axis=1.25)
axis(1, at=1:3, labels=c(1, 0.5, 0), tck=-0.05, las=1, cex.axis=1.25)
mtext("Reliability\nWeight", side=1, line=4, cex=1.25)
dev.off()

#---------------------------------------------------------
# Similar figure as Figure 3, but for individual fishery press perturbations
# Modified 9/14/2017

# matrix of values indicating good - 1, neutral - 2, and bad - 3 objective outcomes
objMat <- matrix(data = c(c(3,2,2,1,1,3,3,1,1,1,1,2), # Pelagic fishery, Compex-Reduced
                          c(3,2,2,1,1,3,3,1,1,1,1,2), # Groundfish fishery, Compex-Reduced
                          c(0,0,0,0,0,0,0,0,0,0,0,0), # Shellfish fishery, Compex-Reduced
                          c(3,2,2,1,1,3,3,1,1,1,1,1), # Pelagic fishery, Compex-Detailed
                          c(3,2,2,1,1,3,3,1,1,1,1,1), # Groundfish fishery, Compex-Detailed
                          c(0,0,0,0,0,0,0,0,0,0,0,0)), # Shellfish fishery, Compex-Detailed
                   # c(3, 3, 3, 3, 3, 3, 3, 3, # seafloor and demersal habitat
                   #        2, 2, 2, 2, 3, 2, 3, 2, # nearshore habitat
                   #        2, 2, 2, 2, 2, 2, 2, 2, # pelagic habitat
                   #        3, 1, 3, 1, 3, 1, 3, 1, # Fished inverts
                   #        1, 1, 1, 1, 1, 1, 1, 1, # Forage fish
                   #        3, 3, 3, 3, 3, 3, 3, 3, # Groundfish
                   #        3, 3, 3, 3, 3, 3, 3, 3, # TEP species
                   #        1, 1, 1, 1, 1, 1, 1, 1, # overall wellbeing (Cultural Practices & Attachments)
                   #        1, 1, 1, 1, 1, 1, 1, 1, # profits
                   #        1, 1, 1, 1, 1, 1, 1, 1, # Employment
                   #        1, 1, 1, 1, 1, 1, 1, 1, # food production
                   #        2, 2, 1, 1, 2, 2, 1, 1), # recreation
                 nrow = 12, ncol = 6, byrow = FALSE)

# matrix of predictability weights, rounded to one decimal place
weightMat <- matrix(data = c(0.5, 0.4, 0.0, 0.5, 0.4, 0.0, # seafloor and demersal habitat
                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, # nearshore habitat
                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, # pelagic habitat
                             0.2, 0.2, 0.0, 0.2, 0.2, 0.0, # Fished inverts
                             0.3, 0.5, 0.0, 0.3, 0.5, 0.0, # Forage fish
                             0.5, 0.4, 0.0, 0.5, 0.4, 0.0, # Groundfish
                             0.7, 0.5, 0.0, 0.7, 0.5, 0.0, # TEP species
                             0.6, 0.5, 0.0, 0.4, 0.4, 0.0, # overall wellbeing (Cultural Practices & Attachments)
                             0.2, 0.3, 0.0, 0.2, 0.3, 0.0, # profits
                             0.2, 0.3, 0.0, 0.2, 0.3, 0.0, # Employment
                             0.2, 0.2, 0.0, 0.2, 0.2, 0.0, # food production
                             1.0, 1.0, 1.0, 0.2, 0.3, 0.0), # recreation
                    nrow = 12, ncol = 6, byrow = TRUE)

# implement colors that are colorblind-friendly and print friendly from http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3
colKey <- c("#fc8d62", "#ffeda0", "#bae4bc")
xMat <- matrix(NA, nrow = 12, ncol = 6, byrow = FALSE)
colMat <- matrix(NA, nrow = 12, ncol = 6, byrow = FALSE)
pchMat <- matrix(NA, nrow = 12, ncol = 6, byrow = FALSE)

# create matrix indicating character types and colors
for(i in 1:ncol(objMat)){
  for(j in 1:nrow(objMat)){
    
    if(objMat[j,i] == 3){
      colMat[j,i] <- colKey[1]
      pchMat[j,i] <- 7
    } else if(objMat[j,i] == 1){
      colMat[j,i] <- colKey[3]
      pchMat[j,i] <- 15
    } else {
      colMat[j,i] <- colKey[2]
      pchMat[j,i] <- 16
    }
    
  }
  xMat[, i] <- i
}

png("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/QualLoopAnalysis/ReviewRevisions/RevisionComments/ObjPlot_Fisheries.png",
    width = 12, height = 8, units = 'in', res = 800)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(6,2))
layout.show(2)
par(mar=c(5,20,5,1)) 

plot(as.vector(t(xMat)), rep(c(1:12), each = 6), pch = as.vector(t(pchMat)), 
     col = as.vector(t(colMat)), bg = "white", cex = (as.vector(t(weightMat))*3)+1.5,
     ylab='', xlab='', xaxt='n', yaxt='n', xlim=c(0.75,6.25))
text(x=seq(1,6)-0.4, y=par("usr")[3]-0.25, adj=0, labels=rep(c('Pelagic','Groundfish', 'Shellfish'), 2), cex=1.1, srt=-45, xpd = TRUE)
# mtext('Complex', SOUTH<-1, line=1, at=2, cex=1.1, col=1)
mtext('Reduced', SOUTH<-1, line=4, at=2, cex=1.1, col=1)
# mtext('Simple', SOUTH<-1, line=1, at=3, cex=1.1, col=1)
# mtext('Complex', SOUTH<-1, line=1, at=4, cex=1.1, col=1)
mtext('Detailed', SOUTH<-1, line=4, at=5, cex=1.1, col=1)
# mtext('Simple', SOUTH<-1, line=1, at=5, cex=1.1, col=1)
# mtext('Complex', SOUTH<-1, line=1, at=6, cex=1.1, col=1)
#mtext('Reduced', SOUTH<-1, line=4, at=5.5, cex=1.1, col=1)
# mtext('Simple', SOUTH<-1, line=1, at=7, cex=1.1, col=1)
# mtext('Complex', SOUTH<-1, line=1, at=8, cex=1.1, col=1)
#mtext('Detailed', SOUTH<-1, line=4, at=7.5, cex=1.1, col=1)
mtext("Fishery:", SOUTH<-1, line=1, at=-0.56, cex=1.1, col=1, adj = 1)
mtext("Trophic interactions:", SOUTH<-1, line=4, at=-0.56, cex=1.1, col=1, adj = 1)

axis(side = 2,
     at = nrow(objMat):1,
     labels = c("Improve recreation",
                "Optimize seafood *",
                "Increase employment *",
                "Increase profits *",
                "Improve human wellbeing",
                "Increase protected species",
                "Increase groundfish stocks",
                "Increase forage fish stocks",
                "Increase fished invert. stocks",
                "Maintain pelagic habitat",
                "Maintain nearshore habitat",
                "Maintain demersal habitat"),
     tck=0, las=2, cex.axis=1.5)
#abline(v=2.5)
abline(v=3.5)
#abline(v=6.5)
#mtext("Fishing", NORTH<-3, line=1, at=2.5, cex=1.5,col=1)
#mtext("Energy Production", NORTH<-3, line=1, at=6.5, cex=1.5,col=1)
mtext("Objectives", NORTH<-3, line=3, at=-0.56, cex=1.5, col=1, adj = 1)
#mtext("Management Strategy", NORTH<-3, line=3, at = 4.5, cex=1.5, col=1)

par(mar=c(15,0.1,15,6.5))
plot(x=rep(c(1,2,3),3), y=rep(c(1.5,2.5,3.5), each=3), pch=rep(c(7,16,15), each=3), 
     col=rep(colKey, each=3), bg=rep(colKey, each=3), cex=(rep(c(1,0.5,0), 3)*3)+1.5, ylab='', xlab='', xaxt='n', yaxt='n', bty="n", 
     xlim = c(0.5,3.5), ylim=c(1,5))
axis(4, at=c(1.5, 2.5, 3.5), labels=c("Negative", "Neutral", "Positive"), tck=0, las=1, cex.axis=1.25)
axis(1, at=1:3, labels=c(1, 0.5, 0), tck=-0.05, las=1, cex.axis=1.25)
mtext("Reliability\nWeight", side=1, line=4, cex=1.25)
dev.off()
#---------------------------------------------------------

graph.cm(commMat, 
         file='C:/Users/rwildermuth/Dropbox/PhD_UMass/WGNARS/2016 docs/GeorgesBankQualitativeAnalysis/GBcommMat.dot',
         color = 'bw')
# !!! RW: Get "Error in validate.cm(CM) : Two or more variables in the system are isolated from one another."
#         --> one or more subsystem is not connected to the rest of the system
graph [bgcolor = "transparent", size = "18!,18!", nodesep="1", ranksep="1", rankdir="LR"];
node [fixedsize=true, fontname="Sans", fontsize="76.85", shape=circle, height="2", width="2", style="setlinewidth(4)"];
edge [style="setlinewidth(3)", arrowsize=3];

grViz(diagram = 'digraph G {
      graph [bgcolor = "transparent", size = "20,20", nodesep="1", ranksep="1", rankdir="LR"];
      node [fixedsize=true, fontname="Sans", fontsize="76.85", shape=circle, height="4", width="4", style="setlinewidth(4)"];
      edge [style="setlinewidth(3)", arrowsize=3];
      W -> W [arrowhead=odot];
      W -> STemp;
      W -> SSal;
      W -> BTemp;
      W -> BSal;
      W -> Pel;
      W -> Rec [arrowhead=odot];
      W -> GF [arrowhead=odot];
      W -> PF [arrowhead=odot];
      W -> SF [arrowhead=odot];
      STemp -> STemp [arrowhead=odot];
      STemp -> Pel;
      STemp -> Near;
      STemp -> Str;
      SSal -> SSal [arrowhead=odot];
      SSal -> Pel ;
      SSal -> Near ;
      SSal -> Str ;
      BTemp -> BTemp [arrowhead=odot];
      BTemp -> Near;
      BTemp -> Str;
      BTemp -> Dem;
      BSal -> BSal [arrowhead=odot];
      BSal -> Near ;
      BSal -> Str ;
      BSal -> Dem ;
      Pel -> Pel [arrowhead=odot];
      Pel -> FF;
      Pel -> PS;
      Pel -> PP;
      Pel -> G;
      Rec -> Rec [arrowhead=odot];
      Rec -> MidG [arrowhead=odot];
      Rec -> C ;
      Rec -> SFood ;
      Rec -> Emp ;
      Rec -> Prof ;
      GF -> GF [arrowhead=odot];
      GF -> Dem [arrowhead=odot];
      GF -> PS [arrowhead=odot];
      GF -> G [arrowhead=odot];
      GF -> C;
      GF -> SFood;
      GF -> Emp;
      GF -> Prof;
      PF -> PF [arrowhead=odot];
      PF -> FF [arrowhead=odot];
      PF -> PS [arrowhead=odot];
      PF -> Emp ;
      PF -> Prof ;
      SF -> SF [arrowhead=odot];
      SF -> Dem [arrowhead=odot];
      SF -> Invt [arrowhead=odot];
      SF -> SFood;
      SF -> Emp;
      SF -> Prof;
      SW -> STemp ;
      SW -> SSal ;
      SW -> BTemp ;
      SW -> BSal ;
      SW -> SW [arrowhead=odot];
      Near -> Near [arrowhead=odot];
      Near -> FF;
      Near -> PS;
      Near -> PP;
      Near -> G;
      Near -> Invt;
      Str -> Pel ;
      Str -> Near ;
      Str -> Str [arrowhead=odot];
      Dem -> Dem [arrowhead=odot];
      Dem -> FF;
      Dem -> PP;
      Dem -> G;
      Dem -> Invt;
      FF -> PF ;
      FF -> FF [arrowhead=odot];
      FF -> PS ;
      FF -> G ;
      PS -> GF [arrowhead=odot];
      PS -> FF [arrowhead=odot];
      PS -> PS [arrowhead=odot];
      PP -> PP [arrowhead=odot];
      PP -> Invt ;
      PP -> CM ;
      PP -> Ben ;
      G -> GF;
      G -> SF;
      G -> FF [arrowhead=odot];
      G -> PS;
      G -> G [arrowhead=odot];
      G -> Invt [arrowhead=odot];
      Invt -> SF ;
      Invt -> G ;
      Invt -> Invt [arrowhead=odot];
      Tide -> Str;
      Tide -> Tide [arrowhead=odot];
      AT -> STemp ;
      AT -> SSal ;
      AT -> BTemp ;
      AT -> BSal ;
      AT -> Rec ;
      AT -> AT [arrowhead=odot];
      Per -> Rec [arrowhead=odot];
      Per -> Per [arrowhead=odot];
      MidG -> Rec ;
      MidG -> MidG [arrowhead=odot];
      C -> Rec;
      C -> GF;
      C -> C [arrowhead=odot];
      SFood -> SFood [arrowhead=odot];
      Emp -> Emp [arrowhead=odot];
      Prof -> Rec ;
      Prof -> GF ;
      Prof -> PF ;
      Prof -> SF ;
      Prof -> Prof [arrowhead=odot];
      CM -> FF;
      CM -> PS;
      CM -> CM [arrowhead=odot];
      Ben -> G ;
      Ben -> Invt [arrowhead=odot];
      Ben -> Ben [arrowhead=odot];
      DB -> Invt;
      DB -> DB [arrowhead=odot];
      GZ -> FF [arrowhead=odot];
      GZ -> GZ [arrowhead=odot];
      }
      ')

