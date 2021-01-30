get.pars <- function(model="prodmodel")
{
  repfile <- paste(model,".rep",sep="")
  parfile <- paste(model,".par",sep="")
  results <- NULL
  # negative log-likelihood
  results$nll <- scan(repfile,n=1)
  # gradient
  pardat <- read.table(parfile,nrows=1,col.names=1:16,header=FALSE,fill=TRUE,
                       comment.char="")
  results$grad <- pardat[1,16]
  #parameter estimates
  results$par <- scan(repfile,n=6,skip=1)
  #biomass estimates
  results$bio <- read.table(repfile,skip=7,header=FALSE)
  return(results)
}

# writes .dat file, may need changes for z=1 part of lab
# modified write.datfile() to accept changes in z values and year-dependent catchabilities
write.datfile <- function(Nsp=1,BioObs=NULL,CatObs=1:10,zVal = 2, #Syear=1,
                          #yrBreak1=length(CatObs),yrBreak2=length(CatObs),
                          # add pieces to scan over theta values and adjust initial values
                          thetaPhase=2, theta=0.5, K=300, r=0.4,
                          outfile="bdnMSY.dat"
)
{
  fyear <- CatObs[1,1]
  lyear <- CatObs[nrow(CatObs), 1]
  
  #outfile <- "bdnMSY.dat"
  write("#Nsp",outfile)
  write(Nsp,outfile,append=TRUE)
  write("# r phase",outfile,append=TRUE)
  write(1,outfile,append=TRUE)
  write("# rinit",outfile,append=TRUE)
  write(r,outfile,append=TRUE)
  write("# k phase",outfile,append=TRUE)
  write(2,outfile,append=TRUE)
  write("# Kinit",outfile,append=TRUE)
  write(K,outfile,append=TRUE)
  write("# z phase",outfile,append=TRUE) # fix z as a constant
  write(-3,outfile,append=TRUE)
  write("# Z init",outfile,append=TRUE) # setting value of z
  write(zVal,outfile,append=TRUE)
  write("# theta phase",outfile,append=TRUE)
  write(thetaPhase,outfile,append=TRUE)
  write("# Theta init",outfile,append=TRUE)
  write(theta,outfile,append=TRUE)
  write("# fyear",outfile,append=TRUE)
  write(fyear,outfile,append=TRUE)
  #write("# Syear",outfile,append=TRUE)
  #write(Syear,outfile,append=TRUE) # first year of survey
  #write("# yrBreak1",outfile,append=TRUE) # first year of change in catchability
  #write(yrBreak1,outfile,append=TRUE) # first year of change in catchability
  #write("# yrBreak2",outfile,append=TRUE) # second year of change in catchability
  #write(yrBreak2,outfile,append=TRUE) # second year of change in catchability
  write("# lyear",outfile,append=TRUE)
  write(lyear,outfile,append=TRUE)
  write("# catches",outfile,append=TRUE)
  write.table(CatObs[,2],outfile,append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
  write("# nbio",outfile,append=TRUE)
  write(nrow(BioObs),outfile,append=TRUE)
  write("# nsurvey",outfile,append=TRUE)
  write(3,outfile,append=TRUE)
  write("# obs bio",outfile,append=TRUE)
  write.table(BioObs,outfile,append=TRUE,col.names=FALSE,
              row.names=FALSE,quote=FALSE)
  #write("# q1 init",outfile,append=TRUE)
  #write(1e-9,outfile,append=TRUE)
  #write("# qAdj phase",outfile,append=TRUE)
  #write(2,outfile,append=TRUE)
  #write("# qAdj1 init",outfile,append=TRUE)
  #write(1,outfile,append=TRUE)
  #write("# qAdj2 init",outfile,append=TRUE)
  #write(1,outfile,append=TRUE)
}


