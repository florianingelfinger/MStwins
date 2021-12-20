# Variance component analysis to estimate the extent to which each population in the reference framework is modulated by genetics, 
# early and late environment as presented in Ingelfinger, Gerdes et al. "Twin study reveals non-heritable immune perturbations in multiple sclerosis"
# Author: Reinhard Furrer

require( umx)
options(digits=5)

# read data
setwd("/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/2020-06-30_variance_component_decomposition")
Data <- read.csv(file="2020-05-13_twin_data_frequencies_combined.csv", stringsAsFactors=F)[,-1]

# define variables for umx
Data <- Data[41:nrow(Data),]
n <- nrow(Data)
names(Data) <- stringr::str_remove_all(names(Data),"[.]") 
origNames <- names(Data)
origVarNames <- origNames[2:26]
jFam <- as.numeric(as.factor(Data$pair))
iInd <- rep(2, n)
for (i in unique(jFam)) iInd[ match(i, jFam)]<- 1
Data$iInd <- iInd
Data$DS <- as.numeric( as.factor( Data$disease_state))-1
table(Data[,c("DS","disease_state")])

cat("n:",n,", number of twins:",length(unique(jFam)), ", number of MS cases:",sum(Data$DS) )
cat("MZ/DZ individuals:",sum(Data$Zygosity=="MZ"), sum(Data$Zygosity=="DZ"))

wideData <- reshape(Data, drop="X", direction="wide", idvar="pair", timevar="iInd", sep="")
wideData$Zygosity2 <- NULL
mzD <- subset(wideData, Zygosity1 == "MZ")
dzD <- subset(wideData, Zygosity1 == "DZ")
names(mzD) # just checking

nVar <- length(origVarNames) 
res <- array(0, c(nVar, 2, 5)) 
sdOfVar <- numeric(nVar)

# run umx
for (va in 1:nVar){
  cat("\n\n\n",origVarNames[va],"\n")
  m0 <- umxACE(selDVs = origVarNames[va], selCovs = "DS", dzData = dzD, mzData = mzD,
               sep = "", equateMeans=TRUE, intervals=TRUE) 
  res[va,1,] <- m0$output$estimate
  res[va,2,] <- m0$output$standardErrors
  sdOfVar[va] <- sd(c(dzD[,paste0(origVarNames[va],"1")],dzD[,paste0(origVarNames[va],"2")],
                      mzD[,paste0(origVarNames[va],"1")],mzD[,paste0(origVarNames[va],"2")]))
}

# export VCs
VC <- res[,1,3:5]^2
VC <- VC/rowSums(VC)

VC_data <- VC
rownames(VC_data) <- origVarNames
colnames(VC_data) <- c("A", "C", "E")

write.csv(VC_data, file = "VC_pops_wo_MS.csv")