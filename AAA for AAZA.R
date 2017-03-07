#### Parameter Initialisation --------

library(xts)
library(rrcov)
library(tawny)
library(urca)
library(vars)
library(quadprog)
library(Rblpapi)
library(rJava)
library(xlsx)
# library(corpcor)

#The following are run with bloomberg 

# library(Rblpapi)
# library(Rbbg)


# setwd("C:/Users/Marco Cora/Dropbox/docs/R/Asset Allocation")
# 
# #Download from BBG
# 
# #Load funcrions to compute NAV and returns
# source("C:/Users/Marco Cora/Dropbox/docs/R/functions v2.3.R")

#Load securities tickers

securities <- c("700 HK","1970 HK","883 HK","1211 HK","2331 HK","1508 HK","2202 HK","2208 HK","941 HK","1448 HK","2186 HK","2357 HK","2318 HK","1288 HK","857 HK","175 HK","327 HK","BSDE IJ","SMRA IJ","JSMR IJ","012330 KS","005930 KS","051900 KS","URC PM","AC PM","TEL PM","CIT SP","KDCREIT SP","FIRT SP","PREIT SP","ART SP","ADVANC TB","KBANK TB","3008 TT","2498 TT","BIDU US","VIPS US","CTRP US","JD US")

NumOfAssets <- length(securities)

# securities.names <- c("RWO",  "VTI",  "VGK",  "VWO",  "IEI" , "TIP",  "LQD" , "IGOV", "BNDX", "EMB")

suffix <- c(rep(" Equity",length(securities)))


#Initialise fields to download
fields <- c("PX_LAST")
fieldsName <- c("Close")


#Initialise start date and end date of download
startdate <- as.Date("2001-06-01")
enddate <- as.Date(Sys.Date())



#### Bloomberg connection ----------

conn <- blpConnect()
BBGData <- bdh(con=conn, 
               securities=paste(securities, suffix), 
               fields=fields, 
               start.date=startdate, 
               end.date=enddate, 
               options= c("nonTradingDayFillOption"="NON_TRADING_WEEKDAYS",
                          "nonTradingDayFillMethod"="PREVIOUS_VALUE",
                          "currency"="USD"), 
               verbose = F)
blpDisconnect(conn)

VarLst <- list()

for(iBBG in 1:length(securities)){
  BBGData2 <- BBGData[[paste(securities[iBBG], suffix[iBBG])]]
  BBGData2 <- as.xts(BBGData2[,-1], order.by =as.Date(BBGData2$date))
  names(BBGData2) <- paste(securities[iBBG],".",fieldsName, sep="")
  VarLst[[length(VarLst)+1]] <- BBGData2
  # names(VarLst[[iBBG]]) <- securities[iBBG]
}

for(jBBG in 1:length(fields)){
  
  for(iBBG in 1:length(securities)){
    dummy <- as.xts(VarLst[[iBBG]])
    dummy <- dummy[,jBBG]
    names(dummy) <- paste(securities[iBBG],".",fieldsName[jBBG], sep="")
    if(iBBG==1){
      assign(paste("PriceLst.",fieldsName[jBBG],sep=""),dummy)
    } else {
      dummy2 <- get(paste("PriceLst.",fieldsName[jBBG],sep=""))
      dummy2 <- cbind(dummy2,dummy)
      assign(paste("PriceLst.",fieldsName[jBBG],sep=""),dummy2)
    }
  }
}





####setting pricelist as a variable----------

PriceLst <- PriceLst.Close

dim(PriceLst)
# head(PriceLst, 20)

PriceLst <- na.omit(PriceLst)


#### Building of Price arrays and SAA Covariance Estimations  ------------

#in order to avoid asyncronous data problems and autocorrelation issues we use weekly data

PriceLst.ToBeUsed <- PriceLst #PriceLst[endpoints(PriceLst[],on="weeks", k=1)]
periodicity <- 1


#calculate historic log returns

PriceLst.LogRet <- na.omit(log(PriceLst.ToBeUsed/lag(PriceLst.ToBeUsed)))

# PriceLst.cov <- getCov(CovOgk(as.matrix(head(PriceLst.LogRet,nrow(PriceLst.LogRet)))))
# PriceLst.cor <- cov2cor(PriceLst.cov)
# 
# ShrinkageIntensity <- estimate.lambda(PriceLst.LogRet, verbose=TRUE)
# ShrinkageIntensity.var <- estimate.lambda.var(PriceLst.LogRet, verbose=TRUE)
# 
# PriceLst.cor.SH <- ShrinkageIntensity * diag(ncol(PriceLst.cov)) + (1-ShrinkageIntensity) * PriceLst.cor
# PriceLst.Sigma.SH <- (ShrinkageIntensity.var * rep(median(diag(PriceLst.cov)),ncol(PriceLst.cov)) + (1-ShrinkageIntensity.var) * diag(PriceLst.cov))^0.5
# PriceLst.cov.SH <- diag(PriceLst.Sigma.SH) %*%  PriceLst.cor.SH %*% diag(PriceLst.Sigma.SH)


#Calculate Exponential factor for variance estimation

ExponentialDecayFactor <- 0.00

Rescale <- mean(rev((1-ExponentialDecayFactor)^seq(0,nrow(PriceLst.ToBeUsed)-2)))
#by multiplicating for exponentially decaying numbers we are artificially reducing the covariance of the sample
#this does not affect the Factor extraction (apart from a scaling parameter) but will affect the implied returns 
#or the Risk Aversion coefficient
#In order to avoid to have to rescale the risk aversion coefficient we need to first multiply the weights
#by the inverse of the average reduction (parameter Rescale) and then by the time frequency to move the matrix from 
#periodal to annual

ExpoWeights <- rep(rev((1-ExponentialDecayFactor)^seq(0,nrow(PriceLst.ToBeUsed)-2)),ncol(PriceLst.LogRet))/Rescale
ExpoWeights <- matrix(ExpoWeights, ncol=NumOfAssets)


#Estimate Shrinked var, covar and correl
#Log returns are pre-multipled by the resclaed exponentially decaying weights

PriceLst.cov.SH <- cov_shrink(PriceLst.LogRet*ExpoWeights)
PriceLst.cor.SH <- cov2cor(PriceLst.cov.SH)
PriceLst.var.SH <- diag(PriceLst.cov.SH)


#Annualise the covariance matrix

PriceLst.cov.SH <- PriceLst.cov.SH * 252/periodicity
PriceLst.cor.SH <- PriceLst.cor.SH
PriceLst.var.SH <- PriceLst.var.SH * 252/periodicity


PriceLst.sigma.SH <- PriceLst.var.SH^0.5



#### Meucci max diversification -------

#Meucci Maximum number of bets

PriceLst.cor.eigen <- eigen(PriceLst.cor.SH)

# PriceLst.sigma <- diag(PriceLst.cov)^0.5

Meucci.Gamma <- diag(PriceLst.cor.eigen$values^0.5)
Meucci.G <- PriceLst.cor.eigen$vectors
Meucci.C <- Meucci.G %*% Meucci.Gamma %*% t(Meucci.G)

Meucci.d <- diag(NumOfAssets)
Meucci.Pi <- diag(NumOfAssets)
Meucci.Counter <- 0

repeat{
  Meucci.Counter <- Meucci.Counter+1
  Meucci.u <- sqrtm(Meucci.d %*% Meucci.C %*% Meucci.C %*% Meucci.d)
  Meucci.q <- solve(Meucci.u) %*% Meucci.d %*% Meucci.C
  diag(Meucci.d) <- diag(Meucci.q %*% Meucci.C)
  Meucci.Pi1 <- Meucci.d %*% Meucci.q
  if (sum(abs(Meucci.Pi1-Meucci.Pi)) < 0.0000000001)
  {
    break
  }
  if (Meucci.Counter > 10000) 
  {
    break
  }
  Meucci.Pi <- Meucci.Pi1
}

#Final torsion matrix
Meucci.i <- t(diag(PriceLst.sigma.SH) %*% Meucci.Pi1 %*% solve(Meucci.C) %*% solve(diag(PriceLst.sigma.SH)))

#rescaled torsion matrix
Meucci.i.norm <- Meucci.i / t(matrix(rep(colSums(Meucci.i*Meucci.i)^0.5,NumOfAssets),nrow=NumOfAssets))

#tests: with EW portfolio
Portfolio.EW <- rep(1/NumOfAssets,NumOfAssets)
Portfolio.EW.FactCont <- (solve(Meucci.i.norm) %*% Portfolio.EW)^2 * diag(t(Meucci.i.norm) %*% PriceLst.cov.SH %*% Meucci.i.norm) / sum((solve(Meucci.i.norm) %*% Portfolio.EW)^2 * diag(t(Meucci.i.norm) %*% PriceLst.cov.SH %*% Meucci.i.norm))

#test that factor correlation is identity
PriceLst.Factcov <- t(Meucci.i.norm) %*% PriceLst.cov.SH %*% Meucci.i.norm
PriceLst.Factcor <- cov2cor(PriceLst.Factcov)
sum(abs(PriceLst.Factcor-diag(1,NumOfAssets))) < 0.0000000001

#define entropy function
Meucci.Entropy <- function(Weights, Covariance, FactWeights) {
  dummy1 <- (solve(FactWeights) %*% Weights)^2 * diag(t(FactWeights) %*% Covariance %*% FactWeights) / sum((solve(FactWeights) %*% Weights)^2 * diag(t(FactWeights) %*% Covariance %*% FactWeights))
  dummy2 <- dummy1*log(dummy1)
  exp(-sum(dummy2))
}

#define function to be optimise (Optim naturally minises)
Meucci.FunctToOptim <- function(Arguments, Covariance, FactWeights) {
  #   Lambda <- tail(Arguments,1)
  #   param.to.optim <- head(Arguments,-1)
  -Meucci.Entropy(Arguments,Covariance, FactWeights) 
}

Portfolio.ERFC.optimisation <- optim(
  par = rep(1/NumOfAssets,NumOfAssets), 
  fn = Meucci.FunctToOptim,
  Covariance = PriceLst.cov.SH,
  FactWeights = Meucci.i.norm,
  method = "L-BFGS-B",
  lower = rep(-0.1,NumOfAssets),
  upper= rep(1,NumOfAssets))

Portfolio.ERFC.optimisation


Portfolio.ERFC <- Portfolio.ERFC.optimisation$par
sum(Portfolio.ERFC)

#to implement constrained optimisation so that parameters sum to 1, currently rescaling param to 1
Portfolio.ERFC <- Portfolio.ERFC / sum(Portfolio.ERFC)
Portfolio.ERFC

#check that portfolio factor contribution is uniform
Portfolio.ERFC.FactCont <- (solve(Meucci.i.norm) %*% Portfolio.ERFC)^2 * diag(t(Meucci.i.norm) %*% PriceLst.cov.SH %*% Meucci.i.norm) / sum((solve(Meucci.i.norm) %*% Portfolio.ERFC)^2 * diag(t(Meucci.i.norm) %*% PriceLst.cov.SH %*% Meucci.i.norm))
Portfolio.ERFC.FactCont

#### Portfolio Optimisation BL -------

#rescaling risk aversion to move covariance matrix from weekly to yearly
RiskAversion <- 6

#Implied returns for prior portfolio
Portfolio.ERFC.EV <- RiskAversion  * PriceLst.cov.SH %*% Portfolio.ERFC
Portfolio.ERFC.EV


####Black Litterman model-------------

BL.tau <- 0.05


#initialise number of views similar to number of parameters with equal confidence (can be changed as necessary)
# BL.P <- matrix(0,nrow=ncol(PriceLst.cov.SH),ncol=10)
BL.P <- diag(rep(1,NumOfAssets))
BL.conf.dummy <- c(0.5,0.7,0.3,0.7,0.3,0.7,0.1,0.5,0.7,0.5,0.5,0.7,0.5,0.9,0.7,0.5,0.7,0.7,0.7,0.3,0.3,0.5,0.7,0.3,0.3,0.5,0.7,0.9,0.7,0.7,0.7,0.5,0.5,0.7,0.5,0.7,0.7,0.7,0.7)
# BL.conf.dummy[BL.conf.dummy==0.7] <- 0.6
# BL.conf.dummy[BL.conf.dummy==0.9] <- 0.7
# BL.conf.dummy[BL.conf.dummy==0.3] <- 0.4
# BL.conf.dummy[BL.conf.dummy==0.1] <- 0.3
BL.conf <- BL.conf.dummy*1

PriceLst.sigma.SH.rescaled <- diag(PriceLst.sigma.SH*(1/(BL.conf+0.00000001)-1))
BL.Omega <- t(PriceLst.sigma.SH.rescaled) %*% PriceLst.cor.SH %*% PriceLst.sigma.SH.rescaled
# #note: this is different from standard BL model where the error of the view are independent

# BL.Omega <- diag(PriceLst.var.SH*BL.tau*(1/(BL.conf+0.00000001)-1))




# #add two views 
# BL.P <- cbind(BL.P,c(0,1,0,0,0,0,0,0,0,0))
# BL.P <- cbind(BL.P,c(0,0,1,0,0,0,0,0,0,0))
# BL.Omega <- diag(c(PriceLst.sigma.SH*BL.tau*(1/(BL.conf+0.00000001)-1),BL.tau*0.085,BL.tau*0.21))


BL.PriceTarget<- c(22.5625950852221,7.09110131249839,1.22482659034063,6.4464557386349,0.515716459090792,0.309429875454475,2.44965318068126,1.67607849204507,14.8268481988603,0.838039246022537,0.773574688636188,0.838039246022537,5.15716459090792,0.309429875454475,0.838039246022537,0.966968360795235,0.902399092444341,0.167708492148193,0.152462265589267,0.4573867967678,244.758097413723,1359.76720785402,1087.81376628321,4.71799270855672,20.373150332404,51.4690113660733,7.43439149505613,0.892126979406736,0.966470894357297,1.85859787376403,0.892126979406736,5.17613227893602,6.90150970524802,108.166576527853,3.31910858226648,200,14,50,28)

BL.LastPrice <- coredata(tail(PriceLst,1))
BL.TimeFrame <- 1


#Views expectations
BL.Q <- matrix((BL.PriceTarget/BL.LastPrice)^(1/BL.TimeFrame)-1,ncol=1)

#In percentage
PBL<- paste(round(100*BL.Q, 2), "%", sep="")

####write to Excel ------------------
write.xlsx(PBL, "expectations.xlsx")


#calculate posterior
BL.Prior <- Portfolio.ERFC.EV
BL.EV <- BL.Prior + BL.tau * PriceLst.cov.SH %*% BL.P %*% 
  solve(BL.tau * t(BL.P) %*% PriceLst.cov.SH %*% BL.P + BL.Omega) %*% 
  (BL.Q - t(BL.P) %*% BL.Prior)
BL.EV

apply(BL.tau * PriceLst.cov.SH %*% BL.P %*% solve(BL.tau * t(BL.P) %*% PriceLst.cov.SH %*% BL.P + BL.Omega),
      MARGIN = 1, FUN="sum")


#calculate implied portfolio and tilts from prior as well as difference in EV in percentage
Portfolio.BL <- solve(RiskAversion*PriceLst.cov.SH) %*% BL.EV
Portfolio.BL
sum(Portfolio.BL)
Portfolio.BL-Portfolio.ERFC


#constrained optimisation close form
Portfolio.Gamma <- (rep(1,NumOfAssets) %*% solve(PriceLst.cov.SH) %*% BL.EV - RiskAversion ) / (rep(1,NumOfAssets) %*% solve(PriceLst.cov.SH)%*% t(t(rep(1,NumOfAssets))))
Portfolio.BL.Constr <- solve(PriceLst.cov.SH) %*% (BL.EV - t(t(rep(Portfolio.Gamma,NumOfAssets)))) / RiskAversion
Portfolio.BL.Constr
sum(Portfolio.BL.Constr)


#constrained doptimisation using optimiser and adding inequality constraints
Dmat <- PriceLst.cov.SH * RiskAversion
dvec <- matrix(BL.EV,ncol=1,nrow=NumOfAssets)
Amat <- matrix(rep(1,NumOfAssets),ncol=1)
Amat <- cbind(Amat,diag(1,NumOfAssets,NumOfAssets))
Amat <- cbind(Amat,diag(-1,NumOfAssets,NumOfAssets))
bvec <- c(1,rep(-0.06,NumOfAssets),rep(-0.06,NumOfAssets))
meq <- 1

Portfolio.BL.Constr2 <- t(t(solve.QP(Dmat,dvec,Amat,bvec=bvec, meq=meq)$solution))
row.names(Portfolio.BL.Constr2) <- row.names(Portfolio.BL)
Portfolio.BL.Constr2*100



####Write BL allocation ------------

write.xlsx(Portfolio.BL.Constr2, "BL Allocation.xlsx")



# barplot(t(Portfolio.BL.Constr2))

barplot(t(cbind(Portfolio.ERFC,Portfolio.BL.Constr2)), 
        beside=TRUE,
        col=c("darkgreen","darkblue"),
        legend.text=c("EFRC", "Optimal"),
        horiz=FALSE,
        names.arg=securities,
        main="Optimal Portfolios")

barplot(t(cbind(Portfolio.ERFC.EV,BL.EV, BL.Q)), 
        beside=TRUE,
        col=c("darkgreen","darkblue","red"),
        legend.text=c("EFRC", "BL","Aldric"),
        horiz=FALSE,
        names.arg=securities,
        main="Expected Returns")

write.xlsx(BL.EV, "BLER.xlsx")
write.xlsx(BL.Q, "AldricER.xlsx")


#### Cointegration ---------------


library(corrplot)

colnames(PriceLst.cor.SH) <- securities
#colnames(PriceLst.cor.SH.TAA) <- securities

corrplot.mixed (PriceLst.cor.SH, 
                lower ="pie",
                upper = "number")

