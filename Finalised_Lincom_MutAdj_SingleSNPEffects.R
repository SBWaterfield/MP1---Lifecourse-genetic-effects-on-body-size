library(R2MLwiN)
library(dplyr)
library(reshape2)
library(lspline)
library(ggplot2)
library(Jmisc)
library(tidyverse)
library(biostat3)
library(dplyr)


# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = F)
}
options(MLwiN_path = mlwin)


# Load in ALSPAC data, and get relevant variables
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP1/ALSPACSingleSNP.Rdata")
ALSPAC <- ALSPACSingleSNP

#total fat
DXATF <- c('f9dx135', 'fedx135', 'fg3254', 'fh2254', 'FJDX135','FKDX1001',
          'Total Fat Mass')

# Trunk Fat
DXATrF <- c('f9dx126', 'fedx126', 'fg3245', 'fh2245', 'FJDX126','FKDX1031',
            'Trunk Fat Mass')

#Lean Mass
DXALM <- c('f9dx136', 'fedx136', 'fg3255', 'fh2255', 'FJDX136','FKDX1002',
          'Total Lean Mass')

#Bone Mass
DXABM <- c('f9dx134', 'fedx134', 'fg3253', 'fh2253', 'FJDX134','FKDX1000',
          'Total Bone Mass')

DXAlist <- list(DXATrF, DXATF, DXALM, DXABM)
DXAlist <- list(DXATF, DXALM)


FemTotFat_SNP <- data.frame()
MalTotFat_SNP <- data.frame()
MalLean_SNP <- data.frame()
FemLean_SNP <- data.frame()

FemTotFatDiff_SNP <- data.frame()
MalTotFatDiff_SNP <- data.frame()
MalLeanDiff_SNP <- data.frame()
FemLeanDiff_SNP <- data.frame()

SNPList <- RichardsonSNP$SNPlist
progcount <- 0

for (SNP in SNPList){
  
  progcount <- progcount + 1
  cat('Progress:', progcount, 'of', length(SNPList), 'begun', '\n')
  
  
  for (yVal in DXAlist){
    # Create dataframe with all variables neccesary
    DX_TotFat <-  ALSPAC[,c('ID', 'kz021',
                            'f9003c', yVal[1], 'f9ms010',
                            'fe003c', yVal[2], 'fems010',
                            'fg0011a', yVal[3], 'fg3100',
                            'fh0011a', yVal[4], 'fh3000',
                            'FJ003a', yVal[5], 'FJMR020',
                            'FKAR0010', yVal[6], 'FKMS1000',
                            SNP)]
    
    names(DX_TotFat) <- c('ID', 'Sex', 'Age1', 'DX1','Height1', 'Age2', 'DX2','Height2',
                          'Age3', 'DX3', 'Height3', 'Age4', 'DX4', 'Height4', 
                          'Age5', 'DX5','Height5', 'Age6', 'DX6', 'Height6', 'SNP')
    
    DX_TotFat[DX_TotFat < 0] <- NA # Remove -10 (and other) NA placeholders
    DX_TotFat <- DX_TotFat[!is.na(DX_TotFat$Sex),] # Remove NAs
    # Keep individuals who have at least one DXA measure
    DX_TotFat <- DX_TotFat[(rowSums(is.na(DX_TotFat[, c('DX1','DX2','DX3','DX4','DX5','DX6')])) !=6) ,] #DXA measures
    
    
    
    # Change sex to dummy variable
    # Males 1, Females 0
    DX_TotFat <- transform(DX_TotFat,Sex=ifelse(Sex ==1, 1, 0))
    
    
    # Change to kilos, meters, and years
    DX_TotFat[c('Age1','Age2','Age3','Age4','Age5','Age6')] <- lapply(DX_TotFat[c('Age1','Age2','Age3','Age4','Age5','Age6')], function(x) x/12)
    DX_TotFat[c('DX1','DX2','DX3','DX4','DX5','DX6')] <- lapply(DX_TotFat[c('DX1','DX2','DX3','DX4','DX5','DX6')], function(x) x/1000)
    #DX_TotFat[c('Height1','Height2','Height3','Height4','Height5')] <- lapply(DX_TotFat[c('Height1','Height2','Height3','Height4','Height5')], function(x) x/100)
    #DX_TotFat[c('Height6')] <- lapply(DX_TotFat[c('Height6')], function(x) x/1000)
    # Height in centimeters
    DX_TotFat[c('Height6')] <- lapply(DX_TotFat[c('Height6')], function(x) x/10)
    
    
    # Adjust for height 
    #height raised to the power of each age-and sex-specific β H^β.
    DX_TotFat<- transform(DX_TotFat, Height1=ifelse(Sex==0, Height1^5.2, Height1^6.6))
    DX_TotFat<- transform(DX_TotFat, Height2=ifelse(Sex==0, Height2^4.2, Height2^5.4))
    DX_TotFat<- transform(DX_TotFat, Height3=ifelse(Sex==0, Height3^3.0, Height3^2.0))
    DX_TotFat<- transform(DX_TotFat, Height4=ifelse(Sex==0, Height4^2.4, Height4^2.4))
    DX_TotFat<- transform(DX_TotFat, Height5=ifelse(Sex==0, Height5^1.8, Height5^1.9))
    DX_TotFat<- transform(DX_TotFat, Height6=ifelse(Sex==0, Height6^2, Height6^2))
    
    
    # melt data to long form
    TotFat <- reshape(DX_TotFat, idvar = "ID", timevar = "occasion", varying = 
                        c('Age1', 'DX1', 'Height1',
                          'Age2', 'DX2', 'Height2',
                          'Age3', 'DX3', 'Height3',
                          'Age4', 'DX4', 'Height4',
                          'Age5', 'DX5', 'Height5',
                          'Age6', 'DX6', 'Height6'),
                      sep = "", direction = "long")
    
    # Center age 
    agecenter <- TotFat[TotFat$occasion==1,] # Age 9
    agecenter <- mean(agecenter$Age, na.rm = T)
    #TotFat$Age_c <- TotFat$Age - agecenter
    TotFat$Age_c <- TotFat$Age - 9
    
    
    
    # Log DXA measures
    TotFat$DX <- log(TotFat$DX)
    

    
    # set knot points for age
    # splines: <13, 13-15, 15-18, 18+
    
    TotFatspline <- as.data.frame(matrix(
      lspline(TotFat$Age_c, c((13-9),(15-9),(18-9)), 
              marginal = F),ncol=4))
    
    TotFat <- cbind(TotFat, TotFatspline)
    
    # Sort data by model hierarchy
    TotFat <-  TotFat[order(TotFat$ID, TotFat$occasion), ]
    
    # Remove non DXA values
    TotFat <- TotFat[!is.na(TotFat$DX),]
    #TotFat <- TotFat[!is.na(TotFat$Height),]
    
    # Sex based data
    TotFat_Fem <- TotFat[TotFat$Sex==0,]
    TotFat_Mal <- TotFat[TotFat$Sex==1,]
    
    # Centre height around mean (Height - mean)
    TotFat_Fem$Height <- scale(TotFat_Fem$Height, scale = F)
    TotFat_Mal$Height <- scale(TotFat_Mal$Height, scale = F)
    
    #------------------------------------------------------------------------------#
    # Sort out covariance matrix
    smat <- matrix(nrow=2, ncol=2)
    smat[1, 1] <- 1 # Level of covariance matrix
    smat[1, 2] <- 2 # Level of covariance matrix
    smat[2,1] <- 0 #Set matrix type to diagonal 
    smat[2,2] <- 1 
    
    #print(yVal[7])
    
    (Mod_TotFat_GRSF <- runMLwiN(DX ~ 1 + V1 + V2 + V3 + V4 +
                                   Height + SNP +
                                   V1*SNP + V2*SNP + V3*SNP + V4*SNP + 
                                   (1 + V1+V2+V3+V4|ID) + (1|occasion),
                                 estoption = list(resi.store = T, debugmode=F,
                                                  reset=c(0, 0), maxiter=150, smat=smat),
                                 data = TotFat_Fem))
    #print(Mod_TotFat_GRSF)
    
    (Mod_TotFat_GRSM <- runMLwiN(DX ~ 1 + V1 + V2 + V3 + V4 +
                                   Height + SNP +
                                   V1*SNP + V2*SNP + V3*SNP + V4*SNP + 
                                   (1+V1+V2+V3+V4|ID) + (1|occasion),
                                 estoption = list(resi.store = T, debugmode=F,
                                                  reset=c(0, 0), maxiter=150, smat=smat),
                                 data = TotFat_Mal))
    
    #print(Mod_TotFat_GRSM)
    #------------------------------------------------------------------------------#
    #Lincom reporting of data
    
    test <- coef(Mod_TotFat_GRSF)
    
    #Mean Child and adult GRS values 
    Fem9  <- lincom(Mod_TotFat_GRSF, c("FP_Intercept"),eform=T)
    Fem13 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + 4*FP_V1"),eform=T)
    Fem15 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + 4*FP_V1 + 2*FP_V2"),eform=T)
    Fem18 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + 4*FP_V1 + 2*FP_V2 + 3*FP_V3"),eform=T)
    Fem25 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + 4*FP_V1 + 2*FP_V2 + 3*FP_V3 + 7*FP_V4"),eform=T)
    
    Mal9  <- lincom(Mod_TotFat_GRSM, c("FP_Intercept"),eform=T)
    Mal13 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + 4*FP_V1"),eform=T)
    Mal15 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + 4*FP_V1 + 2*FP_V2"),eform=T)
    Mal18 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + 4*FP_V1 + 2*FP_V2 + 3*FP_V3"),eform=T)
    Mal25 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + 4*FP_V1 + 2*FP_V2 + 3*FP_V3 + 7*FP_V4"),eform=T)
    

    # +1SD Child GRS calculations
    FemC9 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP"),eform=T)
    MalC9 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP"),eform=T)
    
    FemC13 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP"),eform=T)
    MalC13 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP"),eform=T)
    
    FemC15 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP"),eform=T)
    MalC15 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP"),eform=T)
    
    FemC18 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP"),eform=T)
    MalC18 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP"),eform=T)
    
    FemC25 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP + 
                                       7*FP_V4 + 7*FP_V4:SNP"),eform=T)
    MalC25 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP + 
                                       7*FP_V4 + 7*FP_V4:SNP"),eform=T)
    
    # +1SD Adult GRS calculations
    FemA9 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP"),eform=T)
    MalA9 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP"),eform=T)
    
    FemA13 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP"),eform=T)
    MalA13 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP"),eform=T)
    
    FemA15 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP"),eform=T)
    MalA15 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP"),eform=T)
    
    FemA18 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP+ 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP"),eform=T)
    MalA18 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP+ 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP"),eform=T)
    
    FemA25 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP + 
                                       7*FP_V4 + 7*FP_V4:SNP"),eform=T)
    MalA25 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_SNP + 4*FP_V1 + 4*FP_V1:SNP +
                                       2*FP_V2 + 2*FP_V2:SNP + 
                                       3*FP_V3 + 3*FP_V3:SNP + 
                                       7*FP_V4 + 7*FP_V4:SNP"),eform=T)
    
    # Tables of results
    lincommales <- rbind(Mal9, MalC9, MalA9,
                         Mal13, MalC13, MalA13,
                         Mal15, MalC15, MalA15,
                         Mal18, MalC18, MalA18,
                         Mal25, MalC25, MalA25)
    lincommales$age <- c(9,9,9,13,13,13,15,15,15,18,18,18,25,25,25)
    lincommales$GRS <- c('Mean Child & Adult GRS', '+SNP', '+SNP', 
                         'Mean Child & Adult GRS', '+SNP', '+SNP',
                         'Mean Child & Adult GRS', '+SNP', '+SNP',
                         'Mean Child & Adult GRS', '+SNP', '+SNP',
                         'Mean Child & Adult GRS', '+SNP', '+SNP')
    lincomfemales <- rbind(Fem9, FemC9, FemA9,
                           Fem13, FemC13, FemA13,
                           Fem15, FemC15, FemA15,
                           Fem18, FemC18, FemA18,
                           Fem25, FemC25, FemA25)
    lincomfemales$age <- c(9,9,9,13,13,13,15,15,15,18,18,18,25,25,25)
    lincomfemales$GRS <- c('Mean Child & Adult GRS', '+SNP', '+SNP', 
                           'Mean Child & Adult GRS', '+SNP', '+SNP',
                           'Mean Child & Adult GRS', '+SNP', '+SNP',
                           'Mean Child & Adult GRS', '+SNP', '+SNP',
                           'Mean Child & Adult GRS', '+SNP', '+SNP')
    
    femresults <- data.frame(
      Age = c('9', '13', '15','18', '25'), 
      Mean = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[1],lincomfemales$`2.5 %`[1],lincomfemales$`97.5 %`[1]),
               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[4],lincomfemales$`2.5 %`[4],lincomfemales$`97.5 %`[4]),
               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[7],lincomfemales$`2.5 %`[7],lincomfemales$`97.5 %`[7]),
               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[10],lincomfemales$`2.5 %`[10],lincomfemales$`97.5 %`[10]),
               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[13],lincomfemales$`2.5 %`[13],lincomfemales$`97.5 %`[13])),
      Child = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[2],lincomfemales$`2.5 %`[2],lincomfemales$`97.5 %`[2]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[5],lincomfemales$`2.5 %`[5],lincomfemales$`97.5 %`[5]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[8],lincomfemales$`2.5 %`[8],lincomfemales$`97.5 %`[8]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[11],lincomfemales$`2.5 %`[11],lincomfemales$`97.5 %`[11]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[14],lincomfemales$`2.5 %`[14],lincomfemales$`97.5 %`[14])),
      Adult = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[3],lincomfemales$`2.5 %`[3],lincomfemales$`97.5 %`[3]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[6],lincomfemales$`2.5 %`[6],lincomfemales$`97.5 %`[6]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[9],lincomfemales$`2.5 %`[9],lincomfemales$`97.5 %`[9]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[12],lincomfemales$`2.5 %`[12],lincomfemales$`97.5 %`[12]),
                sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[15],lincomfemales$`2.5 %`[15],lincomfemales$`97.5 %`[15])))
    
    malresults <- data.frame(
      Age = c('9', '13', '15','18', '25'), 
      Mean = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[1],lincommales$`2.5 %`[1],lincommales$`97.5 %`[1]),
               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[4],lincommales$`2.5 %`[4],lincommales$`97.5 %`[4]),
               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[7],lincommales$`2.5 %`[7],lincommales$`97.5 %`[7]),
               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[10],lincommales$`2.5 %`[10],lincommales$`97.5 %`[10]),
               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[13],lincommales$`2.5 %`[13],lincommales$`97.5 %`[13])),
      Child = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[2],lincommales$`2.5 %`[2],lincommales$`97.5 %`[2]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[5],lincommales$`2.5 %`[5],lincommales$`97.5 %`[5]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[8],lincommales$`2.5 %`[8],lincommales$`97.5 %`[8]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[11],lincommales$`2.5 %`[11],lincommales$`97.5 %`[11]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[14],lincommales$`2.5 %`[14],lincommales$`97.5 %`[14])),
      Adult = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[3],lincommales$`2.5 %`[3],lincommales$`97.5 %`[3]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[6],lincommales$`2.5 %`[6],lincommales$`97.5 %`[6]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[9],lincommales$`2.5 %`[9],lincommales$`97.5 %`[9]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[12],lincommales$`2.5 %`[12],lincommales$`97.5 %`[12]),
                sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[15],lincommales$`2.5 %`[15],lincommales$`97.5 %`[15])))
    
    
    if (yVal[7] == 'Total Fat Mass'){
      #write.table(femresults, file = "FemTotFat.tsv", row.names=FALSE, sep="\t")
      #write.table(malresults, file = "MalTotFat.tsv", row.names=FALSE, sep="\t")
      femresults$SNP <- SNP
      malresults$SNP <- SNP
      FemTotFat_SNP <- rbind(FemTotFat_SNP, femresults)
      MalTotFat_SNP <- rbind(MalTotFat_SNP, malresults)
    }
    
    if (yVal[7] == 'Trunk Fat Mass'){
      #write.table(femresults, file = "FemTrunkFat.tsv", row.names=FALSE, sep="\t")
      #write.table(malresults, file = "MalTrunkFat.tsv", row.names=FALSE, sep="\t")
      
    }
    
    if (yVal[7] == 'Total Lean Mass'){
      #write.table(femresults, file = "FemLean.tsv", row.names=FALSE, sep="\t")
      #write.table(malresults, file = "MalLean.tsv", row.names=FALSE, sep="\t")
      femresults$SNP <- SNP
      malresults$SNP <- SNP
      FemLean_SNP <- rbind(FemLean_SNP, femresults)
      MalLean_SNP <- rbind(MalLean_SNP, malresults)
    }
    
    if (yVal[7] == 'Total Bone Mass'){
      #write.table(femresults, file = "FemBone.tsv", row.names=FALSE, sep="\t")
      #write.table(malresults, file = "MalBone.tsv", row.names=FALSE, sep="\t")
    }
    #------------------------------------------------------------------------------#
    # Calculate differences between mean GRS and +1SD GRS
    # +1SD Child GRS calculations
    FemC9 <- lincom(Mod_TotFat_GRSF, c("FP_SNP"),eform=T)
    MalC9 <- lincom(Mod_TotFat_GRSM, c("FP_SNP"),eform=T)
    
    FemC13 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP"),eform=T)
    MalC13 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP"),eform=T)
    
    FemC15 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP"),eform=T)
    MalC15 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP"),eform=T)
    
    FemC18 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP"),eform=T)
    MalC18 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP"),eform=T)
    
    FemC25 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP + 
                                       7*FP_V4:SNP"),eform=T)
    MalC25 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP + 
                                       7*FP_V4:SNP"),eform=T)
    # +1SD Adult GRS calculations
    FemA9 <- lincom(Mod_TotFat_GRSF, c("FP_SNP"),eform=T)
    MalA9 <- lincom(Mod_TotFat_GRSM, c("FP_SNP"),eform=T)
    
    FemA13 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP"),eform=T)
    MalA13 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP"),eform=T)
    
    FemA15 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP"),eform=T)
    MalA15 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP"),eform=T)
    
    FemA18 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP"),eform=T)
    MalA18 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP"),eform=T)
    
    FemA25 <- lincom(Mod_TotFat_GRSF, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP + 
                                       7*FP_V4:SNP"),eform=T)
    MalA25 <- lincom(Mod_TotFat_GRSM, c("FP_SNP +
                                       4*FP_V1:SNP +
                                       2*FP_V2:SNP + 
                                       3*FP_V3:SNP + 
                                       7*FP_V4:SNP"),eform=T)
    
    # Tables of results
    lincommales <- rbind(MalC9, MalA9,
                         MalC13, MalA13,
                         MalC15, MalA15,
                         MalC18, MalA18,
                         MalC25, MalA25)
    lincommales$age <- c(9,9,13,13,15,15,18,18,25,25)
    lincommales$GRS <- c( '+SNP', '+SNP', 
                          '+SNP', '+SNP',
                          '+SNP', '+SNP',
                          '+SNP', '+SNP',
                          '+SNP', '+SNP')
    lincomfemales <- rbind(FemC9, FemA9,
                           FemC13, FemA13,
                           FemC15, FemA15,
                           FemC18, FemA18,
                           FemC25, FemA25)
    lincomfemales$age <- c(9,9,13,13,15,15,18,18,25,25)
    lincomfemales$GRS <- c( '+SNP', '+SNP', 
                            '+SNP', '+SNP',
                            '+SNP', '+SNP',
                            '+SNP', '+SNP',
                            '+SNP', '+SNP')
    
    # Calculate percentage changes
    Pctfun <- function(x) {                                                  
      (x-1)*100
    }
    
    lincomfemales[ , c(1: 3)] <- apply(lincomfemales[ , c(1: 3)], 2, Pctfun)   
    lincommales[ , c(1: 3)] <- apply(lincommales[ , c(1: 3)], 2, Pctfun)    
    
    # # Set levels of GRS values
    # lincomfemales$GRS <- factor(lincomfemales$GRS, 
    #                             levels = c('+SNP', '+SNP'))
    # lincommales$GRS <- factor(lincommales$GRS, 
    #                           levels = c('+SNP', '+SNP'))
  
    # Format tables
    data_reshapemale <- reshape(lincommales,                               
                             idvar = "age",
                             timevar = "GRS",
                             direction = "wide")
    data_reshapefemale <- reshape(lincomfemales,                               
                                idvar = "age",
                                timevar = "GRS",
                                direction = "wide")
    data_reshapemale$`+SNP`<- sprintf("%.2f (%.2f,%.2f)", data_reshapemale$`Estimate.+SNP`, 
                                      data_reshapemale$`2.5 %.+SNP`, data_reshapemale$`97.5 %.+SNP`)
    data_reshapemale$`+SNP`<- sprintf("%.2f (%.2f,%.2f)", data_reshapemale$`Estimate.+SNP`, 
                                                data_reshapemale$`2.5 %.+SNP`, data_reshapemale$`97.5 %.+SNP`)
    data_reshapefemale$`+SNP`<- sprintf("%.2f (%.2f,%.2f)", data_reshapefemale$`Estimate.+SNP`, 
                                                data_reshapefemale$`2.5 %.+SNP`, data_reshapefemale$`97.5 %.+SNP`)
    data_reshapefemale$`+SNP`<- sprintf("%.2f (%.2f,%.2f)", data_reshapefemale$`Estimate.+SNP`, 
                                                data_reshapefemale$`2.5 %.+SNP`, data_reshapefemale$`97.5 %.+SNP`)
    data_reshapemale <- data_reshapemale[,c('age', '+SNP', '+SNP' )]  
    data_reshapefemale <- data_reshapefemale[,c('age', '+SNP', '+SNP' )]  
    
    
    # Set axis values
    MaxVal <- max(c(lincomfemales$`97.5 %`, lincommales$`97.5 %`))
    MinVal <- min(c(lincomfemales$`2.5 %`, lincommales$`2.5 %`))
    
    # Factor age for axis
    lincommales$age <- factor(lincommales$age)
    lincomfemales$age <- factor(lincomfemales$age)
    
    # Plot theme
    My_Theme = theme(plot.title = element_text(size = 26,hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size = 18),
                     axis.text.x = element_text(size = 16),
                     axis.title.y = element_text(size = 18),
                     axis.text.y = element_text(size = 16),
                     legend.title = element_text(16), legend.text = element_text(14))
    
    Overallmax <-16
    Overallmin <- -3
    
    if (yVal[7] == 'Total Fat Mass'){
      test <- MinVal
    }
    
    #Total fat
    if (yVal[7] == 'Total Fat Mass'){
      print(ggplot(lincommales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Total fat mass difference males') +ylab('Percentage difference (%)')+        
          xlab('Age (Years)')+
        geom_hline(yintercept = 0) + ylim(Overallmin, Overallmax) + My_Theme)
      Sys.sleep(0)
      #write.table(data_reshapemale, file = "Results/MalTotFatDiff.tsv", row.names=FALSE, sep="\t")
      #write.table(data_reshapefemale, file = "Results/FemTotFatDiff.tsv", row.names=FALSE, sep="\t")
      data_reshapemale$SNP <- SNP
      data_reshapemale$SNP <- SNP
      
      MalTotFatDiff_SNP <- rbind(MalTotFatDiff_SNP, data_reshapemale)
      FemTotFatDiff_SNP <- rbind(FemTotFatDiff_SNP, data_reshapefemale)
      
      
      
    }
    if (yVal[7] == 'Total Fat Mass'){
      print(ggplot(lincomfemales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Total fat mass difference females') +ylab('Percentage difference (%)')+       
          xlab('Age (Years)')+
        geom_hline(yintercept = 0) +  ylim(Overallmin, Overallmax) + My_Theme)
  
    }
    # Trunk fat
    if (yVal[7] == 'Trunk Fat Mass'){
      print(ggplot(lincommales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Trunk fat mass difference males') +ylab('Percentage difference (%)')+
          xlab('Age (Years)')+
        geom_hline(yintercept = 0) + ylim(Overallmin, Overallmax) +My_Theme)
      Sys.sleep(0)
      #write.table(data_reshapemale, file = "Results/MalTrunkFatDiff.tsv", row.names=FALSE, sep="\t")
      #write.table(data_reshapefemale, file = "Results/FemTrunkFatDiff.tsv", row.names=FALSE, sep="\t")
    }
    if (yVal[7] == 'Trunk Fat Mass'){
      print(ggplot(lincomfemales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Trunk fat mass difference females') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
        geom_hline(yintercept = 0) +  ylim(Overallmin, Overallmax) + My_Theme)
      Sys.sleep(0)
    }
    # Total Lean Mass
    if (yVal[7] == 'Total Lean Mass'){
      print(ggplot(lincommales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Total lean mass difference males') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
        geom_hline(yintercept = 0) + ylim(Overallmin, Overallmax) +My_Theme)
      Sys.sleep(0)
      #write.table(data_reshapemale, file = "Results/MalLeanFatDiff.tsv", row.names=FALSE, sep="\t")
      #write.table(data_reshapefemale, file = "Results/FemLeanFatDiff.tsv", row.names=FALSE, sep="\t")
      data_reshapemale$SNP <- SNP
      data_reshapemale$SNP <- SNP
      MalLeanDiff_SNP <- rbind(MalLeanDiff_SNP, data_reshapemale)
      FemLeanDiff_SNP <- rbind(FemLeanDiff_SNP, data_reshapefemale)
      
      
      
    }
    if (yVal[7] == 'Total Lean Mass'){
      print(ggplot(lincomfemales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Total lean mass difference females') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
        geom_hline(yintercept = 0) +  ylim(Overallmin, Overallmax) +My_Theme)
      Sys.sleep(0)
    }
    # Total Bone Mass
    if (yVal[7] == 'Total Bone Mass'){
      print(ggplot(lincommales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Total bone mass difference males') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
        geom_hline(yintercept = 0) + ylim(MinVal, MaxVal) +My_Theme)
    }
    if (yVal[7] == 'Total Bone Mass'){
      print(ggplot(lincomfemales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
              geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                            position=position_dodge(width=0.3),width=0.2) +
              geom_point(position=position_dodge(width=0.3)) + 
        ggtitle('Total bone mass difference females') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
        geom_hline(yintercept = 0) + ylim(MinVal, MaxVal) +My_Theme)
      Sys.sleep(0)
    }
    #System sleep for plotting
    Sys.sleep(0)
    
    
  # # Rate of change per spline
  #   femmeans0 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept"),eform=T)
  #   malmeans0 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept"),eform=T)
  #   femmeans1 <- lincom(Mod_TotFat_GRSF, c("FP_V1"),eform=T)
  #   malmeans1 <- lincom(Mod_TotFat_GRSM, c("FP_V1"),eform=T)
  #   femmeans2 <- lincom(Mod_TotFat_GRSF, c("FP_V2"),eform=T)
  #   malmeans2 <- lincom(Mod_TotFat_GRSM, c("FP_V2"),eform=T)
  #   femmeans3 <- lincom(Mod_TotFat_GRSF, c("FP_V3"),eform=T)
  #   malmeans3 <- lincom(Mod_TotFat_GRSM, c("FP_V3"),eform=T)
  #   femmeans4 <- lincom(Mod_TotFat_GRSF, c("FP_V4"),eform=T)
  #   malmeans4 <- lincom(Mod_TotFat_GRSM, c("FP_V4"),eform=T)
  #   
  #   # femmeans1C <- lincom(Mod_TotFat_GRSF, c("FP_V1 + FP_V1:GRS_Child"),eform=T)
  #   # malmeans1C <- lincom(Mod_TotFat_GRSM, c("FP_V1 + FP_V1:GRS_Child"),eform=T)
  #   # femmeans2C <- lincom(Mod_TotFat_GRSF, c("FP_V2 + FP_V2:GRS_Child"),eform=T)
  #   # malmeans2C <- lincom(Mod_TotFat_GRSM, c("FP_V2 + FP_V2:GRS_Child"),eform=T)
  #   # femmeans3C <- lincom(Mod_TotFat_GRSF, c("FP_V3 + FP_V3:GRS_Child"),eform=T)
  #   # malmeans3C <- lincom(Mod_TotFat_GRSM, c("FP_V3 + FP_V3:GRS_Child"),eform=T)
  #   # femmeans4C <- lincom(Mod_TotFat_GRSF, c("FP_V4 + FP_V4:GRS_Child"),eform=T)
  #   # malmeans4C <- lincom(Mod_TotFat_GRSM, c("FP_V4 + FP_V4:GRS_Child"),eform=T)
  #   # 
  #   # femmeans1A <- lincom(Mod_TotFat_GRSF, c("FP_V1 + FP_V1:GRS_Adult"),eform=T)
  #   # malmeans1A <- lincom(Mod_TotFat_GRSM, c("FP_V1 + FP_V1:GRS_Adult"),eform=T)
  #   # femmeans2A <- lincom(Mod_TotFat_GRSF, c("FP_V2 + FP_V2:GRS_Adult"),eform=T)
  #   # malmeans2A <- lincom(Mod_TotFat_GRSM, c("FP_V2 + FP_V2:GRS_Adult"),eform=T)
  #   # femmeans3A <- lincom(Mod_TotFat_GRSF, c("FP_V3 + FP_V3:GRS_Adult"),eform=T)
  #   # malmeans3A <- lincom(Mod_TotFat_GRSM, c("FP_V3 + FP_V3:GRS_Adult"),eform=T)
  #   # femmeans4A <- lincom(Mod_TotFat_GRSF, c("FP_V4 + FP_V4:GRS_Adult"),eform=T)
  #   # malmeans4A <- lincom(Mod_TotFat_GRSM, c("FP_V4 + FP_V4:GRS_Adult"),eform=T)
  # 
  #   #femmeans0C <- lincom(Mod_TotFat_GRSF, c("FP_Intercept:GRS_Child"),eform=T)
  #   #malmeans0C <- lincom(Mod_TotFat_GRSM, c("FP_Intercept:GRS_Child"),eform=T)
  #   femmeans1C <- lincom(Mod_TotFat_GRSF, c("FP_V1:GRS_Child"),eform=T)
  #   malmeans1C <- lincom(Mod_TotFat_GRSM, c("FP_V1:GRS_Child"),eform=T)
  #   femmeans2C <- lincom(Mod_TotFat_GRSF, c("FP_V2:GRS_Child"),eform=T)
  #   malmeans2C <- lincom(Mod_TotFat_GRSM, c("FP_V2:GRS_Child"),eform=T)
  #   femmeans3C <- lincom(Mod_TotFat_GRSF, c("FP_V3:GRS_Child"),eform=T)
  #   malmeans3C <- lincom(Mod_TotFat_GRSM, c("FP_V3:GRS_Child"),eform=T)
  #   femmeans4C <- lincom(Mod_TotFat_GRSF, c("FP_V4:GRS_Child"),eform=T)
  #   malmeans4C <- lincom(Mod_TotFat_GRSM, c("FP_V4:GRS_Child"),eform=T)
  #   
  #   #femmeans0A <- lincom(Mod_TotFat_GRSF, c("FP_Intercept:GRS_Adult"),eform=T)
  #   #malmeans0A <- lincom(Mod_TotFat_GRSM, c("FP_Intercept:GRS_Adult"),eform=T)
  #   femmeans1A <- lincom(Mod_TotFat_GRSF, c("FP_V1:GRS_Adult"),eform=T)
  #   malmeans1A <- lincom(Mod_TotFat_GRSM, c("FP_V1:GRS_Adult"),eform=T)
  #   femmeans2A <- lincom(Mod_TotFat_GRSF, c("FP_V2:GRS_Adult"),eform=T)
  #   malmeans2A <- lincom(Mod_TotFat_GRSM, c("FP_V2:GRS_Adult"),eform=T)
  #   femmeans3A <- lincom(Mod_TotFat_GRSF, c("FP_V3:GRS_Adult"),eform=T)
  #   malmeans3A <- lincom(Mod_TotFat_GRSM, c("FP_V3:GRS_Adult"),eform=T)
  #   femmeans4A <- lincom(Mod_TotFat_GRSF, c("FP_V4:GRS_Adult"),eform=T)
  #   malmeans4A <- lincom(Mod_TotFat_GRSM, c("FP_V4:GRS_Adult"),eform=T)
  #   
  #   # Tables of results
  #   lincommales <- rbind(malmeans1, malmeans1C, malmeans1A,
  #                        malmeans2, malmeans2C, malmeans2A,
  #                        malmeans3, malmeans3C, malmeans3A,
  #                        malmeans4, malmeans4C, malmeans4A)
  #   lincommales$age <- c('9-13','9-13','9-13','13-15','13-15','13-15',
  #                        '15-18','15-18','15-18','18-25','18-25','18-25')
  #   lincommales$GRS <- c('Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS', 
  #                        'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
  #                        'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
  #                        'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS')
  #   lincommales$Estimate <- (lincommales$Estimate -1) * 100
  #   lincommales$`2.5 %` <- (lincommales$`2.5 %` -1) * 100
  #   lincommales$`97.5 %` <- (lincommales$`97.5 %` -1) * 100
  #   
  #   
  #   
  #   lincomfemales <- rbind(femmeans1, femmeans1C, femmeans1A,
  #                        femmeans2, femmeans2C, femmeans2A,
  #                        femmeans3, femmeans3C, femmeans3A,
  #                        femmeans4, femmeans4C, femmeans4A)
  #   lincomfemales$age <- c('9-13','9-13','9-13','13-15','13-15','13-15',
  #                        '15-18','15-18','15-18','18-25','18-25','18-25')
  #   lincomfemales$GRS <- c('Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS', 
  #                        'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
  #                        'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
  #                        'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS')
  #   lincomfemales$Estimate <- (lincomfemales$Estimate -1) * 100
  #   lincomfemales$`2.5 %` <- (lincomfemales$`2.5 %` -1) * 100
  #   lincomfemales$`97.5 %` <- (lincomfemales$`97.5 %` -1) * 100
  #   
  #   
  #   femmean <- data.frame(
  #     Age = c('9-13', '13-15', '15-18','18-25'), 
  #     Mean = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[1],lincomfemales$`2.5 %`[1],lincomfemales$`97.5 %`[1]),
  #              sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[4],lincomfemales$`2.5 %`[4],lincomfemales$`97.5 %`[4]),
  #              sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[7],lincomfemales$`2.5 %`[7],lincomfemales$`97.5 %`[7]),
  #              sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[10],lincomfemales$`2.5 %`[10],lincomfemales$`97.5 %`[10])),
  #     Child = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[2],lincomfemales$`2.5 %`[2],lincomfemales$`97.5 %`[2]),
  #               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[5],lincomfemales$`2.5 %`[5],lincomfemales$`97.5 %`[5]),
  #               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[8],lincomfemales$`2.5 %`[8],lincomfemales$`97.5 %`[8]),
  #               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[11],lincomfemales$`2.5 %`[11],lincomfemales$`97.5 %`[11])),
  #     Adult = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[3],lincomfemales$`2.5 %`[3],lincomfemales$`97.5 %`[3]),
  #               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[6],lincomfemales$`2.5 %`[6],lincomfemales$`97.5 %`[6]),
  #               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[9],lincomfemales$`2.5 %`[9],lincomfemales$`97.5 %`[9]),
  #               sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[12],lincomfemales$`2.5 %`[12],lincomfemales$`97.5 %`[12])))
  #   
  #   malmean <- data.frame(
  #     Age = c('9-13', '13-15', '15-18','18-25'), 
  #     Mean = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[1],lincommales$`2.5 %`[1],lincommales$`97.5 %`[1]),
  #              sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[4],lincommales$`2.5 %`[4],lincommales$`97.5 %`[4]),
  #              sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[7],lincommales$`2.5 %`[7],lincommales$`97.5 %`[7]),
  #              sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[10],lincommales$`2.5 %`[10],lincommales$`97.5 %`[10])),
  #     Child = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[2],lincommales$`2.5 %`[2],lincommales$`97.5 %`[2]),
  #               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[5],lincommales$`2.5 %`[5],lincommales$`97.5 %`[5]),
  #               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[8],lincommales$`2.5 %`[8],lincommales$`97.5 %`[8]),
  #               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[11],lincommales$`2.5 %`[11],lincommales$`97.5 %`[11])),
  #     Adult = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[3],lincommales$`2.5 %`[3],lincommales$`97.5 %`[3]),
  #               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[6],lincommales$`2.5 %`[6],lincommales$`97.5 %`[6]),
  #               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[9],lincommales$`2.5 %`[9],lincommales$`97.5 %`[9]),
  #               sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[12],lincommales$`2.5 %`[12],lincommales$`97.5 %`[12])))
  #   
  #   
  #   if (yVal[7] == 'Total Fat Mass'){
  #     write.table(femmean, file = "FemTotFatrate.tsv", row.names=FALSE, sep="\t")
  #     write.table(malmean, file = "MalTotFatrate.tsv", row.names=FALSE, sep="\t")
  #   }
  #   if (yVal[7] == 'Trunk Fat Mass'){
  #     write.table(femmean, file = "FemTrunkFatrate.tsv", row.names=FALSE, sep="\t")
  #     write.table(malmean, file = "MalTrunkFatrate.tsv", row.names=FALSE, sep="\t")
  #   }
  #   if (yVal[7] == 'Total Lean Mass'){
  #     write.table(femmean, file = "FemLeanrate.tsv", row.names=FALSE, sep="\t")
  #     write.table(malmean, file = "MalLeanrate.tsv", row.names=FALSE, sep="\t")
  #   }
  #   if (yVal[7] == 'Total Bone Mass'){
  #     write.table(femmean, file = "FemBonerate.tsv", row.names=FALSE, sep="\t")
  #     write.table(malmean, file = "MalBonerate.tsv", row.names=FALSE, sep="\t")
  #   }
  #   
  #   
    

  } # Model level loop
  

  cat('Progress:', progcount, 'of', length(SNPList), 'finished', '\n')
  
} # SNP level loop

#Write data to file
# write.table(MalTotFat_SNP, file = "Results/MalTotFat_SNP.tsv", row.names=T, sep="\t")
# write.table(FemTotFat_SNP, file = "Results/FemTotFat_SNP.tsv", row.names=T, sep="\t")
# write.table(MalLean_SNP, file = "Results/MalLean_SNP.tsv", row.names=T, sep="\t")
# write.table(FemLean_SNP, file = "Results/FemLean_SNP.tsv", row.names=T, sep="\t")
# 
# 
# write.table(MalTotFatDiff_SNP, file = "Results/MalTotFatDiff_SNP.tsv", row.names=T, sep="\t")
# write.table(FemTotFatDiff_SNP, file = "Results/FemTotFatDiff_SNP.tsv", row.names=T, sep="\t")
# write.table(MalLeanDiff_SNP, file = "Results/MalLeanDiff_SNP.tsv", row.names=T, sep="\t")
# write.table(FemLeanDiff_SNP, file = "Results/FemLeanDiff_SNP.tsv", row.names=T, sep="\t")
# 
 save(MalTotFat_SNP, MalTotFatDiff_SNP, FemTotFat_SNP, FemTotFatDiff_SNP,
      MalLean_SNP, MalLeanDiff_SNP, FemLean_SNP, FemLeanDiff_SNP,
      file = 'SingleSNP_Results.RData')
