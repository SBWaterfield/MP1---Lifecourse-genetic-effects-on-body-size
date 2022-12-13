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
ALSPAC <- read.csv('Data/ALSPAC_data.csv')


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

DXAlist <- list(DXATF, DXATrF, DXALM, DXABM)

for (yVal in DXAlist){
  # Create dataframe with all variables neccesary
  DX_TotFat <-  ALSPAC[,c('cidB3618', 'kz021',
                          'f9003c', yVal[1], 'f9ms010',
                          'fe003c', yVal[2], 'fems010',
                          'fg0011a', yVal[3], 'fg3100',
                          'fh0011a', yVal[4], 'fh3000',
                          'FJ003a', yVal[5], 'FJMR020',
                          'FKAR0010', yVal[6], 'FKMS1000',
                          'grs_childbmi_both_oct20_g1',
                          'grs_adultbmi_both_oct20_g1')]
  
  names(DX_TotFat) <- c('ID', 'Sex', 'Age1', 'DX1','Height1', 'Age2', 'DX2','Height2',
                        'Age3', 'DX3', 'Height3', 'Age4', 'DX4', 'Height4', 
                        'Age5', 'DX5','Height5', 'Age6', 'DX6', 'Height6', 'GRS_Child', 'GRS_Adult')
  
  DX_TotFat[DX_TotFat < 0] <- NA # Remove -10 (and other) NA placeholders
  DX_TotFat <- DX_TotFat[!is.na(DX_TotFat$Sex),] # Remove NAs
  # Keep individuals with GRS values
  DX_TotFat <- DX_TotFat[!is.na(DX_TotFat$GRS_Child),]
  # Keep individuals who have at least one DXA measure
  DX_TotFat <- DX_TotFat[(rowSums(is.na(DX_TotFat[, c('DX1','DX2','DX3','DX4','DX5','DX6')])) !=6) ,] #DXA measures
  
  
  
  # Change sex to dummy variable
  # Males 1, Females 0
  DX_TotFat <- transform(DX_TotFat,Sex=ifelse(Sex ==1, 1, 0))
  
  
  # Scale GRS
  DX_TotFat$GRS_Child <- scale(DX_TotFat$GRS_Child)
  DX_TotFat$GRS_Adult <- scale(DX_TotFat$GRS_Adult)
  
  
  # Change to kilos, meters, and years
  DX_TotFat[c('Age1','Age2','Age3','Age4','Age5','Age6')] <- lapply(DX_TotFat[c('Age1','Age2','Age3','Age4','Age5','Age6')], function(x) x/12)
  DX_TotFat[c('DX1','DX2','DX3','DX4','DX5','DX6')] <- lapply(DX_TotFat[c('DX1','DX2','DX3','DX4','DX5','DX6')], function(x) x/1000)
  #DX_TotFat[c('Height1','Height2','Height3','Height4','Height5')] <- lapply(DX_TotFat[c('Height1','Height2','Height3','Height4','Height5')], function(x) x/100)
  #DX_TotFat[c('Height6')] <- lapply(DX_TotFat[c('Height6')], function(x) x/1000)
  # Height in centimeters
  DX_TotFat[c('Height6')] <- lapply(DX_TotFat[c('Height6')], function(x) x/10)
  
  
  # Adjust for height 
  #height raised to the power of each age-and sex-specific β H^β.
  # Used age 18 values for age 25 height adjustment
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
  
  # Info on model
  print(yVal[7])
  print((Mod_TotFat_GRSF <- runMLwiN(DX ~ 1 + V1 + V2 + V3 + V4 +
                                 Height + GRS_Child + GRS_Adult + 
                                 V1*GRS_Child + V1*GRS_Adult +
                                 V2*GRS_Child + V2*GRS_Adult +
                                 V3*GRS_Child + V3*GRS_Adult +
                                 V4*GRS_Child + V4*GRS_Adult + 
                                 (1 + V1+V2+V3+V4|ID) + (1|occasion),
                               estoption = list(resi.store = T, debugmode=F,
                                                reset=c(0, 0), maxiter=150, smat=smat),
                               data = TotFat_Fem)))
  
  print((Mod_TotFat_GRSM <- runMLwiN(DX ~ 1 + V1 + V2 + V3 + V4 +
                                 Height +  GRS_Child + GRS_Adult + 
                                 V1*GRS_Child + V1*GRS_Adult +
                                 V2*GRS_Child + V2*GRS_Adult +
                                 V3*GRS_Child + V3*GRS_Adult +
                                 V4*GRS_Child + V4*GRS_Adult + 
                                 (1+V1+V2+V3+V4|ID) + (1|occasion),
                               estoption = list(resi.store = T, debugmode=F,
                                                reset=c(0, 0), maxiter=150, smat=smat),
                               data = TotFat_Mal)))
  #------------------------------------------------------------------------------#
  CoefF <- data.frame(Mod_TotFat_GRSF@FP)
  CoefF <- data.frame(t(CoefF))
  CoefM <- data.frame(Mod_TotFat_GRSM@FP)
  CoefM <- data.frame(t(CoefM))
  
  Traj_f_GRS0 <-  CoefF$FP_Intercept   +
    CoefF$FP_V1*Mod_TotFat_GRSF@data$V1 + 
    CoefF$FP_V2*Mod_TotFat_GRSF@data$V2 + 
    CoefF$FP_V3*Mod_TotFat_GRSF@data$V3 + 
    CoefF$FP_V4*Mod_TotFat_GRSF@data$V4 
  
  Traj_m_GRS0 <-  CoefM$FP_Intercept   +
    CoefM$FP_V1*Mod_TotFat_GRSM@data$V1 + 
    CoefM$FP_V2*Mod_TotFat_GRSM@data$V2 + 
    CoefM$FP_V3*Mod_TotFat_GRSM@data$V3 + 
    CoefM$FP_V4*Mod_TotFat_GRSM@data$V4 
  
  Traj_f_GRSC <-  CoefF$FP_Intercept + CoefF$FP_GRS_Child  +
    ((CoefF$FP_V1+CoefF$FP_V1.GRS_Child)*Mod_TotFat_GRSF@data$V1) + 
    ((CoefF$FP_V2+CoefF$FP_V2.GRS_Child)*Mod_TotFat_GRSF@data$V2) + 
    ((CoefF$FP_V3+CoefF$FP_V3.GRS_Child)*Mod_TotFat_GRSF@data$V3) + 
    ((CoefF$FP_V4+CoefF$FP_V4.GRS_Child)*Mod_TotFat_GRSF@data$V4) 
  
  Traj_m_GRSC <-  CoefM$FP_Intercept + CoefM$FP_GRS_Child  +
    ((CoefM$FP_V1+CoefM$FP_V1.GRS_Child)*Mod_TotFat_GRSM@data$V1) + 
    ((CoefM$FP_V2+CoefM$FP_V2.GRS_Child)*Mod_TotFat_GRSM@data$V2) + 
    ((CoefM$FP_V3+CoefM$FP_V3.GRS_Child)*Mod_TotFat_GRSM@data$V3) + 
    ((CoefM$FP_V4+CoefM$FP_V4.GRS_Child)*Mod_TotFat_GRSM@data$V4) 
  
  Traj_f_GRSA <-  CoefF$FP_Intercept + CoefF$FP_GRS_Adult  +
    ((CoefF$FP_V1+CoefF$FP_V1.GRS_Adult)*Mod_TotFat_GRSF@data$V1) + 
    ((CoefF$FP_V2+CoefF$FP_V2.GRS_Adult)*Mod_TotFat_GRSF@data$V2) + 
    ((CoefF$FP_V3+CoefF$FP_V3.GRS_Adult)*Mod_TotFat_GRSF@data$V3) + 
    ((CoefF$FP_V4+CoefF$FP_V4.GRS_Adult)*Mod_TotFat_GRSF@data$V4) 
  
  Traj_m_GRSA <-  CoefM$FP_Intercept + CoefM$FP_GRS_Adult  +
    ((CoefM$FP_V1+CoefM$FP_V1.GRS_Adult)*Mod_TotFat_GRSM@data$V1) + 
    ((CoefM$FP_V2+CoefM$FP_V2.GRS_Adult)*Mod_TotFat_GRSM@data$V2) + 
    ((CoefM$FP_V3+CoefM$FP_V3.GRS_Adult)*Mod_TotFat_GRSM@data$V3) + 
    ((CoefM$FP_V4+CoefM$FP_V4.GRS_Adult)*Mod_TotFat_GRSM@data$V4) 
  #------------------------------------------------------------------------------#
  # Sex specific plots with all trajectory values
  Females <- cbind(TotFat_Fem, Traj_f_GRS0, Traj_f_GRSC, Traj_f_GRSA)
  
  Females <- Females %>%
    dplyr::select(Age, occasion, Traj_f_GRS0,Traj_f_GRSC,Traj_f_GRSA) %>%
    gather(key = "GRS", value = "Traj", -c(Age, occasion))
  
  Females$GRS[Females$GRS == "Traj_f_GRS0"] <- "Mean Child & Adult GRS"  
  Females$GRS[Females$GRS == "Traj_f_GRSC"] <- "+1 SD Child GRS"  
  Females$GRS[Females$GRS == "Traj_f_GRSA"] <- "+1 SD Adult GRS" 
  
  Females$Traj <- exp(Females$Traj)
  
  
  Males <- cbind(TotFat_Mal, Traj_m_GRS0, Traj_m_GRSC, Traj_m_GRSA)
  
  Males <- Males %>%
    dplyr::select(Age, occasion, Traj_m_GRS0,Traj_m_GRSC,Traj_m_GRSA) %>%
    gather(key = "GRS", value = "Traj", -c(Age, occasion))
  
  Males$GRS[Males$GRS == "Traj_m_GRS0"] <- "Mean Child & Adult GRS"  
  Males$GRS[Males$GRS == "Traj_m_GRSC"] <- "+1 SD Child GRS"  
  Males$GRS[Males$GRS == "Traj_m_GRSA"] <- "+1 SD Adult GRS" 
  
  Males$Traj <- exp(Males$Traj)
  
  MaxVal <- max(c(Females$Traj, Males$Traj))
  MinVal <- min(c(Females$Traj, Males$Traj))
  
  if (yVal[7] == 'Total Lean Mass'){
    plottitle_M <- sprintf('Male %s ', yVal[7])
    print(ggplot(data=Males, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
      geom_line() + ylab('Lean Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_M) +
        theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      + ylim(MinVal, MaxVal))
  }
  
  if (yVal[7] == 'Total Lean Mass'){
  plottitle_F <- sprintf('Female %s ', yVal[7])
  print(ggplot(data=Females, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
    geom_line() + ylab('Lean Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_F) +
      theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    + ylim(MinVal, MaxVal))
  }
  
  if (yVal[7] == 'Total Bone Mass'){
    plottitle_M <- sprintf('Male %s ', yVal[7])
    print(ggplot(data=Males, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
      geom_line() + ylab('Bone Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_M) +
        theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      + ylim(MinVal, MaxVal))
  }
  
  if (yVal[7] == 'Total Bone Mass'){
    plottitle_F <- sprintf('Female %s ', yVal[7])
    print(ggplot(data=Females, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
      geom_line() + ylab('Bone Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_F) +
        theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
      + ylim(MinVal, MaxVal))
  }
  
  if (yVal[7] == 'Total Fat Mass' | yVal[7] == 'Trunk Fat Mass'){
  plottitle_M <- sprintf('Male %s ', yVal[7])
  print(ggplot(data=Males, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
    geom_line() + ylab('Fat Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_M) +
      theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    + ylim(MinVal, MaxVal))
  }
  
  if (yVal[7] == 'Total Fat Mass' | yVal[7] == 'Trunk Fat Mass'){
  plottitle_F <- sprintf('Female %s ', yVal[7])
  print(ggplot(data=Females, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
    geom_line() + ylab('Fat Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_F) +
      theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    + ylim(MinVal, MaxVal))
  }
  
  
} # Code loop finished  




