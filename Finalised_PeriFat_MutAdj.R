library(R2MLwiN)
library(dplyr)
library(reshape2)
library(lspline)
library(ggplot2)
library(Jmisc)
library(tidyverse)
library(biostat3)
library(dplyr)
library(tidyr)


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

#Peripheral fat
yVal <- c('f9dx117','f9dx108',
          'fedx117', 'fedx108',
          'fg3236', 'fg3227',
          'fh2236', 'fh2227',
          'FJDX117','FJDX108',
          'FKDX1021','FKDX1011',
         'Peripheral Fat Mass')

# Create dataframe with all variables neccesary
DX_TotFat <-  ALSPAC[,c('cidB3618', 'kz021',
                        'f9003c', yVal[1], yVal[2], 'f9ms010',
                        'fe003c', yVal[3], yVal[4], 'fems010',
                        'fg0011a', yVal[5], yVal[6], 'fg3100',
                        'fh0011a', yVal[7], yVal[8], 'fh3000',
                        'FJ003a', yVal[9], yVal[10],'FJMR020',
                        'FKAR0010', yVal[11], yVal[12], 'FKMS1000',
                        'grs_childbmi_male_oct20_g1', 'grs_childbmi_female_oct20_g1',
                        'grs_adultbmi_male_oct20_g1', 'grs_adultbmi_female_oct20_g1')]

DX_TotFat[DX_TotFat < 0] <- NA # Remove -10 (and other) NA placeholders

DX_TotFat$DX1 <- rowSums(DX_TotFat[,c('f9dx117', 'f9dx108')], na.rm=F)
DX_TotFat$DX2 <- rowSums(DX_TotFat[,c('fedx117', 'fedx108')], na.rm=F)
DX_TotFat$DX3 <- rowSums(DX_TotFat[,c('fg3236', 'fg3227')], na.rm=F)
DX_TotFat$DX4 <- rowSums(DX_TotFat[,c('fh2236', 'fh2227')], na.rm=F)
DX_TotFat$DX5 <- rowSums(DX_TotFat[,c('FJDX117', 'FJDX108')], na.rm=F)
DX_TotFat$DX6 <- rowSums(DX_TotFat[,c('FKDX1021', 'FKDX1011')], na.rm=F)

DX_TotFat <- DX_TotFat[,!names(DX_TotFat) %in% c("f9dx117", "f9dx108",'fedx117', 'fedx108',
                                                    'fg3236', 'fg3227','fh2236', 'fh2227',
                                                    'FJDX117', 'FJDX108','FKDX1021', 'FKDX1011')]

names(DX_TotFat) <- c('ID', 'Sex', 'Age1','Height1', 'Age2','Height2',
                      'Age3', 'Height3', 'Age4', 'Height4', 'Age5', 'Height5',
                      'Age6',  'Height6', 'GRS_Child_M','GRS_Child_F', 'GRS_Adult_M', 'GRS_Adult_F',
                      'DX1', 'DX2', 'DX3', 'DX4', 'DX5', 'DX6')


DX_TotFat <- DX_TotFat[!is.na(DX_TotFat$Sex),] # Remove NAs
# Keep individuals with GRS values
DX_TotFat <- DX_TotFat[!is.na(DX_TotFat$GRS_Child_M),]
# Keep individuals who have at least one DXA measure
DX_TotFat <- DX_TotFat[(rowSums(is.na(DX_TotFat[, c('DX1','DX2','DX3','DX4','DX5','DX6')])) !=6) ,] #DXA measures



# Change sex to dummy variable
# Males 1, Females 0
DX_TotFat <- transform(DX_TotFat,Sex=ifelse(Sex ==1, 1, 0))

# Remove Female GRS values from Males, and vice versa
DX_TotFat<- transform(DX_TotFat, GRS_Child_M=ifelse(Sex==1, GRS_Child_M, NA))
DX_TotFat<- transform(DX_TotFat, GRS_Child_F=ifelse(Sex==0, GRS_Child_F, NA))
DX_TotFat<- transform(DX_TotFat, GRS_Adult_M=ifelse(Sex==1, GRS_Adult_M, NA))
DX_TotFat<- transform(DX_TotFat, GRS_Adult_F=ifelse(Sex==0, GRS_Adult_F, NA))


# Scale GRS
DX_TotFat$GRS_Child_M <- scale(DX_TotFat$GRS_Child_M)
DX_TotFat$GRS_Child_F <- scale(DX_TotFat$GRS_Child_F)
DX_TotFat$GRS_Adult_M <- scale(DX_TotFat$GRS_Adult_M)
DX_TotFat$GRS_Adult_F <- scale(DX_TotFat$GRS_Adult_F)

# Alter GRS into single columns
DX_TotFat<- transform(DX_TotFat, GRS_Child=ifelse(Sex==0, GRS_Child_F, GRS_Child_M ))
DX_TotFat<- transform(DX_TotFat, GRS_Adult=ifelse(Sex==0, GRS_Adult_F, GRS_Adult_M ))


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


(Mod_TotFat_GRSF <- runMLwiN(DX ~ 1 + V1 + V2 + V3 + V4 +
                               Height + GRS_Child + GRS_Adult + 
                               V1*GRS_Child + V1*GRS_Adult +
                               V2*GRS_Child + V2*GRS_Adult +
                               V3*GRS_Child + V3*GRS_Adult +
                               V4*GRS_Child + V4*GRS_Adult + 
                               (1 + V1+V2+V3+V4|ID) + (1|occasion),
                             estoption = list(resi.store = T, debugmode=F,
                                              reset=c(0, 0), maxiter=150, smat=smat),
                             data = TotFat_Fem))

(Mod_TotFat_GRSM <- runMLwiN(DX ~ 1 + V1 + V2 + V3 + V4 +
                               Height +  GRS_Child + GRS_Adult + 
                               V1*GRS_Child + V1*GRS_Adult +
                               V2*GRS_Child + V2*GRS_Adult +
                               V3*GRS_Child + V3*GRS_Adult +
                               V4*GRS_Child + V4*GRS_Adult + 
                               (1+V1+V2+V3+V4|ID) + (1|occasion),
                             estoption = list(resi.store = T, debugmode=F,
                                              reset=c(0, 0), maxiter=150, smat=smat),
                             data = TotFat_Mal))

#------------------------------------------------------------------------------#
# Trajectories
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


plottitle_M <- sprintf('Male %s ', yVal[13])
ggplot(data=Males, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
          geom_line() + ylab('Fat Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_M) +
  theme(plot.title = element_text(size = 15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
       + ylim(MinVal, MaxVal)


plottitle_F <- sprintf('Female %s ', yVal[13])
ggplot(data=Females, aes(x=Age, y=Traj, group=GRS, color = GRS)) +
          geom_line() + ylab('Fat Mass (Kg)') + xlab('Age (Years)') + ggtitle(plottitle_F) +
          theme(plot.title = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
          + ylim(MinVal, MaxVal)


#------------------------------------------------------------------------------#
#Lincom reporting of data

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
FemC9 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Child"),eform=T)
MalC9 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Child"),eform=T)

FemC13 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child"),eform=T)
MalC13 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child"),eform=T)

FemC15 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child +
                                   2*FP_V2 + 2*FP_V2:GRS_Child"),eform=T)
MalC15 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child +
                                   2*FP_V2 + 2*FP_V2:GRS_Child"),eform=T)

FemC18 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child +
                                   2*FP_V2 + 2*FP_V2:GRS_Child + 
                                   3*FP_V3 + 3*FP_V3:GRS_Child"),eform=T)
MalC18 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child +
                                   2*FP_V2 + 2*FP_V2:GRS_Child + 
                                   3*FP_V3 + 3*FP_V3:GRS_Child"),eform=T)

FemC25 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child +
                                   2*FP_V2 + 2*FP_V2:GRS_Child + 
                                   3*FP_V3 + 3*FP_V3:GRS_Child + 
                                   7*FP_V4 + 7*FP_V4:GRS_Child"),eform=T)
MalC25 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Child + 4*FP_V1 + 4*FP_V1:GRS_Child +
                                   2*FP_V2 + 2*FP_V2:GRS_Child + 
                                   3*FP_V3 + 3*FP_V3:GRS_Child + 
                                   7*FP_V4 + 7*FP_V4:GRS_Child"),eform=T)

# +1SD Adult GRS calculations
FemA9 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Adult"),eform=T)
MalA9 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Adult"),eform=T)

FemA13 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Adult + 4*FP_V1 + 4*FP_V1:GRS_Adult"),eform=T)
MalA13 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Adult + 4*FP_V1 + 4*FP_V1:GRS_Adult"),eform=T)

FemA15 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Adult + 4*FP_V1 + 4*FP_V1:GRS_Adult +
                                   2*FP_V2 + 2*FP_V2:GRS_Adult"),eform=T)
MalA15 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Adult + 4*FP_V1 + 4*FP_V1:GRS_Adult +
                                   2*FP_V2 + 2*FP_V2:GRS_Adult"),eform=T)

FemA18 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Adult+ 4*FP_V1 + 4*FP_V1:GRS_Adult +
                                   2*FP_V2 + 2*FP_V2:GRS_Adult + 
                                   3*FP_V3 + 3*FP_V3:GRS_Adult"),eform=T)
MalA18 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Adult+ 4*FP_V1 + 4*FP_V1:GRS_Adult +
                                   2*FP_V2 + 2*FP_V2:GRS_Adult + 
                                   3*FP_V3 + 3*FP_V3:GRS_Adult"),eform=T)

FemA25 <- lincom(Mod_TotFat_GRSF, c("FP_Intercept + FP_GRS_Adult + 4*FP_V1 + 4*FP_V1:GRS_Adult +
                                   2*FP_V2 + 2*FP_V2:GRS_Adult + 
                                   3*FP_V3 + 3*FP_V3:GRS_Adult + 
                                   7*FP_V4 + 7*FP_V4:GRS_Adult"),eform=T)
MalA25 <- lincom(Mod_TotFat_GRSM, c("FP_Intercept + FP_GRS_Adult + 4*FP_V1 + 4*FP_V1:GRS_Adult +
                                   2*FP_V2 + 2*FP_V2:GRS_Adult + 
                                   3*FP_V3 + 3*FP_V3:GRS_Adult + 
                                   7*FP_V4 + 7*FP_V4:GRS_Adult"),eform=T)

# Tables of results
lincommales <- rbind(Mal9, MalC9, MalA9,
                     Mal13, MalC13, MalA13,
                     Mal15, MalC15, MalA15,
                     Mal18, MalC18, MalA18,
                     Mal25, MalC25, MalA25)
lincommales$age <- c(9,9,9,13,13,13,15,15,15,18,18,18,25,25,25)
lincommales$GRS <- c('Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS', 
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS')
lincomfemales <- rbind(Fem9, FemC9, FemA9,
                       Fem13, FemC13, FemA13,
                       Fem15, FemC15, FemA15,
                       Fem18, FemC18, FemA18,
                       Fem25, FemC25, FemA25)
lincomfemales$age <- c(9,9,9,13,13,13,15,15,15,18,18,18,25,25,25)
lincomfemales$GRS <- c('Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS', 
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS')

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



write.table(femresults, file = "FemPerFat.tsv", row.names=FALSE, sep="\t")
write.table(malresults, file = "MalPerFat.tsv", row.names=FALSE, sep="\t")

#------------------------------------------------------------------------------#
# Calculate differences between mean GRS and +1SD GRS
# +1SD Child GRS calculations
FemC9 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Child"),eform=T)
MalC9 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Child"),eform=T)

FemC13 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child"),eform=T)
MalC13 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child"),eform=T)

FemC15 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child +
                                   2*FP_V2:GRS_Child"),eform=T)
MalC15 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child +
                                   2*FP_V2:GRS_Child"),eform=T)

FemC18 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child +
                                   2*FP_V2:GRS_Child + 
                                   3*FP_V3:GRS_Child"),eform=T)
MalC18 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child +
                                   2*FP_V2:GRS_Child + 
                                   3*FP_V3:GRS_Child"),eform=T)

FemC25 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child +
                                   2*FP_V2:GRS_Child + 
                                   3*FP_V3:GRS_Child + 
                                   7*FP_V4:GRS_Child"),eform=T)
MalC25 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Child +
                                   4*FP_V1:GRS_Child +
                                   2*FP_V2:GRS_Child + 
                                   3*FP_V3:GRS_Child + 
                                   7*FP_V4:GRS_Child"),eform=T)
# +1SD Adult GRS calculations
FemA9 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Adult"),eform=T)
MalA9 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Adult"),eform=T)

FemA13 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult"),eform=T)
MalA13 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult"),eform=T)

FemA15 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult +
                                   2*FP_V2:GRS_Adult"),eform=T)
MalA15 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult +
                                   2*FP_V2:GRS_Adult"),eform=T)

FemA18 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult +
                                   2*FP_V2:GRS_Adult + 
                                   3*FP_V3:GRS_Adult"),eform=T)
MalA18 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult +
                                   2*FP_V2:GRS_Adult + 
                                   3*FP_V3:GRS_Adult"),eform=T)

FemA25 <- lincom(Mod_TotFat_GRSF, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult +
                                   2*FP_V2:GRS_Adult + 
                                   3*FP_V3:GRS_Adult + 
                                   7*FP_V4:GRS_Adult"),eform=T)
MalA25 <- lincom(Mod_TotFat_GRSM, c("FP_GRS_Adult +
                                   4*FP_V1:GRS_Adult +
                                   2*FP_V2:GRS_Adult + 
                                   3*FP_V3:GRS_Adult + 
                                   7*FP_V4:GRS_Adult"),eform=T)

# Tables of results
lincommales <- rbind(MalC9, MalA9,
                     MalC13, MalA13,
                     MalC15, MalA15,
                     MalC18, MalA18,
                     MalC25, MalA25)
lincommales$age <- c(9,9,13,13,15,15,18,18,25,25)
lincommales$GRS <- c( '+1SD Child GRS', '+1SD Adult GRS', 
                      '+1SD Child GRS', '+1SD Adult GRS',
                      '+1SD Child GRS', '+1SD Adult GRS',
                      '+1SD Child GRS', '+1SD Adult GRS',
                      '+1SD Child GRS', '+1SD Adult GRS')
lincomfemales <- rbind(FemC9, FemA9,
                       FemC13, FemA13,
                       FemC15, FemA15,
                       FemC18, FemA18,
                       FemC25, FemA25)
lincomfemales$age <- c(9,9,13,13,15,15,18,18,25,25)
lincomfemales$GRS <- c( '+1SD Child GRS', '+1SD Adult GRS', 
                        '+1SD Child GRS', '+1SD Adult GRS',
                        '+1SD Child GRS', '+1SD Adult GRS',
                        '+1SD Child GRS', '+1SD Adult GRS',
                        '+1SD Child GRS', '+1SD Adult GRS')

# Calculate percentage changes
Pctfun <- function(x) {                                                  
  (x-1)*100
}

lincomfemales[ , c(1: 3)] <- apply(lincomfemales[ , c(1: 3)], 2, Pctfun)   
lincommales[ , c(1: 3)] <- apply(lincommales[ , c(1: 3)], 2, Pctfun)    

# Set levels of GRS values
lincomfemales$GRS <- factor(lincomfemales$GRS, 
                            levels = c('+1SD Child GRS', '+1SD Adult GRS'))
lincommales$GRS <- factor(lincommales$GRS, 
                          levels = c('+1SD Child GRS', '+1SD Adult GRS'))

# Format tables
data_reshapemale <- reshape(lincommales,                               
                            idvar = "age",
                            timevar = "GRS",
                            direction = "wide")
data_reshapefemale <- reshape(lincomfemales,                               
                              idvar = "age",
                              timevar = "GRS",
                              direction = "wide")
data_reshapemale$`+1SD Child GRS`<- sprintf("%.2f (%.2f,%.2f)", data_reshapemale$`Estimate.+1SD Child GRS`, 
                                            data_reshapemale$`2.5 %.+1SD Child GRS`, data_reshapemale$`97.5 %.+1SD Child GRS`)
data_reshapemale$`+1SD Adult GRS`<- sprintf("%.2f (%.2f,%.2f)", data_reshapemale$`Estimate.+1SD Adult GRS`, 
                                            data_reshapemale$`2.5 %.+1SD Adult GRS`, data_reshapemale$`97.5 %.+1SD Adult GRS`)
data_reshapefemale$`+1SD Child GRS`<- sprintf("%.2f (%.2f,%.2f)", data_reshapefemale$`Estimate.+1SD Child GRS`, 
                                              data_reshapefemale$`2.5 %.+1SD Child GRS`, data_reshapefemale$`97.5 %.+1SD Child GRS`)
data_reshapefemale$`+1SD Adult GRS`<- sprintf("%.2f (%.2f,%.2f)", data_reshapefemale$`Estimate.+1SD Adult GRS`, 
                                              data_reshapefemale$`2.5 %.+1SD Adult GRS`, data_reshapefemale$`97.5 %.+1SD Adult GRS`)
data_reshapemale <- data_reshapemale[,c('age', '+1SD Child GRS', '+1SD Adult GRS' )]  
data_reshapefemale <- data_reshapefemale[,c('age', '+1SD Child GRS', '+1SD Adult GRS' )]  

# Set axis values
#MaxVal <- max(c(lincomfemales$`97.5 %`, lincommales$`97.5 %`))
#MinVal <- min(c(lincomfemales$`2.5 %`, lincommales$`2.5 %`))
# need to be taken from the other set of results due to effect sizes 

Overallmax <-16
Overallmin <- -3

# Plot theme
My_Theme = theme(plot.title = element_text(size = 26,hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 axis.title.x = element_text(size = 18),
                 axis.text.x = element_text(size = 16),
                 axis.title.y = element_text(size = 18),
                 axis.text.y = element_text(size = 16),)

#Peripheral fat
lincommales$age <- factor(lincommales$age)
lincomfemales$age <- factor(lincomfemales$age)


ggplot(lincommales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
          geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), 
                        position=position_dodge(width=0.3),width=0.2) +
          geom_point(position=position_dodge(width=0.3)) + 
          ggtitle('Peripheral fat mass difference males') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
          geom_hline(yintercept = 0) + ylim(Overallmin, Overallmax) + My_Theme

ggplot(lincomfemales, aes(x=age, y=Estimate, group=GRS, color=GRS)) + 
          geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`),
                        position=position_dodge(width=0.3), width=0.2) +
          geom_point(position=position_dodge(width=0.3))+
          ggtitle('Peripheral fat mass difference females') +ylab('Percentage difference (%)')+xlab('Age (Years)')+
          geom_hline(yintercept = 0) + ylim(Overallmin, Overallmax) + My_Theme

write.table(data_reshapemale, file = "Results/MalPeriFatDiff.tsv", row.names=FALSE, sep="\t")
write.table(data_reshapefemale, file = "Results/FemPeriFatDiff.tsv", row.names=FALSE, sep="\t")

# Rate of change per spline
femmeans1 <- lincom(Mod_TotFat_GRSF, c("FP_V1"),eform=T)
malmeans1 <- lincom(Mod_TotFat_GRSM, c("FP_V1"),eform=T)
femmeans2 <- lincom(Mod_TotFat_GRSF, c("FP_V2"),eform=T)
malmeans2 <- lincom(Mod_TotFat_GRSM, c("FP_V2"),eform=T)
femmeans3 <- lincom(Mod_TotFat_GRSF, c("FP_V3"),eform=T)
malmeans3 <- lincom(Mod_TotFat_GRSM, c("FP_V3"),eform=T)
femmeans4 <- lincom(Mod_TotFat_GRSF, c("FP_V4"),eform=T)
malmeans4 <- lincom(Mod_TotFat_GRSM, c("FP_V4"),eform=T)

femmeans1C <- lincom(Mod_TotFat_GRSF, c("FP_V1 + FP_V1:GRS_Child"),eform=T)
malmeans1C <- lincom(Mod_TotFat_GRSM, c("FP_V1 + FP_V1:GRS_Child"),eform=T)
femmeans2C <- lincom(Mod_TotFat_GRSF, c("FP_V2 + FP_V2:GRS_Child"),eform=T)
malmeans2C <- lincom(Mod_TotFat_GRSM, c("FP_V2 + FP_V2:GRS_Child"),eform=T)
femmeans3C <- lincom(Mod_TotFat_GRSF, c("FP_V3 + FP_V3:GRS_Child"),eform=T)
malmeans3C <- lincom(Mod_TotFat_GRSM, c("FP_V3 + FP_V3:GRS_Child"),eform=T)
femmeans4C <- lincom(Mod_TotFat_GRSF, c("FP_V4 + FP_V4:GRS_Child"),eform=T)
malmeans4C <- lincom(Mod_TotFat_GRSM, c("FP_V4 + FP_V4:GRS_Child"),eform=T)

femmeans1A <- lincom(Mod_TotFat_GRSF, c("FP_V1 + FP_V1:GRS_Adult"),eform=T)
malmeans1A <- lincom(Mod_TotFat_GRSM, c("FP_V1 + FP_V1:GRS_Adult"),eform=T)
femmeans2A <- lincom(Mod_TotFat_GRSF, c("FP_V2 + FP_V2:GRS_Adult"),eform=T)
malmeans2A <- lincom(Mod_TotFat_GRSM, c("FP_V2 + FP_V2:GRS_Adult"),eform=T)
femmeans3A <- lincom(Mod_TotFat_GRSF, c("FP_V3 + FP_V3:GRS_Adult"),eform=T)
malmeans3A <- lincom(Mod_TotFat_GRSM, c("FP_V3 + FP_V3:GRS_Adult"),eform=T)
femmeans4A <- lincom(Mod_TotFat_GRSF, c("FP_V4 + FP_V4:GRS_Adult"),eform=T)
malmeans4A <- lincom(Mod_TotFat_GRSM, c("FP_V4 + FP_V4:GRS_Adult"),eform=T)

# Tables of results
lincommales <- rbind(malmeans1, malmeans1C, malmeans1A,
                     malmeans2, malmeans2C, malmeans2A,
                     malmeans3, malmeans3C, malmeans3A,
                     malmeans4, malmeans4C, malmeans4A)
lincommales$age <- c('9-13','9-13','9-13','13-15','13-15','13-15',
                     '15-18','15-18','15-18','18-25','18-25','18-25')
lincommales$GRS <- c('Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS', 
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                     'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS')
lincommales$Estimate <- (lincommales$Estimate -1) * 100
lincommales$`2.5 %` <- (lincommales$`2.5 %` -1) * 100
lincommales$`97.5 %` <- (lincommales$`97.5 %` -1) * 100



lincomfemales <- rbind(femmeans1, femmeans1C, femmeans1A,
                       femmeans2, femmeans2C, femmeans2A,
                       femmeans3, femmeans3C, femmeans3A,
                       femmeans4, femmeans4C, femmeans4A)
lincomfemales$age <- c('9-13','9-13','9-13','13-15','13-15','13-15',
                       '15-18','15-18','15-18','18-25','18-25','18-25')
lincomfemales$GRS <- c('Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS', 
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS',
                       'Mean Child & Adult GRS', '+1SD Child GRS', '+1SD Adult GRS')
lincomfemales$Estimate <- (lincomfemales$Estimate -1) * 100
lincomfemales$`2.5 %` <- (lincomfemales$`2.5 %` -1) * 100
lincomfemales$`97.5 %` <- (lincomfemales$`97.5 %` -1) * 100

femmean <- data.frame(
  Age = c('9-13', '13-15', '15-18','18-25'), 
  Mean = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[1],lincomfemales$`2.5 %`[1],lincomfemales$`97.5 %`[1]),
           sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[4],lincomfemales$`2.5 %`[4],lincomfemales$`97.5 %`[4]),
           sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[7],lincomfemales$`2.5 %`[7],lincomfemales$`97.5 %`[7]),
           sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[10],lincomfemales$`2.5 %`[10],lincomfemales$`97.5 %`[10])),
  Child = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[2],lincomfemales$`2.5 %`[2],lincomfemales$`97.5 %`[2]),
            sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[5],lincomfemales$`2.5 %`[5],lincomfemales$`97.5 %`[5]),
            sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[8],lincomfemales$`2.5 %`[8],lincomfemales$`97.5 %`[8]),
            sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[11],lincomfemales$`2.5 %`[11],lincomfemales$`97.5 %`[11])),
  Adult = c(sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[3],lincomfemales$`2.5 %`[3],lincomfemales$`97.5 %`[3]),
            sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[6],lincomfemales$`2.5 %`[6],lincomfemales$`97.5 %`[6]),
            sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[9],lincomfemales$`2.5 %`[9],lincomfemales$`97.5 %`[9]),
            sprintf("%.2f (%.2f,%.2f)", lincomfemales$Estimate[12],lincomfemales$`2.5 %`[12],lincomfemales$`97.5 %`[12])))

malmean <- data.frame(
  Age = c('9-13', '13-15', '15-18','18-25'), 
  Mean = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[1],lincommales$`2.5 %`[1],lincommales$`97.5 %`[1]),
           sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[4],lincommales$`2.5 %`[4],lincommales$`97.5 %`[4]),
           sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[7],lincommales$`2.5 %`[7],lincommales$`97.5 %`[7]),
           sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[10],lincommales$`2.5 %`[10],lincommales$`97.5 %`[10])),
  Child = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[2],lincommales$`2.5 %`[2],lincommales$`97.5 %`[2]),
            sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[5],lincommales$`2.5 %`[5],lincommales$`97.5 %`[5]),
            sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[8],lincommales$`2.5 %`[8],lincommales$`97.5 %`[8]),
            sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[11],lincommales$`2.5 %`[11],lincommales$`97.5 %`[11])),
  Adult = c(sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[3],lincommales$`2.5 %`[3],lincommales$`97.5 %`[3]),
            sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[6],lincommales$`2.5 %`[6],lincommales$`97.5 %`[6]),
            sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[9],lincommales$`2.5 %`[9],lincommales$`97.5 %`[9]),
            sprintf("%.2f (%.2f,%.2f)", lincommales$Estimate[12],lincommales$`2.5 %`[12],lincommales$`97.5 %`[12])))


if (yVal[7] == 'Total Fat Mass'){
  write.table(femmean, file = "FemTotFatrate.tsv", row.names=FALSE, sep="\t")
  write.table(malmean, file = "MalTotFatrate.tsv", row.names=FALSE, sep="\t")
}

write.table(femmean, file = "FemPeriFatrate.tsv", row.names=FALSE, sep="\t")
write.table(malmean, file = "MalPeriFatrate.tsv", row.names=FALSE, sep="\t")







