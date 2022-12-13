library(haven)
library(matrixStats)
library(tidyr)


# Load in ALSPAC data, and get relevant variables
ALSPAC <- read.csv('Data/ALSPAC_data.csv')
ALSPACd <- ALSPAC[!is.na(ALSPAC$kz021),] # Remove NAs

# Only individuals with GRS and DXA measures
ALSPACd <- ALSPACd[!is.na(ALSPACd$grs_childbmi_both_oct20_g1),]
ALSPACd[ALSPACd < 0] <- NA # Remove -10 (and other) NA placeholders
ALSPACd <- ALSPACd[(rowSums(is.na(ALSPACd[, c('f9dx135','fedx135','fg3254','fh2254','FJDX135','FKDX1001')])) !=6) ,]
# Individuals with GRS scores
ALSPACd <- ALSPACd[!is.na(ALSPACd$grs_childbmi_male_oct20_g1),]

# Only individuals without GRS and DXA measures
#ALSPACd <- ALSPACd[is.na(ALSPACd$grs_childbmi_both_oct20_g1),] 
#ALSPACd[ALSPACd < 0] <- NA # Remove -10 (and other) NA placeholders
#ALSPACd <- ALSPACd[(rowSums(is.na(ALSPACd[, c('f9dx135','fedx135','fg3254','fh2254','FJDX135','FKDX1001')])) ==6) ,]




fem <- ALSPACd[ALSPACd$kz021==2,]
mal <- ALSPACd[ALSPACd$kz021==1,]
fem$grs_childbmi_both_oct20_g1
#Variance explained by each GRS

FemFatGRS <- data.frame()
for (GRS in c('grs_childbmi_female_oct20_g1', 'grs_adultbmi_female_oct20_g1',
              'grs_childbmi_both_oct20_g1', 'grs_adultbmi_both_oct20_g1')){
  vals <- c(
    summary(lm(fem[[GRS]] ~ fem$f9dx135))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$fedx135))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$fg3254))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$fh2254))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$FJDX135))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$FKDX1001))[['r.squared']])
  
  FemFatGRS <- rbind(FemFatGRS, vals)

}

colnames(FemFatGRS) <- c('Age9', 'Age11', 'Age13', 'Age15', 'Age18', 'Age25')
rownames(FemFatGRS) <- c('grs_childbmi_female', 'grs_adultbmi_female',
                         'grs_childbmi_both', 'grs_adultbmi_both')
FemFatGRS <- signif(FemFatGRS, 2)
#------------------------------------------------------------------------------#

MalFatGRS <- data.frame()
for (GRS in c('grs_childbmi_male_oct20_g1', 'grs_adultbmi_male_oct20_g1',
              'grs_childbmi_both_oct20_g1', 'grs_adultbmi_both_oct20_g1')){
  vals <- c(
    summary(lm(mal[[GRS]] ~ mal$f9dx135))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$fedx135))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$fg3254))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$fh2254))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$FJDX135))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$FKDX1001))[['r.squared']])
  
  MalFatGRS <- rbind(MalFatGRS, vals)
  
}

colnames(MalFatGRS) <- c('Age9', 'Age11', 'Age13', 'Age15', 'Age18', 'Age25')
rownames(MalFatGRS) <- c('grs_childbmi_male', 'grs_adultbmi_male',
                         'grs_childbmi_both', 'grs_adultbmi_both')
MalFatGRS <- signif(MalFatGRS, 2)
#------------------------------------------------------------------------------#
FemLeanGRS <- data.frame()
for (GRS in c('grs_childbmi_female_oct20_g1', 'grs_adultbmi_female_oct20_g1',
              'grs_childbmi_both_oct20_g1', 'grs_adultbmi_both_oct20_g1')){
  vals <- c(
    summary(lm(fem[[GRS]] ~ fem$f9dx136))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$fedx136))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$fg3255))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$fh2255))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$FJDX136))[['r.squared']],
    summary(lm(fem[[GRS]] ~ fem$FKDX1002))[['r.squared']])
  
  FemLeanGRS <- rbind(FemLeanGRS, vals)
  
}

colnames(FemLeanGRS) <- c('Age9', 'Age11', 'Age13', 'Age15', 'Age18', 'Age25')
rownames(FemLeanGRS) <- c('grs_childbmi_female', 'grs_adultbmi_female',
                         'grs_childbmi_both', 'grs_adultbmi_both')
FemLeanGRS <- signif(FemLeanGRS, 2)
#------------------------------------------------------------------------------#

MalLeanGRS <- data.frame()
for (GRS in c('grs_childbmi_male_oct20_g1', 'grs_adultbmi_male_oct20_g1',
              'grs_childbmi_both_oct20_g1', 'grs_adultbmi_both_oct20_g1')){
  vals <- c(
    summary(lm(mal[[GRS]] ~ mal$f9dx136))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$fedx136))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$fg3255))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$fh2255))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$FJDX136))[['r.squared']],
    summary(lm(mal[[GRS]] ~ mal$FKDX1002))[['r.squared']])
  
  MalLeanGRS <- rbind(MalLeanGRS, vals)
  
}

colnames(MalLeanGRS) <- c('Age9', 'Age11', 'Age13', 'Age15', 'Age18', 'Age25')
rownames(MalLeanGRS) <- c('grs_childbmi_male', 'grs_adultbmi_male',
                         'grs_childbmi_both', 'grs_adultbmi_both')
MalLeanGRS <- signif(MalLeanGRS, 2)


MalFatGRS <- MalFatGRS*100
MalLeanGRS <- MalLeanGRS*100
FemFatGRS <- FemFatGRS*100
FemLeanGRS <- FemLeanGRS*100

write.csv(MalFatGRS, file = 'MalFatVariance.csv')
write.csv(MalLeanGRS, file = 'MalLeanVariance.csv')
write.csv(FemFatGRS, file = 'FemFatVariance.csv')
write.csv(FemLeanGRS, file = 'FemLeanVariance.csv')
#------------------------------------------------------------------------------#


# Grs values
# femgrschild <- scale(fem$grs_childbmi_female_oct20_g1)
# malgrschild <- scale(mal$grs_childbmi_male_oct20_g1)
# femgrsadult <- scale(fem$grs_adultbmi_female_oct20_g1)
# malgrsadult <- scale(mal$grs_adultbmi_male_oct20_g1)

mean(fem$grs_childbmi_female_oct20_g1) #0.0086
sd(fem$grs_childbmi_female_oct20_g1) #0.0005

mean(fem$grs_adultbmi_female_oct20_g1) #0.0079
sd(fem$grs_adultbmi_female_oct20_g1) #0.0004

mean(mal$grs_childbmi_male_oct20_g1) #0.0098
sd(mal$grs_childbmi_male_oct20_g1) #0.0008

mean(mal$grs_adultbmi_male_oct20_g1)# 0.0076
sd(mal$grs_adultbmi_male_oct20_g1) #0.0004

# sex combined GRS
mean(mal$grs_childbmi_both_oct20_g1)# 0.006
sd(mal$grs_childbmi_both_oct20_g1) #0.0003

mean(fem$grs_childbmi_both_oct20_g1)# 0.006
sd(fem$grs_childbmi_both_oct20_g1) #0.0002

mean(mal$grs_adultbmi_both_oct20_g1)# 0.006
sd(mal$grs_adultbmi_both_oct20_g1) #0.0002

mean(fem$grs_adultbmi_both_oct20_g1)# 0.006
sd(fem$grs_adultbmi_both_oct20_g1) #0.0002

# Sex combined GRS in both sexes
mean(ALSPACd$grs_childbmi_both_oct20_g1)# 0.006
sd(ALSPACd$grs_childbmi_both_oct20_g1) #0.0003

mean(ALSPACd$grs_adultbmi_both_oct20_g1)# 0.006
sd(ALSPACd$grs_adultbmi_both_oct20_g1) #0.0002


# age 9 DXA clinic measures # Total fat 
mean(fem$f9dx135, na.rm = T)
sd(fem$f9dx135, na.rm = T)
mean(mal$f9dx135, na.rm = T)
sd(mal$f9dx135, na.rm = T)
colSums(!is.na(fem['f9dx135']))
colSums(!is.na(mal['f9dx135']))


# age 25 DXA clinic measures # Total fat 
mean(fem$FKDX1001, na.rm = T)
sd(fem$FKDX1001, na.rm = T)
mean(mal$FKDX1001, na.rm = T)
sd(mal$FKDX1001, na.rm = T)
colSums(!is.na(fem['FKDX1001']))
colSums(!is.na(mal['FKDX1001']))

# age 9 DXA clinic measures # trunk fat
mean(fem$f9dx126, na.rm = T)
sd(fem$f9dx126, na.rm = T)
mean(mal$f9dx126, na.rm = T)
sd(mal$f9dx126, na.rm = T)
colSums(!is.na(fem['f9dx126']))
colSums(!is.na(mal['f9dx126']))

# age 25 DXA clinic measures # trunk fat
mean(fem$FKDX1031, na.rm = T)
sd(fem$FKDX1031, na.rm = T)
mean(mal$FKDX1031, na.rm = T)
sd(mal$FKDX1031, na.rm = T)
colSums(!is.na(fem['FKDX1031']))
colSums(!is.na(mal['FKDX1031']))

fem$DX1 <- rowSums(fem[,c('f9dx117', 'f9dx108')], na.rm=F)
fem$DX6 <- rowSums(fem[,c('FKDX1021', 'FKDX1011')], na.rm=F)
mal$DX1 <- rowSums(mal[,c('f9dx117', 'f9dx108')], na.rm=F)
mal$DX6 <- rowSums(mal[,c('FKDX1021', 'FKDX1011')], na.rm=F)

# age 9 DXA clinic measures # peripheral fat
mean(fem$DX1, na.rm = T)
sd(fem$DX1, na.rm = T)
mean(mal$DX1, na.rm = T)
sd(mal$DX1, na.rm = T)
colSums(!is.na(fem['DX1']))
colSums(!is.na(mal['DX1']))

# age 25 DXA clinic measures # peripheral fat
mean(fem$DX6, na.rm = T)
sd(fem$DX6, na.rm = T)
mean(mal$DX6, na.rm = T)
sd(mal$DX6, na.rm = T)
colSums(!is.na(fem['DX6']))
colSums(!is.na(mal['DX6']))


# age 9 DXA clinic measures # total lean
mean(fem$f9dx136, na.rm = T)
sd(fem$f9dx136, na.rm = T)
mean(mal$f9dx136, na.rm = T)
sd(mal$f9dx136, na.rm = T)
colSums(!is.na(fem['f9dx136']))
colSums(!is.na(mal['f9dx136']))


# age 25 DXA clinic measures # total lean
mean(fem$FKDX1002, na.rm = T)
sd(fem$FKDX1002, na.rm = T)
mean(mal$FKDX1002, na.rm = T)
sd(mal$FKDX1002, na.rm = T)
colSums(!is.na(fem['FKDX1002']))
colSums(!is.na(mal['FKDX1002']))


# Social class
femclass <- fem[,c('c755', 'c765')]
malclass <- mal[,c('c755', 'c765')]

# Make NA values a number (use value higher than max social class (6))
femclass[is.na(femclass)] <- 7
malclass[is.na(malclass)] <- 7

femclass$row_minimum = rowMins(as.matrix(femclass),na.rm = F)
malclass$row_minimum = rowMins(as.matrix(malclass),na.rm = F)

#remove rows where 7 is minimum (no data)
femclass <- femclass[femclass$row_minimum<7,]
femtot <- nrow(femclass)
malclass <- malclass[malclass$row_minimum<7,]
maltot <- nrow(malclass)


# summarise class 
femclass <- data.frame(table(femclass$row_minimum))
malclass <- data.frame(table(malclass$row_minimum))

femclass$percent <- (femclass$Freq/femtot)*100
malclass$percent <- (malclass$Freq/maltot)*100
#------------------------------------------------------------------------------#
# maternal education
femedclass <- fem['c645a']
maledclass <- mal['c645a']

femedclass <- femedclass %>% drop_na()
maledclass <- maledclass %>% drop_na()

femtot <- nrow(femedclass)
maltot <- nrow(maledclass)

femedclass <- data.frame(table(femedclass))
maledclass <- data.frame(table(maledclass))
femedclass$percent <- (femedclass$Freq/femtot)*100
maledclass$percent <- (maledclass$Freq/maltot)*100

#------------------------------------------------------------------------------#
# paternal education
femedclass <- fem['c666a']
maledclass <- mal['c666a']

femedclass <- femedclass %>% drop_na()
maledclass <- maledclass %>% drop_na()

femtot <- nrow(femedclass)
maltot <- nrow(maledclass)

femedclass <- data.frame(table(femedclass))
maledclass <- data.frame(table(maledclass))
femedclass$percent <- (femedclass$Freq/femtot)*100
maledclass$percent <- (maledclass$Freq/maltot)*100

#------------------------------------------------------------------------------#
# Smoking during pregnancy
femsmok <- fem['b665']
malsmok <- mal['b665']

femsmok <- femsmok %>% drop_na()
malsmok <- malsmok %>% drop_na()

femtot <- nrow(femsmok)
maltot <- nrow(malsmok)

femsmok <- data.frame(table(femsmok))
malsmok <- data.frame(table(malsmok))
femsmok$percent <- (femsmok$Freq/femtot)*100
malsmok$percent <- (malsmok$Freq/maltot)*100
#------------------------------------------------------------------------------#
# Peak height velocity
femapv <- fem['apv']
malapv <- mal['apv']
mean(femapv$apv, na.rm = T) #11.75
sd(femapv$apv, na.rm = T) #0.82
mean(malapv$apv, na.rm = T) #13.58
sd(malapv$apv, na.rm = T) #0.92

#------------------------------------------------------------------------------#
# Ethnicity 
femtot <- nrow(fem)
maltot <- nrow(mal)

maleth <- data.frame(table(mal$c804, useNA = 'ifany')/maltot)
femeth <- data.frame(table(fem$c804, useNA = 'ifany')/femtot)
#------------------------------------------------------------------------------#
#Birth weight
mean(fem$kz030, na.rm = T)
sd(fem$kz030, na.rm = T)
mean(mal$kz030, na.rm = T)
sd(mal$kz030, na.rm = T)




