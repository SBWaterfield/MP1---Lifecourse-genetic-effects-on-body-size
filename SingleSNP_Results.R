library(readxl)
library(haven)
library(data.table)
library(tidyverse)
library(psych)
library(ggpmisc)
library(ggplot2)

# Load in Examplar single SNP results: Not actally factor, just incorrect name saved
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP1/SingleSNP_Results.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP1/SingleSNP_ResultsSexCombo.RData")


FemTotFatDiff_SNP$SNP <- FemLean_SNP$SNP 
FemLeanDiff_SNP$SNP <- FemLean_SNP$SNP
rm(FemLean_SNP, FemTotFat_SNP, MalLean_SNP, MalTotFat_SNP)

# Get results in integer form
FemLeanDiff_SNP[c('PctChange', 'LI', 'UI')] <- t(sapply(strsplit(FemLeanDiff_SNP$`+SNP`, split = '[(,)]'), function(x) x))
MalLeanDiff_SNP[c('PctChange', 'LI', 'UI')] <- t(sapply(strsplit(MalLeanDiff_SNP$`+SNP`, split = '[(,)]'), function(x) x))
FemTotFatDiff_SNP[c('PctChange', 'LI', 'UI')] <- t(sapply(strsplit(FemTotFatDiff_SNP$`+SNP`, split = '[(,)]'), function(x) x))
MalTotFatDiff_SNP[c('PctChange', 'LI', 'UI')] <- t(sapply(strsplit(MalTotFatDiff_SNP$`+SNP`, split = '[(,)]'), function(x) x))


# Make values numeric
FemLeanDiff_SNP[c(5,6,7)] <- sapply(FemLeanDiff_SNP[c(5,6,7)], as.numeric)
MalLeanDiff_SNP[c(5,6,7)] <- sapply(MalLeanDiff_SNP[c(5,6,7)], as.numeric)
FemTotFatDiff_SNP[c(5,6,7)] <- sapply(FemTotFatDiff_SNP[c(5,6,7)], as.numeric)
MalTotFatDiff_SNP[c(5,6,7)] <- sapply(MalTotFatDiff_SNP[c(5,6,7)], as.numeric)

# Define Number of SNPS with 95% confidence level effects
FemLeanDiff_SNP$TrueEffect <- ifelse( (FemLeanDiff_SNP$LI <0 & FemLeanDiff_SNP$UI <0) | 
                                        (FemLeanDiff_SNP$LI >0 & FemLeanDiff_SNP$UI >0),
                                       T, F)
FemTotFatDiff_SNP$TrueEffect <- ifelse( (FemTotFatDiff_SNP$LI <0 & FemTotFatDiff_SNP$UI <0) | 
                                        (FemTotFatDiff_SNP$LI >0 & FemTotFatDiff_SNP$UI >0),
                                      T, F)
MalLeanDiff_SNP$TrueEffect <- ifelse( (MalLeanDiff_SNP$LI <0 & MalLeanDiff_SNP$UI <0) | 
                                        (MalLeanDiff_SNP$LI >0 & MalLeanDiff_SNP$UI >0),
                                      T, F)
MalTotFatDiff_SNP$TrueEffect <- ifelse( (MalTotFatDiff_SNP$LI <0 & MalTotFatDiff_SNP$UI <0) | 
                                          (MalTotFatDiff_SNP$LI >0 & MalTotFatDiff_SNP$UI >0),
                                        T, F)

# Load in Ricardson SNP details
SNPs_FemChild <- read_excel("SNPs_FemChild.xlsx")
SNPs_FemAdult <- read_excel("SNPs_FemAdult.xlsx")
SNPs_MalChild <- read_excel("SNPs_MalChild.xlsx")
SNPs_MalAdult <- read_excel("SNPs_MalAdult.xlsx")

FemIntersect <- intersect(SNPs_FemChild$SNP, SNPs_FemAdult$SNP)
MalIntersect <- intersect(SNPs_MalChild$SNP, SNPs_MalAdult$SNP)
AllIntersect <- intersect(FemIntersect, MalIntersect)

# SNPs in both GRS with effects at multiple time points
MalFatSNPMulitEffect <- MalTotFatDiff_SNP[(MalTotFatDiff_SNP$TrueEffect==T & MalTotFatDiff_SNP$SNP %in% MalIntersect),]
FemFatSNPMulitEffect <- FemTotFatDiff_SNP[(FemTotFatDiff_SNP$TrueEffect==T & FemTotFatDiff_SNP$SNP %in% FemIntersect),]
MalLeanSNPMulitEffect <- MalLeanDiff_SNP[(MalLeanDiff_SNP$TrueEffect==T & MalLeanDiff_SNP$SNP %in% MalIntersect),]
FemLeanSNPMulitEffect <- FemLeanDiff_SNP[(FemLeanDiff_SNP$TrueEffect==T & FemLeanDiff_SNP$SNP %in% FemIntersect),]


# Get GRS specific results 
FemChildLeanDiff_SNP <- FemLeanDiff_SNP[FemLeanDiff_SNP$SNP %in% SNPs_FemChild$SNP,]
FemAdultLeanDiff_SNP <- FemLeanDiff_SNP[FemLeanDiff_SNP$SNP %in% SNPs_FemAdult$SNP,]
FemChildTotFatDiff_SNP <- FemTotFatDiff_SNP[FemTotFatDiff_SNP$SNP %in% SNPs_FemChild$SNP,]
FemAdultTotFatDiff_SNP <- FemTotFatDiff_SNP[FemTotFatDiff_SNP$SNP %in% SNPs_FemAdult$SNP,]

MalChildLeanDiff_SNP <- MalLeanDiff_SNP[MalLeanDiff_SNP$SNP %in% SNPs_MalChild$SNP,]
MalAdultLeanDiff_SNP <- MalLeanDiff_SNP[MalLeanDiff_SNP$SNP %in% SNPs_MalAdult$SNP,]
MalChildTotFatDiff_SNP <- MalTotFatDiff_SNP[MalTotFatDiff_SNP$SNP %in% SNPs_MalChild$SNP,]
MalAdultTotFatDiff_SNP <- MalTotFatDiff_SNP[MalTotFatDiff_SNP$SNP %in% SNPs_MalAdult$SNP,]

# SNP longitudinal effects
ggplot(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$SNP %in% (FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_point() +  geom_line() + geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Female child SNPs on total fat mass') + 
  stat_poly_eq(aes(x=age, y=PctChange), inherit.aes = F)

ggplot(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$SNP %in% (FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Female child SNPs on total fat mass') 

ggplot(FemChildTotFatDiff_SNP,
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Female child SNPs on total fat mass') 




ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$SNP %in% (FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_point() +  geom_line() + geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Female adult SNPs on total fat mass') + 
  stat_poly_eq(aes(x=age, y=PctChange), inherit.aes = F)

ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$SNP %in% (FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Female adult SNPs on total fat mass') 

ggplot(FemAdultTotFatDiff_SNP,
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Female adult SNPs on total fat mass') 



ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$SNP %in% (MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_point() +  geom_line() + geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) +  
  theme(legend.position = 'none') + ggtitle('Male child SNPs on total fat mass') + 
  stat_poly_eq(aes(x=age, y=PctChange), inherit.aes = F)

ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$SNP %in% (MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Male child SNPs on total fat mass') 

ggplot(MalChildTotFatDiff_SNP,
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Male child SNPs on total fat mass') 



ggplot(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$SNP %in% (MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_point() +  geom_line() + geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) +  
  theme(legend.position = 'none') + ggtitle('Male adult SNPs on total fat mass') + 
  stat_poly_eq(aes(x=age, y=PctChange), inherit.aes = F)

ggplot(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$SNP %in% (MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$TrueEffect==T,]$SNP),],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Male adult SNPs on total fat mass') 

ggplot(MalAdultTotFatDiff_SNP,
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_smooth(aes(x=age, y=PctChange), method = 'loess', inherit.aes = F) + 
  theme(legend.position = 'none') + ggtitle('Male adult SNPs on total fat mass') 
# # Format to wide format
# FemChildLeanDiff_SNP <- FemChildLeanDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# FemAdultLeanDiff_SNP <- FemAdultLeanDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# FemChildTotFatDiff_SNP <- FemChildTotFatDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# FemAdultTotFatDiff_SNP <- FemAdultTotFatDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# MalChildLeanDiff_SNP <- MalChildLeanDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# MalAdultLeanDiff_SNP <- MalAdultLeanDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# MalChildTotFatDiff_SNP <- MalChildTotFatDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# MalAdultTotFatDiff_SNP <- MalAdultTotFatDiff_SNP %>% select('age', 'SNP', '+SNP') %>% reshape(idvar = "SNP", timevar = "age", sep = "", direction = "wide", varying = c('Age9', 'Age13', 'Age15', 'Age18', 'Age25'))
# 
# #Add in Annotation
# FemChildLeanDiff_SNP <- merge(FemChildLeanDiff_SNP, SNPs_FemChild[, c("SNP", "Closest gene")], by="SNP")
# FemChildTotFatDiff_SNP <- merge(FemChildTotFatDiff_SNP, SNPs_FemChild[, c("SNP", "Closest gene")], by="SNP")
# MalChildLeanDiff_SNP <- merge(MalChildLeanDiff_SNP, SNPs_MalChild[, c("SNP", "Closest gene")], by="SNP")
# MalChildTotFatDiff_SNP <- merge(MalChildTotFatDiff_SNP, SNPs_MalChild[, c("SNP", "Closest gene")], by="SNP")
# FemAdultLeanDiff_SNP <- merge(FemAdultLeanDiff_SNP, SNPs_FemAdult[, c("SNP", "Closest gene")], by="SNP")
# FemAdultTotFatDiff_SNP <- merge(FemAdultTotFatDiff_SNP, SNPs_FemAdult[, c("SNP", "Closest gene")], by="SNP")
# MalAdultLeanDiff_SNP <- merge(MalAdultLeanDiff_SNP, SNPs_MalAdult[, c("SNP", "Closest gene")], by="SNP")
# MalAdultTotFatDiff_SNP <- merge(MalAdultTotFatDiff_SNP, SNPs_MalAdult[, c("SNP", "Closest gene")], by="SNP")

#EXport to CSV
# write.csv(FemChildTotFatDiff_SNP, row.names = F, file = 'FemchildTotFatDiff_SNP.csv')
# write.csv(FemAdultTotFatDiff_SNP, row.names = F, file = 'FemAdultTotFatDiff_SNP.csv')
# write.csv(FemChildLeanDiff_SNP, row.names = F, file = 'FemchildLeanDiff_SNP.csv')
# write.csv(FemAdultLeanDiff_SNP, row.names = F, file = 'FemAdultLeanDiff_SNP.csv')
# write.csv(MalChildTotFatDiff_SNP, row.names = F, file = 'MalchildTotFatDiff_SNP.csv')
# write.csv(MalAdultTotFatDiff_SNP, row.names = F, file = 'MalAdultTotFatDiff_SNP.csv')
# write.csv(MalChildLeanDiff_SNP, row.names = F, file = 'MalchildLeanDiff_SNP.csv')
# write.csv(MalAdultLeanDiff_SNP, row.names = F, file = 'MalAdultLeanDiff_SNP.csv')

# Clean up workspace
rm(SNPs_MalAdult, SNPs_MalChild, SNPs_FemAdult, SNPs_FemChild, 
   FemLeanDiff_SNP, FemTotFatDiff_SNP, MalLeanDiff_SNP, MalTotFatDiff_SNP,
   MalFatSNPMulitEffect, MalLeanSNPMulitEffect, FemFatSNPMulitEffect, FemLeanSNPMulitEffect)

# Summarise Single SNP effects
FemChildTotFatDiff_SNP_Summary <- FemChildTotFatDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

FemAdultTotFatDiff_SNP_Summary <- FemAdultTotFatDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

MalChildTotFatDiff_SNP_Summary <- MalChildTotFatDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

MalAdultTotFatDiff_SNP_Summary <- MalAdultTotFatDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

FemChildLeanDiff_SNP_Summary <- FemChildLeanDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

FemAdultLeanDiff_SNP_Summary <- FemAdultLeanDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

MalChildLeanDiff_SNP_Summary <- MalChildLeanDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

MalAdultLeanDiff_SNP_Summary <- MalAdultLeanDiff_SNP %>%
  group_by(SNP) %>%
  summarise(Freq = sum(TrueEffect))

sum(FemChildTotFatDiff_SNP$TrueEffect) #76
sum(FemAdultTotFatDiff_SNP$TrueEffect) #67
sum(MalChildTotFatDiff_SNP$TrueEffect) #58
sum(MalAdultTotFatDiff_SNP$TrueEffect) #62

length(unique(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$TrueEffect==T,]$SNP)) #38
length(unique(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$TrueEffect==T,]$SNP)) #37
length(unique(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$TrueEffect==T,]$SNP)) #24
length(unique(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$TrueEffect==T,]$SNP)) #30



#------------------------------------------------------------------------------#
limits <- ylim(-40,40) 

#female child total fat
nrow(FemChildTotFatDiff_SNP)/5
FemChildTotFatDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Female Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Female adult total fat
nrow(FemAdultTotFatDiff_SNP)/5
FemAdultTotFatDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Female Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Male child total fat
nrow(MalChildTotFatDiff_SNP)/5
MalChildTotFatDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Male Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Male adult total fat
nrow(MalAdultTotFatDiff_SNP)/5
MalAdultTotFatDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Male Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Female child lean
FemChildLeanDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(FemChildLeanDiff_SNP[FemChildLeanDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total lean mass (Female Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Female adult lean
FemAdultLeanDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(FemAdultLeanDiff_SNP[FemAdultLeanDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total lean mass (Female Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Male child lean
MalChildLeanDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(MalChildLeanDiff_SNP[MalChildLeanDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total lean mass (Male Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

# Male adult lean
MalAdultLeanDiff_SNP %>%
  group_by(age) %>%
  summarise(Freq = sum(TrueEffect))

ggplot(MalAdultLeanDiff_SNP[MalAdultLeanDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total lean mass (Male Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')




# Mess around with data presemtation
ggplot(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Female Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Male Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Female Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')


ggplot(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$TrueEffect==T,],
       aes(x=age, y=PctChange, group=SNP, color=SNP)) + 
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Male Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$SNP=='rs12364470',],
       aes(x=age, y=PctChange)) + geom_line() +
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Female Adult GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$SNP=='rs115597956',],
       aes(x=age, y=PctChange)) + geom_line() +
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total fat mass (Male Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')

ggplot(MalChildLeanDiff_SNP[MalChildLeanDiff_SNP$SNP=='rs115597956',],
       aes(x=age, y=PctChange)) + geom_line() +
  geom_errorbar(aes(ymin=LI, ymax=UI), 
                position=position_dodge(width=1),width=0.2) +
  geom_point(position=position_dodge(width=1)) + 
  geom_hline(yintercept = 0) +   theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Single SNP effect on total lean mass (Male Child GRS)') + limits +
  xlab('Age (Years)') + ylab('Change from mean model (%)')


# Fat Mass at age 9 
ggplot(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$age==9,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black") + ylim(0,25) + xlim(-18,18) +
  ggtitle('Female Child SNPs effect on total fat mass (Age 9)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$age==9,]) * 0.5)), inherit.aes = F)

ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$age==9,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black")+  ylim(0,25) + xlim(-18,18) +
  ggtitle('Female Adult SNPs effect on total fat mass (Age 9)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$age==9,]) * 0.5)), inherit.aes = F)


ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$age==9,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black")+ ylim(0,25) + xlim(-18,18) +
  ggtitle('Male Child SNPs effect on total fat mass (Age 9)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$age==9,]) * 0.5)), inherit.aes = F)

ggplot(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$age==9,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black")+ ylim(0,25) + xlim(-18,18) +
  ggtitle('Male Adult SNPs effect on total fat mass (Age 9)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$age==9,]) * 0.5)), inherit.aes = F)

# Fat Mass at age 25
ggplot(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$age==25,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black") + ylim(0,25) + xlim(-18,18) +
  ggtitle('Female Child SNPs effect on total fat mass (Age 25)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(FemChildTotFatDiff_SNP[FemChildTotFatDiff_SNP$age==25,]) * 0.5)), inherit.aes = F)

ggplot(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$age==25,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black")+  ylim(0,25) + xlim(-18,18) +
  ggtitle('Female Adult SNPs effect on total fat mass (Age 25)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(FemAdultTotFatDiff_SNP[FemAdultTotFatDiff_SNP$age==25,]) * 0.5)), inherit.aes = F)


ggplot(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$age==25,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black")+ ylim(0,25) + xlim(-18,18) +
  ggtitle('Male Child SNPs effect on total fat mass (Age 25)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(MalChildTotFatDiff_SNP[MalChildTotFatDiff_SNP$age==25,]) * 0.5)), inherit.aes = F)

ggplot(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$age==25,], aes(x=PctChange, fill= TrueEffect)) + 
  geom_histogram(aes(), binwidth=.5,colour="black")+ ylim(0,25) + xlim(-18,18) +
  ggtitle('Male Adult SNPs effect on total fat mass (Age 25)') + xlab('Percentage Change')+ylab('count') +
  geom_vline(xintercept = 0) +   theme(plot.title = element_text(hjust = 0.5))+
  geom_density(aes(x=PctChange, y = ..density.. * (nrow(MalAdultTotFatDiff_SNP[MalAdultTotFatDiff_SNP$age==25,]) * 0.5)), inherit.aes = F)


reshape(test, idvar = "SNP", timevar = "age", sep = "", direction = "wide")
