library(readxl)
library(haven)
library(data.table)
library(tidyverse)

SNPs_FemChild <- read_excel("SNPs_FemChild.xlsx")
SNPs_FemAdult <- read_excel("SNPs_FemAdult.xlsx")
SNPs_MalChild <- read_excel("SNPs_MalChild.xlsx")
SNPs_MalAdult <- read_excel("SNPs_MalAdult.xlsx")
SNPs_ComboChild <- read_excel("SNPs_ComboChild.xlsx")
SNPs_ComboAdult <- read_excel("SNPs_ComboAdult.xlsx")

FemIntersect <- intersect(SNPs_FemChild$SNP, SNPs_FemAdult$SNP)
MalIntersect <- intersect(SNPs_MalChild$SNP, SNPs_MalAdult$SNP)

ComboSNPlist <- c(SNPs_ComboChild$SNP, SNPs_ComboAdult$SNP)
SNPlist <- c(SNPs_FemChild$SNP, SNPs_FemAdult$SNP, SNPs_MalChild$SNP, SNPs_MalAdult$SNP)
SNPlistCombo <- c(SNPs_ComboChild$SNP, SNPs_ComboAdult$SNP)
Betalist <- c(SNPs_FemChild$Beta, SNPs_FemAdult$Beta, SNPs_MalChild$Beta, SNPs_MalAdult$Beta)
BetalistCombo <- c(SNPs_ComboChild$Beta, SNPs_ComboAdult$Beta)
Effectlist <- c(SNPs_FemChild$`Effect allele`, SNPs_FemAdult$`Effect allele`, SNPs_MalChild$`Effect allele`, SNPs_MalAdult$`Effect allele`)
EffectlistCombo <- c(SNPs_ComboChild$`Effect allele`, SNPs_ComboAdult$`Effect allele`)
Otherlist <- c(SNPs_FemChild$`Other allele`, SNPs_FemAdult$`Other allele`, SNPs_MalChild$`Other allele`, SNPs_MalAdult$`Other allele`)
OtherlistCombo <- c(SNPs_ComboChild$`Other allele`, SNPs_ComboAdult$`Other allele`)

RichardsonSNP <- data.frame(SNPlist, Betalist, Effectlist, Otherlist)
RichardsonSNP <- RichardsonSNP[!duplicated(RichardsonSNP$SNPlist),]
RichardsonSNP$PosSNP <- ifelse(RichardsonSNP$Betalist<0, RichardsonSNP$Otherlist, RichardsonSNP$Effectlist)

RichardsonSNPCombo <- data.frame(SNPlistCombo, BetalistCombo, EffectlistCombo, OtherlistCombo)
RichardsonSNPCombo <- RichardsonSNPCombo[!duplicated(RichardsonSNPCombo$SNPlist),]
RichardsonSNPCombo$PosSNP <- ifelse(RichardsonSNPCombo$Betalist<0, RichardsonSNPCombo$Otherlist, RichardsonSNPCombo$Effectlist)


# Get SNP names for external SNP extraction
SNPs <- c(SNPs_FemChild$SNP, SNPs_FemAdult$SNP, SNPs_MalChild$SNP, SNPs_MalAdult$SNP)
SNPs <- unique(SNPs) #545
SNPs <- as.matrix(SNPs)
#write.csv(RichardsonSNP$SNPlist, file='snp-names_HRC.txt', row.names = F)
#write.csv(RichardsonSNP[c('SNPlist', 'PosSNP')], file = 'snp-names_HRCRef.txt', row.names = F)
#ComboSNPs <- unique(ComboSNPlist) 
#write.csv(RichardsonSNPCombo$SNPlist, file = 'snp-names_SexCombo.txt', row.names = F)
#write.csv(RichardsonSNPCombo[c('SNPlistCombo', 'PosSNP')], file = 'snp-names_SexComboRef.txt', row.names = F)

# Function to deal with extracted genotype data
meffil.extract.genotypes <- function(filenames, verbose=F) {
  stopifnot(all(sapply(filenames, file.exists)))
  
  table.list <- lapply(filenames, function(filename) {
    read.table(filename, header=T)
  })
  
  sample.names <- table.list[[1]][,1] ## family id
  genotypes <- lapply(table.list, function(genotype.table) {
    genotype.table[match(sample.names, genotype.table[,1]), ## match family ids
                   -(1:6),
                   drop=F]
  })
  genotypes <- do.call(cbind, genotypes)
  colnames(genotypes) <- sub("_.*", "", colnames(genotypes)) 
  rownames(genotypes) <- sample.names
  
  t(as.matrix(genotypes))
}

filenames <- list.files('Data/SNPs_HRC', full.names = T)
filenames <- list.files('Data/SNPs_HRCSexCombo', full.names = T)
SNPs <- meffil.extract.genotypes(filenames)
SNPs <- t(SNPs)
SNPs <- data.frame(SNPs)
SNPs$ID <- rownames(SNPs)
SNPs$ID <- gsub('X', '', SNPs$ID)
# Get IDs that end in a A or B (kids only) 
SNPs <- subset(SNPs,grepl("^.+(A|B)$",ID))



# Load in ALSPAC data
cp_3a <- read_dta("Data/cp_3a.dta")
F9 <- read_dta('Data/f09.dta')
F11 <- read_dta('Data/f11.dta')
FTF2 <- read_dta('Data/fTF2.dta')
FTF3 <- read_dta('Data/fTF3.dta')
FTF4 <- read_dta('Data/fTF4.dta')
F24 <- read_dta('Data/f24.dta')

ALSPACSingleSNP <- merge(F9, F11, by='aln_qlet')
ALSPACSingleSNP <- merge(ALSPACSingleSNP, FTF2, by='aln_qlet')
ALSPACSingleSNP <- merge(ALSPACSingleSNP, FTF3, by='aln_qlet')
ALSPACSingleSNP <- merge(ALSPACSingleSNP, FTF4, by='aln_qlet')
ALSPACSingleSNP <- merge(ALSPACSingleSNP, F24, by='aln_qlet')
ALSPACSingleSNP <- merge(ALSPACSingleSNP, cp_3a, by='aln_qlet')

ALSPACSingleSNP$ID <- gsub('_', '', ALSPACSingleSNP$aln_qlet)


# Merge and save data
ALSPACSingleSNP <- merge(ALSPACSingleSNP, SNPs, by='ID')

#Reduce reference down to only Available SNPs
RichardsonSNP <- RichardsonSNP[RichardsonSNP$SNPlist %in% colnames(SNPs),]
RichardsonSNPCombo <- RichardsonSNPCombo[RichardsonSNPCombo$SNPlistCombo %in% colnames(SNPs),]


save(ALSPACSingleSNP, RichardsonSNP, file = 'ALSPACSingleSNP.Rdata')
save(ALSPACSingleSNP, RichardsonSNPCombo, file = 'ALSPACSingleSNPSexCombo.Rdata')




#------------------------------------------------------------------------------#
# Correlations of BMI and fat mass

ALSPACSingleSNP[ALSPACSingleSNP < 0] <- NA 


# Age 9
cor.test(ALSPACSingleSNP$f9ms026a, ALSPACSingleSNP$f9dx125) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$f9ms026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$f9dx125) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$f9ms026a,
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$f9dx125) #

# Age 11
cor.test(ALSPACSingleSNP$fems026a, ALSPACSingleSNP$fedx135) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fems026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fedx135) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fems026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fedx135) #

# Age 13
cor.test(ALSPACSingleSNP$fg3139, ALSPACSingleSNP$fg3254) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fg3139, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fg3254) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fg3139, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fg3254) #

# Age 15
cor.test(ALSPACSingleSNP$fh3019, ALSPACSingleSNP$fh2254) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fh3019, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fh2254) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fh3019, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fh2254) #

# Age 18
cor.test(ALSPACSingleSNP$FJMR022a, ALSPACSingleSNP$FJDX135) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FJMR022a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FJDX135) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FJMR022a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FJDX135) #

# Age 25
cor.test(ALSPACSingleSNP$FKMS1040, ALSPACSingleSNP$FKDX1001) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FKMS1040, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FKDX1001) #
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FKMS1040, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FKDX1001) #


# Correlations of BMI and lean mass

# Age 9
cor.test(ALSPACSingleSNP$f9ms026a, ALSPACSingleSNP$f9dx136) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$f9ms026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$f9dx136) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$f9ms026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$f9dx136) 

# Age 11
cor.test(ALSPACSingleSNP$fems026a, ALSPACSingleSNP$fedx136) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fems026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fedx136) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fems026a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fedx136) 

# Age 13
cor.test(ALSPACSingleSNP$fg3139, ALSPACSingleSNP$fg3255) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fg3139, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fg3255) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fg3139, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fg3255) 

# Age 15
cor.test(ALSPACSingleSNP$fh3019, ALSPACSingleSNP$fh2255) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fh3019, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fh2255) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fh3019, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fh2255) 

# Age 18
cor.test(ALSPACSingleSNP$FJMR022a, ALSPACSingleSNP$FJDX136) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FJMR022a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FJDX136) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FJMR022a, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FJDX136) 

# Age 25
cor.test(ALSPACSingleSNP$FKMS1040, ALSPACSingleSNP$FKDX1002) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FKMS1040, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$FKDX1002) 
cor.test(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FKMS1040, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$FKDX1002) 



plot(ALSPACSingleSNP$fg3139, ALSPACSingleSNP$fg3255) 
plot(ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fg3139, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==1,]$fg3255) 
plot(ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fg3139, 
         ALSPACSingleSNP[ALSPACSingleSNP$kz021==2,]$fg3255) 





test <- ALSPACSingleSNP[c('kz021', 'FJDX136', 'FJMR022a')]
colnames(test) <- c('Sex', 'Lean Mass', 'BMI')
test$Sex <- as.factor(ifelse(test$Sex==2, 'Female', 'Male'))
ggplot(test, aes(x = BMI, y = `Lean Mass`)) +
  geom_point(aes(color=Sex)) + 
  geom_smooth(method = "lm", aes(color=Sex), fullrange=F, level=0.95) + 
  geom_smooth(method='lm',level=0.95) +
  xlab('BMI') + ylab('Total Lean Mass (g)') + ggtitle('Correlation between BMI and lean mass')











