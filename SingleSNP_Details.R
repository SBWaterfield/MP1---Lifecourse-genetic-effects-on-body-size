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

filenames <- list.files('Data/SNPs_HRCSexCombo', full.names = T)
SNPs <- meffil.extract.genotypes(filenames)
SNPs <- t(SNPs)
SNPs <- data.frame(SNPs)

SNPs_ComboChild$Present <- ifelse(SNPs_ComboChild$SNP %in% colnames(SNPs), T, F)
SNPs_ComboAdult$Present <- ifelse(SNPs_ComboAdult$SNP %in% colnames(SNPs), T, F)

sum(!duplicated(SNPs_ComboChild$SNP)) #298
sum(!duplicated(SNPs_ComboAdult$SNP)) #557
sum(SNPs_ComboChild$Present) # 279
sum(SNPs_ComboAdult$Present) # 477

#------------------------------------------------------------------------------#
filenames <- list.files('Data/SNPs_HRC', full.names = T)
SNPs <- meffil.extract.genotypes(filenames)
SNPs <- t(SNPs)
SNPs <- data.frame(SNPs)

SNPs_FemChild$Present <- ifelse(SNPs_FemChild$SNP %in% colnames(SNPs), T, F)
SNPs_FemAdult$Present <- ifelse(SNPs_FemAdult$SNP %in% colnames(SNPs), T, F)
SNPs_MalChild$Present <- ifelse(SNPs_MalChild$SNP %in% MalChildTotFatDiff_SNP$SNP, T, F)
SNPs_MalAdult$Present <- ifelse(SNPs_MalAdult$SNP %in% colnames(SNPs), T, F)

sum(!duplicated(SNPs_MalChild$SNP)) #68
sum(!duplicated(SNPs_FemChild$SNP)) #138
sum(!duplicated(SNPs_MalAdult$SNP)) #159
sum(!duplicated(SNPs_FemAdult$SNP)) #215


sum(SNPs_FemChild$Present) # 134
sum(SNPs_FemAdult$Present) # 187
sum(SNPs_MalChild$Present) # 60
sum(SNPs_MalAdult$Present) # 133

#------------------------------------------------------------------------------#
# Export data to XLSX
SNPs_FemChild$PositiveEffectSNP <- ifelse(SNPs_FemChild$Beta<0, SNPs_FemChild$`Other allele`, SNPs_FemChild$`Effect allele`)
SNPs_FemAdult$PositiveEffectSNP <- ifelse(SNPs_FemAdult$Beta<0, SNPs_FemAdult$`Other allele`, SNPs_FemAdult$`Effect allele`)
SNPs_MalChild$PositiveEffectSNP <- ifelse(SNPs_MalChild$Beta<0, SNPs_MalChild$`Other allele`, SNPs_MalChild$`Effect allele`)
SNPs_MalAdult$PositiveEffectSNP <- ifelse(SNPs_MalAdult$Beta<0, SNPs_MalAdult$`Other allele`, SNPs_MalAdult$`Effect allele`)
SNPs_ComboChild$PositiveEffectSNP <- ifelse(SNPs_ComboChild$Beta<0, SNPs_ComboChild$`Other allele`, SNPs_ComboChild$`Effect allele`)
SNPs_ComboAdult$PositiveEffectSNP <- ifelse(SNPs_ComboAdult$Beta<0, SNPs_ComboAdult$`Other allele`, SNPs_ComboAdult$`Effect allele`)

SNPs_FemChild <- SNPs_FemChild[c('SNP', 'PositiveEffectSNP', 'Present')]
SNPs_FemAdult <- SNPs_FemAdult[c('SNP', 'PositiveEffectSNP', 'Present')]
SNPs_MalChild <- SNPs_MalChild[c('SNP', 'PositiveEffectSNP', 'Present')]
SNPs_MalAdult <- SNPs_MalAdult[c('SNP', 'PositiveEffectSNP', 'Present')]
SNPs_ComboChild <- SNPs_ComboChild[c('SNP', 'PositiveEffectSNP', 'Present')]
SNPs_ComboAdult <- SNPs_ComboAdult[c('SNP', 'PositiveEffectSNP', 'Present')]


library(openxlsx)

dataset_names <- list('SNPs_FemChild' = SNPs_FemChild, 'SNPs_FemAdult' = SNPs_FemAdult,
                      'SNPs_MalChild' = SNPs_MalChild, 'SNPs_MalAdult' = SNPs_MalChild,
                      'SNPs_ComboChild' = SNPs_ComboChild, 'SNPs_ComboAdult' = SNPs_ComboAdult)
write.xlsx(dataset_names, file = 'MP1_SI1.xlsx')





























