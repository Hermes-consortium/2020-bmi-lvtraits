library(devtools)
require(data.table)
library(TwoSampleMR)
library(readr)
library(dplyr)
library(MRMix)

# devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")
# install_github("mrcieu/ieugwasr")

memory.limit(size = 56000)

# setwd("data/BMI_WHR_GIANT/")
BMI <- fread("data/BMI_WHR_GIANT/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt",key = "SNP")

# setwd("~/MSc dissertation/Data/LVtraits_AungN2019")
lvmvr<- fread("data/LVtraits_AungN2019/LVMVR/AungN_31554410_LVMVR.txt.gz",key = "SNP")
#lvef<- fread("LVEF-LV ejection fraction.txt",key = "SNP")
#lvesv<- fread("LVESV-LV end-systolic volume.txt",key = "SNP")
#lvm<- fread("LVM-LV mass.txt",key = "SNP")
#lvedv<- fread("LVEDV-LV end-diastolic volume.txt",key = "SNP")

# Get instruments (worked)
exposure_dat <- read_exposure_data("data/BMI_WHR_GIANT/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt", 
                                   snp_col = "SNP",beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "Tested_Allele",
                                   other_allele_col = "Other_Allele",
                                   samplesize_col = "N", 
                                   eaf_col = "Freq_Tested_Allele")

#write.table(bmi, file = "bmi.txt",row.names = FALSE, col.names = TRUE)
BMI <- fread("bmi.txt",key = "SNP")

# Get effects of instruments on outcome (does not work) 
outcome_dat <- read_outcome_data("data/LVtraits_AungN2019/LVEF/AungN_31554410_LVEF.txt", 
                                 sep="\t",
                                 snp_col = "SNP",beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "EA",
                                 other_allele_col = "NEA",
                                 eaf_col = "EAF",
                                 pval_col = "P")


make_DT.MR <- function(DT.exposure, DT.outcome){
  DT.x <- DT.exposure[,.(SNP, Tested_Allele, BETA, SE, P)]
  DT.y <- DT.outcome[,.(SNP, EA, BETA, SE, P)]
  DT.MR <- merge(DT.x, DT.y, by="SNP")
  DT.MR[, BETA.y.harmonised := ifelse( Tested_Allele== EA, BETA.y, -BETA.y)]
  DT.MR
}

DT.MR_FA <- make_DT.MR(BMI, lvmvr)

write.table(DT.MR_FA, file = "bmi_lvmvr_h.txt",row.names = FALSE, col.names = TRUE)

run_MR <- function(DT.MR, ...){
  MR_res <- mr_allmethods(mr_input(bx = DT.MR$BETA.x,
                                   bxse = DT.MR$SE.x,
                                   by = DT.MR$BETA.y.harmonised,
                                   byse = DT.MR$SE.y,
                                   ...))
  MR_res
}

MR_lvmvr <- run_MR(DT.MR_FA)
