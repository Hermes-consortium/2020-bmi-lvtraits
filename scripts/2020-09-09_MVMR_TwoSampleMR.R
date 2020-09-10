# Example script for MVMR

library(data.table)
library(TwoSampleMR)

# Helper function to read, prune, and format exposure data
read_expo_mv <- function(
  data_file,  # path to data file
  ins_files,   # path to instrument files (vector of instrument files used for MVMR
  trait_name, # trait name
  key_data = "SNP",   # key column for data file
  key_ins = "SNP",    # key column for instrument file
  ... # additional arguments to format_data() function
){
  data <- fread(data_file, key = key_data)
  ins <- rbindlist(lapply(ins_files, fread, key = key_ins))
  data <- unique(data[ins[,..key_ins], nomatch=0L])
  
  data_mr_expo <- format_data(data, ...)
  data_mr_expo$exposure <- trait_name
  return(data_mr_expo)
}

data_bmi <- read_expo_mv(
  data_file = "data/BMI_WHR_GIANT/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt",
  ins_files = c("data/BMI_GIANT_UKB_2018_instrument",
               "data/SBP_UKB_noDrug/SBP_UKB-no-BP-drug_instrument"),
  trait_name = "BMI",
  type = "exposure",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  samplesize_col = "N",
  eaf_col = "Freq_Tested_Allele"
)

data_sbp <- read_expo_mv(
  data_file = "data/SBP_UKB_noDrug/GWAS_SBP_UKB-no-BP-drug.txt",
  ins_files = c("data/BMI_GIANT_UKB_2018_instrument",
                "data/SBP_UKB_noDrug/SBP_UKB-no-BP-drug_instrument"),
  trait_name = "SBP",
  type = "exposure",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  pval_col = "P",
  samplesize_col = "N",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq"
)

# combine exposure data
data_expo_mv <- rbind(data_bmi, data_sbp)


# Read outcome data
data_lvedv <- read_outcome_data(
  filename = "data/LVtraits_Pirruccello2020/MRI_lvedv_filtered.tsv",
  sep = "\t",
  snps = data_bmi$SNP,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM"
)
data_lvedv$outcome <- "LVEDV"

# Run MVMR
data_mvmr <- mv_harmonise_data(data_expo_mv, data_lvedv)

res_mvmr <- mv_multiple(data_mvmr)

res_mvmr

