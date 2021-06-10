library(devtools)
require(data.table)
library(TwoSampleMR)

# setwd("~/BMI_LV/BMI")
WHRadjBMI <- read_exposure_data(
  "data/BMI_WHR_GIANT/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt",
  snp_col = "SNP", beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  phenotype_col = "WHRadjBMI",
  samplesize_col = "N",
  eaf_col = "Freq_Tested_Allele"
)

# setwd("~/BMI_LV")
bmi <- fread("data/WHRadjBMI_GIANT_UKB_2018_instrument", key = "SNP")
exposure_dat <- merge(WHRadjBMI, bmi, by = "SNP")

# setwd("~/BMI_LV/LVtraits_Pirruccello2020")
# Get effects of instruments on outcome
outcome_LVEDV <- read_outcome_data("data/LVtraits_Pirruccello2020/MRI_lvedv_filtered.tsv",
  sep = "\t",
  snp_col = "SNP", beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_BOLT_LMM" # previously P_BOLT_LMM_INF
)


# Harmonise the exposure and outcome data
hardata_LVEDV <- harmonise_data(exposure_dat, outcome_LVEDV)
write.table(hardata_LVEDV, file = "hardata_lvedv.txt", row.names = FALSE, col.names = TRUE)

# Perform MR
mr_wwhrbmi_lvedv <- mr(hardata_LVEDV)
generate_odds_ratios(mr_wwhrbmi_lvedv)

# Results
# id.exposure id.outcome outcome exposure                    method   nsnp     b         se         pval       lo_ci
# 1      ecwNrQ     eeDtek outcome exposure                  MR Egger  577 -0.2242149 0.05511819 5.406501e-05 -0.3322465
# 2      ecwNrQ     eeDtek outcome exposure           Weighted median  577 -0.2095517 0.03106072 1.514368e-11 -0.2704307
# 3      ecwNrQ     eeDtek outcome exposure Inverse variance weighted  577 -0.2272837 0.02319925 1.159485e-22 -0.2727542
# 4      ecwNrQ     eeDtek outcome exposure               Simple mode  577 -0.2421231 0.10524201 2.176882e-02 -0.4483975
# 5      ecwNrQ     eeDtek outcome exposure             Weighted mode  577 -0.2044936 0.06194424 1.022097e-03 -0.3259044
# up_ci        or  or_lci95  or_uci95
# 1 -0.11618322 0.7991434 0.7173105 0.8903121
# 2 -0.14867267 0.8109477 0.7630508 0.8618512
# 3 -0.18181314 0.7966948 0.7612799 0.8337571
# 4 -0.03584878 0.7849595 0.6386508 0.9647862
# 5 -0.08308293 0.8150599 0.7218742 0.9202748

#######################################################################################

clump_data(bmi)

p1 <- mr_scatter_plot(mr_wwhrbmi_lvedv, hardata_LVEDV)

p1[[1]]

# Heterogeneity statistics
mr_heterogeneity(hardata_LVEDV)

# Horizontal pleiotropy
mr_pleiotropy_test(hardata_LVEDV)

# Single SNP analysis
res_single <- mr_singlesnp(hardata_LVEDV)

# Leave-one-out analysis
res_loo <- mr_leaveoneout(hardata_LVEDV)

# Forest plot
res_single <- mr_singlesnp(hardata_LVEDV)
p2 <- mr_forest_plot(res_single)
p2[[1]]
