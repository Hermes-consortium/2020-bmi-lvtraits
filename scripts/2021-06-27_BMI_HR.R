# scripts to run MR analysis of BMI on LV traits (data from Nay Aung 2021)

pacman::p_load(MendelianRandomization, tidyverse, data.table, fs, glue)

DT_assoc_BMI_ins <- fread("data/BMI_GIANT+UKB2018_instrument_rsq0.01_assoc.tsv", key="varID")
DT_assoc_HR_ins <- fread("data/HeartRate_Eppinga2016_BMI_instrument_rsq0.01_assoc.tsv", key="varID")

file_gwas_BMI <- "data/GWAS/BMI_GIANT+UKB2018.CLEAN.tsv.gz"
file_ins <- "data/instruments/BMI_GIANT_UKB_2018_instrument_rsq0.01.txt"

cmd <- glue("tabix -R <(cut -f1,3 {file_ins}) {file_gwas_BMI}")

fread(cmd = glue("tabix -R <(awk -v OFS='\\t' 'NR>1 {{print $1, $3}}' {file_ins}) {file_gwas_BMI}"))

DT_ins_BMI <- fread("data/instruments/BMI_GIANT_UKB_2018_instrument_rsq0.01.txt",
  key = "SNP")

DT_gwas_ins_BMI <- DT_gwas_BMI[DT_ins_BMI[, .(SNP)], nomatch = 0L]

DT_gwas_HR <- fread("data/GWAS/HeartRate_Eppinga2016.CLEAN.tsv.gz", key = "varID")


merge_DT <- function(DT_gwas_x, DT_gwas_y) {
  merge(DT_gwas_x, DT_gwas_y[, .(varID, beta, se, P)], by = "varID")
}

list_DT_mr <- map(list_DT_gwas, ~merge_DT(DT_gwas_ins_BMI, .x))

run_MR <- function(DT.MR, ...) {
  MR_res <- mr_allmethods(mr_input(bx = DT.MR$beta.x, bxse = DT.MR$se.x, by = DT.MR$beta.y,
    byse = DT.MR$se.y, ...))
  MR_res
}

list_res_MR <- map(list_DT_mr, run_MR)

# concatenate MR results
format_res_MR <- function(res_MR) {
  res_MR@Values %>%
    as_tibble() %>%
    mutate(Method = ifelse(Method == "(intercept)",
                           paste(lag(Method), Method),
                           Method)) %>%
    select(Method, beta = Estimate, SE = `Std Error`, beta_LCI = `95% CI `, beta_UCI = ` `, P = `P-value`)
}

df_res <- map_df(list_res_MR, format_res_MR, .id = "Outcome") %>%
  mutate(Exposure = "BMI") %>%
  relocate(Exposure)

# write results
dir_create("results")
write_tsv(df_res, glue("results/{Sys.Date()}_BMI_HR.tsv"))
