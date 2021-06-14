# scripts to run MR analysis of BMI on LV traits (data from Nay Aung 2021)

pacman::p_load(MendelianRandomization, tidyverse, data.table, fs, glue)

list_DT_gwas <- dir_ls("data/GWAS/LVtraits-subset_AungN2021", glob = "*.tsv.gz") %>%
  set_names(~basename(.) %>%
    str_extract(".*(?=-subset)")) %>%
  map(fread, key = "varID")

DT_gwas_BMI <- fread("data/GWAS/BMI_GIANT+UKB2018.CLEAN.tsv.gz", key = "varID")
DT_ins_BMI <- fread("data/instruments/BMI_GIANT_UKB_2018_instrument_rsq0.01.txt",
  key = "SNP")

merge_DT <- function(DT_gwas_x, DT_gwas_y) {
  merge(DT_gwas_x, DT_gwas_y[, .(varID, beta, se, P)], by = "varID")
}

list_DT_mr <- map(list_DT_gwas, ~merge_DT(DT_gwas_BMI[DT_ins_BMI[, .(SNP)], nomatch = 0L],
  .x))

run_MR <- function(DT.MR, ...) {
  MR_res <- mr_allmethods(mr_input(bx = DT.MR$beta.x, bxse = DT.MR$se.x, by = DT.MR$beta.y,
    byse = DT.MR$se.y, ...))
  MR_res
}

list_res_MR <- map(list_DT_mr, run_MR)

# concatenate MR results
format_res_MR <- function(res_MR) {
  MR_res@Values %>%
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
write_tsv(df_res, glue("results/{Sys.Date()}_BMI_LVtrait.tsv"))
