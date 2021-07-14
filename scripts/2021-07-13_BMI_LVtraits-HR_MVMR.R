# R script to run MVMR analysis with BMI and SBP exposure
# using LVtraits outcome data from Nay Aung 2021 (personal communication)
# as per 2021-07-13 we only have access to BMI-SBP data

if (!require(pacman)) install.packages("pacman")
pacman::p_load(MendelianRandomization, tidyverse, data.table, fs, glue)

expos <- c("BMI-SBP") %>% set_names(.)

# GWAS data files
expo_files <- c(BMI = "data/GWAS/BMI_GIANT+UKB2018.CLEAN.tsv.gz",
                SBP = "data/GWAS/BloodPressure/SBP_UKB-no-BP-drug.CLEAN.tsv.gz")

outcome_files <- dir_ls("data/GWAS/LVtraits-subset_AungN2021/", glob="*.CLEAN.tsv.gz") %>%
  set_names(~basename(.x) %>% str_remove_all("-subset.*")) %>%
  c(HR = "data/GWAS/HeartRate_Eppinga2016.CLEAN.tsv.gz")

# write instrument file for tabix query
ins_files <- c(`BMI-SBP` = "data/instruments/BMI-SBP_instrument_rsq0.01.tsv")
query_files <- map_chr(ins_files, path_ext_set, ".pos.tsv")

walk2(ins_files, query_files, ~fread(.x) %>%
  .[order(CHR, BP), 1:2] %>%
  fwrite(.y, sep = "\t", col.names=FALSE))

read_tabix <- function(file_input, file_instrument){
  fread(cmd = glue("tabix -hR {file_instrument} {file_input}")) %>%
    unique() %>%
    rename(uniqID = `00,000000000000,0:0`)
}

trait_files <- c(expo_files, outcome_files)
list.DT_trait <- map2(trait_files, query_files, read_tabix)

walk(traits, function(x)list.DT_trait[[x]][,Trait := x])

read_format_ins <- function(DT_ins){
  fread(DT_ins) %>%
    mutate(A1 = pmin(EA, OA),
           A2 = pmax(EA, OA),
           fCHR = formatC(CHR, width = 2, format = "d", flag = "0"),
           fBP = formatC(BP, width = 12, format = "d", flag = "0"),
           uniqID = glue("{fCHR},{fBP},{A1}:{A2}"))
}

list.DT_ins <- map(ins_files, read_format_ins)

expo_order <- c("BMI", "SBP")

# List all traits estimates for BMI SNP
combine_df <- function(expos, outcome, filter_NA=F){
  expos <- factor(expos, levels = expo_order)
  ins <- expos %>% sort %>% paste(collapse = "-")

  DT.traits <- list.DT_trait[c(as.character(expos), outcome)] %>%
    bind_rows(.id = "Trait")

  df.MR <- DT.traits[uniqID %in% list.DT_ins[[ins]][,uniqID]] %>%
    pivot_wider(id_cols = uniqID, names_from = Trait,
                values_from = c("beta", "se"))
  if (filter_NA){
    df.MR <- df.MR %>%
      filter_all(all_vars(!is.na(.)))
  }
  df.MR
}

# Run Multi-variable MR
run_MRMV <- function(expos, outcome, filter_NA=FALSE){
  df.MR <- combine_df(expos, outcome = outcome, filter_NA = filter_NA)
  input.MRMV <- mr_mvinput(
    bx = as.matrix(select(df.MR, one_of(paste0("beta_", expos)))),
    bxse = as.matrix(select(df.MR, one_of(paste0("se_", expos)))),
    by = df.MR[[paste0("beta_", outcome)]],
    byse = df.MR[[paste0("se_", outcome)]],
    snps = df.MR$uniqID,
    exposure = expos,
    outcome = outcome
  )

  mr_mvivw(input.MRMV)
}

list.expos <- list(
  c("BMI", "SBP")
)

list.outcomes <- names(outcome_files)


format_res <- function(output.MRMV, expo = "BMI"){
  expo_pos <- which(output.MRMV@Exposure == expo)

  data.table(
    Exposure = expo,
    Outcome = output.MRMV@Outcome,
    Covars = paste(output.MRMV@Exposure[-expo_pos], collapse = ","),
    beta = output.MRMV@Estimate[expo_pos],
    SE = output.MRMV@StdError[expo_pos],
    beta_LCI = output.MRMV@CILower[expo_pos],
    beta_UCI = output.MRMV@CIUpper[expo_pos],
    P = output.MRMV@Pvalue[expo_pos],
    N_snps = output.MRMV@SNPs
  )
}

list.res_DT_MVMR <- map2(list.expos, list.outcomes, run_MRMV)

DT.res <- map(list.res_DT_MVMR, format_res) %>% bind_rows

# write results
write_tsv(DT.res, glue("results/{Sys.Date()}_MRMV_BMI_LVtraits.tsv"))
