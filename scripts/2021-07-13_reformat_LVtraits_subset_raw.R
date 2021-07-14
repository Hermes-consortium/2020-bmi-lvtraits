# reformat data to standardized format

pacman::p_load(data.table, readxl, tidyverse, glue)

files <- list(
  SBP = "data/2021-07-01_lookup/SBP_variants_lookup_in_CMR_INT_no_height_adj.xlsx",
  `SBP_heightadj` = "data/2021-07-01_lookup/SBP_variants_lookup_in_CMR_INT_with_height_adj.xlsx",
  BMI = "data/2021-07-01_lookup/BMI_VEGFA_variants_lookup_in_CMR_INT_no_height_adj.xlsx",
  `BMI_heightadj` = "data/2021-07-01_lookup/BMI_VEGFA_variants_lookup_in_CMR_INT_with_height_adj.xlsx"
)

read_format <- function(file){
  col <- excel_sheets(file) %>% set_names(~str_remove_all(., "_INT$"))
  map(col, ~read_excel(file, sheet = .x)) %>% bind_rows(.id = "LVtrait")
}

df_SBP <- map(files[1:2], read_format) %>% bind_rows(.id = "SNP") %>%
  separate(SNP, c("locus", "adjustment"))

df_BMI <- map(files[3:4], read_format) %>% bind_rows(.id = "SNP") %>%
  separate(SNP, c("locus", "adjustment")) %>% 
  mutate(locus = Source) %>% 
  select(locus, adjustment, CHR, BP = pos_b37, everything(), -Source)

df <- bind_rows(df_BMI, df_SBP) %>% 
  group_by(LVtrait, adjustment) %>% 
  nest()

write_format <- function(df_subset, LVtrait, adjustment){
  outfile <- ifelse(!is.na(adjustment),
                    glue("data/2021-07-01_lookup/{LVtrait}-{adjustment}-subset_AungN2021.tsv.gz"),
                    glue("data/2021-07-01_lookup/{LVtrait}-subset_AungN2021.tsv.gz"))
  write_tsv(df_subset, outfile)
}

with(df, pwalk(list(data, LVtrait, adjustment), write_format))
