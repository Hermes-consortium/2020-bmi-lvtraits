# produce forest plot of estimates (using new LV traits - Aung N 2021)

pacman::p_load(tidyverse, patchwork, fs, glue)

outcome_order <- c("LVEF", "LVEDV", "LVESV", "LVSV", "LVM", "LVMVR", "HeartRate")

df <- read_tsv("results/2021-06-14_BMI_LVtrait.tsv") %>%
  separate(Outcome, c("Outcome", "Adjustment")) %>%
  mutate(Adjustment = replace_na(Adjustment, "unadjusted"),
         Outcome = factor(Outcome, levels = rev(outcome_order)))

df_plot <- df %>% 
  filter(Method %in% c("IVW", "MR-Egger", "Weighted median"))

ggplot(df_plot, aes(y = beta, ymin = beta_LCI, ymax = beta_UCI, x = Outcome)) +
  theme_bw() +
  facet_grid(cols = vars(Adjustment)) +
  geom_point(aes(color = Method), position = position_dodge(0.5)) +
  geom_errorbar(aes(color = Method), position = position_dodge(0.5), width = 0.3) +
  geom_hline(yintercept = 0, color = "gray21", linetype = "dashed") +
  labs(y = "SD change per additional SD BMI") +
  coord_flip()

dir_create("figures/")
ggsave(glue("figures/{Sys.Date()}_ForestPlot_BMI_LVtraits.png"))
