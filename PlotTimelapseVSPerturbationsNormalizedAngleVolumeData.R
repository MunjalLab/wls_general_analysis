# For plotting angles vs volumes measured from pooled timecourse data
# Written by Kira L. Heikes, with help from Jiacheng Wang and Copilot GPT-5, 7/19/2026
# For use in manuscript: Hydrostatic pressure shapes and canalizes semicircular canal morphology to ensure vestibular function.
# doi.org/10.64898/2026.07.22.740136

# Plots and correlation calculations for WT timelapse and perturbation before and after angle and volume data

# pull needed libraries
library(tidyverse)
library(dplyr)

setwd("E:/ScriptsForGithub")

# -----------------------------------------------
# # 1. Read data from csv
# -----------------------------------------------
# USER INPUT: set csv file name
df <- read_csv("./TimelapseDataForPlotting.csv")

# -----------------------------------------------
# 2. Plot WT Change in Angle vs Change in Volume with linear regression
# -----------------------------------------------
# build dataframe of only WT timelapse data from initial dataframe
df_subset <- df %>%
  filter(
    (Group == "WT")
  ) %>%
  # group by timelapse
  group_by(Sample) %>%
  # for each timelapse, compute nomralized Angle and normalized Volume values, based on initial ('first') value for each timelapse
  mutate(
    dAdApot = (Angle - first(Angle)) / (180-first(Angle)),
    dVd0 = (Volume - first(Volume)) / first(Volume)
) %>%
  ungroup()

# compute spearman correlation for the full set of timelapses, use exact = FALSE because there will be tied ranks in the pooled data
spearman_WT <- cor.test(df_subset$dVd0, df_subset$dAdApot, method = "spearman", exact = FALSE)

# compute spearman correlation for each individual timelapse, use exact = TRUE because there will be no tied ranks in the individual timelapse data
spearman_WT1 <-cor.test(~ dVd0 + dAdApot, data = df_subset, subset = Sample == 1, method = "spearman", exact = TRUE)
spearman_WT2 <-cor.test(~ dVd0 + dAdApot, data = df_subset, subset = Sample == 2, method = "spearman", exact = TRUE)
spearman_WT3 <-cor.test(~ dVd0 + dAdApot, data = df_subset, subset = Sample == 3, method = "spearman", exact = TRUE)
spearman_WT4 <-cor.test(~ dVd0 + dAdApot, data = df_subset, subset = Sample == 4, method = "spearman", exact = TRUE)
spearman_WT5 <-cor.test(~ dVd0 + dAdApot, data = df_subset, subset = Sample == 5, method = "spearman", exact = TRUE)

# plot pooled timelapse normalized data in scatterplot
p_WT <- ggplot(df_subset,
                      aes(x = dVd0, y = dAdApot, color = factor(Sample))) +
  geom_point(size = 4, alpha = 1) +
  #  Set colors by timelapse ID (sample)
  scale_color_manual(values = c(
    "1" = "#009E73",
    "2" = "#E69F00",
    "3" = "#56B4E9",
    "4" = "#cc79a7",
    "5" = "#D55E00"
  )) +
  labs(
    title = "WT dVd0 dAdApotential",
    y = "dA/dApotential",
    x = "dV/d0",
    color = "Embryo"
  ) +
  theme_classic(base_size = 14)+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30),
    legend.text = element_text(size = 30),
    axis.ticks.length = unit(0.4, "cm"))+
  ylim(NA, 1)

p_WT

ggsave("./WT_dVd0_dAdApot.pdf", p_WT, width = 9, height = 6, units="in")

# -----------------------------------------------
# 3. Plot WT Change in Angle vs Change in Volume with linear regression
# -----------------------------------------------
df_all <- df %>%
  # group by data type (WT, forskolin perturbation, ablation A, ablation P) and by sample (timelapse or each indiviudal perturbation replicate)
  group_by(Group, Sample) %>%
  # for each timelapse or each individual perturbed embryo, compute nomralized Angle and normalized Volume values, based on initial ('first') value for each timelapse or sample
  mutate(
    dAdApot = (Angle - first(Angle)) / (180-first(Angle)),
    dVd0 = (Volume - first(Volume)) / first(Volume)
) %>%
  ungroup()

# compute spearman correlation for the full set of pooled WT and perturbation data, use exact = FALSE because there will be tied ranks in the pooled data
spearman_All <- cor.test(df_all$dVd0, df_all$dAdApot, method = "spearman", exact = FALSE)

# plot pooled timelapse normalized data in scatterplot
p_All <- ggplot(df_all,
                      aes(x = dVd0, y = dAdApot, color = factor(Group))) +
  geom_point(size = 4, alpha = 1) +
  #  Set colors by Group (WT timelapse or perturbation type)
  scale_color_manual(values = c(
    "WT" = "#A9A9A9",
    "ForskolinP" = "#E69F00",
    "AblationA" = "#56B4E9",
    "AblationP" = "#cc79a7",
    "ForskolinA" = "#D55E00"
  )) +
  labs(
    title = "All dVd0 dAdApotential",
    y = "dA/dApotential",
    x = "dV/d0",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14)+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30),
    legend.text = element_text(size = 30),
    axis.ticks.length = unit(0.4, "cm"))


p_All

ggsave("./ALL_dVd0_dAdApot.pdf", p_All, width = 9, height = 6, units="in")