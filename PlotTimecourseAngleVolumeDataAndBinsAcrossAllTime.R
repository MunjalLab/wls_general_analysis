# For plotting angles vs volumes measured from pooled timecourse data
# Written by Kira L. Heikes, with help from Jiacheng Wang and Copilot GPT-5, 6/09/2026
# For use in manuscript: Hydrostatic pressure shapes and canalizes semicircular canal morphology to ensure vestibular function.
# doi.org/10.64898/2026.07.22.740136

# Plot for dataset '2025-11-4_mNG_TimeCourse'

# pull needed libraries
library(tidyverse)
library(dplyr)

setwd("E:/ScriptsForGithub")

# -----------------------------------------------
# 1. Read data and extract timepoint IDs from entries
# -----------------------------------------------
# USER INPUT: set csv file name
df <- read_csv("./WT_AllAspectRatiosandAngles.csv")

# USER INPUT: set binning size
# set size for binning data by OV volume (e.g. every 90 picoliters (pL))
binsize = 90 # optimized for good spread of data representation across bins and no bins with only 2 values that would be a consistent bin size across all pillars


# Extract full string "xxhpf" as timepoint (no numeric conversion)
df <- df %>%
  mutate(Timepoint = str_extract(Embryo, "\\d+hpf"))

# -----------------------------------------------
# 2A. Plot A-Pillar Angle vs OV Volume for 54hpf, 60hpf, 66hpf, & 72hpf timepoints combined
# -----------------------------------------------
df_A <- df %>% filter(Timepoint %in% c("54hpf", "60hpf", "66hpf", "72hpf"))

# plot all Anterior data
pA <- ggplot(df_A, aes(x = `Volume (pl)`, y = `A-Angle`, color = Timepoint)) +
  geom_point(size = 4, alpha = 1) +
  labs(
    title = "A-Angle vs Volume",
    color = "Timepoint"
  ) +
  theme_classic(base_size = 14)+
  ylim(90,190)+
  xlim(330,1350)+
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9", "#cc79a7"))+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30),
    legend.text = element_text(size = 30),
    axis.ticks.length = unit(0.4, "cm"))

pA

ggsave("./A_plot_alltime.pdf", pA, width = 11.328, height = 6, units="in")

# -----------------------------------------------
# 2B. Bin Anterior data, calculate CV, and plot in ggplot2
# -----------------------------------------------
# make binned dataframe for computing CV and plotting
binned_dataAnew <- df %>%
  filter(Timepoint %in% c("54hpf", "60hpf", "66hpf", "72hpf")) %>%
  filter(!is.na(`A-Angle`)) # only keep entries in dataframe that have a value for Anterior Angle
binned_dataAnew <- binned_dataAnew %>%
  mutate(
    # Define custom bins using cut
    bin = cut(`Volume (pl)`, 
              breaks = seq(from=315,to=1350,by=binsize),
              labels = round(head(seq(from=315,to=1350,by=binsize), -1)+(binsize/2),digits=2),
              include.lowest = TRUE, 
              dig.lab = 4)
  ) %>%
  filter(!is.na(bin)) %>%  # or drop_na(bin)
  group_by(bin) %>%
  summarise(
    mean_val = mean(`A-Angle`, na.rm = TRUE),
    sd_val = sd(`A-Angle`, na.rm = TRUE),
    n = n(),
    # Calculate Coefficient of Variation (CV)
    cv = sd_val / mean_val, 
    .groups = "drop"
  ) %>%
  # Remove bins with no data or 0 standard deviation (to avoid Inf/NA in CV)
  filter(n > 0 & sd_val > 0) %>%
  mutate(label_text = paste0("n = ", n))

# Plot in ggplot2
pAbincv <- ggplot(binned_dataAnew, aes(x = bin, y = cv)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = label_text), vjust = -1.5, size = 4) +
  annotate(geom = "text", x = Inf, y = Inf, label = sum(binned_dataAnew$n, na.rm = TRUE), vjust=3, hjust=3, color = "blue", size = 5, fontface = "bold")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = paste("A-Angle CV by", binsize, "pL bins"),
    x = "Bin by Volume (pL)",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_classic(base_size = 14)+
  ylim(0,0.15)+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 20),
    axis.ticks.length = unit(0.4, "cm"))

pAbincv
ggsave("./A_cvbinplot_alltime.pdf", pAbincv, width = 9, height = 6, units="in")

# -----------------------------------------------
# 3A. Plot P-Pillar Angle vs OV Volume for 54hpf, 60hpf, 66hpf, & 72hpf timepoints combined
# -----------------------------------------------
df_P <- df %>% filter(Timepoint %in% c("54hpf", "60hpf", "66hpf", "72hpf"))

# plot all Posterior data
pP <- ggplot(df_P, aes(x = `Volume (pl)`, y = `P-Angle`, color = Timepoint)) +
  geom_point(size = 4, alpha = 1) +
  labs(
    title = "P-Angle vs Volume",
    color = "Timepoint"
  ) +
  theme_classic(base_size = 14)+
  ylim(90,190)+
  xlim(330,1350)+
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9", "#cc79a7"))+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30),
    legend.text = element_text(size = 30),
    axis.ticks.length = unit(0.4, "cm"))

pP

ggsave("./P_plot_alltime.pdf", pP, width = 11.328, height = 6, units="in")

# -----------------------------------------------
# 3B. Bin Posterior data, calculate CV, and plot in ggplot2
# -----------------------------------------------
# make binned dataframe for computing CV and plotting
binned_dataPnew <- df %>%
  filter(Timepoint %in% c("54hpf", "60hpf", "66hpf", "72hpf")) %>%
  filter(!is.na(`P-Angle`)) # only keep entries in dataframe that have a value for Posterior Angle
binned_dataPnew <- binned_dataPnew %>%
  mutate(
    # Define custom bins using cut
    bin = cut(`Volume (pl)`, 
              breaks = seq(from=315,to=1350,by=binsize),
              labels = round(head(seq(from=315,to=1350,by=binsize), -1)+(binsize/2),digits=2),
              include.lowest = TRUE, 
              dig.lab = 4)
  ) %>%
  filter(!is.na(bin)) %>%  # or drop_na(bin)
  group_by(bin) %>%
  summarise(
    mean_val = mean(`P-Angle`, na.rm = TRUE),
    sd_val = sd(`P-Angle`, na.rm = TRUE),
    n = n(),
    # Calculate Coefficient of Variation (CV)
    cv = sd_val / mean_val, 
    .groups = "drop"
  ) %>%
  # Remove bins with no data or 0 standard deviation (to avoid Inf/NA in CV)
  filter(n > 0 & sd_val > 0) %>%
  mutate(label_text = paste0("n = ", n))

# Plot in ggplot2
pPbincv <- ggplot(binned_dataPnew, aes(x = bin, y = cv)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = label_text), vjust = -1.5, size = 4) +
  annotate(geom = "text", x = Inf, y = Inf, label = sum(binned_dataPnew$n, na.rm = TRUE), vjust=3, hjust=3, color = "blue", size = 5, fontface = "bold")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = paste("P-Angle CV by", binsize, "pL bins"),
    x = "Bin by Volume (pL)",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_classic(base_size = 14)+
  ylim(0,0.15)+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 20),
    axis.ticks.length = unit(0.4, "cm"))

pPbincv
ggsave("./P_cvbinplot_alltime.pdf", pPbincv, width = 9, height = 6, units="in")

# -----------------------------------------------
# 4A. Plot V-Pillar Angle vs OV Volume for 54hpf, 60hpf, 66hpf, & 72hpf timepoints combined
# -----------------------------------------------
df_V <- df %>% filter(Timepoint %in% c("54hpf", "60hpf", "66hpf", "72hpf"))

# plot all Ventral data
pV <- ggplot(df_V, aes(x = `Volume (pl)`, y = `V-Angle`, color = Timepoint)) +
  geom_point(size = 4, alpha = 1) +
  labs(
    title = "V-Angle vs Volume",
    color = "Timepoint"
  ) +
  theme_classic(base_size = 14)+
  ylim(90,190)+
  xlim(330,1350)+
  scale_color_manual(values = c("#009E73", "#E69F00", "#56B4E9", "#cc79a7"))+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30),
    legend.text = element_text(size = 30),
    axis.ticks.length = unit(0.4, "cm"))

pV

ggsave("./V_plot_alltime.pdf", pV, width = 11.328, height = 6, units="in")

# -----------------------------------------------
# 4B. Bin Ventral data, calculate CV, and plot in ggplot2
# -----------------------------------------------
# make binned dataframe for computing CV and plotting
binned_dataVnew <- df %>%
  filter(Timepoint %in% c("54hpf", "60hpf", "66hpf", "72hpf")) %>%
  filter(!is.na(`V-Angle`)) # only keep entries in dataframe that have a value for Ventral Angle
binned_dataVnew <- binned_dataVnew %>%
  mutate(
    # Define custom bins using cut
    bin = cut(`Volume (pl)`, 
              breaks = seq(from=315,to=1350,by=binsize),
              labels = round(head(seq(from=315,to=1350,by=binsize), -1)+(binsize/2),digits=2),
              include.lowest = TRUE, 
              dig.lab = 4)
  ) %>%
  filter(!is.na(bin)) %>%  # or drop_na(bin)
  group_by(bin) %>%
  summarise(
    mean_val = mean(`V-Angle`, na.rm = TRUE),
    sd_val = sd(`V-Angle`, na.rm = TRUE),
    n = n(),
    # Calculate Coefficient of Variation (CV)
    cv = sd_val / mean_val, 
    .groups = "drop"
  ) %>%
  # Remove bins with no data or 0 standard deviation (to avoid Inf/NA in CV)
  filter(n > 0 & sd_val > 0) %>%
  mutate(label_text = paste0("n = ", n))

# Plot in ggplot2
pVbincv <- ggplot(binned_dataVnew, aes(x = bin, y = cv)) +
  geom_point(size = 4, alpha = 1) +
  geom_text(aes(label = label_text), vjust = -1.5, size = 4) +
  annotate(geom = "text", x = Inf, y = Inf, label = sum(binned_dataVnew$n, na.rm = TRUE), vjust=3, hjust=3, color = "blue", size = 5, fontface = "bold")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = paste("V-Angle CV by", binsize, "pL bins"),
    x = "Bin by Volume (pL)",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_classic(base_size = 14)+
  ylim(0,0.15)+
  theme(
    title = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 20),
    axis.ticks.length = unit(0.4, "cm"))

pVbincv
ggsave("./V_cvbinplot_alltime.pdf", pVbincv, width = 9, height = 6, units="in")