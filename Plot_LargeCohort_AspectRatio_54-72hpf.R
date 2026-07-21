# Plot for dataset '2025-11-4_mNG_TimeCourse'
library(tidyverse)
setwd("/Users/jiacheng/Documents/plot-Inflation-deflation-WT")
# 1. Read data
df <- read_csv("./LumenVolume-combined.csv")

# 2. Extract full string "xxhpf" as timepoint (no numeric conversion)
df <- df %>%
  mutate(Timepoint = str_extract(Embryo, "\\d+hpf"))

# -----------------------------------------------
# 3. Plot 54 hpf A + 60 hpf P + 72 hpf V
# -----------------------------------------------
df_subset <- df %>%
  filter(
    (Timepoint == "54hpf" & !is.na(`A-aspect ratio`)) |
      (Timepoint == "60hpf" & !is.na(`P-aspect ratio`)) |
      (Timepoint == "72hpf" & !is.na(`V-aspect ratio`))
  )

df_long <- df_subset %>%
  mutate(
    AspectType = case_when(
      Timepoint == "54hpf" ~ "Anterior",
      Timepoint == "60hpf" ~ "Posterior",
      Timepoint == "72hpf" ~ "Ventral"
    ),
    AspectValue = case_when(
      Timepoint == "54hpf" ~ `A-aspect ratio`,
      Timepoint == "60hpf" ~ `P-aspect ratio`,
      Timepoint == "72hpf" ~ `V-aspect ratio`
    )
  )


model_combined <- lm(AspectValue ~ `Volume (pl)`, data = df_long)

R2_combined <- summary(model_combined)$r.squared
R2_combined

p_combined <- ggplot(df_long,
                      aes(x = `Volume (pl)`, y = AspectValue, color = AspectType)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.3) +
  labs(
    title = "Combined regression: Anterior (54hpf) + Posterior (60hpf) + Ventral (72hpf)",
    y = "Aspect Ratio",
    x = "Volume (pl)",
    color = "Pillar Type"
  ) +
  theme_classic(base_size = 14)

p_combined

ggsave("./combined_plot.pdf", p_combined, width = 7, height = 5)


# -----------------------------------------------
# 4. Plot deflated and inflated conditions (Pillar Aspect ratio) over plot 3 (WT)
# -----------------------------------------------
# Deflated conditions: ablation experiments (before and after ablation)
# Inflated conditions: forskolin treatments (before and after treatment)
library(readxl)
# Reuse the dataframe preprocessing section from plot 6.
df_subset <- df %>%
  filter(
    (Timepoint == "54hpf" & !is.na(`A-aspect ratio`)) |
      (Timepoint == "60hpf" & !is.na(`P-aspect ratio`)) |
      (Timepoint == "72hpf" & !is.na(`V-aspect ratio`))
  )

df_long <- df_subset %>%
  mutate(
    AspectType = case_when(
      Timepoint == "54hpf" ~ "Anterior",
      Timepoint == "60hpf" ~ "Posterior",
      Timepoint == "72hpf" ~ "Ventral"
    ),
    AspectValue = case_when(
      Timepoint == "54hpf" ~ `A-aspect ratio`,
      Timepoint == "60hpf" ~ `P-aspect ratio`,
      Timepoint == "72hpf" ~ `V-aspect ratio`
    )
  )

# Create a new dataframe
df_long_new <- data.frame(
  LumenVolume = df_long$`Volume (pl)`,
  PillarAR = df_long$AspectValue,
  Age = df_long$Timepoint
)

df_long_new$Condition <- 'Control'

df_treated <- read_excel('Inflation-Deflation-summary-forPlot.xlsx')
df_treated$`Age (hpf)` <- paste0(df_treated$`Age (hpf)`, "hpf")

df_combined <- data.frame(
  Age = c(df_long_new$Age, df_treated$`Age (hpf)`),
  LumenVolume = c(df_long_new$LumenVolume, df_treated$`LumenVolume (pl)`),
  PillarAR = c(df_long_new$PillarAR, df_treated$PillarAR),
  Condition = c(df_long_new$Condition, df_treated$Condition)
)

model_control <- lm(PillarAR ~ LumenVolume,
                    data = df_combined %>% filter(Condition == "Control"))

x_seq <- seq(
  min(df_combined$LumenVolume, na.rm = TRUE),
  max(df_combined$LumenVolume, na.rm = TRUE),
  length.out = 100
)

df_line <- data.frame(LumenVolume = x_seq)

df_line$predicted <- predict(model_control, newdata = df_line)


p <- ggplot(df_combined,
            aes(x = LumenVolume, y = PillarAR, color = Condition)) +
  
  # All points (color-coded)
  geom_point(size = 3, alpha = 0.8) +
  
  # Manual regression line (FULL range)
  geom_line(data = df_line,
            aes(x = LumenVolume, y = predicted),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 1.3) +
  
  # Force control to be gray
  scale_color_manual(values = c(
    "Control" = "gray",
    "Before ablation" = "#1b9e77",
    "Post ablation" = "#d95f02",
    "Before forskolin treatment" = "#7570b3",
    "Post forskolin treatment" = "#e7298a"
  )) +
  
  labs(
    x = "Lumen Volume (pl)",
    y = "Pillar Aspect Ratio",
    color = "Condition",
    title = "Aspect Ratio vs Volume with Control Fit"
  ) +
  
  theme_classic(base_size = 14)

p

ggsave("./combined_inflation_deflation.pdf", p, width = 8, height = 5)


