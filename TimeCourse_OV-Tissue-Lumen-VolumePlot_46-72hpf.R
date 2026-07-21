library(tidyverse)
library(readxl)
library(extrafont)

fonts()
font_import(pattern = "arial", prompt = FALSE)
#loadfonts(device = 'win')
fonts()

# Load data
setwd("/Users/Kira/OneDrive - Duke University/Jiacheng Wang's files - Wls paper/Fig1")

summary_df <- read_csv('summary_file_combined.csv')
plot_df <- summary_df %>% select(-tissue_fraction_mean, -tissue_fraction_sd)
plot_df <- plot_df %>%
  pivot_longer(
    cols = -Age,
    names_to = c("type", ".value"),
    names_pattern = "(.*)_(mean|sd)"
  )
#%>%
# group_by(plot_df$experiment)

# ---- Step 4: Plot ----
p <- ggplot(plot_df, aes(x = Age, y = mean, color = type, fill = type)) +
#p <- ggplot(plot_df, aes(x = cut(Age,breaks=c(0,59,120)), y = mean, color = type, fill = type)) +

  # Shaded SD area
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  
  # Mean line + points
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +

  scale_x_continuous(
    breaks = unique(plot_df$Age),
    labels = round(unique(plot_df$Age))
    ) +
  
  # Colors to match your style
  scale_color_manual(values = c(
    "lumen" = "#1f77b4",
    "ov" = "#ff7f0e",
    "tissue" = "#2ca02c"
  )) +
  scale_fill_manual(values = c(
    "lumen" = "#1f77b4",
    "ov" = "#ff7f0e",
    "tissue" = "#2ca02c"
  )) +
  
  labs(
    x = "Time (hpf)",
    y = "Volume (pL)",
    color = NULL,
    fill = NULL
  ) +
  ylim(0, NA) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 30.48),
    axis.text  = element_text(size = 19.05),
    legend.text = element_text(size = 18),
    # Set tick length to 0.5 centimeters
    axis.ticks.length = unit(0.4, "cm")
  )
p
ggsave(
  "46-72hpf_volume_plot_new.pdf",
  plot = p,
  width = 12,        # inches
  height = 5,
  units = "in",
)


# ---- Step 4: Plot ----
p <- ggplot(summary_df, aes(x = Age, y = tissue_fraction_mean)) +
  
  # Shaded SD area
  geom_ribbon(aes(ymin = tissue_fraction_mean - tissue_fraction_sd, ymax = tissue_fraction_mean + tissue_fraction_sd), alpha = 0.2, color = NA) +
  
  # Mean line + points
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  
  scale_x_continuous(
    breaks = unique(summary_df$Age),
    labels = round(unique(summary_df$Age))
  ) +
  
  labs(
    x = "Time (hpf)",
    y = "Tissue fraction",
    color = NULL,
    fill = NULL
  ) +
  #ylim(0, NA) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 30.48),
    axis.text  = element_text(size = 19.05),
    legend.text = element_text(size = 18),
    # Set tick length to 0.5 centimeters
    axis.ticks.length = unit(0.4, "cm")
  )


p
ggsave(
  "46-72hpf_tissue_fraction_new.pdf",
  plot = p,
  width = 12,        # inches
  height = 5,
  units = "in",
)
