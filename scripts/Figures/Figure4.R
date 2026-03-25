rm(list = ls())
source("R/utilities.R")
source("R/ThoraxDataParameters.R")


results <- compute_metrics(
  "data/Thorax_X_101.xlsx",
  "data/Thorax_Y_101.xlsx",
  "data/Thorax_Z_101.xlsx"
)

mean_X_male <- apply(results$row_means$X_male, 1, mean)
mean_X_female <- apply(results$row_means$X_female, 1, mean)
mean_Y_male <- apply(results$row_means$Y_male, 1, mean)
mean_Y_female <- apply(results$row_means$Y_female, 1, mean)
mean_Z_male <- apply(results$row_means$Z_male, 1, mean)
mean_Z_female <- apply(results$row_means$Z_female, 1, mean)
fwhm_X <- results$fwhm$X
fwhm_Y <- results$fwhm$Y
fwhm_Z <- results$fwhm$Z
sd_X <- results$sd$X
sd_Y <- results$sd$Y
sd_Z <- results$sd$Z

set.seed(12345)
number_functions <- 20
con_size <- length(sd_X)

noise_X_male <- Noise_generator(Sample_size = number_functions, Continuum_size = con_size,
                                Noise_mu = 0, Noise_sig = sd_X , Noise_fwhm = fwhm_X)$noise
noise_X_female <- Noise_generator(Sample_size = number_functions, Continuum_size = con_size,
                                Noise_mu = 0, Noise_sig = sd_X , Noise_fwhm = fwhm_X)$noise
noise_Y_male <- Noise_generator(Sample_size = number_functions, Continuum_size = con_size,
                                Noise_mu = 0, Noise_sig = sd_Y , Noise_fwhm = fwhm_Y)$noise
noise_Y_female <- Noise_generator(Sample_size = number_functions, Continuum_size = con_size,
                                  Noise_mu = 0, Noise_sig = sd_Y , Noise_fwhm = fwhm_Y)$noise
noise_Z_male <- Noise_generator(Sample_size = number_functions, Continuum_size = con_size,
                                Noise_mu = 0, Noise_sig = sd_Z , Noise_fwhm = fwhm_Z)$noise
noise_Z_female <- Noise_generator(Sample_size = number_functions, Continuum_size = con_size,
                                  Noise_mu = 0, Noise_sig = sd_Z , Noise_fwhm = fwhm_Z)$noise


X_male <- data_generator(signal = mean_X_male, noise = noise_X_male)
X_female <- data_generator(signal = mean_X_female, noise = noise_X_female)
Y_male <- data_generator(signal = mean_Y_male, noise = noise_Y_male)
Y_female <- data_generator(signal = mean_Y_female, noise = noise_Y_female)
Z_male <- data_generator(signal = mean_Z_male, noise = noise_Z_male)
Z_female <- data_generator(signal = mean_Z_female, noise = noise_Z_female)

#####Plotting####
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale) #Allows multiple color scales in the same plot, which 
#is essential for our layered approach (mean lines vs. tau thresholds). Instead 
#of hardcoding 10 tau colors
library(ggh4x) #Y-Axis labels: I used breaks = NULL for the D1prime row. Note: 
#This requires the ggh4x library (install via install.packages("ggh4x")).
#It is the most reliable way to hide labels on specific facets.

Domain <- 0:100

build_complete_df <- function(complete_data, group, axis_name) {
  
  complete_data %>%
    mutate(time = Domain) %>%              # add domain
    pivot_longer(
      cols = -time,
      names_to = "id",
      values_to = "value"
    ) %>%
    mutate(
      Group = group,
      Axis = axis_name,
      Panel = "GeneratedData"
    )
}


df_complete <- bind_rows(
  build_complete_df(as_tibble(X_male),   "Male",   "X"),
  build_complete_df(as_tibble(X_female), "Female", "X"),
  build_complete_df(as_tibble(Y_male),   "Male",   "Y"),
  build_complete_df(as_tibble(Y_female), "Female", "Y"),
  build_complete_df(as_tibble(Z_male),   "Male",   "Z"),
  build_complete_df(as_tibble(Z_female), "Female", "Z")
)





df_plot_X <- data.frame(
  time = Domain,
  Male = mean_X_male,
  Female = mean_X_female
)

df_plot_Y <- data.frame(
  time = Domain,
  Male = mean_Y_male,
  Female = mean_Y_female
)

df_plot_Z <- data.frame(
  time = Domain,
  Male = mean_Z_male,
  Female = mean_Z_female
)



# Combine mean data
df_means <- bind_rows(
  df_plot_X  %>% mutate(Axis = "X"),
  df_plot_Y %>% mutate(Axis = "Y"),
  df_plot_Z %>% mutate(Axis = "Z")
) %>%
  pivot_longer(cols = c(Male, Female),
               names_to = "Group",
               values_to = "Value") %>%
  mutate(Row = "Mean", Panel = "Mean")


build_noise_df <- function(noise_data, group, axis_name) {
  
  noise_data %>%
    mutate(id = row_number()) %>%
    pivot_longer(
      cols = -id,
      values_to = "value"
    ) %>%
    group_by(id) %>%
    mutate(time = Domain) %>%   # directly assign Domain vector
    ungroup() %>%
    mutate(
      Group = group,
      Axis = axis_name,
      Panel = "Noise"
    ) %>%
    dplyr::select(time, value, id, Group, Axis, Panel)
}

df_noise <- bind_rows(
  build_noise_df(as_tibble(t(noise_X_male)),   "Male",   "X"),
  build_noise_df(as_tibble(t(noise_X_female)), "Female", "X"),
  build_noise_df(as_tibble(t(noise_Y_male)),   "Male",   "Y"),
  build_noise_df(as_tibble(t(noise_Y_female)), "Female", "Y"),
  build_noise_df(as_tibble(t(noise_Z_male)),   "Male",   "Z"),
  build_noise_df(as_tibble(t(noise_Z_female)), "Female", "Z")
)


axis_labels <- c(X = "Sagittal",
                 Y = "Frontal",
                 Z = "Transverse")

panel_labels <- c(
  Noise = "Noise",
  Mean = "Mean",
  GeneratedData = "Generated Data"
)

df_noise$Axis <- factor(df_noise$Axis, levels = names(axis_labels), labels = axis_labels)
df_means$Axis  <- factor(df_means$Axis,  levels = names(axis_labels), labels = axis_labels)
df_complete$Axis <- factor(df_complete$Axis, levels = names(axis_labels), labels = axis_labels)

df_noise$Panel <- factor(df_noise$Panel, levels = names(panel_labels), labels = panel_labels)
df_means$Panel  <- factor(df_means$Panel,  levels = names(panel_labels), labels = panel_labels)
df_complete$Panel <- factor(df_complete$Panel, levels = names(panel_labels), labels = panel_labels)


####################Y_SCALE AND LIMITS####################
# Disable lims and y_scale if one wants to use free_y in facet_grid.
#Note: free_y will make the y-axis limits different across panels, 
#which can be useful for visualizing patterns within each panel but may hinder 
#direct comparisons across panels.
lims <- bind_rows(
  df_complete %>% dplyr::select(Axis, value),
  df_noise %>% dplyr::select(Axis, value),
) %>%
  group_by(Axis) %>%
  summarise(ymin = min(value), ymax = max(value))


y_scales <- list(
  scale_y_continuous(limits = lims[lims$Axis == "Sagittal",  c("ymin","ymax")] |> unlist()),
  scale_y_continuous(limits = lims[lims$Axis == "Frontal",   c("ymin","ymax")] |> unlist()),
  scale_y_continuous(limits = lims[lims$Axis == "Transverse",c("ymin","ymax")] |> unlist())
)

#######################
plot4 <- ggplot() +
  
  # --- Noise---
  geom_line(data = df_noise,
            aes(time, value,
                group = interaction(id, Group),
                color = Group),
            alpha = 0.3) +
  
  # --- Complete ---
  geom_line(
    data = df_complete,
    aes(time, value,
        group = interaction(id, Group),
        color = Group)
  ) +
  
  # --- MEAN ---
  
  geom_line(data = df_means,
            aes(time, Value, color = Group),
            linewidth = 1.5) +
  

  
  scale_color_manual(
    name = "Groups:",
    values = c("Male" = "cadetblue",
               "Female" = "tomato"),
    labels = c("Male" = "Group 2",
               "Female" = "Group 1"),
    guide = guide_legend(nrow = 1)
  ) +
  
  #facet_grid(Panel ~ Axis, scales = "free_y") + #try "fixed"
  facet_grid2(
    Panel ~ Axis,
    scales = "free_y",
    independent = "y"
  ) +
  #-----------
#For free_y, use facetted_pos_scales to set y-axis limits for each facet separately.
facetted_pos_scales(
  y = list(
    Axis == "Sagittal"   ~ scale_y_continuous(limits = lims[lims$Axis=="Sagittal",c("ymin","ymax")] |> unlist()),
    Axis == "Frontal"    ~ scale_y_continuous(limits = lims[lims$Axis=="Frontal",c("ymin","ymax")] |> unlist()),
    Axis == "Transverse" ~ scale_y_continuous(limits = lims[lims$Axis=="Transverse",c("ymin","ymax")] |> unlist())
  )
) +
  #-----------
theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    strip.text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Domain",
    y = "Value"
  )

print(plot4)

ggsave("outputs/Figure4.jpeg", plot4, width = 190, height = 180, units = "mm", dpi = 1200)


