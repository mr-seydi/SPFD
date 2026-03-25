rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale) #Allows multiple color scales in the same plot, which 
#is essential for our layered approach (mean lines vs. tau thresholds). Instead 
#of hardcoding 10 tau colors
library(ggh4x) #Y-Axis labels: I used breaks = NULL for the D1prime row. Note: 
#This requires the ggh4x library (install via install.packages("ggh4x")).
#It is the most reliable way to hide labels on specific facets.

path_X <- "data/Thorax_X_101.xlsx"
path_Y <- "data/Thorax_Y_101.xlsx"
path_Z <- "data/Thorax_Z_101.xlsx"
SheetsNames <- excel_sheets(path_X)
X_list <- lapply(SheetsNames, function(sheet) read_excel(path_X, sheet = sheet))
Y_list <- lapply(SheetsNames, function(sheet) read_excel(path_Y, sheet = sheet))
Z_list <- lapply(SheetsNames, function(sheet) read_excel(path_Z, sheet = sheet))
names(X_list) <- SheetsNames
names(Y_list) <- SheetsNames
names(Z_list) <- SheetsNames

row_means_X_male <- X_list[["Male"]]
row_means_X_female <- X_list[["Female"]]

row_means_Y_male <- Y_list[["Male"]]
row_means_Y_female <- Y_list[["Female"]]

row_means_Z_male <- Z_list[["Male"]]
row_means_Z_female <- Z_list[["Female"]]

dim(row_means_X_male)
dim(row_means_X_female)
dim(row_means_Y_male)
dim(row_means_Y_female)
dim(row_means_Z_male)
dim(row_means_Z_female)

Domain <- 0:100

build_raw_df <- function(raw_data, group, axis_name) {
  
  raw_data %>%
    mutate(time = Domain) %>%              # add domain
    pivot_longer(
      cols = -time,
      names_to = "id",
      values_to = "value"
    ) %>%
    mutate(
      Group = group,
      Axis = axis_name,
      Panel = "Raw"
    )
}


df_raw <- bind_rows(
  build_raw_df(row_means_X_male,   "Male",   "X"),
  build_raw_df(row_means_X_female, "Female", "X"),
  build_raw_df(row_means_Y_male,   "Male",   "Y"),
  build_raw_df(row_means_Y_female, "Female", "Y"),
  build_raw_df(row_means_Z_male,   "Male",   "Z"),
  build_raw_df(row_means_Z_female, "Female", "Z")
)


mean_male_X <- rowMeans(row_means_X_male, na.rm = TRUE)
mean_female_X <- rowMeans(row_means_X_female, na.rm = TRUE)

mean_male_Y <- rowMeans(row_means_Y_male, na.rm = TRUE)
mean_female_Y <- rowMeans(row_means_Y_female, na.rm = TRUE)

mean_male_Z <- rowMeans(row_means_Z_male, na.rm = TRUE)
mean_female_Z <- rowMeans(row_means_Z_female, na.rm = TRUE)



df_plot_X <- data.frame(
  time = Domain,
  Male = mean_male_X,
  Female = mean_female_X
)

df_plot_Y <- data.frame(
  time = Domain,
  Male = mean_male_Y,
  Female = mean_female_Y
)

df_plot_Z <- data.frame(
  time = Domain,
  Male = mean_male_Z,
  Female = mean_female_Z
)



# Combine mean data
df_means <- bind_rows(
  df_plot_X  %>% mutate(Axis = "X"),
  df_plot_Y %>% mutate(Axis = "Y"),
  df_plot_Z %>% mutate(Axis = "Z")
) %>%
  pivot_longer(cols = c(Male, Female),
               names_to = "Sex",
               values_to = "Value") %>%
  mutate(Row = "Mean")

sd_male_X <- apply(row_means_X_male, 1, sd)
mean(sd_male_X)
sd_female_X <- apply(row_means_X_female, 1, sd)
mean(sd_female_X)

mean(apply(cbind(row_means_X_male,row_means_X_female), 1, sd))

sd_male_Y <- apply(row_means_Y_male, 1, sd)
mean(sd_male_Y)
sd_female_Y <- apply(row_means_Y_female, 1, sd)
mean(sd_female_Y)

mean(apply(cbind(row_means_Y_male,row_means_Y_female), 1, sd))

sd_male_Z <- apply(row_means_Z_male, 1, sd)
mean(sd_male_Z)
sd_female_Z <- apply(row_means_Z_female, 1, sd)
mean(sd_female_Z)

mean(apply(cbind(row_means_Z_male,row_means_Z_female), 1, sd))

residuals_data <- function(data) {
  # Subtract the column-wise means of each group
  r <- data - matrix(colMeans(data), nrow = nrow(data), ncol = ncol(data), byrow = TRUE)
  return(as.data.frame(r))
}


r_male_X <- residuals_data(t(row_means_X_male))
r_female_X <- residuals_data(t(row_means_X_female))

r_male_Y <- residuals_data(t(row_means_Y_male))
r_female_Y <- residuals_data(t(row_means_Y_female))

r_male_Z <- residuals_data(t(row_means_Z_male))
r_female_Z <- residuals_data(t(row_means_Z_female))


df_sd <- bind_rows(
  data.frame(time = Domain,
             mean = mean_male_X,
             sd = sd_male_X,
             Group = "Male",
             Axis = "X"),
  data.frame(time = Domain,
             mean = mean_female_X,
             sd = sd_female_X,
             Group = "Female",
             Axis = "X"),
  data.frame(time = Domain,
             mean = mean_male_Y,
             sd = sd_male_Y,
             Group = "Male",
             Axis = "Y"),
  data.frame(time = Domain,
             mean = mean_female_Y,
             sd = sd_female_Y,
             Group = "Female",
             Axis = "Y"),
  data.frame(time = Domain,
             mean = mean_male_Z,
             sd = sd_male_Z,
             Group = "Male",
             Axis = "Z"),
  data.frame(time = Domain,
             mean = mean_female_Z,
             sd = sd_female_Z,
             Group = "Female",
             Axis = "Z")
) %>%
  mutate(Panel = "MeanSD")


build_residual_df <- function(r_data, group, axis_name) {
  
  r_data %>%
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
      Panel = "Residuals"
    ) %>%
    select(time, value, id, Group, Axis, Panel)
}

df_res <- bind_rows(
  build_residual_df(r_male_X, "Male", "X"),
  build_residual_df(r_female_X, "Female", "X"),
  build_residual_df(r_male_Y, "Male", "Y"),
  build_residual_df(r_female_Y, "Female", "Y"),
  build_residual_df(r_male_Z, "Male", "Z"),
  build_residual_df(r_female_Z, "Female", "Z")
)



axis_labels <- c(X = "Sagittal",
                 Y = "Frontal",
                 Z = "Transverse")

panel_labels <- c(
  Raw = "Data",
  MeanSD = "Mean ± SD",
  Residuals = "Residuals"
)

df_raw$Axis <- factor(df_raw$Axis, levels = names(axis_labels), labels = axis_labels)
df_sd$Axis  <- factor(df_sd$Axis,  levels = names(axis_labels), labels = axis_labels)
df_res$Axis <- factor(df_res$Axis, levels = names(axis_labels), labels = axis_labels)

df_raw$Panel <- factor(df_raw$Panel, levels = names(panel_labels), labels = panel_labels)
df_sd$Panel  <- factor(df_sd$Panel,  levels = names(panel_labels), labels = panel_labels)
df_res$Panel <- factor(df_res$Panel, levels = names(panel_labels), labels = panel_labels)

####################Y_SCALE AND LIMITS####################
# Disable lims and y_scale if one wants to use free_y in facet_grid.
#Note: free_y will make the y-axis limits different across panels, 
#which can be useful for visualizing patterns within each panel but may hinder 
#direct comparisons across panels.
lims <- bind_rows(
  df_raw %>% select(Axis, value),
  df_res %>% select(Axis, value),
  df_sd %>% transmute(Axis, value = mean + sd),
  df_sd %>% transmute(Axis, value = mean - sd)
) %>%
  group_by(Axis) %>%
  summarise(ymin = min(value), ymax = max(value))


y_scales <- list(
  scale_y_continuous(limits = lims[lims$Axis == "Sagittal",  c("ymin","ymax")] |> unlist()),
  scale_y_continuous(limits = lims[lims$Axis == "Frontal",   c("ymin","ymax")] |> unlist()),
  scale_y_continuous(limits = lims[lims$Axis == "Transverse",c("ymin","ymax")] |> unlist())
)
#######################
plot2 <- ggplot() +
  
  # --- RAW ---
  geom_line(
    data = df_raw,
    aes(time, value,
        group = interaction(id, Group),
        color = Group)
  ) +
  
  # --- MEAN + SD ---
  geom_ribbon(data = df_sd,
              aes(time,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  fill = Group),
              alpha = 0.1) +
  
  geom_line(data = df_sd,
            aes(time, mean, color = Group),
            linewidth = 1.5) +
  
  # --- RESIDUALS ---
  geom_line(data = df_res,
            aes(time, value,
                group = interaction(id, Group),
                color = Group),
            alpha = 0.3) +
  
  scale_color_manual(
    name = "Groups:",
    values = c("Male" = "cadetblue",
               "Female" = "tomato"),
    labels = c("Male" = "Group 2",
               "Female" = "Group 1"),
    guide = guide_legend(nrow = 1)
  ) +
  
  scale_fill_manual(
    name = "Groups:",
    values = c("Male" = "cadetblue",
               "Female" = "tomato"),
    labels = c("Male" = "Group 2",
               "Female" = "Group 1")
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

print(plot2)



ggsave("outputs/Figure2.jpeg", plot2, width = 190, height = 180, units = "mm", dpi = 1200)











