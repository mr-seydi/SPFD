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

n_X_male <- ncol(row_means_X_male)
n_X_female <- ncol(row_means_X_female)

n_Y_male <- ncol(row_means_Y_male)
n_Y_female <- ncol(row_means_Y_female)

n_Z_male <- ncol(row_means_Z_male)
n_Z_female <- ncol(row_means_Z_female)



mean_male_X <- rowMeans(row_means_X_male, na.rm = TRUE)
mean_female_X <- rowMeans(row_means_X_female, na.rm = TRUE)

mean_male_Y <- rowMeans(row_means_Y_male, na.rm = TRUE)
mean_female_Y <- rowMeans(row_means_Y_female, na.rm = TRUE)

mean_male_Z <- rowMeans(row_means_Z_male, na.rm = TRUE)
mean_female_Z <- rowMeans(row_means_Z_female, na.rm = TRUE)


sd_male_X <- apply(row_means_X_male, 1, sd)
sd_female_X <- apply(row_means_X_female, 1, sd)

sd_male_Y <- apply(row_means_Y_male, 1, sd)
sd_female_Y <- apply(row_means_Y_female, 1, sd)

sd_male_Z <- apply(row_means_Z_male, 1, sd)
sd_female_Z <- apply(row_means_Z_female, 1, sd)


Domain <- 0:100


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

df_sd <- bind_rows(
  data.frame(time = Domain,
             mean = mean_male_X,
             sd = sd_male_X,
             Sex = "Male",
             Axis = "X"),
  data.frame(time = Domain,
             mean = mean_female_X,
             sd = sd_female_X,
             Sex = "Female",
             Axis = "X"),
  data.frame(time = Domain,
             mean = mean_male_Y,
             sd = sd_male_Y,
             Sex = "Male",
             Axis = "Y"),
  data.frame(time = Domain,
             mean = mean_female_Y,
             sd = sd_female_Y,
             Sex = "Female",
             Axis = "Y"),
  data.frame(time = Domain,
             mean = mean_male_Z,
             sd = sd_male_Z,
             Sex = "Male",
             Axis = "Z"),
  data.frame(time = Domain,
             mean = mean_female_Z,
             sd = sd_female_Z,
             Sex = "Female",
             Axis = "Z")
) %>%
  mutate(Panel = "Mean")

compute_tau_df <- function(m1, m2, sd1, sd2, n1, n2, axis_name,
                           taus = seq(0.1, 1, by = 0.1)) {
  
  #pooled_sd <- sqrt((sd1^2 + sd2^2) / 2) # We need to assume equal sample sizes 
  pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  #for the Absolute Standardized Mean Difference (ASMD) , so we can use the average of the two SDs
  effect <- abs(m1 - m2) / pooled_sd
  
  df_effect <- data.frame(
    time = seq_along(effect)-1, # Adjust to match Domain (0-100)
    value = effect,
    Axis = axis_name,
    Row = "SAE"
  )
  
  df_tau <- data.frame(
    tau = taus,
    Axis = axis_name,
    Row = "SAE"
  )
  
  list(effect = df_effect,
       tau = df_tau)
}

tau_X <- compute_tau_df(mean_male_X, mean_female_X, sd_male_X, sd_female_X,
                        n_X_male, n_X_female, "X")
tau_Y <- compute_tau_df(mean_male_Y, mean_female_Y, sd_male_Y, sd_female_Y,
                        n_Y_male, n_Y_female, "Y")
tau_Z <- compute_tau_df(mean_male_Z, mean_female_Z, sd_male_Z, sd_female_Z,
                        n_Z_male, n_Z_female, "Z")

df_tau_effect <- bind_rows(
  tau_X$effect,
  tau_Y$effect,
  tau_Z$effect
)

df_tau_lines <- bind_rows(
  tau_X$tau,
  tau_Y$tau,
  tau_Z$tau
) %>%
  mutate(Panel = "SAE")

df_tau_lines <- df_tau_lines %>%
  group_by(Axis) %>%
  mutate(tau_id = row_number()) %>%
  ungroup()

compute_D1_df <- function(m1, m2, sd1, sd2, n1, n2, axis_name,
                          taus = seq(0.1, 1, by = 0.1)) {
  
  #pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  effect <- abs(m1 - m2) / pooled_sd
  
  results <- list()
  
  for (i in seq_along(taus)) {
    tau <- taus[i]
    idx <- which(effect > tau)
    
    if (length(idx) > 0) {
      runs <- split(idx, cumsum(c(1, diff(idx) != 1)))
      
      for (r in runs) {
        results[[length(results)+1]] <-
          data.frame(
            xmin = min(r)-1,
            xmax = max(r)-1,
            tau_id = i,
            Axis = axis_name,
            Row = "D1prime"
          )
      }
    }
  }
  
  bind_rows(results)
}

df_D1 <- bind_rows(
  compute_D1_df(mean_male_X, mean_female_X, sd_male_X, sd_female_X,
                n_X_male, n_X_female, "X"),
  compute_D1_df(mean_male_Y, mean_female_Y, sd_male_Y, sd_female_Y,
                n_Y_male, n_Y_female, "Y"),
  compute_D1_df(mean_male_Z, mean_female_Z, sd_male_Z, sd_female_Z,
                n_Z_male, n_Z_female, "Z")
)




df_means2 <- df_means %>%
  rename(value = Value) %>%
  mutate(
    Panel = "Mean",
    time = time
  ) %>%
  select(time, value, Axis, Panel, Sex)


df_tau_effect2 <- df_tau_effect %>%
  mutate(
    Panel = "SAE",
    Sex = NA
  ) %>%
  select(time, value, Axis, Panel, Sex)

df_all_lines <- bind_rows(df_means2, df_tau_effect2)
df_D1 <- df_D1 %>% mutate(Panel = "D1prime")

df_all_lines$Panel <- factor(
  df_all_lines$Panel,
  levels = c("Mean", "SAE", "D1prime"),
  labels = c("Mean ~ '\u00B1' ~ SD", "SAE", "D[1]*\"'\"")
)

df_D1$Panel <- factor(
  df_D1$Panel,
  levels = c("Mean", "SAE", "D1prime"),
  labels = c("Mean ~ '\u00B1' ~ SD", "SAE", "D[1]*\"'\"")
)

df_tau_lines$Panel <- factor(
  df_tau_lines$Panel,
  levels = c("Mean", "SAE", "D1prime"),
  labels = c("Mean ~ '\u00B1' ~ SD", "SAE", "D[1]*\"'\"")
)

df_sd$Panel <- factor(
  df_sd$Panel,
  levels = c("Mean", "SAE", "D1prime"),
  labels = c("Mean ~ '\u00B1' ~ SD", "SAE", "D[1]*\"'\"")
)



df_all_lines$Axis <- factor(
  df_all_lines$Axis,
  levels = c("X", "Y", "Z"),
  labels = c("Sagittal", "Frontal", "Transverse")
)

df_sd$Axis <- factor(
  df_sd$Axis,
  levels = c("X", "Y", "Z"),
  labels = c("Sagittal", "Frontal", "Transverse")
)

df_D1$Axis <- factor(
  df_D1$Axis,
  levels = c("X", "Y", "Z"),
  labels = c("Sagittal", "Frontal", "Transverse")
)

df_tau_lines$Axis <- factor(
  df_tau_lines$Axis,
  levels = c("X", "Y", "Z"),
  labels = c("Sagittal", "Frontal", "Transverse")
)


paired_10 <- c(
  "#332288", # dark blue
  "#B2DF8A", # light green
  "#A6761D", # beige
  "#CAB2D6", # light purple
  "#C61019", # dark red
  "#88CCEE", # light blue
  "#B200ED", # dark purple
  "#FFFF99",  # light yellow
  "#FF9DA7", # light pink
  "#33A02C" # dark green
)



plot <- ggplot() +
  
  # ---- MEAN ----
  geom_ribbon(data = df_sd,
              aes(time,
                  ymin = mean - sd,
                  ymax = mean + sd,
                  fill = Sex),
              alpha = 0.1) +
  geom_line(data = df_sd,
            aes(x = time, y = mean, color = Sex),
            linewidth = 1.5) +
  scale_fill_manual(
    name = "Groups:",
    values = c("Male" = "cadetblue",
               "Female" = "tomato"),
    labels = c("Male" = "Group 2",
               "Female" = "Group 1"),
    guide = guide_legend(nrow = 2)
  ) +
  scale_color_manual(
    name = "Groups:",
    values = c("Male" = "cadetblue",
               "Female" = "tomato"),
    labels = c("Male" = "Group 2",
               "Female" = "Group 1"),
    guide = guide_legend(nrow = 2)
  ) +
  new_scale_color() +
  
  # ---- SAE ----
geom_line(
  data = subset(df_all_lines, Panel == "SAE"),
  aes(x = time, y = value),
  linewidth = 1, color = "black"
) +
  
  # ---- TAU + D1 MATCHED ----
geom_hline(
  data = df_tau_lines,
  aes(yintercept = tau, color = factor(tau_id)),
  linetype = 2, linewidth = 1
) +
  geom_segment(
    data = df_D1,
    aes(x = xmin, xend = xmax, y = tau_id, yend = tau_id,
        color = factor(tau_id)),
    linewidth = 2
  ) +
  
  scale_color_manual(
    name = expression(paste("Thresholds", " ", delta, ":")),
    values = paired_10,
    labels = seq(0.1, 1, by = 0.1)
  ) +
  
  facet_grid(
    Panel ~ Axis,
    scales = "free_y",
    labeller = labeller(Panel = label_parsed,
    Axis = c(
      X = "Sagittal",
      Y = "Frontal",
      Z = "Transverse"
      )
    )
  ) +
  
  guides(color = guide_legend(override.aes = list(linewidth = 1, linetype = 1))) +
  
  theme_minimal() +
  force_panelsizes(rows = c(1, 1, 0.3)) + #This is a ggh4x function that allows us to set relative heights for the rows of facets. The D1prime row is set to be shorter than the others.
  theme(
    strip.text = element_text(size = 14),
    strip.placement = "outside",
    panel.spacing = unit(1, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.title.y = element_text(hjust = .35)
  ) +
  
  facetted_pos_scales(
    y = list(
      Panel == "D[1]*\"'\"" ~ scale_y_continuous(breaks = NULL)
    )
  ) +
  labs(
    x = "Domain",
    y = "Value"
  )
plot

ggsave("outputs/Figure3.jpeg", plot, width = 190, height = 150, units = "mm", dpi = 1200)




