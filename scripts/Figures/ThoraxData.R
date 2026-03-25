rm(list=ls())
library(readxl)
path_X <- "data/Thorax_X.xlsx"
path_Y <- "data/Thorax_Y.xlsx"
path_Z <- "data/Thorax_Z.xlsx"
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

mean_male_X <- rowMeans(row_means_X_male, na.rm = TRUE)
mean_female_X <- rowMeans(row_means_X_female, na.rm = TRUE)

mean_male_Y <- rowMeans(row_means_Y_male, na.rm = TRUE)
mean_female_Y <- rowMeans(row_means_Y_female, na.rm = TRUE)

mean_male_Z <- rowMeans(row_means_Z_male, na.rm = TRUE)
mean_female_Z <- rowMeans(row_means_Z_female, na.rm = TRUE)


df_plot_X <- data.frame(
  time = seq_len(length(mean_male_X)),
  Male = mean_male_X,
  Female = mean_female_X
)

df_plot_Y <- data.frame(
  time = seq_len(length(mean_male_Y)),
  Male = mean_male_Y,
  Female = mean_female_Y
)

df_plot_Z <- data.frame(
  time = seq_len(length(mean_male_Z)),
  Male = mean_male_Z,
  Female = mean_female_Z
)

library(ggplot2)

ggplot(df_plot_X, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax X Over Time",
    x = "Time",
    y = "Mean Thorax X",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-4, 16)

# Note: -Male and -Female
ggplot(df_plot_Y, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax Y Over Time",
    x = "Time",
    y = "Mean Thorax Y",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-5, 4)


ggplot(df_plot_Z, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax Z Over Time",
    x = "Time",
    y = "Mean Thorax Z",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-8, 6)

set.seed(123)
#install and load fda package
#install.packages("fda")
library(fda)

data_to_101evalpoints <- function(data){
  # data: matrix with rows = original time/eval points, columns = functions
  data_to_fd <- Data2fd(argvals = 0:(nrow(data) - 1), y = data)
  rangeval <- data_to_fd$basis$rangeval
  eval.points <- seq(rangeval[1], rangeval[2], length.out = 101)
  final_data <- eval.fd(fdobj = data_to_fd, evalarg = eval.points)
  return(final_data)
}


row_means_X_male_101 <- data_to_101evalpoints(as.matrix(row_means_X_male))
row_means_X_female_101 <- data_to_101evalpoints(as.matrix(row_means_X_female))

row_means_Y_male_101 <- data_to_101evalpoints(as.matrix(row_means_Y_male))
row_means_Y_female_101 <- data_to_101evalpoints(as.matrix(row_means_Y_female))

row_means_Z_male_101 <- data_to_101evalpoints(as.matrix(row_means_Z_male))
row_means_Z_female_101 <- data_to_101evalpoints(as.matrix(row_means_Z_female))

mean_male_X_101 <- rowMeans(row_means_X_male_101, na.rm = TRUE)
mean_female_X_101 <- rowMeans(row_means_X_female_101, na.rm = TRUE)

mean_male_Y_101 <- rowMeans(row_means_Y_male_101, na.rm = TRUE)
mean_female_Y_101 <- rowMeans(row_means_Y_female_101, na.rm = TRUE)

mean_male_Z_101 <- rowMeans(row_means_Z_male_101, na.rm = TRUE)
mean_female_Z_101 <- rowMeans(row_means_Z_female_101, na.rm = TRUE)



df_plot_X_101 <- data.frame(
  time = seq_len(length(mean_male_X_101)),
  Male = mean_male_X_101,
  Female = mean_female_X_101
)

df_plot_Y_101 <- data.frame(
  time = seq_len(length(mean_male_Y_101)),
  Male = mean_male_Y_101,
  Female = mean_female_Y_101
)

df_plot_Z_101 <- data.frame(
  time = seq_len(length(mean_male_Z_101)),
  Male = mean_male_Z_101,
  Female = mean_female_Z_101
)

library(ggplot2)

ggplot(df_plot_X_101, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax X Over Time",
    x = "Time",
    y = "Mean Thorax X",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-4, 16)

# Note: -Male and -Female
ggplot(df_plot_Y_101, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax Y Over Time",
    x = "Time",
    y = "Mean Thorax Y",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-5, 4)


ggplot(df_plot_Z_101, aes(x = time)) +
  geom_line(aes(y = Male, color = "Male"), linewidth = 1) +
  geom_line(aes(y = Female, color = "Female"), linewidth = 1) +
  labs(
    title = "Mean Thorax Z Over Time",
    x = "Time",
    y = "Mean Thorax Z",
    color = "Sex"
  ) +
  theme_minimal() + ylim(-8, 6)
























