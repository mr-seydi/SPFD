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


set.seed(123)
# min_max <- function(y1, y2){
#   min_y1 <- min(y1)
#   max_y1 <- max(y1)
#   min_y2 <- min(y2)
#   max_y2 <- max(y2)
#   
# }



y1=mean_male_X_101
y2=mean_female_X_101


# x : numeric vector of domain points
# y1, y2 : numeric vectors of mean function values
# tau : threshold in (0,1)

stopifnot(length(y1) == length(y2))
n <- length(y1)

x=seq(0, 1, length.out = n)

# grid spacing (assumed uniform)
dx <- if (n > 1) mean(diff(x)) else 1


rel_diff <- abs(y1 - y2) / (abs(y1) + abs(y2))

# avoid 0/0
rel_diff[is.nan(rel_diff)] <- 0

flag_rel <- rel_diff > tau

plot(x, rel_diff, type = "l", main = "Relative Difference", ylab = "Relative Difference",
     xlab = "x", ylim = c(0, max(rel_diff)*1.1))


scale_Linf <- max(max(abs(y1)), max(abs(y2)))

tau <- 0.1

flag_Linf <- abs(y1 - y2) > tau * scale_Linf

plot(x, abs(y1 - y2), type = "l", main = "Relative Difference", ylab = "Relative Difference",
     xlab = "x", ylim = c(0, max(abs(y1 - y2))*1.1))
abline(h=tau*scale_Linf)


L2_norm <- function(y, dx) sqrt(sum(y^2) * dx)

scale_L2 <- max(L2_norm(y1, dx), L2_norm(y2, dx))

flag_L2 <- abs(y1 - y2) > tau * scale_L2
plot(x, abs(y1 - y2), type = "l", main = "Relative Difference", ylab = "Relative Difference",
     xlab = "x", ylim = c(0, max(abs(y1 - y2))*1.1))
abline(h=tau*scale_L2)


## Good

diff_energy <- (y1 - y2)^2
total_energy <- sum(diff_energy) * dx
cum_energy <- cumsum(diff_energy) * dx
flag_energy <- diff_energy > tau * total_energy
plot(x, diff_energy, type = "l", main = "Energy Difference", ylab = "Energy Difference",
     xlab = "x", ylim = c(0, max(diff_energy)*1.1))
tau=0.5
abline(h=tau*total_energy)







global_L2 <- sqrt(sum((y1 - y2)^2))
global_L2
dif = (y2-y1)^2
plot(dif,type="l",ylim = c(0,max(dif)))
abline(h=max(dif)*0.6)

max(y1)
max(y2)

global_L2 <- sqrt(mean((y1 - y2)^2))
global_L2
max(dif)*0.6
abline(h=max(dif)*0.6)
total=dif*global_L2

plot(total,type="l",ylim = c(0,max(total)))
max(total)*0.6
abline(h=max(total)*0.6)

sorted_total=sort(total,decreasing = TRUE)
ss=sum(total)
ss95=0.95*ss
cum_sum=cumsum(sorted_total)
cum_sum/cum_sum[length(cum_sum)]*100
abline(h=sorted_total[5])



finder2 <- function(y1, y2, threshold = 1){
  
  n <- length(y1)
  delta <- y1 - y2
  e <- delta^2
  
  e_bar <- mean(e)
  s_e <- sd(e)
  
  N <- n * (n + 1) / 2
  start <- integer(N)
  end   <- integer(N)
  score <- numeric(N)
  
  idx <- 1
  for (i in 1:n) {
    for (j in i:n) {
      local_mean <- mean(e[i:j])
      score[idx] <- (local_mean - e_bar) / s_e
      start[idx] <- i
      end[idx] <- j
      idx <- idx + 1
    }
  }
  
  res <- data.frame(start, end, score)
  
  if (!is.null(threshold)) {
    res <- res[res$score > threshold, ]
  }
  
  sel <- res[, c("start","end")]
  
  
  n <- length(y1)
  coverage <- integer(n)
  
  for (k in 1:nrow(sel)) {
    coverage[ sel$start[k] : sel$end[k] ] <-
      coverage[ sel$start[k] : sel$end[k] ] + 1
  }
  
  
  return(coverage)
}

y1=mean_male_Y_101
y2=mean_female_Y_101

finder2(y1, y2, threshold = 0.42)




finder <- function(y1, y2, threshold = 1){
  
  n <- length(y1)
  diff_full <- y1 - y2
  
  global_rms <- sqrt(mean(diff_full^2))
  
  N <- n * (n + 1) / 2
  start <- integer(N)
  end   <- integer(N)
  score <- numeric(N)
  
  idx <- 1
  for (i in 1:n) {
    for (j in i:n) {
      diff <- diff_full[i:j]
      local_rms <- sqrt(mean(diff^2))
      
      start[idx] <- i
      end[idx]   <- j
      score[idx] <- local_rms / global_rms
      idx <- idx + 1
    }
  }
  
  threshold=0.2
  
  res <- data.frame(start, end, score)
  
  if (!is.null(threshold)) {
    res <- res[res$score > threshold, ]
  }
  
  sel <- res[, c("start","end")]
  
  
  n <- length(y1)
  coverage <- integer(n)
  
  for (k in 1:nrow(sel)) {
    coverage[ sel$start[k] : sel$end[k] ] <-
      coverage[ sel$start[k] : sel$end[k] ] + 1
  }
  
  coverage
  
  return(res)
}




finder <- function(y1, y2, threshold = NULL){
  
  n <- length(y1)
  N <- n * (n + 1) / 2
  
  start <- integer(N)
  end   <- integer(N)
  L2    <- numeric(N)
  
  idx <- 1
  for (i in 1:n) {
    for (j in i:n) {
      diff <- y1[i:j] - y2[i:j]
      start[idx] <- i
      end[idx]   <- j
      L2[idx]    <- sqrt(sum(diff^2)/length(diff))
      idx <- idx + 1
    }
  }
  
  res <- data.frame(start, end, L2)
  # 
  # dim(res)
  # sd(res$L2)
  # mean(res$L2)
  # max(res$L2)
  # min(res$L2)
  # which.max(res$L2)
  # res[which.max(res$L2), ]
  
  
  
  #hist(res$L2, breaks=50, main="Histogram of L2 Norms", xlab="L2 Norm")
  
  sel <- res[res$L2 > mean(res$L2) + 1.5*sd(res$L2), c("start","end")]
  
  global_L2 <- sqrt(mean((y1 - y2)^2))
  
  res$score <- res$L2 / global_L2
  
  tau <- quantile(res$score, 0.50)
  sel <- res[res$score > tau, c("start","end")]
  
  
  
  n <- length(y1)
  coverage <- integer(n)
  
  for (k in 1:nrow(sel)) {
    coverage[ sel$start[k] : sel$end[k] ] <-
      coverage[ sel$start[k] : sel$end[k] ] + 1
  }
  
  coverage
  
  support <- coverage > 0
  support <- coverage >= quantile(coverage[coverage>0], 0.05)
  support <- coverage > 0.9 * max(coverage)
  
  idx <- which(support)
  
  regions <- split(idx, cumsum(c(1, diff(idx) != 1)))
  regions <- lapply(regions, range)
  
  do.call(rbind, regions)
  
}

finder(y1, y2)
