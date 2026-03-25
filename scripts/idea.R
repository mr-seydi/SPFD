y1=c(rep(0,45),1,2,1.1,1.2,2,rep(0,45),1,2,0.1,0.2,0,0)
length(y1)
y2=c(rep(0,45),0,0,0.2,0.3,1.5,rep(0,45),1,1,1.1,1.2,0,1)
length(y2)

plot(y1,type="l",col="blue",ylim=c(0,2.5),xlab="Time",ylab="Value",main="Line Plot of y1 and y2")
lines(y2,col="red")

diff=y1-y2
plot(diff^2,type="l",col="green",ylim=c(-2,2))

L2_norm=sqrt(sum(diff^2))

y1_norm=sqrt(sum(y1^2))
y2_norm=sqrt(sum(y2^2))



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

sd(res$L2)
mean(res$L2)
max(res$L2)
min(res$L2)
which.max(res$L2)
res[which.max(res$L2), ]

hist(res$L2, breaks=50, main="Histogram of L2 Norms", xlab="L2 Norm")

sel <- res[res$L2 > mean(res$L2) + 2.5*sd(res$L2), c("start","end")]
n <- length(y1)
coverage <- integer(n)

for (k in 1:nrow(sel)) {
  coverage[ sel$start[k] : sel$end[k] ] <-
    coverage[ sel$start[k] : sel$end[k] ] + 1
}

coverage

support <- coverage > 0
support <- coverage > quantile(coverage[coverage>0], 0.5)
support <- coverage > 0.9 * max(coverage)


idx <- which(support)

regions <- split(idx, cumsum(c(1, diff(idx) != 1)))
regions <- lapply(regions, range)

do.call(rbind, regions)
plot(coverage, type="l", ylab="Interval coverage", xlab="Time")
abline(h=0, col="gray")



res$len <- res$end - res$start + 1
res$RMS <- res$L2 / sqrt(res$len)


res$signal_rms <- sapply(1:nrow(res), function(k) {
  i <- res$start[k]
  j <- res$end[k]
  sqrt(mean(y1[i:j]^2 + y2[i:j]^2))
})



theta <- 0.5 * L2_norm / sqrt(n)
which(res$RMS > theta)

res$rel_diff <- res$RMS / res$signal_rms

substantial <- res[res$rel_diff > 0.5, ]


score <- numeric(n)

for (k in which(res$rel_diff > 0.5)) {
  score[res$start[k]:res$end[k]] <-
    score[res$start[k]:res$end[k]] + 1
}

plot(score, type="l", ylab="Substantial difference support")

