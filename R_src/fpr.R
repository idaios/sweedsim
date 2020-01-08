a <- read.table("likelihood_neutral58.txt")
b <- read.table("likelihood_mssel58.txt")
sum(b[,1] > quantile(a[,1], probs=0.95))
