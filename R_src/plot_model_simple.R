rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

if(!is.null(dev.list())){
    dev.off()
}

## cmd=paste("./SweeD -name m", migrate, " -grid 10 -input mssel_mig2.out -length 100000 -strictPolymorphic -mssel 20 1000 0 20 ./trajectory_continent_island.txt 500 -I 2 0 20 0 0 ", migrate, " -t 2000 -r 500 1000 -ej 2 2 1", sep="")
## -en 0.02 2 0.01 -en 0.1 1 0.1"
##system(cmd)

suffix="continent_island_twoway_bottleneck3"

if(length(args) > 1){
    suffix <- args[1]
}

name= paste("SweeD_SimulationModel.", suffix, sep="")
name
a <- read.table(name, header=F)[,-1]
pos <- as.numeric(a[1,])
maxvalue <- max(a[3:nrow(a),])

b <- t(apply(a[3:nrow(a),], 1, log))
maxvalue <- max(b)
minvalue <- min(b)
pdf(paste("SweeD_plotModel.", suffix, ".pdf", sep=""), height=40, width = 10 )
plot(pos, as.numeric(b[1,]), type='l', ylim=c(minvalue, maxvalue), xlim = c(0, 13))
lastindex <- length(b[1,])
text(x=13, y=b[1,lastindex], labels=1, col="red")



for(i in 2:(nrow(b))){
    lastindex <- length(b[i,])
    points(pos, as.numeric(b[i,]), type='l', col="black")
    text(x=13, y=b[i,lastindex], labels=i, col="red")
}
dev.off()





