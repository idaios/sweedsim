rm(list=ls())
dev.off()

migrates=c(0.01, 0.1, 1, 2, 10, 40, 80, 200, 1000)

for( migrate in migrates){
    print(migrate)
    cmd=paste("./SweeD -name m", migrate, " -grid 10 -input mssel_mig2.out -length 100000 -strictPolymorphic -mssel 20 1000 0 20 ./trajectory_continent_island.txt 500 -I 2 0 20 0 0 ", migrate, " -t 2000 -r 500 1000 -ej 2 2 1", sep="")
    ## -en 0.02 2 0.01 -en 0.1 1 0.1"
    ##system(cmd)
    name= paste("SweeD_SimulationModel.m", migrate, sep="")
        a <- read.table(name, header=F)[,-1]
    pos <- as.numeric(a[1,])
    maxvalue <- max(a[3:20,])
    pdf(paste("SweeD_plotModel_m", migrate, ".pdf", sep="") )
    plot(pos, as.numeric(a[3,]), type='l', ylim=c(0, maxvalue), xlim = c(0, 13))
    lastindex <- length(a[3,])
    text(x=13, y=a[3,lastindex], labels=1, col="red")
    for(i in 4:20){
        lastindex <- length(a[i,])
        points(pos, as.numeric(a[i,]), type='l', col="black")
        text(x=13, y=a[i,lastindex], labels=i-2, col="red")
    }
    dev.off()
}



