rm(list=ls())

numTransmAll <- c(0 ,1,2,3,5,13)
frqTransmAll <- c(46,6,1,1,1,1)
transmAll <- rep(numTransmAll, frqTransmAll)

numTransmTrav <- c(0,2,5,13)
frqTransmTrav <- c(4,1,1,1)
transmTrav <- rep(numTransmTrav, frqTransmTrav)

numTransmEvac <- c(0 ,1)
frqTransmEvac <- c(19,1)
transmEvac <- rep(numTransmEvac, frqTransmEvac)

numTransmLoc <- c(0 ,1,3)
frqTransmLoc <- c(23,5,1)
transmLoc <- rep(numTransmLoc, frqTransmLoc)

getR0andk <- function(transm){
  R0 <- mean(transm)
  k <- R0^2/(var(transm)-R0)
  c(R0=R0,k=k)
}

getMeanTableRow <- function(transm){
  c(num = length(transm), transm = sum(transm), getR0andk(transm))
}

meanTable <- rbind(getMeanTableRow(transmAll),
                   getMeanTableRow(transmTrav),
                   getMeanTableRow(transmEvac),
                   getMeanTableRow(transmLoc))
rownames(meanTable) <- c("All","Traveler","Evacuated","LocallyAcquired")

print(meanTable)