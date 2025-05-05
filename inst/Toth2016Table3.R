Rest <- 0.8652294
kest <- 0.03539390
R0est1 <- 5.2123787
kest1 <- 0.06778359
Rcest1 <- 0.18875788 
R0estGeom1 <- 5.4998656
k0estGeom1 <- 0.06113228 
RcestGeom1 <- 0.14344697
R0estPois1 <- 5.5349927 
k0estPois1 <- 0.06104179
RcestPois1 <- 0.13782021
R0est2 <- 2.0233019 
kest2 <- 0.07788500
Rcest2 <- 0.06386585 
R0estGeom2 <- 2.1945074
k0estGeom2 <- 0.07574887
RcestGeom2 <- 0.06033546
R0estPois2 <- 2.1966033
k0estPois2 <- 0.07574568
RcestPois2 <- 0.05991646
 
dist0 <- pFinalSize(1,1:1000,Rest,kest)
dist1a <- Vectorize(pFinalSizeSwitch1)(1,1:1000,R0est1,kest1,Rcest1,kest1)
dist1b <- Vectorize(pFinalSizeSwitch1)(1,1:1000,R0estGeom1,k0estGeom1,RcestGeom1,1)
dist1c <- Vectorize(pFinalSizeSwitch1)(1,1:1000,R0estPois1,k0estPois1,RcestPois1,Inf)
dist2a <- Vectorize(pFinalSizeSwitch2)(1,1:1000,R0est2,kest2,Rcest2,kest2)
dist2b <- Vectorize(pFinalSizeSwitch2)(1,1:1000,R0estGeom2,k0estGeom2,RcestGeom2,1)
dist2c <- Vectorize(pFinalSizeSwitch2)(1,1:1000,R0estPois2,k0estPois2,RcestPois2,Inf)

tableTotVals <- c(10,100,500,1000)
tableLeft <- rbind((1-cumsum(dist0))[tableTotVals]*100,
			(1-cumsum(dist1a))[tableTotVals]*100,
			(1-cumsum(dist1b))[tableTotVals]*100,
			(1-cumsum(dist1c))[tableTotVals]*100,
			(1-cumsum(dist2a))[tableTotVals]*100,
			(1-cumsum(dist2b))[tableTotVals]*100,
			(1-cumsum(dist2c))[tableTotVals]*100)

maxGen <- 5
pgl <- pGen(maxGen,Rest,kest)
pgl1a <- pGenSwitch1(maxGen,R0est1,kest1,Rcest1,kest1)
pgl1b <- pGenSwitch1(maxGen,R0estGeom1,k0estGeom1,RcestGeom1,1)
pgl1c <- pGenSwitch1(maxGen,R0estPois1,k0estPois1,RcestPois1,Inf)
pgl2a <- pGenSwitch2(maxGen,R0est2,kest2,Rcest2,kest2)
pgl2b <- pGenSwitch2(maxGen,R0estGeom2,k0estGeom2,RcestGeom2,1)
pgl2c <- pGenSwitch2(maxGen,R0estPois2,k0estPois2,RcestPois2,Inf)

#pgg[h] is the probability of observing at least h generations of transmission		
pgg <- 1-pgl
pgg1a <- 1-pgl1a
pgg1b <- 1-pgl1b
pgg1c <- 1-pgl1c

pgg2a <- 1-pgl2a
pgg2b <- 1-pgl2b
pgg2c <- 1-pgl2c

tableGenVals <- 1:5
tableRight <- rbind(pgg[tableGenVals],pgg1a[tableGenVals],pgg1b[tableGenVals],pgg1c[tableGenVals],
			 			 pgg2a[tableGenVals],pgg2b[tableGenVals],pgg2c[tableGenVals])*100

Table3 <- cbind(tableLeft,tableRight)

rownames(Table3) <- c('Model0','Model1a','Model1b','Model1c','Model2a','Model2b','Model2c')
colnames(Table3) <- c(paste0('tr',tableTotVals), paste0('gen',tableGenVals))

print(Table3)
