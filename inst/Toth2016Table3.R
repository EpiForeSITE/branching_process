rm(list=ls())

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
 
#p is the probability of y total transmission from x independent patients
lp <- function(x,y,R,k) lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k)
p <- function(x,y,R,k) exp(lp(x,y,R,k))

lpPois <- function(x,y,R) y*(log(x)+log(R))-x*R-lfactorial(y)
pPois <- function(x,y,R) exp(lpPois(x,y,R))

qAny <- function(n,j,R,k) p(j,j-n,R,k)*n/j
qAnyPois <- function(n,j,R) pPois(j,j-n,R)*n/j

qAnyGen1 <- function(n,j,R0,k0,Rc,kc=k0) ifelse(j==n,p(n,0,R0,k0),sum(p(n,1:(j-n),R0,k0)*qAny(1:(j-n),j-n,Rc,kc)))
qAnyPoisGen1 <- function(n,j,R0,k0,Rc) ifelse(j==n,p(n,0,R0,k0),sum(p(n,1:(j-n),R0,k0)*qAnyPois(1:(j-n),j-n,Rc)))

qAnyGen2 <- function(n,j,R0,k0,Rc,kc=k0){
	f <- function(x) p(n,x,R0,k0) * sum(p(x,1:(j-n-x),R0,k0) * qAny(1:(j-n-x),j-n-x,Rc,kc))
	ifelse(j==n,p(n,0,R0,k0),ifelse(j == n+1,p(n,1,R0,k0) * p(1,0,R0,k0), 
		 p(n,j-n,R0,k0) * p(j-n,0,R0,k0) + sum(Vectorize(f)(1:(j-n-1)))))
}
qAnyPoisGen2 <- function(n,j,R0,k0,Rc){
	f <- function(x) p(n,x,R0,k0) * sum(p(x,1:(j-n-x),R0,k0) * qAnyPois(1:(j-n-x),j-n-x,Rc))
	ifelse(j==n,p(n,0,R0,k0),ifelse(j == n+1,p(n,1,R0,k0) * p(1,0,R0,k0), 
		 p(n,j-n,R0,k0) * p(j-n,0,R0,k0) + sum(Vectorize(f)(1:(j-n-1)))))
}


dist0 <- qAny(1,1:20000,Rest,kest)
dist1a <- Vectorize(qAnyGen1)(1,1:1000,R0est1,kest1,Rcest1)
dist1b <- Vectorize(qAnyGen1)(1,1:1000,R0estGeom1,k0estGeom1,RcestGeom1,1)
dist1c <- Vectorize(qAnyPoisGen1)(1,1:1000,R0estPois1,k0estPois1,RcestPois1)

dist2a <- Vectorize(qAnyGen2)(1,1:1000,R0est2,kest2,Rcest2)
dist2b <- Vectorize(qAnyGen2)(1,1:1000,R0estGeom2,k0estGeom2,RcestGeom2,1)
dist2c <- Vectorize(qAnyPoisGen2)(1,1:1000,R0estPois2,k0estPois2,RcestPois2)

tableTotVals <- c(10,100,500,1000)
tableLeft <- rbind((1-cumsum(dist0))[tableTotVals]*100,
			(1-cumsum(dist1a))[tableTotVals]*100,
			(1-cumsum(dist1b))[tableTotVals]*100,
			(1-cumsum(dist1c))[tableTotVals]*100,
			(1-cumsum(dist2a))[tableTotVals]*100,
			(1-cumsum(dist2b))[tableTotVals]*100,
			(1-cumsum(dist2c))[tableTotVals]*100)

#pgl[h] is the probability of observing less than h generations of transmission
getPgl <- function(R0,Rc,k0,kc,genSwitch=Inf){
	maxGen <- 100
	R <- rep(R0,maxGen); k <- rep(k0,maxGen)
	if(genSwitch < maxGen){
		R[(genSwitch+1):maxGen] <- Rc
		k[(genSwitch+1):maxGen] <- kc
	}
	pgl <- (1+R[1]/k[1])^(-k[1])
	for(h in 2:maxGen) pgl[h] <- (1+R[h]/k[h]*(1-pgl[h-1]))^(-k[h])	
	pgl
}
getPglPois <- function(R0,Rc,k0,genSwitch=Inf){
	maxGen <- 100
	pgl <- (1+R0/k0)^(-k0)
	for(h in 2:maxGen) pgl[h] <- ifelse(h <= genSwitch,(1+R0/k0*(1-pgl[h-1]))^(-k0),
									   exp(-Rc*(1-pgl[h-1])))	
	pgl
}

pgl <- getPgl(Rest,Rest,kest,kest)
pgl1a <- getPgl(R0est1,Rcest1,kest1,kest1,1)
pgl1b <- getPgl(R0estGeom1,RcestGeom1,k0estGeom1,1,1)
pgl1c <- getPglPois(R0estPois1,RcestPois1,k0estPois1,1)

pgl2a <- getPgl(R0est2,Rcest2,kest2,kest2,2)
pgl2b <- getPgl(R0estGeom2,RcestGeom2,k0estGeom2,1,2)
pgl2c <- getPglPois(R0estPois2,RcestPois2,k0estPois2,2)

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
