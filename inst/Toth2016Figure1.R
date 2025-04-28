rm(list=ls())

#p is the probability of y total transmission from x independent patients
lp <- function(x,y,R,k) lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k)
p <- function(x,y,R,k) exp(lp(x,y,R,k))

lpPois = function(x,y,R) y*(log(x)+log(R))-x*R-lfactorial(y)
pPois <- function(x,y,R) exp(lpPois(x,y,R))

qAny = function(n,j,R,k) p(j,j-n,R,k)*n/j
qAnyPois = function(n,j,R) pPois(j,j-n,R)*n/j

qAnyGen1 = function(n,j,R0,k0,Rc,kc=k0) ifelse(j==n,p(n,0,R0,k0),sum(p(n,1:(j-n),R0,k0)*qAny(1:(j-n),j-n,Rc,kc)))
qAnyPoisGen1 = function(n,j,R0,k0,Rc) ifelse(j==n,p(n,0,R0,k0),sum(p(n,1:(j-n),R0,k0)*qAnyPois(1:(j-n),j-n,Rc)))

qAnyGen2 = function(n,j,R0,k0,Rc,kc=k0){
	f = function(x) p(n,x,R0,k0) * sum(p(x,1:(j-n-x),R0,k0) * qAny(1:(j-n-x),j-n-x,Rc,kc))
	ifelse(j==n,p(n,0,R0,k0),ifelse(j == n+1,p(n,1,R0,k0) * p(1,0,R0,k0), 
		 p(n,j-n,R0,k0) * p(j-n,0,R0,k0) + sum(Vectorize(f)(1:(j-n-1)))))
}
qAnyPoisGen2 = function(n,j,R0,k0,Rc){
	f = function(x) p(n,x,R0,k0) * sum(p(x,1:(j-n-x),R0,k0) * qAnyPois(1:(j-n-x),j-n-x,Rc))
	ifelse(j==n,p(n,0,R0,k0),ifelse(j == n+1,p(n,1,R0,k0) * p(1,0,R0,k0), 
		 p(n,j-n,R0,k0) * p(j-n,0,R0,k0) + sum(Vectorize(f)(1:(j-n-1)))))
}

#pgl[h] is the probability of observing less than h generations of transmission
getPgl <- function(R0,Rc,k0,kc,genSwitch=Inf){
	maxGen = 100
	R = rep(R0,maxGen); k = rep(k0,maxGen)
	if(genSwitch < maxGen){
		R[(genSwitch+1):maxGen] = Rc
		k[(genSwitch+1):maxGen] = kc
	}
	pgl <- (1+R[1]/k[1])^(-k[1])
	for(h in 2:maxGen) pgl[h] <- (1+R[h]/k[h]*(1-pgl[h-1]))^(-k[h])	
	pgl
}

pClusts <- c(1,.5,.25)
titles <- c('A','B','C','D','E','F')

for(i in 1:length(pClusts)){

	pClust <- pClusts[i]
	if(pClust == 1){
		Rest <- 0.8652294
		kest <- 0.03539390
		R0estGeom1 <- 5.4998656
		k0estGeom1 <- 0.06113228 
		RcestGeom1 <- 0.14344697
		R0estGeom2 <- 2.1945074
		k0estGeom2 <- 0.07574887
		RcestGeom2 <- 0.06033546
	}
	if(pClust == .5){
		Rest <- 0.7631224 
		kest <- 0.02814292
		R0estGeom1 <- 2.7178632 
		k0estGeom1 <- 0.03189325
		RcestGeom1 <- 0.14414662
		R0estGeom2 <- 1.4687318
		k0estGeom2 <- 0.04242669
		RcestGeom2 <- 0.06326285
		x11()
	}
	if(pClust == .25){
		Rest <- 0.6224364 
		kest <- 0.02245021
		R0estGeom1 <- 1.4025546
		k0estGeom1 <- 0.01921805
		RcestGeom1 <- 0.14568040
		R0estGeom2 <- 1.4265324
		k0estGeom2 <- 0.02554593
		RcestGeom2 <- 0.06494508
		x11()
	}
	dist0 = qAny(1,1:20000,Rest,kest)
	dist1b = Vectorize(qAnyGen1)(1,1:1000,R0estGeom1,k0estGeom1,RcestGeom1,1)
	dist2b = Vectorize(qAnyGen2)(1,1:1000,R0estGeom2,k0estGeom2,RcestGeom2,1)
	
	smallticksX = c(2:9,seq(20,90,10),seq(200,900,100),seq(2000,9000,1000))
	smallticksp = c((2:9)*1e-4,(2:9)*1e-3,(2:9)*1e-2,(2:9)*1e-1)
	freqlabels = c('0.01%','0.1%','1%','10%','50%')
	labelticks = c(1,10,100,1000,10000)

	plot(1-cumsum(dist0),log='xy',type='l',ylim=c(1e-4,1),lty=2,lwd=2,
		xlab = 'Total transmissions X',
		ylab = 'Prob. X or more transmissions',yaxt='n',xaxt='n')

	axis(side=4, at=10^seq(-4,0), labels=rep('',5),tck=1,lty=2)
	axis(side=3, at=10^seq(0,4),labels=rep('',5),tck=1,lty=2)
	axis(side=1, at=labelticks, labels=c('1','10','100','1,000','10,000'))
	axis(side=1, at=smallticksX,tck=-0.016,labels=rep('',length(smallticksX)))
	axis(side=2, at=c(10^(-4:-1),0.5), labels=freqlabels, par(mgp = c(3,.5,0)),las=1)
	axis(side=2, at=smallticksp,tck=-0.016,labels=rep('',length(smallticksp)))

	loc <- par('usr')
	text(exp(loc[1]),exp(loc[4]),titles[i*2-1],xpd=TRUE,adj=c(2.5,-4),font=2)

	lines(1-cumsum(dist1b),lty=3,lwd=2)
	lines(1-cumsum(dist2b),lty=1,lwd=2)

	pgl <- getPgl(Rest,Rest,kest,kest)
	pgl1b <- getPgl(R0estGeom1,RcestGeom1,k0estGeom1,1,1)
	pgl2b <- getPgl(R0estGeom2,RcestGeom2,k0estGeom2,1,2)

	#pgg[h] is the probability of observing at least h generations of transmission		
	pgg <- 1-pgl
	pgg1b <- 1-pgl1b
	pgg2b <- 1-pgl2b

	x11()

	plot(pgg,log='y',pch=19,ylim=c(1e-4,1),xlim=c(1,7),
		xlab = 'Total transmission generations G',
		ylab = 'Prob. G or more generations',yaxt='n')
	axis(side=2, at=c(10^(-4:-1),0.5), labels=freqlabels, par(mgp = c(3,.5,0)),las=1)
	axis(side=2, at=smallticksp,tck=-0.016,labels=rep('',length(smallticksp)))
	axis(side=4, at=10^seq(-4,0), labels=rep('',5),tck=1,lty=2)
	axis(side=3, at=1:7,labels=rep('',7),tck=1,lty=2)

	loc <- par('usr')
	text(loc[1],exp(loc[4]),titles[i*2],xpd=TRUE,adj=c(2.5,-4),font=2)

	lines(pgg,lty=2,lwd=2)
	points(pgg1b,pch=19)
	lines(pgg1b,lty=3,lwd=2)
	points(pgg2b,pch=19)
	lines(pgg2b,lty=1,lwd=2)
}
