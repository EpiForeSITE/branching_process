rm(list=ls())

#p is the probability of y total transmission from x independent patients
#lp <- function(x,y,R,k) lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k)
#p <- function(x,y,R,k) exp(lp(x,y,R,k))

#probability of no transmissions from n initial cases:
#q0 <- function(n,R0,k0) p(n,0,R0,k0)

#probability of j total cases in an outbreak ending after exactly one generation of transmission:
#q1 <- function(n,j,R0,Rc,k0,kc,genSwitch=Inf) p(n,j-n,R0,k0)*p(j-n,0,ifelse(genSwitch>1,R0,Rc),ifelse(genSwitch>1,k0,kc))

#probability of j total cases in an outbreak ending after exactly two generations of transmission:
#q2 <- function(n,j,R0,Rc,k0,kc,genSwitch=Inf){
#	x <- 1:(j-n-1)
#	sum(p(n,x,R0,k0) * p(x,j-n-x,ifelse(genSwitch>1,R0,Rc),ifelse(genSwitch>1,k0,kc)) * p(j-n-x,0,ifelse(genSwitch>2,R0,Rc),ifelse(genSwitch>2,k0,kc)))
#}

#probability of j total cases in an outbreak ending after exactly three generations of transmission:
#q3 <- function(n,j,R0,Rc,k0,kc,genSwitch=Inf){
#	f <- function(x) {
#		y <- 1:(j-n-x-1)
#		p(n,x,R0,k0) * sum(p(x,y,ifelse(genSwitch>1,R0,Rc),ifelse(genSwitch>1,k0,kc)) * p(y,j-n-x-y,ifelse(genSwitch>2,R0,Rc),ifelse(genSwitch>2,k0,kc)) * p(j-n-x-y,0,ifelse(genSwitch>3,R0,Rc),ifelse(genSwitch>3,k0,kc)))
#	}
#	sum(Vectorize(f)(1:(j-n-2)))
#}

#probability of j total cases in an outbreak (any number of generations)
#qAny <- function(n,j,R,k) p(j,j-n,R,k)*n/j
#qAnyGen1 <- function(n,j,R0,k0,Rc,kc=k0) ifelse(j==n,p(n,0,R0,k0),sum(p(n,1:(j-n),R0,k0)*qAny(1:(j-n),j-n,Rc,kc)))
#qAnyGen2 <- function(n,j,R0,k0,Rc,kc=k0){
#	f <- function(x) p(n,x,R0,k0) * sum(p(x,1:(j-n-x),R0,k0) * qAny(1:(j-n-x),j-n-x,Rc,kc))
#	ifelse(j==n,p(n,0,R0,k0),ifelse(j == n+1,p(n,1,R0,k0) * p(1,0,R0,k0), 
#		 p(n,j-n,R0,k0) * p(j-n,0,R0,k0) + sum(Vectorize(f)(1:(j-n-1)))))
#}

#probability of a size-j outbreak that goes entirely undetected
#qAnyMissed <- function(n,j,R,k,probDetect) qAny(n,j,R,k)*(1-probDetect)^j

#Assuming that an undetected outbreak would not switch to Rc/kc (no control if not detected)
#probUndetectedCluster <- function(R0,k0,pDet) sum(Vectorize(qAnyMissed)(1,1:40,R0,k0,pDet))

pFinalSizeUndetected <- function(n,j,R,k,pDet) pFinalSize(n,j,R,k)*(1-pDet)^j
pUndetectedCluster <- function(R,k,pDet){
  maxj <- ceiling(log(1e-7)/log(1-pDet))
  sum(Vectorize(pFinalSizeUndetected)(1,1:maxj,R,k,pDet))
}


LL <- function(numUndetectedClusters,R,k,pDet){
  pDetCluster = 1-(1-pDet)^clusterSize
  sum(clusterNum * log(Vectorize(pFinalSizeAndGen)(clusterGen,1,clusterSize,R,k) * pDetCluster)) + 
        numUndetectedClusters * log(pUndetectedCluster(R,k,pDet))
}

LLswitch1 <- function(numUndetectedClusters,R0,k0,Rc,kc,pDet){
  pDetCluster = 1-(1-pDet)^clusterSize
  sum(clusterNum * log(Vectorize(pFinalSizeAndGenSwitch1)(clusterGen,1,clusterSize,R0,k0,Rc,kc) * pDetCluster)) + 
    numUndetectedClusters * log(pUndetectedCluster(R0,k0,pDet))
}

LLswitch2 <- function(numUndetectedClusters,R0,k0,Rc,kc,pDet){
  pDetCluster = 1-(1-pDet)^clusterSize
  sum(clusterNum * log(Vectorize(pFinalSizeAndGenSwitch2)(clusterGen,1,clusterSize,R0,k0,Rc,kc) * pDetCluster)) + 
    numUndetectedClusters * log(pUndetectedCluster(R0,k0,pDet))
}

#LLtest <- function(numUndetectedClusters,R0,Rc,k0,kc,pDet,genSwitch=Inf){
# LL <- 0
#  for(i in 1:length(clusterSize)){
#    pDetCluster = 1-(1-pDet)^clusterSize[i]
#    if(clusterGen[i]==0){
#      LL <- LL + clusterNum[i] * log(q0(1,R0,k0) * pDetCluster)
#    }else if(clusterGen[i]==1){
#      LL <- LL + clusterNum[i] * log(q1(1,clusterSize[i],R0,Rc,k0,kc,genSwitch) * pDetCluster)
#    }else if(clusterGen[i]==2){
#      LL <- LL + clusterNum[i] * log(q2(1,clusterSize[i],R0,Rc,k0,kc,genSwitch) * pDetCluster)
#    }else if(clusterGen[i]==3){
#      LL <- LL + clusterNum[i] * log(q3(1,clusterSize[i],R0,Rc,k0,kc,genSwitch) * pDetCluster)
#    }
#  }
#  LL + numUndetectedClusters * log(probUndetectedCluster(R0,k0,pDet))
#}


getAIC <- function(LL, params) 2*params - 2*LL

#Toth31
clusterSize <- c( 1,2,3,7,186)
clusterGen <-  c( 0,1,1,3,  3)
clusterR0 <-   c( 0,1,2,1, 31)
clusterNum <-  c(23,4,2,1,  1)

Table4 <- {}
undetClusterFrac <- c(0.5,0.75)

for(ucf in undetClusterFrac){

	numUndetClus <- sum(clusterNum) * ucf/(1-ucf)	

	#optRk <- optim(f = function(x) ifelse(all(x>0), -LLtest(numUndetClus,x[1],x[1],x[2],x[2],x[3]), Inf), par = c(.8,.2,.5))
	optRk <- optim(f = function(x) ifelse(all(x>0), -LL(numUndetClus,x[1],x[2],x[3]), Inf), par = c(.8,.2,.5))
	R_mle <- optRk$par[1]
	k_mle <- optRk$par[2]
	pDet_mle <- optRk$par[3]
	LLhatRk <- -optRk$value
	AIC_Rk <- getAIC(LLhatRk,3)

	#optR0Rck1 <- optim(f = function(x) ifelse(all(x>0), -LLtest(numUndetClus,x[1],x[2],x[3],x[3],x[4],genSwitch=1), Inf), par = c(.8,.8,.2,.5))
	optR0Rck1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(numUndetClus,x[1],x[3],x[2],x[3],x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mle1 <- optR0Rck1$par[1]
	Rc_mle1 <- optR0Rck1$par[2]
	k_mle1 <- optR0Rck1$par[3]
	pDet_mle1 <- optR0Rck1$par[4]
	LLhatR0Rck1 <- -optR0Rck1$value
	AIC_R0Rck1 <- getAIC(LLhatR0Rck1,4)

	#optR0Rck2 <- optim(f = function(x) ifelse(all(x>0), -LLtest(numUndetClus,x[1],x[2],x[3],x[3],x[4],genSwitch=2), Inf), par = c(.8,.8,.2,.5))
	optR0Rck2 <- optim(f = function(x) ifelse(all(x>0), -LLswitch2(numUndetClus,x[1],x[3],x[2],x[3],x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mle2 <- optR0Rck2$par[1]
	Rc_mle2 <- optR0Rck2$par[2]
	k_mle2 <- optR0Rck2$par[3]
	pDet_mle2 <- optR0Rck2$par[4]
	LLhatR0Rck2 <- -optR0Rck2$value
	AIC_R0Rck2 <- getAIC(LLhatR0Rck2,4)

	#optR0RckGeomPost1 <- optim(f = function(x) ifelse(all(x>0), -LLtest(numUndetClus,x[1],x[2],x[3],1,x[4],genSwitch=1), Inf), par = c(.8,.8,.2,.5))
	optR0RckGeomPost1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(numUndetClus,x[1],x[3],x[2],1,x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mleGeomPost1 <- optR0RckGeomPost1$par[1]
	Rc_mleGeomPost1 <- optR0RckGeomPost1$par[2]
	k_mleGeomPost1 <- optR0RckGeomPost1$par[3]
	pDet_mleGeomPost1 <- optR0RckGeomPost1$par[4]
	LL_R0RckGeomPost1 <- -optR0RckGeomPost1$value
	AIC_R0RckGeomPost1 <- getAIC(LL_R0RckGeomPost1,4)

	#optR0RckGeomPost2 <- optim(f = function(x) ifelse(all(x>0), -LLtest(numUndetClus,x[1],x[2],x[3],1,x[4],genSwitch=2), Inf), par = c(.8,.8,.2,.5))
	optR0RckGeomPost2 <- optim(f = function(x) ifelse(all(x>0), -LLswitch2(numUndetClus,x[1],x[3],x[2],1,x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mleGeomPost2 <- optR0RckGeomPost2$par[1]
	Rc_mleGeomPost2 <- optR0RckGeomPost2$par[2]
	k_mleGeomPost2 <- optR0RckGeomPost2$par[3]
	pDet_mleGeomPost2 <- optR0RckGeomPost2$par[4]
	LL_R0RckGeomPost2 <- -optR0RckGeomPost2$value
	AIC_R0RckGeomPost2 <- getAIC(LL_R0RckGeomPost2,4)

	R0 <- c(R_mle,R0_mle1,R0_mleGeomPost1,R0_mle2,R0_mleGeomPost2)
	k0 <- c(k_mle,k_mle1 ,k_mleGeomPost1 ,k_mle2 ,k_mleGeomPost2)
	Rc <- c(R_mle,Rc_mle1,Rc_mleGeomPost1,Rc_mle2,Rc_mleGeomPost2)
	kc <- c(k_mle,k_mle1 ,1              ,k_mle2 ,1)
	u <- c(pDet_mle,pDet_mle1,pDet_mleGeomPost1,pDet_mle2,pDet_mleGeomPost2)
	optLL <- c(LLhatRk, LLhatR0Rck1, LL_R0RckGeomPost1, LLhatR0Rck2, LL_R0RckGeomPost2)
	AIC <- c(AIC_Rk,AIC_R0Rck1,AIC_R0RckGeomPost1,AIC_R0Rck2,AIC_R0RckGeomPost2)

	Table4 <- rbind(Table4,cbind(undetFrac=ucf,R0,k0,Rc,kc,u,optLL,AIC))
}
rownames(Table4) <- rep(c('Model0','Model1a','Model1b','Model2a','Model2b'),2)
print(Table4)
