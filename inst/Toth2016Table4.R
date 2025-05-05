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

getAIC <- function(LL, params) 2*params - 2*LL

clusterSize <- c( 1,2,3,7,186)
clusterGen <-  c( 0,1,1,3,  3)
clusterR0 <-   c( 0,1,2,1, 31)
clusterNum <-  c(23,4,2,1,  1)

Table4 <- {}
undetClusterFrac <- c(0.5,0.75)

for(ucf in undetClusterFrac){

	numUndetClus <- sum(clusterNum) * ucf/(1-ucf)	

	optRk <- optim(f = function(x) ifelse(all(x>0), -LL(numUndetClus,x[1],x[2],x[3]), Inf), par = c(.8,.2,.5))
	R_mle <- optRk$par[1]
	k_mle <- optRk$par[2]
	pDet_mle <- optRk$par[3]
	LLhatRk <- -optRk$value
	AIC_Rk <- getAIC(LLhatRk,3)

	optR0Rck1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(numUndetClus,x[1],x[3],x[2],x[3],x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mle1 <- optR0Rck1$par[1]
	Rc_mle1 <- optR0Rck1$par[2]
	k_mle1 <- optR0Rck1$par[3]
	pDet_mle1 <- optR0Rck1$par[4]
	LLhatR0Rck1 <- -optR0Rck1$value
	AIC_R0Rck1 <- getAIC(LLhatR0Rck1,4)

	optR0Rck2 <- optim(f = function(x) ifelse(all(x>0), -LLswitch2(numUndetClus,x[1],x[3],x[2],x[3],x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mle2 <- optR0Rck2$par[1]
	Rc_mle2 <- optR0Rck2$par[2]
	k_mle2 <- optR0Rck2$par[3]
	pDet_mle2 <- optR0Rck2$par[4]
	LLhatR0Rck2 <- -optR0Rck2$value
	AIC_R0Rck2 <- getAIC(LLhatR0Rck2,4)

	optR0RckGeomPost1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(numUndetClus,x[1],x[3],x[2],1,x[4]), Inf), par = c(.8,.8,.2,.5))
	R0_mleGeomPost1 <- optR0RckGeomPost1$par[1]
	Rc_mleGeomPost1 <- optR0RckGeomPost1$par[2]
	k_mleGeomPost1 <- optR0RckGeomPost1$par[3]
	pDet_mleGeomPost1 <- optR0RckGeomPost1$par[4]
	LL_R0RckGeomPost1 <- -optR0RckGeomPost1$value
	AIC_R0RckGeomPost1 <- getAIC(LL_R0RckGeomPost1,4)

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
