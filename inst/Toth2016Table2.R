rm(list=ls())

#p is the probability of y total transmission from x independent patients
lp <- function(x,y,R,k) lgamma(k*x+y)-lgamma(k*x)-lgamma(y+1)+y*log(R/k)-(k*x+y)*log(1+R/k)
p <- function(x,y,R,k) exp(lp(x,y,R,k))

lpPois <- function(x,y,R) y*(log(x)+log(R))-x*R-lfactorial(y)
pPois <- function(x,y,R) exp(lpPois(x,y,R))

#probability of no transmissions from n initial cases:
#q0 <- function(n,R0,k0) p(n,0,R0,k0)
q0 <- function(n,R0,k0) pNextGenSize(n,0,R0,k0)

#probability of j total cases in an outbreak ending after exactly one generation of transmission:
q1 <- function(n,j,R0,Rc,k0,kc,genSwitch=Inf) p(n,j-n,R0,k0)*p(j-n,0,ifelse(genSwitch>1,R0,Rc),ifelse(genSwitch>1,k0,kc))
q1PoisPost <- function(n,j,R0,Rc,k0,genSwitch=Inf) p(n,j-n,R0,k0)*ifelse(genSwitch>1,p(j-n,0,R0,k0),pPois(j-n,0,Rc))

#probability of j total cases in an outbreak ending after exactly two generations of transmission:
q2 <- function(n,j,R0,Rc,k0,kc,genSwitch=Inf){
	x = 1:(j-n-1)
	sum(p(n,x,R0,k0) * p(x,j-n-x,ifelse(genSwitch>1,R0,Rc),ifelse(genSwitch>1,k0,kc)) * p(j-n-x,0,ifelse(genSwitch>2,R0,Rc),ifelse(genSwitch>2,k0,kc)))
}
q2PoisPost <- function(n,j,R0,Rc,k0,genSwitch=Inf){
	x = 1:(j-n-1)
	if(genSwitch>1){
		v1 = p(x,j-n-x,R0,k0)
	}else{
		v1 = pPois(x,j-n-x,Rc)
	}
	if(genSwitch>2){
		v2 = p(j-n-x,0,R0,k0)
	}else{
		v2 = pPois(j-n-x,0,Rc)
	}
	sum(p(n,x,R0,k0) * v1 * v2)
}

#probability of j total cases in an outbreak ending after exactly three generations of transmission:
q3 <- function(n,j,R0,Rc,k0,kc,genSwitch=Inf){
	f = function(x) {
		y = 1:(j-n-x-1)
		p(n,x,R0,k0) * sum(p(x,y,ifelse(genSwitch>1,R0,Rc),ifelse(genSwitch>1,k0,kc)) * p(y,j-n-x-y,ifelse(genSwitch>2,R0,Rc),ifelse(genSwitch>2,k0,kc)) * p(j-n-x-y,0,ifelse(genSwitch>3,R0,Rc),ifelse(genSwitch>3,k0,kc)))
	}
	sum(Vectorize(f)(1:(j-n-2)))
}
q3PoisPost <- function(n,j,R0,Rc,k0,genSwitch=Inf){
	f = function(x) {
		y = 1:(j-n-x-1)
		if(genSwitch>1){
			v1 = p(x,y,R0,k0)
		}else{
			v1 = pPois(x,y,Rc)
		}
		if(genSwitch>2){
			v2 = p(y,j-n-x-y,R0,k0)
		}else{
			v2 = pPois(y,j-n-x-y,Rc)
		}
		if(genSwitch>3){
			v3 = p(j-n-x-y,0,R0,k0)
		}else{
			v3 = pPois(j-n-x-y,0,Rc)
		}
		
		p(n,x,R0,k0) * sum(v1 * v2 * v3)
	}
	sum(Vectorize(f)(1:(j-n-2)))
}
#Toth31
clusterSize <- c( 1,2,3,7,186)
clusterGen <-  c( 0,1,1,3,  3)
clusterR0 <-   c( 0,1,2,1, 31)
clusterNum <-  c(23,4,2,1,  1)

LLnoSwitch <- function(R,k){
	LL <- 0
	for(i in 1:length(clusterSize)){
	  LL <- LL + clusterNum[i] * log(pFinalSizeAndGen(clusterGen[i],1,clusterSize[i],R,k))
	}
	LL
}

LLswitch1 <- function(R0,Rc,k0,kc){
  LL <- 0
  for(i in 1:length(clusterSize)){
    LL <- LL + clusterNum[i] * log(pFinalSizeAndGenSwitch1(clusterGen[i],1,clusterSize[i],R0,k0,Rc,kc))
  }
  LL
}

LLnew <- function(R0,Rc,k0,kc,genSwitch=Inf){
  LL <- 0
  for(i in 1:length(clusterSize)){
    if(clusterGen[i]==0){
      LL <- LL + clusterNum[i] * log(q0(1,R0,k0))
    }else if(clusterGen[i]==1){
    	LL <- LL + clusterNum[i] * log(q1(1,clusterSize[i],R0,Rc,k0,kc,genSwitch))
    }else if(clusterGen[i]==2){
    	LL <- LL + clusterNum[i] * log(q2(1,clusterSize[i],R0,Rc,k0,kc,genSwitch))
    }else if(clusterGen[i]==3){
    	LL <- LL + clusterNum[i] * log(q3(1,clusterSize[i],R0,Rc,k0,kc,genSwitch))
    }
  }
  LL
}

LLpoisPost <- function(R0,Rc,k0,genSwitch=Inf){
	LL <- 0
	for(i in 1:length(clusterSize)){
		if(clusterGen[i]==0){
			LL <- LL + clusterNum[i] * log(q0(1,R0,k0))
		}else if(clusterGen[i]==1){
			LL <- LL + clusterNum[i] * log(q1PoisPost(1,clusterSize[i],R0,Rc,k0,genSwitch))
		}else if(clusterGen[i]==2){
			LL <- LL + clusterNum[i] * log(q2PoisPost(1,clusterSize[i],R0,Rc,k0,genSwitch))
		}else if(clusterGen[i]==3){
			LL <- LL + clusterNum[i] * log(q3PoisPost(1,clusterSize[i],R0,Rc,k0,genSwitch))
		}
	}
	LL
}

AIC <- function(LL, params) 2*params - 2*LL

#optRk <- optim(f = function(x) ifelse(all(x>0), -LLnew(x[1],x[1],x[2],x[2]), Inf), par = c(.8,.2))

optRk <- optim(f = function(x) ifelse(all(x>0), -LLnoSwitch(x[1],x[2]), Inf), par = c(.8,.2))

R_mle <- optRk$par[1]
k_mle <- optRk$par[2]
LL_Rk <- -optRk$value
AIC_Rk <- AIC(LL_Rk,2)

#optR0Rck1 <- optim(f = function(x) ifelse(all(x>0), -LLnew(x[1],x[2],x[3],x[3],genSwitch=1), Inf), par = c(.8,.8,.2))
optR0Rck1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(x[1],x[2],x[3],x[3]), Inf), par = c(.8,.8,.2))
R0_mle1 <- optR0Rck1$par[1]
Rc_mle1 <- optR0Rck1$par[2]
k_mle1 <- optR0Rck1$par[3]
LL_R0Rck1 <- -optR0Rck1$value
AIC_R0Rck1 <- AIC(LL_R0Rck1,3)

optR0Rck2 <- optim(f = function(x) ifelse(all(x>0), -LLnew(x[1],x[2],x[3],x[3],genSwitch=2), Inf), par = c(.8,.8,.2))
R0_mle2 <- optR0Rck2$par[1]
Rc_mle2 <- optR0Rck2$par[2]
k_mle2 <- optR0Rck2$par[3]
LL_R0Rck2 <- -optR0Rck2$value
AIC_R0Rck2 <- AIC(LL_R0Rck2,3)

optR0RckPoisPost1 <- optim(f = function(x) ifelse(all(x>0), -LLpoisPost(x[1],x[2],x[3],genSwitch=1), Inf), par = c(.8,.8,.2))
R0_mlePoisPost1 <- optR0RckPoisPost1$par[1]
Rc_mlePoisPost1 <- optR0RckPoisPost1$par[2]
k_mlePoisPost1 <- optR0RckPoisPost1$par[3]
LL_R0RckPoisPost1 <- -optR0RckPoisPost1$value
AIC_R0RckPoisPost1 <- AIC(LL_R0RckPoisPost1,3)

optR0RckPoisPost2 <- optim(f = function(x) ifelse(all(x>0), -LLpoisPost(x[1],x[2],x[3],genSwitch=2), Inf), par = c(.8,.8,.2))
R0_mlePoisPost2 <- optR0RckPoisPost2$par[1]
Rc_mlePoisPost2 <- optR0RckPoisPost2$par[2]
k_mlePoisPost2 <- optR0RckPoisPost2$par[3]
LL_R0RckPoisPost2 <- -optR0RckPoisPost2$value
AIC_R0RckPoisPost2 <- AIC(LL_R0RckPoisPost2,3)

optR0RckGeomPost1 <- optim(f = function(x) ifelse(all(x>0), -LLnew(x[1],x[2],x[3],1,genSwitch=1), Inf), par = c(.8,.8,.2))
R0_mleGeomPost1 <- optR0RckGeomPost1$par[1]
Rc_mleGeomPost1 <- optR0RckGeomPost1$par[2]
k_mleGeomPost1 <- optR0RckGeomPost1$par[3]
LL_R0RckGeomPost1 <- -optR0RckGeomPost1$value
AIC_R0RckGeomPost1 <- AIC(LL_R0RckGeomPost1,3)

optR0RckGeomPost2 <- optim(f = function(x) ifelse(all(x>0), -LLnew(x[1],x[2],x[3],1,genSwitch=2), Inf), par = c(.8,.8,.2))
R0_mleGeomPost2 <- optR0RckGeomPost2$par[1]
Rc_mleGeomPost2 <- optR0RckGeomPost2$par[2]
k_mleGeomPost2 <- optR0RckGeomPost2$par[3]
LL_R0RckGeomPost2 <- -optR0RckGeomPost2$value
AIC_R0RckGeomPost2 <- AIC(LL_R0RckGeomPost2,3)

R0 <- c(R_mle,R0_mle1,R0_mleGeomPost1,R0_mlePoisPost1,R0_mle2,R0_mleGeomPost2,R0_mlePoisPost2)
k0 <- c(k_mle,k_mle1 ,k_mleGeomPost1 ,k_mlePoisPost1 ,k_mle2 ,k_mleGeomPost2 ,k_mlePoisPost2)
Rc <- c(R_mle,Rc_mle1,Rc_mleGeomPost1,Rc_mlePoisPost1,Rc_mle2,Rc_mleGeomPost2,Rc_mlePoisPost2)
kc <- c(k_mle,k_mle1 ,1              ,Inf            ,k_mle2 ,1              ,Inf)
LL <- c(LL_Rk, LL_R0Rck1, LL_R0RckGeomPost1, LL_R0RckPoisPost1, LL_R0Rck2, LL_R0RckGeomPost2, LL_R0RckPoisPost2)
AIC <- c(AIC_Rk,AIC_R0Rck1,AIC_R0RckGeomPost1,AIC_R0RckPoisPost1,AIC_R0Rck2,AIC_R0RckGeomPost2,AIC_R0RckPoisPost2)

Table2 <- cbind(R0,k0,Rc,kc,LL,AIC)
rownames(Table2) <- c('Model0','Model1a','Model1b','Model1c','Model2a','Model2b','Model2c')

print(Table2)
