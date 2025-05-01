clusterSize <- c( 1,2,3,7,186)
clusterGen <-  c( 0,1,1,3,  3)
clusterNum <-  c(23,4,2,1,  1)

LLnoSwitch <- function(R,k){
  sum(clusterNum * log(Vectorize(pFinalSizeAndGen)(clusterGen,1,clusterSize,R,k)))
}

LLswitch1 <- function(R0,Rc,k0,kc){
  sum(clusterNum * log(Vectorize(pFinalSizeAndGenSwitch1)(clusterGen,1,clusterSize,R0,k0,Rc,kc)))
}

LLswitch2 <- function(R0,Rc,k0,kc){
  sum(clusterNum * log(Vectorize(pFinalSizeAndGenSwitch2)(clusterGen,1,clusterSize,R0,k0,Rc,kc)))
}

AIC <- function(LL, params) 2*params - 2*LL

optRk <- optim(f = function(x) ifelse(all(x>0), -LLnoSwitch(x[1],x[2]), Inf), par = c(.8,.2))
R_mle <- optRk$par[1]
k_mle <- optRk$par[2]
LL_Rk <- -optRk$value
AIC_Rk <- AIC(LL_Rk,2)

optR0Rck1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(x[1],x[2],x[3],x[3]), Inf), par = c(.8,.8,.2))
R0_mle1 <- optR0Rck1$par[1]
Rc_mle1 <- optR0Rck1$par[2]
k_mle1 <- optR0Rck1$par[3]
LL_R0Rck1 <- -optR0Rck1$value
AIC_R0Rck1 <- AIC(LL_R0Rck1,3)

optR0Rck2 <- optim(f = function(x) ifelse(all(x>0), -LLswitch2(x[1],x[2],x[3],x[3]), Inf), par = c(.8,.8,.2))
R0_mle2 <- optR0Rck2$par[1]
Rc_mle2 <- optR0Rck2$par[2]
k_mle2 <- optR0Rck2$par[3]
LL_R0Rck2 <- -optR0Rck2$value
AIC_R0Rck2 <- AIC(LL_R0Rck2,3)

optR0RckPoisPost1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(x[1],x[2],x[3],Inf), Inf), par = c(.8,.8,.2))
R0_mlePoisPost1 <- optR0RckPoisPost1$par[1]
Rc_mlePoisPost1 <- optR0RckPoisPost1$par[2]
k_mlePoisPost1 <- optR0RckPoisPost1$par[3]
LL_R0RckPoisPost1 <- -optR0RckPoisPost1$value
AIC_R0RckPoisPost1 <- AIC(LL_R0RckPoisPost1,3)

optR0RckPoisPost2 <- optim(f = function(x) ifelse(all(x>0), -LLswitch2(x[1],x[2],x[3],Inf), Inf), par = c(.8,.8,.2))
R0_mlePoisPost2 <- optR0RckPoisPost2$par[1]
Rc_mlePoisPost2 <- optR0RckPoisPost2$par[2]
k_mlePoisPost2 <- optR0RckPoisPost2$par[3]
LL_R0RckPoisPost2 <- -optR0RckPoisPost2$value
AIC_R0RckPoisPost2 <- AIC(LL_R0RckPoisPost2,3)

optR0RckGeomPost1 <- optim(f = function(x) ifelse(all(x>0), -LLswitch1(x[1],x[2],x[3],1), Inf), par = c(.8,.8,.2))
R0_mleGeomPost1 <- optR0RckGeomPost1$par[1]
Rc_mleGeomPost1 <- optR0RckGeomPost1$par[2]
k_mleGeomPost1 <- optR0RckGeomPost1$par[3]
LL_R0RckGeomPost1 <- -optR0RckGeomPost1$value
AIC_R0RckGeomPost1 <- AIC(LL_R0RckGeomPost1,3)

optR0RckGeomPost2 <- optim(f = function(x) ifelse(all(x>0), -LLswitch2(x[1],x[2],x[3],1), Inf), par = c(.8,.8,.2))
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
