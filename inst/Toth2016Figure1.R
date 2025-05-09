plotPair <- function(R,k,R01,k01,Rc1,R02,k02,Rc2,titles){
  dist0 <- pFinalSize(1,1:20000,R,k)
  dist1b <- Vectorize(pFinalSizeSwitch1)(1,1:1000,R01,k01,Rc1,1)
  dist2b <- Vectorize(pFinalSizeSwitch2)(1,1:700,R02,k02,Rc2,1)
  
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
  text(exp(loc[1]),exp(loc[4]),titles[1],xpd=TRUE,adj=c(2.5,-4),font=2)
  
  lines(1-cumsum(dist1b),lty=3,lwd=2)
  lines(1-cumsum(dist2b),lty=1,lwd=2)
  
  pgl <- pGen(gMax=10,R,k)
  pgl1b <- pGenSwitch1(gMax=10,R01,k01,Rc1,1)
  pgl2b <- pGenSwitch2(gMax=10,R02,k02,Rc2,1)
  
  #pgg[h] is the probability of observing at least h generations of transmission		
  pgg <- 1-pgl
  pgg1b <- 1-pgl1b
  pgg2b <- 1-pgl2b
  
  plot(pgg,log='y',pch=19,ylim=c(1e-4,1),xlim=c(1,7),
       xlab = 'Total transmission generations G',
       ylab = 'Prob. G or more generations',yaxt='n')
  axis(side=2, at=c(10^(-4:-1),0.5), labels=freqlabels, par(mgp = c(3,.5,0)),las=1)
  axis(side=2, at=smallticksp,tck=-0.016,labels=rep('',length(smallticksp)))
  axis(side=4, at=10^seq(-4,0), labels=rep('',5),tck=1,lty=2)
  axis(side=3, at=1:7,labels=rep('',7),tck=1,lty=2)
  
  loc <- par('usr')
  text(loc[1],exp(loc[4]),titles[2],xpd=TRUE,adj=c(2.5,-4),font=2)
  
  lines(pgg,lty=2,lwd=2)
  points(pgg1b,pch=19)
  lines(pgg1b,lty=3,lwd=2)
  points(pgg2b,pch=19)
  lines(pgg2b,lty=1,lwd=2)
  
}

dev.new(width=7,height=10,unit="in",noRStudioGD = TRUE)

oldpar <- par(mfrow=c(3,2))

plotPair(R = 0.8652294,
         k = 0.03539390,
         R01 = 5.4998656,
         k01 = 0.06113228, 
         Rc1 = 0.14344697,
         R02 = 2.1945074,
         k02 = 0.07574887,
         Rc2 = 0.06033546,
         titles = c("A","B"))

plotPair(R = 0.7631224, 
         k = 0.02814292,
         R01 = 2.7178632,
         k01 = 0.03189325,
         Rc1 = 0.14414662,
         R02 = 1.4687318,
         k02 = 0.04242669,
         Rc2 = 0.06326285,
         titles = c("C","D"))

plotPair(R = 0.6224364, 
         k = 0.02245021,
         R01 = 1.4025546,
         k01 = 0.01921805,
         Rc1 = 0.14568040,
         R02 = 1.4265324,
         k02 = 0.02554593,
         Rc2 = 0.06494508,
         titles = c("E","F"))

par(oldpar)