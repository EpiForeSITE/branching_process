init <- 1

Rvals <- seq(0,1,len=100)
log10kvals <- seq(-2,2,len=100)
kvals <- 10^log10kvals 
kats <- c(.01,.1,1,10,100)
klabs <- c('0.01','0.1','1','10','100')

pAtLeastOne <- function(R,k) 1 - (dnbinom(0,mu=R,size=k))^init
pAtLeastOneVals <- outer(Rvals,kvals,pAtLeastOne)

oldpar <- par(mfrow=c(2,2))

contour(Rvals,kvals,pAtLeastOneVals,log='y',axes=FALSE,
        main=expression(paste(''>=1,' transmission')),xlab='R',ylab='k')
axis(1)
axis(2, at=kats, labels=klabs,las=2)

pAtLeast5 <- function(R,k) 1-sum(pFinalSize(init,init:(init+4),R,k))
pAtLeast5 <- Vectorize(pAtLeast5)
pAtLeast5Vals <- outer(Rvals,kvals,pAtLeast5)

contour(Rvals,kvals,pAtLeast5Vals,log='y',axes=FALSE,
        main=expression(paste(''>=5,' transmissions')),xlab='R',ylab='k',
	levels = c(1e-4,1e-3,1e-2,0.03,0.05,0.1,0.2,0.3,0.4))
axis(1)
axis(2, at=kats, labels=klabs, las=2)

pAtLeast10 <- function(R,k) 1-sum(pFinalSize(init,init:(init+9),R,k))
pAtLeast10 <- Vectorize(pAtLeast10)
pAtLeast10Vals <- outer(Rvals,kvals,pAtLeast10)

contour(Rvals,kvals,pAtLeast10Vals,log='y',axes=FALSE,
        main=expression(paste(''>=10,' transmissions')),xlab='R',ylab='k',
        levels = c(1e-4,1e-3,1e-2,0.02,0.05,0.1,0.2))
axis(1)
axis(2, at=kats, labels=klabs, las=2)

pAtLeast100 <- function(R,k) 1-sum(pFinalSize(init,init:(init+99),R,k))
pAtLeast100 <- Vectorize(pAtLeast100)
pAtLeast100Vals <- outer(Rvals,kvals,pAtLeast100)

contour(Rvals,kvals,pAtLeast100Vals,log='y',axes=FALSE,
        main=expression(paste(''>=100,' transmissions')),xlab='R',ylab='k',
        levels = c(1e-4,1e-3,1e-2,.05))
axis(1)
axis(2, at=kats, labels=klabs, las=2)

par(oldpar)