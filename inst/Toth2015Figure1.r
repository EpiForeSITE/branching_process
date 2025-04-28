R0vals <- c(3,3,NA,NA)
Rcvals <- c(.3,.1,.3,.1)

kvals <- c(1,.1,10)
ltys <- c(1,2,4)

xrange <- 1:1000
labelticks <- c(1,10,100,1000,10000)
smallticksX <- c(2:9,seq(20,90,10),seq(200,900,100),seq(2000,9000,1000))
smallticksp <- c((2:9)*1e-4,(2:9)*1e-3,(2:9)*1e-2,(2:9)*1e-1)
freqlabels <- c('0.01','0.1','1','10','50')

titles <- c('A','B','C','D')

oldpar <- par(mfrow=c(2,2))

for(j in 1:length(Rcvals)){
	R0 <- R0vals[j]; Rc <- Rcvals[j]

	if(is.na(R0)){
		ttl <- bquote(paste(italic(R),' = ') ~.(Rc))
	}else{
	  ttl <- bquote(paste(italic(R)[0],' = ') ~.(R0) *paste(', ',italic(R)[c],' = ') ~.(Rc))
	}

	for(jj in 1:length(kvals)){
		k <- kvals[jj]
		if(!is.na(R0)){
		  pge_plot <- 1-cumsum(Vectorize(pFinalSizeSwitch1,'j')(1,xrange,R0,k,Rc,k))
		}else{
			pge_plot <- 1-cumsum(pFinalSize(1,xrange,Rc,k))	
		}
	
		i <- which(pge_plot > 1e-6)

		if(jj==1){
			plot(xrange[i],pge_plot[i],type='l',log='xy',xlim=c(1,1000),ylim=c(1e-4,1),lwd=2,
				xlab='Total no. transmissions X',
				ylab=bquote(paste('% Probability',"" >= "",'X transmissions')),xaxt='n',yaxt='n',
				bty='l',main=ttl)
			
			axis(side=1, at=labelticks, labels=c('1','10','100','1,000','10,000'))
			axis(side=1, at=smallticksX,tck=-0.016,labels=rep('',length(smallticksX)))
			axis(side=2, at=c(10^(-4:-1),0.5), labels=freqlabels, par(mgp = c(3,.5,0)),las=1)
			axis(side=2, at=smallticksp,tck=-0.016,labels=rep('',length(smallticksp)))
			loc <- par('usr')
			text(exp(loc[1]),exp(loc[4]),titles[j],xpd=TRUE,adj=c(2.5,-4),font=2)
		}else{
			lines(xrange[i],pge_plot[i],lty=ltys[jj],lwd=2,col='black')
		}
	}
}
par(oldpar)