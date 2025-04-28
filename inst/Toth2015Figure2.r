R0vals <- c(3,3,0.3,0.1)
Rcvals <- c(0.3,0.1,0.3,0.1)

kvals <- c(10,1,0.1)

pgeTally <- function(n,maxX,R0,Rc,k){
  if(R0==Rc){
    out <- 1-cumsum(pFinalSize(n,n:(n+maxX-1),Rc,k))
  }else{
    out <- 1-cumsum(Vectorize(pFinalSizeSwitch1,'j')(n,n:(n+maxX-1),R0,k,Rc,k))
  }
  out
}

pLevel <- function(n,maxX,R0,Rc,k,pLevel){
	pges <- pgeTally(n,maxX,R0,Rc,k)
	max(which(pges>pLevel))
}

oneIn100outbreak <- function(n,R0,Rc,k) pLevel(n,500,R0,Rc,k,0.01)
oneIn100outbreak <- Vectorize(oneIn100outbreak)

outA <- matrix(0,length(R0vals),length(kvals))
for(i in 1:length(R0vals)){
	outA[i,] <- oneIn100outbreak(1,R0vals[i],Rcvals[i],kvals)
}

oneIn10000outbreak <- function(n,R0,Rc,k) pLevel(n,500,R0,Rc,k,0.0001)
oneIn10000outbreak <- Vectorize(oneIn10000outbreak)

outB <- matrix(0,length(R0vals),length(kvals))
for(i in 1:length(R0vals)){
	outB[i,] <- oneIn10000outbreak(1,R0vals[i],Rcvals[i],kvals)
}

barlab <- as.expression(c(bquote(paste(italic(R)[0],"="*.(R0vals[1])*", ",italic(R)[c],"="*.(Rcvals[1]))),
            bquote(paste(italic(R)[0],"="*.(R0vals[2])*", ",italic(R)[c],"="*.(Rcvals[2]))),
            bquote(paste(italic(R),"="*.(R0vals[3]))),
            bquote(paste(italic(R),"="*.(R0vals[4])))))

oldpar <- par(mfrow=c(1,2))
barplot(t(outA),beside=T,names.arg=barlab,
        col=c("black","white","gray"),ylim=c(0,250),
	  ylab = 'one-in-100 outbreak size',yaxt='n')

axis(side=2,at=c(0,50,100,150,200,250))
legend("topright",  fill=c("black","white","gray"),
      legend=c(bquote(paste(italic(k),"="*.(kvals[1]))),
               bquote(paste(italic(k),"="*.(kvals[2]))),
               bquote(paste(italic(k),"="*.(kvals[3])))), bty="n")

barplot(t(outB),beside=T,names.arg=barlab,
        col=c("black","white","gray"),ylim=c(0,250),
	  ylab = 'one-in-10,000 outbreak size')

par(oldpar)