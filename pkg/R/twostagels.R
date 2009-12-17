tsls <- function(y,yend,X,Zinst,end, reg, yor, modified=FALSE, HAC=FALSE, distance=distance,  type=c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH"), bandwidth=bandwidth) {
### estimation engine that deals with the various cases
if(modified){
	H <- cbind(reg, Zinst)
	Z <- cbind(yend, X)
	HH <- crossprod(H,H)
	Hye <- crossprod(H,Z)
	bz <- solve(HH,Hye)
	Zp <- H %*% bz
	ZpZp <- crossprod(Zp,Z)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	delta <- crossprod(ZpZpi,Zpy)
	rownames(delta)<-c(colnames(yend), colnames(X))
	Z1<-cbind(end, reg)
	yp <- Z1 %*% delta
    e <- yor - yp
     result <- list(coefficients=delta, yhat=yp, residuals=e)
	}
else{
		K<-ifelse(colnames(X)[1] == "(Intercept)" || all(X[,1]==1), 2, 1)
	#print(K)
	#if x does not have an intercept, one should add it to the instruments anyway
if (K==1)	H <- cbind(1,X,Zinst)
else	H <- cbind(X,Zinst)
#print(H[1,])
	Z <- cbind(yend,X)
	df <- nrow(Z) - ncol(Z)
	HH <- crossprod(H,H)
	Hye <- crossprod(H,yend)
	bz <- solve(HH,Hye)
	yendp <- H %*% bz
	Zp <- cbind(yendp,X)
#	print(Zp[1,])
	ZpZp <- crossprod(Zp)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	delta <- crossprod(ZpZpi,Zpy)
	rownames(delta)<-c(colnames(yend), colnames(X))
	yp <- Z %*% delta
	e <- y - yp
if(HAC){
	n<-nrow(X)
	#print(n)
		Ker<-vector(mode='list',length=n)
	#	print(Ker)
	#	print(match.arg(type))
		ker.fun<-switch(match.arg(type), Triangular={
			triangular
			}, Epanechnikov={
				epan
				}, Bisquare = {
					bisq
					}, TH={
						th
						}, QS = {
							qs
							}, Parzen = {
								parzen
								})
#					print(bandwith)
#print(ker.fun)
#print(is.numeric(bandwith))
Ker<-lapply(distance$weights,ker.fun, bandwidth=bandwidth)
#print(Ker)
#		if (type=='Triangular') Ker<-lapply(dist$weights,triangular)
#		if (type=="Epanechnikov") Ker<-lapply(dist$weights,epan)
#		if (type=='Bisquare') Ker<-lapply(dist$weights,bisq)
#		else stop("Kernel not implemented yet")
He<-matrix(,dim(H)[1],dim(H)[2])
KHpe<-matrix(,dim(H)[1],dim(H)[2])
for (i in 1:dim(H)[2]) He[,i]<- H[,i] * e
for(j in 1:dim(H)[2]){
	for (i in 1:n) KHpe[i,j]<- sum(Ker[[i]]*He[distance$neigh[[i]],j])  + He[i,j]
	 }
KHeHe<-(t(He) %*% KHpe)
#KHeHe<-crossprod(He, KHpe)
HHp<-solve(HH)
ZpH<-crossprod(Z,H)
HpZ<-crossprod(H,Z)
fp<-ZpZpi%*%ZpH%*%HHp
vardelta<- (fp%*%KHeHe%*%t(fp))
s2 <- crossprod(e,e) / df
#sedelta<-sqrt(diag(vardelta))
#tdelta<-delta/sedelta
#pdelta<-pnorm(abs(tdelta), lower.tail=FALSE )*2
	}

else{
    	s2 <- crossprod(e,e) / df
	    vardelta <- ZpZpi * as.numeric(s2)
#	    sedelta <- sqrt(diag(vardelta))
#	    tdelta <- delta / sedelta
#	    pdelta <- pnorm(abs(tdelta),lower.tail=FALSE) * 2
}

	    result <- list(coefficients=delta,vcmat=vardelta,s2=s2,
	          residuals=as.numeric(e),yhat=yp)
}	          
	result
}
