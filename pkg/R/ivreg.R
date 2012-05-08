spatial.ivreg<-function(y, Zmat, Hmat, HAC=FALSE, distance=distance,type=c("Epanechnikov","Triangular","Bisquare","Parzen", "QS","TH"), bandwidth=bandwidth){
	df <- nrow(Zmat) - ncol(Zmat)
	HH <- crossprod(Hmat,Hmat)
	Hye <- crossprod(Hmat, Zmat)
	bz <- solve(HH,Hye)
	Zp <- Hmat %*% bz
	ZpZp <- crossprod(Zp)
	ZpZpi <- solve(ZpZp)
	Zpy <- crossprod(Zp,y)
	delta <- crossprod(ZpZpi,Zpy)
	yp <- Zmat %*% delta
	e <- y - yp
if(HAC){
	n<-nrow(Zmat)
		Ker<-vector(mode='list',length=n)
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
								
if(is.null(attributes(distance)$GeoDa$dist)){
	Ker<-lapply(distance$weights,ker.fun, bandwidth=bandwidth)
	Kern<-nb2listw(distance$neighbours,style="B", glist=Ker)
	} 
else{
	Ker<-lapply(attributes(distance)$GeoDa$dist,ker.fun, bandwidth=bandwidth)
	Kern<-nb2listw(distance,style="B", glist=Ker)
	} 
He<-matrix(,dim(Hmat)[1],dim(Hmat)[2])
for (i in 1:dim(Hmat)[2]) He[,i]<- Hmat[,i] * e

KHpe<-lag.listw(Kern,He) +He
KHeHe<-(t(He) %*% KHpe)
HHp<-solve(HH)
ZpH<-crossprod(Zmat,Hmat)
fp<-ZpZpi%*%ZpH%*%HHp
vardelta<- (fp%*%KHeHe%*%t(fp))
s2 <- crossprod(e,e) / df
	}
else{
	s2 <- crossprod(e,e) / df
	vardelta <- ZpZpi * as.numeric(s2)
	}
	result <- list(coefficients=delta, var=vardelta, s2=s2, residuals=as.numeric(e), yhat=yp)
	return(result)
}