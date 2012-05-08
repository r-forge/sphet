spreg<-function(formula, data=list(), listw, listw2=NULL, endog = NULL, instruments= NULL, initial.value=0.2, abs.tol=1e-20,rel.tol=1e-10,eps=1e-5, sarar=TRUE,  het = FALSE, verbose=FALSE, na.action = na.fail ){
         	
         	
# source("listw2dgCMatrix.R")
# source("sp_ivreg.R")
# source("gm_functions.R")
	
#extract model objects	
	mt<-terms(formula,data=data)
	mf<-lm(formula, data, na.action=na.action, method="model.frame")
	na.act<-attr(mf,'na.action')

# record call
cl<- match.call()


#generates x and y 
	y<-model.extract(mf,"response")
	x<-model.matrix(mt,mf)

#checks on teh dimensions of x and y 	
if (length(y)!=nrow(x)) 
	stop("x and y have different length")

#check that X and y does not have missing values	
if (any(is.na(y))) 
        stop("NAs in dependent variable")
if (any(is.na(x))) 
        stop("NAs in independent variable")

if (!is.null(endog) && is.null(instruments)) stop("No instruments specified")
	
#fix the dimensions of the problem
	n<-nrow(x)
	k<-ncol(x)	
	xcolnames<-colnames(x)

	K<-ifelse(xcolnames[1] == "(Intercept)" || all(x[ ,1]==1), 2, 1)
	
	
#check that W is an object of class listw or a Matrix 
if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
if(inherits(listw,"listw"))  Ws<-listw2dgCMatrix(listw)	
if(inherits(listw,"matrix"))  Ws<-Matrix(listw)	

#check on the dimensions of x and W	
if (nrow(x) != nrow(Ws))
	stop("Input data and weights have different dimension")

if (k > 1) {
        wx <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          Wx <- Ws %*% x[, i]
            wx[, (i - (K - 1))] <- as.matrix(Wx)
				         }
            wwx <- Ws %*% wx                    					         
    }


if(!is.null(listw2) && sarar) stop("listw2 can be specified only with sarar")

if(is.null(listw2)) {
	twow<-FALSE		
	Ws2 <- Ws
	}
	
else{

twow<-TRUE	

if(!inherits(listw2,c("listw", "Matrix", "matrix"))) stop("listw2 format unknown")
if(inherits(listw2,"listw"))  Ws2<-listw2dgCMatrix(listw2)	
if(inherits(listw2,"matrix"))  Ws2<-Matrix(listw2)	
	
if (k > 1) {
        w2x <- matrix(nrow = n, ncol = (k  - (K - 1)))
        for (i in K:k) {
          W2x <- Ws2 %*% x[, i]          
            w2x[, (i - (K - 1))] <- as.matrix(W2x)
				         }
		  w2wx <- Ws2 %*% wx                   	
          w2wwx <- Ws2 %*% wwx                    	          

    }
	}

 wy<-Ws %*% y	
 colnames(wy)<-"lambda"


if(sarar){
 
if (!is.null(endog)) {

        if (is.character(endog)) {
        	endognames <- endog  
            xend <- match(endog, colnames(data))
            endog <- data[, xend]
            endog <- as.matrix(endog) 
            colnames(endog) <- endognames
        }
        
endog <- as.matrix(endog)             
if(is.null(colnames(endog)))  endognames <- rep("endogenous", ncol(endog))  

         if (is.character(instruments)) {
            inst <- match(instruments, colnames(data))
            instruments <- data[, inst]
        }
            instruments <- as.matrix(instruments)            
            winst <- Ws %*% instruments
            wwinst<- Ws %*% winst

if(twow){
	w2i <- Ws2 %*% instruments 
	w2wi <- Ws2 %*% winst 
	w2wwi <- Ws2 %*% wwinst 	
	AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst), as.matrix(w2i), as.matrix(w2wi),as.matrix(w2wwi))        
}

else  AddH <- cbind(instruments, as.matrix(winst), as.matrix(wwinst))        

if (K==2) {
	if(twow) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), AddH)
	else  Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), AddH)
		}
else {
	if(twow) Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx), AddH)
	else  Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx), AddH)
	}
Zmat<- cbind(x, endog, as.matrix(wy))            
colnames(Zmat) <- c(colnames(x), endognames, colnames(wy))               
}
 else {
	
	Zmat<- cbind(x, as.matrix(wy))                    

if (K==2){
	if(twow) Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx))
	else  Hmat <- cbind(x, as.matrix(wx), as.matrix(wwx)) 
}

else {
	if(twow) Hmat <- cbind(1,x, as.matrix(wx), as.matrix(wwx), as.matrix(w2x), as.matrix(w2wx), as.matrix(w2wwx))
	else  Hmat <- cbind(1, x, as.matrix(wx), as.matrix(wwx)) 

	}
    }
    }


else{
if (!is.null(endog)) {
	
        if (is.character(endog)) {
        	endognames<- endog        	
            xend <- match(endog, colnames(data))
            endog <- data[, xend]
            endog <- as.matrix(endog)   
            colnames(endog) <- endognames
        }

endog <- as.matrix(endog)                   
if(is.null(colnames(endog)))  endognames <- rep("endogenous", ncol(endog)) 


	 if (is.character(instruments)) {
            inst <- match(instruments, colnames(data))
            instruments <- data[, inst]
            instruments <- as.matrix(instruments)   
           
        }
 AddH<- cbind(instruments)
 Zmat<- cbind(x, endog)            
colnames(Zmat) <- c(colnames(x), endognames) 

if (K==2) Hmat<-cbind(x,wx, AddH) 
else Hmat<-cbind(1, x, wx, AddH)

 }
 
else {
if (K==2) Hmat<-cbind(x,wx)
else Hmat<-cbind(1, x,wx)	
	Zmat<- x
	}

}
	

firststep<-spatial.ivreg(y =y , Zmat = Zmat, Hmat = Hmat)
ubase<-residuals(firststep)
 # print(firststep$coefficients)

if (initial.value=="SAR"){
		Wubase<-Ws2 %*% ubase
		pars<-coefficients(lm(ubase~Wubase-1))
		}
else pars<-initial.value


if(het){
	
Ggmat<-gg_het(Ws2, ubase, n)
optres <-nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol),v=Ggmat, verbose = verbose)
rhotilde<-optres$par

}

else{
	
Ggmat<-gg_hom(Ws2, ubase, n)
optres <- nlminb(pars, optimfunct, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol), v = Ggmat, verbose = verbose)
rhotilde<-optres$par
	
	}
	

# print(rhotilde)
yt  <- y - rhotilde * Ws2 %*% y
wZmat <- Ws2 %*% Zmat
Zt <- Zmat - rhotilde * wZmat

# if(!sarar && is.matrix(endog)) Hmat <- cbind(x, wx, instruments)	
# else Hmat<- cbind(x,wx)	


secondstep<-spatial.ivreg(y =yt , Zmat = Zt, Hmat = Hmat)
delta <- coefficients(secondstep)
utildeb <- y - Zmat %*% delta
# print(delta)

if(het){
	
Ggmat<-gg_het(Ws2, utildeb, n)

gmm.weghts<-psirhorho_het(rhotilde, utildeb, Hmat, Zmat, Ws2)

optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol))	

rhofin<-optres$par
gmm.weghts<-psirhorho_het(rhofin, utildeb, Hmat, Zmat, Ws2)

vcmat <- Omega_het(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)

}

else{
	
Ggmat<-gg_hom(Ws2, utildeb, n)

gmm.weghts<-psirhorho_hom(rhotilde, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )

optres <- nlminb(rhotilde, optimfunct_eff, v = Ggmat, vcmat= gmm.weghts$Phiinv, verbose = verbose, lower= -0.9 + .Machine$double.eps , upper= 0.9 -  .Machine$double.eps,control=list(abs.tol=abs.tol,rel.tol=rel.tol))	

rhofin<-optres$par
gmm.weghts<-psirhorho_hom(rhofin, utildeb, Hmat, Zmat, Ws2, Ggmat$d, Ggmat$v.vec )

vcmat <- Omega_hom(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2, gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat, Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar)



	}





coeff <- as.matrix(c(as.numeric(delta), rhofin))
rownames(coeff)<-c(colnames(Zmat), 'rho')
s2<-crossprod(utildeb)/(n-k)


model.data<-data.frame(cbind(y,x[,-1]))

method<-"gm spatial"

k<-nrow(coeff)
R<-matrix(0,1,k)
R[,((k-1):k)]<-1
Rbeta<-R%*%coeff
Rvar<-R%*% vcmat$Omega %*%t(R)
stat<-as.numeric(t(Rbeta)%*% Rbeta/Rvar)
pval <- pchisq(stat,df=1,lower.tail=FALSE)
W<-list(stat=stat,pval=pval)



 results<-list(coefficients=coeff,var=vcmat$Omega, s2=s2, call=cl, residuals=as.numeric(utildeb), model=model.data,method=method,W=W)
 
 class(results)<-c("sphet", "gstsls")
 
 return(results)
}


impacts.gstsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
    tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    coefs <- drop(obj$coefficients)
    p2 <- length(coefs)
    rho <- coefs[(p2-1)]
    beta <- coefs[1:(p2-2)]
    lambda <- coefs[p2]
# rho is lag coef., lambda is error coef (reversed from function use)
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(beta, rho, lambda)
    Sigma <- obj$var
    irho <- p2-1
    drop2beta <- c(p2-1, p2)
    res <- spdep:::intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu,
        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        interval=NULL, type="lag", tr=tr, R=R, listw=listw, tol=tol,
        empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}
