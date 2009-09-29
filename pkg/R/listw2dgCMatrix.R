listw2dgCMatrix<-function (listw) 
{
    if (!inherits(listw, "listw")) 
        stop("not a listw object")
    n <- length(listw$neighbours)
    cardw <- card(listw$neighbours)
    p0 <- as.integer(c(0, cumsum(cardw)))
    scard <- sum(cardw)
    t<-unlist(listw$neighbours)
    t<-t-1
    res <- new("dgCMatrix", i = as.integer(t), p = p0,  Dim = as.integer(c(n,n)), x = unlist(listw$weights))
    res<-t(res)
}
