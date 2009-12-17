\name{listw2dgCMatrix}
\alias{listw2dgCMatrix}
\title{Interface between Matrix class objects and weights list }

\description{Interface between Matrix class objects and weights list 
}

\usage{
listw2dgCMatrix(listw) 
}

\arguments{
\item{listw}{a \code{listw} object created for example by \code{nb2listw} }
}




\value{
Matrix class object: a sparse Matrix}


\author{ Gianfranco Piras \email{gpiras@mac.com}}

\examples{
library(spdep)
data(columbus)
listw<-nb2listw(col.gal.nb)
spW<-listw2dgCMatrix(listw)
}
\keyword{spatial}
