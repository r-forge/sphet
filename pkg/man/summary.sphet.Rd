\name{summary.sphet}
\alias{summary.sphet}
\alias{sumres}

\title{print method for class sphet}

\description{
Method used to print objects of class \code{'summary.sphet'} and \code{'sphet'}
}

\usage{
\method{summary}{sphet}(object,...)
}

\arguments{
\item{object}{an object of class 'sphet'}
\item{...}{additional arguments to be passed}
}


\details{
The summary function summary.sphet 
returns an objects of class 'sphet' 
organized in a coefficient matrix.

}

\seealso{
  \code{\link{gstslshet}}, \code{\link{stslshac}}}

\author{Gianfranco Piras\email{gpiras@mac.com}}

\examples{
library(spdep)
data(columbus)
listw<-nb2listw(col.gal.nb)
res<-gstslshet(CRIME~HOVAL + INC, data=columbus, listw=listw)
summary(res)

}
