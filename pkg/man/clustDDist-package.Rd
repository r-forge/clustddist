\name{clustDDist-package}
\alias{clustDDist-package}
\alias{clustDDist}
\docType{package}
\title{
Clustering Discrete Distributions}
\description{
Implements functions for clustering units that are described as discrete distributions.  
There are two methods that are implemented: adapted leaders method and adapted hierarchical clustering method. }
\details{
\tabular{ll}{
Package: \tab clustDDist \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-02-25\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Functions and data in the package:
\itemize{
\item{\code{\link{leaders}}}{computes clustering with adapted leaders method}
\item{\code{\link{hierarch}}}{computes clustering with adapted hierarchical method}
\item{\code{\link{compute.plot.order}}}{computes order for plotting dendrogram of hierarchical method}
\item{\code{\link{patents}}}{data of citation frequency for US patents, years 1980-1999}
}
}
\author{
Kejzar, N., Korenjak-Cerne, S. and Batagelj, V.
Maintainer: Natasa Kejzar <natasa.kejzar@fdv.uni-lj.si>
}
\references{
Kejzar, N., Korenjak-Cerne, S. and Batagelj, V. (2009). Clustering discrete distributions: The case of patent citations (Submitted)
}
\keyword{ package }
\seealso{
\code{\link[pkg]{cluster}}
}
