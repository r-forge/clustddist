\name{clustDDist-package}
\alias{clustDDist-package}
\alias{clustDDist}
\docType{package}
\title{
Clustering (Discrete) Distributions}
\description{
Implements functions for clustering units that are described with distributions.  
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
\item{\code{\link{from.Leaders}}}{computes list that allows to start hierarchical method with leaders from leaders method (takes size and units of clusters from leaders method into consideration)}
\item{\code{\link{compute.plot.order}}}{computes order for plotting dendrogram of hierarchical method}
\item{\code{\link{patents}}}{data of citation frequency for US patents, years 1980-1999}
}
}
\author{
Kejzar, N., Korenjak-Cerne, S. and Batagelj, V.
Maintainer: Natasa Kejzar <natasa.kejzar@mf.uni-lj.si>
}
\references{
N. Kejzar, S. Korenjak-Cerne, and V. Batagelj: Clustering of distributions : A case of patent citations.J. classif., 2011, 28, doi: 10.1007/s00357-011-9084-x.}
\keyword{ package }
\seealso{
\code{\link[pkg]{cluster}}
}
