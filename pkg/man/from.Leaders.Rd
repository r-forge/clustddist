\name{from.Leaders}
\alias{from.Leaders}
\title{Returns List to Help Chaining Leaders and Hierarchical Clustering Methods}
\description{
Computes list that allows to start hierarchical method with leaders from leaders method (takes size and units of clusters from leaders method into consideration)}
\usage{
from.Leaders(data, clustering, as.probs = FALSE, penalty = 1e+06)
}
\arguments{
  \item{data}{Dataframe or matrix of units (in rows). Rows can represent probabilities or frequencies.}
  \item{clustering}{Clustering (vector) that one get from \code{\link{leaders}} method - vector of units with cluster numbers.}
  \item{as.probs}{\code{TRUE} if rows represent probabilities and \code{FALSE} otherwise.}
  \item{penalty}{Parameter to specify how should division with zero be treated.}
}
\details{
Function computes a list that can be used as input into \code{hierarch} method when one chains it after \code{\link{leaders}} method. It takes care of sizes for each cluster (from leaders method), for sum of units and for harmonic sum of units (for each cluster).}
\value{
The result is a list of
  \item{size}{Vector of sizes for each cluster that one gets from \code{\link{leaders}} method.}
  \item{sum}{Matrix of summed units for every cluster that one gets from \code{\link{leaders}} method.}
  \item{harmonic}{Matrix of harmonic sum of units that one gets from \code{\link{leaders}} method.}
}
\references{
Kejzar, N., Korenjak-Cerne, S. and Batagelj, V. (2009). Clustering discrete distributions: The case of patent citations (Submitted)
}
\author{Kejzar, N., Korenjak-Cerne, S. and Batagelj, V.}
\seealso{\code{\link{hierarch}},\code{\link{leaders}},\code{\link{kmeans}} from standard package, \code{\link{pam}} from package cluster}
\examples{
data(patents)
centers <- 10
leaders_clust <- leaders(patents[1:40,], centers=centers)
# get the list
fL <- from.Leaders(patents[1:40,],leaders_clust$clustering)
# input the list into hierarch
hier_clust <- hierarch(leaders_clust$centers,fromLeaders=fL)
plot(hier_clust$object)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
