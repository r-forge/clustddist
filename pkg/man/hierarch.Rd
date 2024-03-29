\name{hierarch}
\alias{hierarch}
\title{Clustering of Units with Adapted Hierarchical Clustering Method for Different Error Functions}
\description{
Compute clustering with adapted hierarchical clustering method, adapted stands for different error functions that are computed between units and current leaders
}
\usage{
hierarch(x, as.probs = TRUE, err.measure = "d4", penalty = 1e+06, nleaders = NULL, fromLeaders = NULL)
}
\arguments{
  \item{x}{Dataframe or matrix of units (in rows). Rows can represent probabilities or frequencies.}
  \item{as.probs}{\code{TRUE} if rows represent probabilities and \code{FALSE} otherwise.}
  \item{err.measure}{String for error function that will be used to calculate the error between current leader and a unit. Possibilities are from d1 to d7 (but not d2, see reference).}
  \item{penalty}{Parameter to specify how should division by zero be treated.}
  \item{nleaders}{How many leaders should be saved (\code{NULL} or between 1 and n (number of units)). Note: leaders cannot be always computed after the whole procedure.}
  \item{fromLeaders}{\code{NULL} or list that one gets from function \code{from.Leaders} in order to take size and units from leaders method into consideration.}

}
\details{
Function generates dissimilarity matrix (according to error function) and then in each step agglomerates the nearest two units, later clusters.  This stops when all units are in the same cluster.  Every step dissimilarity matrix is recomputed according to error functions (see reference for details).
}
\value{
The result is a list of
  \item{object}{S3 object used for hierarchical clustering procedure from base package \code{stats}. See \code{\link{hclust}}. This allows to use methods \code{plot} for dendrogram plot etc.}
  \item{saved_leaders}{Matrix of \code{nleaders} leaders.}
  ...
}
\references{
N. Kejzar, S. Korenjak-Cerne, and V. Batagelj: Clustering of distributions : A case of patent citations.J. classif., 2011, 28, doi: 10.1007/s00357-011-9084-x.}
\author{Kejzar, N., Korenjak-Cerne, S. and Batagelj, V.}
\seealso{\code{\link{leaders}} from package \code{clustDDist},\code{\link{hclust}} from standard package, \code{\link{agnes}} from package \code{cluster}}
\examples{
data(patents)
patents_short <- patents[1:20,] # to get shorter version of data
clust <- hierarch(patents_short,err.measure="d4",as.probs=FALSE)
plot(clust$object)
}
\keyword{cluster}