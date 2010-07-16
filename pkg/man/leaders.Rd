\name{leaders}
\alias{leaders}
\title{Clustering of Units with Adapted Leaders Method for Different Error Functions}
\description{
Compute clustering with adapted leaders method, adapted stands for different error functions that are computed between units and current leaders
}
\usage{
leaders(x, centers, err.measure = "d4", as.probs = FALSE, iter.max = 10, stabil = 1e-06, penalty = 1e+06, echo = FALSE)
}
\arguments{
  \item{x}{Dataframe or matrix of units (in rows). Rows can represent probabilities or frequencies.}
  \item{centers}{Either matrix of initial leaders for the function or number of clusters that a user would like as a result. In the second case initial leaders are determined within the algorithm (randomly selected units). There should be at least 1 and at most n centers (n meaning the number of units).}
  \item{err.measure}{String for error function that will be used to calculate the error between current leader and a unit. Possibilities are from d1 to d7 (see reference).}
  \item{as.probs}{TRUE if rows represent probabilities and FALSE otherwise.}
  \item{iter.max}{Maximal number of iterations of the algorithm (because leaders algorithm does local optimization). If centers represents a matrix, the argument is ignored.}
  \item{stabil}{Stability parameter - used in order to achieve convergence if too many clustering iterations.}
  \item{penalty}{Parameter to specify how should division with zero be treated.}
  \item{echo}{if \code{TRUE} function prints error of each iteration}

}
\details{
Function iterates two steps of computation: (1) computation of the new leaders of clusters (according to error function) and (2) recomputation of new clusters (for each unit decide the cluster according to error function).  The iteration is stopped when clustering of the last two steps is the same or when stability parameter is larger than the difference between clustering errors of the last two steps (iterations).
}
\value{
The result is a list of
  \item{clustering}{Vector of units (its names or consecutive numbers) with cluster numbers.}
  \item{error}{Clustering error: sum of errors for each cluster which gives an idea of the distance to optimal cluster. Computed according to error function.}
  \item{leaders}{Final leaders of clusters. Note: cluster leaders can be computed from units in the cluster.}
}
\references{
Kejzar, N., Korenjak-Cerne, S. and Batagelj, V. (2009). Clustering discrete distributions: The case of patent citations (Submitted)
}
\author{Kejzar, N., Korenjak-Cerne, S. and Batagelj, V.}
\seealso{\code{\link{hierarch}},\code{\link{kmeans}} from standard package, \code{\link{pam}} from package cluster}
\examples{
data(patents)
centers <- 3
tt <- system.time(clust <- leaders(patents[1:40,], centers=centers))
# prints out current iteration and current clustering error
print(tt) # should be done in about 10 sec
print(clust) # print result
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
