\name{compute.plot.order}
\alias{compute.plot.order}
\title{Compute Order for Units where Hierarchical Clustering was Performed}
\description{
Compute order for units where hierarchical clustering was performed and where \code{merge} matrix was obtained
}
\usage{
compute.plot.order(merge, n)
}
\arguments{
  \item{merge}{Matrix of merging for hierarchical clustering. See S3 object \code{\link{hclust}}.}
  \item{n}{Number of mergings (usually n-1 where n is number of units).}
}
\details{
Order is computed in order to be able to plot a nice dendrogram (with no crossings).
}
\value{
Vector of the ordering.
}
\references{See \code{\link{hclust}}.}
\author{Kejzar, N.}
\seealso{\code{\link{hclust}},\code{\link{hierarch}}}
\keyword{cluster}