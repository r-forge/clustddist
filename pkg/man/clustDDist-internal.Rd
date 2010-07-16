\name{clustDDist-internal}
\alias{compute.clusters}
\alias{compute.D}
\alias{compute.hleader}
\alias{compute.leaders}
\alias{d1.cleader}
\alias{d1.D}
\alias{d1.hleader}
\alias{d1}
\alias{diss.matrix}
\alias{harmonic.sum}
\alias{hclustering}
\alias{leaders.algorithm}
\alias{sort.matrix.rows}
\alias{update.diss}
\alias{update.leader}
\alias{update.matrix}
\alias{update.vector}
\title{
Internal Functions in the clustDDist package}
\description{Internal functions of the package.}
\usage{
compute.clusters(x, centers, err.measure, penalty = 1e+06)
compute.D(i, j, leader, err.measure, size, sum, harmonic, penalty, ...)
compute.hleader(i, j, leader, err.measure, size, sum, harmonic, penalty, ...)
compute.leaders(x, k, curr_clusters, err.measure, penalty = 1e+06)
d1.cleader(cluster_units, ...)
d1.D(size1, size2, error, ...)
d1.hleader(size1, size2, leader1, leader2, ...)
d1(unit, centers, penalty = 1e+06)
diss.matrix(leader, err.measure, size, sum, harmonic, penalty, ...)
harmonic.sum(x, penalty = 1e+06)
hclustering(x, err.measure, penalty = 1e+06, nleaders, fromLeaders)
leaders.algorithm(x, centers, err.measure = "d4", stabil = 1e-06, penalty = 1e+06, echo = FALSE)
sort.matrix.rows(m)
update.diss(index_min, diss, leader, err.measure, size, sum, harmonic, penalty, ...)
update.leader(mini, maxi, leader, err.measure, size, sum, harmonic, penalty)
update.matrix(mini, maxi, mat)
update.vector(mini, maxi, vec)

}
\arguments{
  \item{x}{Dataframe or matrix of units (in rows). Rows can represent probabilities or frequencies.}
  \item{centers}{matrix of current leaders}
  \item{err.measure}{String for error function that will be used to calculate the error between current leader and a unit. Possibilities are from d1 to d7 (see reference).}
  \item{penalty}{Parameter to specify how should division by zero be treated.}
  \item{i}{i-th cluster}
  \item{j}{j-th cluster}
  \item{leader}{current leader matrix}
  \item{size}{vector - current size of clusters}
  \item{sum}{vector - current sum of cluster units}
  \item{harmonic}{vector - current harmonic sum of cluster units}
  \item{k}{number of clusters}
  \item{curr_clusters}{vector of cluster numbers for each unit}
  \item{cluster_units}{matrix of units in current cluster}
  \item{size1}{size of first cluster}
  \item{size2}{size of second cluster}
  \item{leader1}{vector - leader of first cluster}
  \item{leader2}{vector - leader of second cluster}
  \item{error}{error among two clusters (computed acc. to error function)}
  \item{unit}{vector - a unit}
  \item{nleaders}{How many leaders should be saved (\code{NULL} or between 1 and n (number of units)). Note: leaders cannot be always computed after the whole procedure.}
  \item{fromLeaders}{\code{NULL} or list that one gets from function \code{from.Leaders} in order to take size and units from leaders method into consideration.}
  \item{stabil}{Stability parameter - used in order to achieve convergence if too many clustering iterations.}
  \item{m}{numeric matrix to be sorted}
  \item{index_min}{vector of 2 indices - two cluster numbers to be merged}
  \item{diss}{current dissimilarity matrix}
  \item{mini}{lowest cluster number to be merged}
  \item{maxi}{highest cluster number to be merged}
  \item{mat}{numeric matrix}
  \item{vec}{numeric vector}
  \item{echo}{if \code{TRUE} function prints error of each iteration}
  \item{\ldots}{other possible parameters}
}
\details{
The functions are written just for error function \code{d1}.  Similarly there are functions for other error functions (\code{d2} to \code{d7}).  Internal functions are used to break the code into a more easily readable functions.
}