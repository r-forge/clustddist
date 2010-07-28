
compute.centers <- function (x, k, method = "rnd_cls", err.measure = "d4", penalty=1e6){
  # determine initial centers (random)
  # x data - matrix
  # k number of centers to determine
  # method for determining centers (units, rnd_cls)
  # err.measure to compute "rnd_cls" centers 
  #
  if (method == "rnd_cls"){
    lenx <- dim(x)[1]
    empty_clust <- TRUE
    while(empty_clust) {
      clust <- sample(1:k,lenx,replace=TRUE)
      empty_clust <- any(!(tabulate(clust)))  # 0s become TRUE
    }
    centers <- compute.leaders(x,k,clust,err.measure,penalty)
  }
  else{
    ux <- unique(x)
    ur <- nrow(ux)
    centers <- ux[sample(1:ur, k), , drop = FALSE]
  }
  return(centers)
}
