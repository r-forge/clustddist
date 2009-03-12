# fromLeaders
# calculate for each cluster:
## sum of units
## harmonic sum of units
## size (number of units)

#source("hierarch_new.R")

from.Leaders <- function(data,clustering,as.probs=FALSE,penalty=1e6){
	if(!as.probs){
		data <- data/apply(data, 1, sum) 	
	}
	size <- tabulate(clustering)
	k <- length(size)
	data_length <- dim(data)[2]
	suma <- matrix(nrow=k,ncol=data_length)
	harmonic <- matrix(nrow=k,ncol=data_length)
	for (i in 1:k){
		cluster_data <- data[which(clustering==i),]
		suma[i,] <- apply(cluster_data,2,sum)
		cluster_matrix <- as.matrix(cluster_data)
		cluster_harmonic <- harmonic.sum(cluster_matrix,penalty)
		# internal function from "hierarch.R"
		harmonic[i,] <- apply(cluster_harmonic,2,sum)
	}
	res <- list(size=size,sum=suma,harmonic=harmonic)
	return(res)
}

#harmonic.sum <- function(x,penalty=1e6){
#	res <- 1/x
#	print(res)
#	cat("**********\n")
#	print(which(is.infinite(res)))
#	res[which(is.infinite(res))] <- penalty
#	return(res)	
#}

