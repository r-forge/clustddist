# setwd("Users/natasa/projekti/R_package/paket/")
#source("leaders.R")

hierarch <- function(x, as.probs = TRUE, err.measure = "d4", penalty=1e6, nleaders=NULL, fromLeaders=NULL){
	#
	# basic function, compute hierarchical clustering
	# data: x (units with variable values), 
	# error measure:err.measure
	# penalty: number that substitutes x/0
	# nleaders: number of leaders, that should be returned
	# fromLeaders: object that one gets when combining leaders and hierarch methods
	# 
    x <- data.matrix(x)
    if (any(is.na(x)))
        stop("Missing values not allowed.")
    if (!is.numeric(x)) 
        stop("x is not a numeric dataframe or matrix.")
    x2 <- if (!as.probs) 
          	x/apply(x, 1, sum)    # transform rows to probabilities
          else x
    n <- dim(x)[1]
    if ((!is.null(nleaders)) && ((nleaders<1) || (nleaders>n)))
    	stop("nleaders should be NULL or in [1,n].")
    # compute clustering
	clu <- hclustering(x2,err.measure,penalty,nleaders,fromLeaders)
	# compute order (to be able to plot dendrogram)
	order <- compute.plot.order(clu$merge,dim(clu$merge)[1])
    res <- list(merge=clu$merge,height=clu$height,order=order,
    		labels=attr(x, "dimnames")[[1]],call=match.call(),method=err.measure) 
    class(res) <- "hclust"
    return(list(object=res,saved_leaders=clu$saved_leaders))
}


hclustering <- function(x,err.measure,penalty=1e6,nleaders,fromLeaders){
	n <- dim(x)[1]
	merge <- matrix(nrow=n-1,ncol=2)
	# temporary vector to save current cluster number - in which row it is written
	# in order to construct "merge"
	cluster_number <- seq(-1,-(n),-1)
	height <- vector(mode="numeric",length=n-1)
	leader <- as.matrix(x)
	if (is.null(fromLeaders)){
		# start clustering from one unit
		sum <- as.matrix(x)
		harmonic <- harmonic.sum(x,penalty)
		size <- rep(1,n)
	}
	else{
		# start clustering from clustering leaders from "leaders" method
		sum <- fromLeaders$sum
		harmonic <- fromLeaders$harmonic
		size <- fromLeaders$size	
	}
	# compute dissimilarity matrix
	diss <- diss.matrix(leader,err.measure,size,sum,harmonic,penalty)
	# (get upper right triangle of dissimilarities)
	saved_leaders <- NULL
	
	for (i in 1:(n-1)){
		# compute minimal merging (diss)
		min_dist <- min(diss[which(!is.na(diss))])
		# update height
		height[i] <- min_dist
		index_min <- which(diss == min_dist,arr.ind = TRUE)
		mini <- min(index_min[1,])
		maxi <- max(index_min[1,])
		# take first one that comes (*** according to rows!)
		# UPDATE MERGE - with new cluster numbers!
		merge[i,] <- c(cluster_number[mini],cluster_number[maxi])
		cluster_number[mini] <- i
		cluster_number <- cluster_number[-maxi] 
		# END UPDATE MERGE
		#	leader (compute new, add it on the first free place)
		leader <- update.leader(mini,maxi,leader,err.measure,size,sum,harmonic,penalty)
		if ((!is.null(nleaders)) && nleaders == (n-i))
			saved_leaders <- leader
		# update size,sum,harmonic --- does not depend on err.measure
		sum <- update.matrix(mini,maxi,sum)
		size <- update.vector(mini,maxi,size)
		harmonic <- update.matrix(mini,maxi,harmonic)
		#	diss (compute new diss for new cluster)
		diss <- update.diss(index_min[1,],diss,leader,err.measure,size,sum,harmonic,penalty)
	}
	rez <- list(merge=merge,height=height,saved_leaders=saved_leaders)
    return(rez)
}


compute.plot.order <- function(merge,n){
	#
	# compute order of units in order to plot them
	# in dendrogram (not having crossings on branches)
	#
	temp <- NULL
	#print(merge[n,])	
	if (min(merge[n,]) < 0){
		temp <- abs(min(merge[n,]))
		if (max(merge[n,]) < 0)
			temp <- c(temp,abs(max(merge[n,])))
		else
			temp <- c(temp,compute.plot.order(merge,max(merge[n,])))
		#print(temp)
		return(temp)
	}
	else{
		temp <- c(temp,compute.plot.order(merge,min(merge[n,])))
		temp <- c(temp,compute.plot.order(merge,max(merge[n,])))
		#print(temp)
		return(temp)
	}
}

update.vector <- function(mini,maxi,vec){
	#
	# subroutine to help merging 2 components (vectors) 
	# (sum them and delete the latest)
	#
	vec[mini] <- vec[mini] + vec[maxi]
	vec <- vec[-maxi]
	return(vec)
}

update.matrix <- function(mini,maxi,mat){
	#
	# subroutine to help merging 2 components (row vectors inside matrix) 
	# (sum them and delete the latest)
	#
	mat[mini,] <- mat[mini,] + mat[maxi,]
	mat <- mat[-maxi,]
	return(mat)
}

update.leader <- function(mini,maxi,leader,err.measure,size,sum,harmonic,penalty){
	# compute leader when merging mini and maxi cluster
	new_leader <- compute.hleader(mini,maxi,leader,err.measure,size,sum,harmonic,penalty)
	# put it in the right position inside all leaders
	dimension <- dim(leader)[1] - 1
	leader <- leader[-maxi,]
	if (dimension > 1)
		leader[mini,] <- new_leader
	else
		leader <- t(as.matrix(new_leader))
	return(leader)
}

update.diss <- function(index_min,diss,leader,err.measure,size,sum,harmonic,penalty,...){
	#
	# update according to merging i and j cluster (index_min)
	# ... parameters: optional (depends upon err.measure)
	# i < j: compute new D[i, ...], move all D[k, ...] one up (k > j)
	# 		 delete D[... ,j]
	# 
	# diminish diss for 1 row and column (maxi)
	mini <- min(index_min)
	maxi <- max(index_min)
	diss <- diss[-maxi,-maxi] 
	if (mini != 1){
		for (j in 1:(mini-1)){
			i <- mini
			d <- compute.D(i,j,leader,err.measure,size,sum,harmonic,penalty,...)
			diss[i,j] <- d
		}
	}
	if (!is.vector(diss) && (mini+1) < dim(diss)[1]){
		for (i in (mini+1):dim(diss)[1]){
			j <- mini
			d <- compute.D(i,j,leader,err.measure,size,sum,harmonic,penalty,...)
			diss[i,j] <- d
		}
	}
	return(diss)
}

diss.matrix <- function(leader,err.measure,size,sum,harmonic,penalty,...){
	# size (vector) - from 1s only, sum, harmonic (matrices)
	n <- dim(leader)[1]
	diss <- matrix(nrow=n,ncol=n)
	for (i in 2:n){
		for (j in 1:(i-1)){
			d <- compute.D(i,j,leader,err.measure,size,sum,harmonic,penalty,...)
			diss[i,j] <- d
		}	
	}
	#print(diss)
	return(diss)
}

compute.hleader <- function(i,j,leader,err.measure,size,sum,harmonic,penalty,...){
	#
	# subroutine that helps computing the one leader according to 
	# err.measure and i and j clusters
	#
        temp_leader <- eval(call(paste(err.measure,"hleader",sep="."),leader[i,],leader[j,],
                                 size[i],size[j],harmonic[i,],harmonic[j,],
                                 sum[i,],sum[j,],penalty))
	return(temp_leader)
}

compute.D <- function(i,j,leader,err.measure,size,sum,harmonic,penalty,...){
	#
	# subroutine that helps computing dissimilarity of i, j clusters
	# according to err.measure
	#	
	if (err.measure == "d1"){
		error <- eval(call(err.measure,leader[i,],t(as.matrix(leader[j,]))))
		error <- error$error
		d <- d1.D(size[i],size[j],error)
	}else{
	    # computation depends on error measure
	    # compute leader (before computing d[i,j])
	    new_leader <- compute.hleader(i,j,leader,err.measure,size,sum,harmonic,penalty,...)
	    error_u <- eval(call(err.measure,leader[i,],t(as.matrix(new_leader))))
	    error_u <- error_u$error
	    error_v <- eval(call(err.measure,leader[j,],t(as.matrix(new_leader))))
	    error_v <- error_v$error
            d <- eval(call(paste(err.measure,".D",sep=""),size[i],size[j],error_u,error_v,
                               leader[i,],leader[j,],sum[i,],sum[j,],harmonic[i,],harmonic[j,],penalty))
	}
	return(d)
}

# new LEADER when two clusters merged

d1.hleader <- function(leader1,leader2,size1,size2,...){
	leader <- (size1*leader1 + size2*leader2)/(size1+size2)
	return(leader)
}

d3.hleader <- function(leader1,leader2,size1=NULL,size2=NULL,
                       harmonic1=NULL,harmonic2=NULL,sum1,sum2,...){
	if (sum(sum1+sum2) == 0){
		leader <- rep(0,length(leader1))	
	}else{
		leader <- (sum1*leader1 + sum2*leader2)/(sum1+sum2)
		leader[which(is.nan(leader))] <- 0
	}
	return(leader)
}

d4.hleader <- function(leader1,leader2,size1,size2,...){
	leader <- sqrt((size1*(leader1)^2 + size2*(leader2)^2)/(size1+size2))
	return(leader)
}

d5.hleader <- function(leader1,leader2,size1=NULL,size2=NULL,harmonic1,harmonic2,
                       sum1=NULL,sum2=NULL,penalty=1e6,...){
	temp1 <- harmonic1/leader1
	temp1[which(is.infinite(temp1))] <- penalty
	temp2 <- harmonic2/leader2
	temp2[which(is.infinite(temp2))] <- penalty
	leader <- (harmonic1+harmonic2)/(temp1 + temp2)
	return(leader)
}

d6.hleader <- function(leader1=NULL,leader2=NULL,size1,size2,
                       harmonic1,harmonic2,sum1=NULL,sum2=NULL,...){
	leader <- (size1+size2)/(harmonic1+harmonic2)
	return(leader)
}

# penalty
d7.hleader <- function(leader1,leader2,size1=NULL,size2=NULL,
                       harmonic1=NULL,harmonic2=NULL,sum1,sum2,penalty=1e6,...){
	if (sum(sum1+sum2) == 0){
		leader <- rep(0,length(leader1))
	}else{
		# DODATI
                temp1 <- sum1/(leader1)^2
                temp1[which(is.infinite(temp1))] <- penalty
                temp2 <- sum2/(leader2)^2
                temp2[which(is.infinite(temp2))] <- penalty
		leader <- sqrt((sum1+sum2)/(temp1+temp2))
                leader[which(is.nan(leader))] <- 0
	}
	return(leader)
}

# DISSIMILARITY among two clusters

d1.D <- function(size1,size2,error,...){
	return(size1*size2/(size1+size2)*error)
}

d3.D <- function(size1=NULL,size2=NULL,error1,error2,leader1,leader2,
                 sum1,sum2,harmonic1=NULL,harmonic2=NULL,penalty=1e6,...){
	temp1 <- sum1/leader1*error1
	temp1[which(is.infinite(temp1))] <- penalty
	temp1[which(is.nan(temp1))] <- 0
	temp2 <- sum2/leader2*error2
	temp2[which(is.infinite(temp2))] <- penalty
	temp2[which(is.nan(temp2))] <- 0
	return(sum(temp1 + temp2))
}

d4.D <- function(size1,size2,error1,error2,...){
	return(size1*error1 + size2*error2)
}

d5.D <- function(size1=NULL,size2=NULL,error1,error2,leader1,leader2,
                 sum1=NULL,sum2=NULL,harmonic1,harmonic2,...){
	return(sum(harmonic1*leader1*error1 + harmonic2*leader2*error2))
}

d6.D <- function(size1,size2,error1,error2,...){
	return(d4.D(size1,size2,error1,error2,...))
}

d7.D <- function(size1=NULL,size2=NULL,error1,error2,leader1,leader2,sum1,sum2,
                 harmonic1=NULL,harmonic2=NULL,penalty=1e6,...){
	return(d3.D(error1=error1,error2=error2,leader1=leader1,leader2=leader2,sum1=sum1,sum2=sum2,penalty=penalty,...))
}

harmonic.sum <- function(x,penalty=1e6){
	res <- 1/x
	res[which(is.infinite(res))] <- penalty
	return(res)	
}


# GLAVNI PROGRAM - TESTIRANJE
#x <- read.table("test_short.txt",sep=" ",row.names = 1)
#x <- as.matrix(x)
#y <- hierarch(x,err.measure="d4")
#plot(y$object)
