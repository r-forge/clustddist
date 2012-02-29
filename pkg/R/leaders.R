# setwd("Users/natasa/projekti/R_package/paket/")

leaders <- function (x, centers, err.measure = "d4", as.probs = FALSE,
					iter.max = 10, method = "rnd_cls", stabil=1e-6, penalty=1e6,echo=FALSE){
	#
	# basic function, compute clusters
	# data: x, number of clusters (or centers already): centers, error measure: err.measure
	# if x is probabilities: as.probs, 
	# number of iterations for method: iter.max
        # method: for determining centers (units, rnd_cls)
	# (only considered if centers means "number of centers)
        # echo: whether iterations are written
	#
    x <- data.matrix(x)
    if (any(is.na(x)))
        stop("Missing values not allowed.")
    if (!is.numeric(x)) 
        stop("x is not a numeric dataframe or matrix.")
    x2 <- if (!as.probs) 
              x/apply(x, 1, sum)    # transform rows to probabilities
          else x

    if (iter.max < 1) 
        stop("'iter.max' must be positive.")
		
	if (any(x == 0) && err.measure == "d7") 
	    warning("x contains zeros, error measure d7 could produce empty clusters.")
	
    if (missing(centers)) 
        stop("'centers' must be a number or a matrix.")
    if (length(centers) == 1) {
        k <- centers
        iter <- TRUE
        if (nrow(unique(x2)) < k || k < 1)
            stop("Number of clusters must be in {1,2, .., n}.")
    }
    else {
    	iter <- FALSE
    	iter.max <- 1
    	centers2 <- as.matrix(centers)
    	centers <- if (!as.probs)
    		centers2/apply(centers2,1,sum)
    			else centers2
        if (any(duplicated(centers))) 
            stop("Initial centers are not distinct.")
        if (!(ncol(x2) == ncol(centers))) 
            stop("Must have same number of columns in 'x' and 'centers'.")

        k <- nrow(centers)
        if (nrow(unique(x2)) < k) 
            stop("More cluster centers than data units.")
    }
	
    res <- NULL
	for (i in 1:iter.max){
		if (iter){
		    # compute initial "leaders"
                    centers <- compute.centers(x2, k, method=method,
                                   err.measure=err.measure, penalty=penalty)
		}
        # compute clustering
        curr_res <- leaders.algorithm(x2, centers, err.measure=err.measure,
        		stabil=stabil,penalty=penalty,echo=echo)
        if (is.null(res) || res$error > curr_res$error)
        	res <- curr_res
	}

    out <- list(clustering=res$clu_vector,error=res$error,centers=res$leaders)
    #class(out) <- c("leaders","partition")
    return(out)
}

leaders.algorithm <- function(x, centers, err.measure="d4",stabil=1e-6,
					penalty=1e6,echo=FALSE){
	#
	# algorithm to get k=dim(centers)[1] clusters
	# error measure between unit and center (leader): err.measure
	# when leaders stabilize: stabil
	#
	if (is.null(rownames(x))){
		temp <- 1:dim(x)[1]
		rownames(x) <- as.character(temp)	# names for each unit
	}
	centers <- sort.matrix.rows(centers)
	k <- dim(centers)[1]
	# compute initial clusters (according to err.measure)
	curr_result <- compute.clusters(x,centers,err.measure,penalty)
	iter = 0
	stabilize = FALSE
	while (!stabilize){
		# compute new leaders
		leaders <- compute.leaders(x,k,curr_result$clu_vector,err.measure,penalty)
		leaders <- sort.matrix.rows(leaders)
                if (is.na(max(centers - leaders))){
                  stop("consider less centers! (not all clusters have units)")
                }
		# compute clustering function and check if leaders stabilized
		if (max(centers - leaders) < stabil)
			stabilize = TRUE
		else{
			# compute new clusters
			curr_result <- compute.clusters(x,leaders,err.measure,penalty)
			centers <- leaders
			iter <- iter+1
                        if (echo) {
                          print(iter)
                          print(curr_result$error)
                        }
		}
	}
	return(list(clu_vector=curr_result$clu_vector,error=curr_result$error,leaders=leaders))
}

### SUBROUTINES TO BREAK THE CODE ###

compute.clusters <- function(x, centers, err.measure, penalty=1e6){
	#
	# compute clustering vector: names: units, values: cluster numbers (1:k)
	# k: number of centers
	# calls "err.measure" function for each unit
	#
	clu_vector <- vector(mode="numeric",length=dim(x)[1])
	names(clu_vector) <- rownames(x)
	error <- 0
	for (i in 1:dim(x)[1]){
		temp_e <- eval(call(err.measure,x[i,],centers,penalty))
		clu_vector[i] <- temp_e$clu
		error <- error + temp_e$error
	}
	return(list(clu_vector=clu_vector,error=error))
}

compute.leaders <- function(x, k, curr_clusters,err.measure,penalty=1e6){
	#
	# data: x, number of clusters: k, curr_clusters: units with current clusters
	# calls "error.measure.leader" function for each cluster
	#
	leaders <- matrix(nrow=k,ncol=dim(x)[2])
	for (i in 1:k){
		cluster_units <- x[which(curr_clusters == i),]
		if (!is.matrix(cluster_units))	# only one row
			cluster_units <- t(as.matrix(cluster_units))
		leaders[i,] <- eval(call(paste(err.measure,"cleader",sep="."),cluster_units,penalty))
	}
	return(leaders)	
}

sort.matrix.rows <- function(m){
	
	df <- as.data.frame(m)
	df_sorted <- df[do.call(order,df),]
	m_sorted <- as.matrix(df_sorted)
	dimnames(m_sorted) <- dimnames(m)
	return(m_sorted)	
}

### ERROR MEASURES ###

d1 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- (unit-centers[i,])^2
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d1.cleader <- function(cluster_units,...){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	leader <- apply(cluster_units,2,mean)
	return(leader)
}

d2 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- abs(unit-centers[i,])
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d2.cleader <- function(cluster_units,...){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	leader <- apply(cluster_units,2,median)
	return(leader)
}

d3 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- ((unit-centers[i,])/centers[i,])^2
		# solving problems of zero values (unit or centers component)
		# both 0 -> NaN becomes 0
		temp[which(is.nan(temp))] <- 0
		# centers 0 -> turn to penalty (large number)
		temp[which(is.infinite(temp))] <- penalty
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d3.cleader <- function(cluster_units,...){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	leader <- apply(cluster_units^2,2,mean)/(apply(cluster_units,2,mean))
	leader[which(is.nan(leader))] <- 0 # define 0/0 = 0
	return(leader)
}

d4 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- (unit-centers[i,])^2/centers[i,]
		# solving problems of zero values (unit or centers component)
		# both 0 -> NaN becomes 0
		temp[which(is.nan(temp))] <- 0
		# centers 0 -> turn to penalty (large number)
		temp[which(is.infinite(temp))] <- penalty
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d4.cleader <- function(cluster_units,...){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	leader <- sqrt(apply(cluster_units^2,2,mean))
	return(leader)
}

d5 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- ((unit-centers[i,])/unit)^2
		# solving problems of zero values (unit or centers component)
		# both 0 -> NaN becomes 0
		temp[which(is.nan(temp))] <- 0
		# units 0 -> turn to penalty (large number)
		temp[which(is.infinite(temp))] <- penalty
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d5.cleader <- function(cluster_units,penalty=1e6){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	#penalty <- 1e6
	leader <- vector(mode="numeric",length=dim(cluster_units)[2])
	x <- 1/cluster_units
	for (i in 1:dim(cluster_units)[2]){
		# compute every column (leader variable) separately
		infinite <- which(is.infinite(x[,i]))
		if (length(infinite) >=1) ## CHANGE (see d6.cleader)== dim(cluster_units)[1])
			leader[i] <- 0	# all zeros
		else{
			## x[infinite,i] <- penalty # zeros get penalty
			x2 <- x[,i]^2
			leader[i] <- sum(x[,i])/sum(x2)
		}	
	}
	return(leader)
}

d6 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- (unit-centers[i,])^2/unit
		# solving problems of zero values (unit or centers component)
		# both 0 -> NaN becomes 0
		temp[which(is.nan(temp))] <- 0
		# units 0 -> turn to penalty (large number)
		temp[which(is.infinite(temp))] <- penalty
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d6.cleader <- function(cluster_units,penalty=1e6){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	#penalty <- 1e6
	leader <- vector(mode="numeric",length=dim(cluster_units)[2])
	x <- 1/cluster_units
	for (i in 1:dim(cluster_units)[2]){
		# compute every column (leader variable) separately
		infinite <- which(is.infinite(x[,i]))
		if (length(infinite) >= 1) ## CHANGE (the case of 1 unit zero, leader should be 0)== dim(cluster_units)[1])
			leader[i] <- 0	# all zeros
		else{
			## x[infinite,i] <- penalty # zeros get penalty
			leader[i] <- 1/mean(x[,i])
		}	
	}
	return(leader)
}

d7 <- function(unit, centers,penalty=1e6){
	#
	# compute error term for one unit
	# return number of center and error there
	#
	possible_d <- vector(mode="numeric",length=dim(centers)[1])
	for (i in 1:dim(centers)[1]){
		temp <- (unit-centers[i,])^2/(unit*centers[i,])
		# solving problems of zero values (unit or centers component)
		# both 0 -> NaN becomes 0
		temp[which(is.nan(temp))] <- 0
		# units 0 -> turn to penalty (large number)
		temp[which(is.infinite(temp))] <- penalty
		possible_d[i] <- sum(temp) 
	}
	clu <- which.min(possible_d)
	error <- possible_d[clu]
	return(list(error=error,clu=clu))
}

d7.cleader <- function(cluster_units,penalty=1e6){
	#
	# compute new leader for a cluster of units
	# return leader
	#
	#penalty <- 1e6
	leader <- vector(mode="numeric",length=dim(cluster_units)[2])
	x <- 1/cluster_units
	for (i in 1:dim(cluster_units)[2]){
		# compute every column (leader variable) separately
		infinite <- which(is.infinite(x[,i]))
		if (length(infinite) >= 1) ## CHANGE (see d6.cleader) == dim(cluster_units)[1])
			leader[i] <- 0	# all zeros
		else{
			## x[infinite,i] <- penalty # zeros get penalty
			leader[i] <- sqrt(mean(cluster_units[,i])/mean(x[,i]))
		}	
	}
	return(leader)
}

# GLAVNI PROGRAM - TESTIRANJE
#x <- read.table("test.txt",sep=" ",row.names = 1)
#centers <- x[1:3,]
#centers <- 3
#tt <- system.time(c <- leaders(x, centers=centers)) # izpisujejo se iteracije in trenutna napaka
#print(tt)
#print(c)

#[1] 41 - iterations
#[1] 2313.275
#   user  system elapsed 
# 31.432   0.299  31.808 

#[1] 32 - iterations
#[1] 2301.801
#   user  system elapsed 
# 24.638   0.217  24.873 
