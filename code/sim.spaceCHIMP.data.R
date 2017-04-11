################################################################
################################################################
#	Simulate sample spaceCHIMP datasets
#		one "layer" with nearest-neighbor migration
#		one or more daughter populations that expand
#			from a point in space on a parental layer
################################################################
################################################################


################################
#	Function Library
################################

# spatial covariance function 
#	(3 params, exponential decay with distance D)
spatial.cov.func <- function(a0,aD,a2,D){
	cov <- a0 * exp((-aD*D)^a2))
	return(cov)
}

sim.sample.coords <- function(daughter.source,nSamps){
	coords <- vector(mode="list",length=2)
		names(coords) <- c("parental","daughter")
	coords[[1]] <- replicate(2,runif(nSamps[[1]]))
	coords[[2]] <- MASS::mvrnorm(nSamps[[2]],mu=daughter.source,Sigma=matrix(c(0.03,0,0,0.03),nrow=2,ncol=2))
	return(coords)
}

get.samp.dists <- function(coords,daughter.source){
	samp.dists <- vector(mode="list",length=4)
		names(samp.dists) <- c("p-p","p-ds","d-ds","dd")
	samp.dists[["p-p"]] <- fields::rdist(coords[[1]])
	samp.dists[["p-ds"]] <- fields::rdist(coords[[1]],daughter.source)
	samp.dists[["d-ds"]] <- fields::rdist(coords[[2]],daughter.source)
	samp.dists[["d-d"]] <- fields::rdist(coords[[2]])
	return(samp.dists)
}

layer.cov <- function(samp.dists,par.list){
	
}

sim.layer.cov.pars <- function(parent,nSamps,admix){
	a0 <- abs(rnorm(1,0,0.1))
	aD <- abs(rnorm(1,0,0.1))
	a2 <- runif(1,0,2)
	mu <- abs(rnorm(1,0,0.1))
	admix.props <- rep(0,nSamps)
	if(!parent){
		mu <- abs(rnorm(1,0,0.1))
	}
	if(admix){
		admix.props <- rbeta(nSamps,5,3)
	}
	return(list("a0" = a0,
				"aD" = aD,
				"a2" = a2,
				"mu" = mu,
				"admix.props" = admix.props))
}

sim.par.list <- function(admix,nSamps){
	par.list <- vector(mode="list",length=2)
		names(par.list) <- c("parent","daughter")
	par.list[[1]] <- sim.layer.cov.pars(parent=TRUE,nSamps[[1]],admix)
	par.list[[2]] <- sim.layer.cov.pars(parent=FALSE,nSamps[[2]],admix)
	return(par.list)
}

sim.spaceCHIMP.data <- function(nSamps,nLoci,sampleSize,admix=FALSE){
	daughter.source <- matrix(rbeta(2,5,2),nrow=1,ncol=2)
	coords <- sim.sample.coords(daughter.source,nSamps)
	samp.dists <- get.samp.dists(coords,daughter.source)
	par.list <- sim.par.list(admix=FALSE,nSamps)
	layer.cov <- 
}


	plot(coords[[1]],col="blue",pch=19,ylim=range(coords),xlim=range(coords)) ; points(daughter.source,col=2,pch=8) ; points(coords[[2]],pch=20,col=2)