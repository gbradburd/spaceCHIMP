##################################################################
##################################################################
##	run ms on a lattice
##################################################################
##################################################################

require(geoStructure)

within.cluster.m <- function(sample.inds,migration.index,pairwise.migration.matrix){
	#recover()
	migration.rate.vector <- c()
	for(i in 1:nrow(migration.index)){
		migration.rate.vector <- c(migration.rate.vector,
									sprintf("-m %s %s %s",
										sample.inds[migration.index[i,1]],
										sample.inds[migration.index[i,2]],
										pairwise.migration.matrix[migration.index[i,1],migration.index[i,2]]))
	}
	return(migration.rate.vector)
}

write.admixture.event.calls <- function(admix.list,n.pops){
# recover()
	admixture.call <- c()
	for(i in 1:length(admix.list$sources)){
		admixture.call <- c(admixture.call,
							sprintf("-es %s %s %s",
										admix.list$time.point,
										admix.list$targets[i],
										1 - admix.list$admixture.proportions[i]),
							sprintf("-ej %s %s %s",
										admix.list$time.point + 0.000001,
										n.pops+1+(i-1),
										admix.list$sources[i]))							
	}
	return(admixture.call)
}

get.migration.dist <- function(pop.dist){
	pop.dist[which(pop.dist > sqrt(2) | pop.dist < 1)] <- Inf
	return(pop.dist)
}

write.cluster.merges <- function(sample.inds,split.times,daughter.source){
	merge.vector <- unlist(lapply(sample.inds[[1]][2]:sample.inds[[1]][length(sample.inds[[1]])],
							function(i){
								sprintf("-ej %s %s %s",
										split.times[1],
										i,
										sample.inds[[1]][1])}))
	merge.vector <- c(merge.vector,
						sprintf("-ej %s %s %s",
										split.times[2],
										sample.inds[[2]][1],
										daughter.source))
	merge.vector <- c(merge.vector,
						unlist(lapply(sample.inds[[2]][2]:sample.inds[[2]][length(sample.inds[[2]])],
							function(i){
								sprintf("-ej %s %s %s",
										split.times[3],
										i,
										sample.inds[[2]][1])})))
	return(merge.vector)
}

write.migration.rates <- function(sampling.coords,n.sampled.pops,migration.rate,split.times,admix.list=NULL,daughter.source){
	#recover()
	pop.dist <- lapply(sampling.coords,function(x){get.migration.dist(fields::rdist(x))})
	pairwise.migration.matrix <- lapply(1:2,function(i){migration.rate[i]/pop.dist[[i]]})
	migration.index <- lapply(pairwise.migration.matrix,function(x){which(x != 0,arr.ind=TRUE)})
	sample.inds <- list("parental" = 1:nrow(sampling.coords[[1]]),"daughter" = 1:nrow(sampling.coords[[2]]) + nrow(sampling.coords[[1]]))
	n.sampled.pops <- sum(unlist(lapply(sampling.coords,nrow)))
	migration.rate.vector <- c()
	for(i in 1:2){
		migration.rate.vector <- c(migration.rate.vector,
									within.cluster.m(sample.inds[[i]],migration.index[[i]],pairwise.migration.matrix[[i]]))
	}
		migration.rate.vector <- c(migration.rate.vector,
									write.cluster.merges(sample.inds,split.times,daughter.source))
	if(!is.null(admix.list)){
		migration.rate.vector <- c(migration.rate.vector,
									write.admixture.event.calls(admix.list,n.sampled.pops))
	}
	return(migration.rate.vector)
}

# code cannibalized from Dan Denison that reads ms output into R
read.ms.haplotype.matrices <- function(nsam, ndraws, ms.output.file) {
    txt <- scan(file=ms.output.file, what=character(0), sep="\n", quiet=TRUE)
    h <- list()
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    marker <- grep("segsites", txt)
    stopifnot(length(marker) == ndraws)
    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=":"), function(vec) as.integer(vec[2]))
    for(draw in seq(along=marker)) {
        if(!(draw %% 100)) cat(draw, " ")
        if(segsites[draw] > 0) {
            haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            haplotypes <- strsplit(haplotypes, split="")
            h[[draw]] <- sapply(haplotypes, function(el) c(as.integer(el)))
            ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
            if(segsites[draw] == 1) h[[draw]] <- as.matrix(h[[draw]])
            ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
            else h[[draw]] <- t(h[[draw]])
        }
        else h[[draw]] <- matrix(nrow=nsam, ncol=0)
        stopifnot(all(dim(h[[draw]]) == c(nsam, segsites[draw])))  
    }
    cat("\n")
    h
}

# function that actually calls ms with the call put together by all these R functions
ms <- function(sampling.coords,n.chromo,theta,migration.rate,split.times,admix.list=NULL,daughter.source){
		ms.output.file <- "ms_output"
		random.seeds <- c(sample(1:100000,3,replace=TRUE))
		n.sampled.pops <- sum(unlist(lapply(sampling.coords,nrow)))
		call <- paste(
					sprintf(
						"/Applications/ms.folder/msdir/ms %s 1 -t %s -s 1 -I %s %s 0.0-m %s -seeds %s %s %s",
							n.sampled.pops*n.chromo,
							theta,
							n.sampled.pops,
							paste0(rep(n.chromo,n.sampled.pops),collapse=" "),
							paste0(write.migration.rates(sampling.coords,
														 n.sampled.pops,
														 migration.rate,
														 split.times,
														 admix.list,
														 daughter.source),
														 collapse=" "),
							random.seeds[1],
							random.seeds[2],
							random.seeds[3]
					),
					"-T", ">", ms.output.file
				)
		cat(call,file="ms_call.txt")
		system(call)
		read.ms.haplotype.matrices(nsam=n.sampled.pops*n.chromo,ndraws=1,ms.output.file=ms.output.file)
}

# generate times and migration rates in ms units
generate.ms.command.line.values <- function(diploid.population.size,locus.size,per.bp.mu,migration.fraction,split.times){
	ms.command.line.values <- vector("list",length=4)
		names(ms.command.line.values) <- c("theta","m","split.times")
			ms.command.line.values$theta <- 4*diploid.population.size*per.bp.mu*locus.size
			ms.command.line.values$m <- 4*diploid.population.size*migration.fraction
			ms.command.line.values$split.times <- split.times/(4*diploid.population.size)
			# ms.command.line.values$time <- 4*diploid.population.size*generations.ago
			#admixture on time scale more recent than 1/(4Nm k) won't spread
	return(ms.command.line.values)
}

check.migration.rate <- function(migration.rate){
	if(length(migration.rate != 2)){
		param <- rep(migration.rate,2)
	}
	return(param)
}

check.split.times <- function(split.times){
	if(length(split.times) != 3){
		stop("must specify 3 separate split times\n\n")
	}
	if((!split.times[3] <= split.times[2]) |
	   (!split.times[2] <= split.times[1])){
		stop("split times must be specified in order from oldest to most recent\n\n")	
	}
	return(split.times)
}

sim.spCor.admix.props <- function(N,K,geoDist){
	#recover()
	spCor <- 5*exp(-geoDist/0.7)
	t.ad.props <- matrix(NA,nrow=N,ncol=K)
	t.ad.props <- apply(t.ad.props,2,function(x){MASS::mvrnorm(1,rep(0,N),Sigma=spCor)})
	t.ad.props <- apply(t.ad.props,2,function(x){x + abs(min(x))})
	lt.ad.props <- exp(t.ad.props)
	sum.ad.props <- rowSums(lt.ad.props)
	ad.props <- apply(lt.ad.props,2,function(x){x/sum.ad.props})
	return(ad.props)
}


make.par.list <- function(sampling.coords,sampled.pops,admix.list,drop,allele.counts,sample.sizes,daughter.source){
	#recover()
	layer.indices <- c(rep("parental",nrow(sampling.coords[[1]])),
					   rep("daughter",nrow(sampling.coords[[2]])))
	layer.indices[admix.list$targets] <- "admixed"
	#simulation parameters
	admix.props  <- matrix(NA,nrow=sampled.pops,ncol=2)
		if(!is.null(admix.list)){
			for(i in 1:sampled.pops){
				admix.props[i,] <- switch(layer.indices[i],
											"parental" = c(1,0),
											"daughter" = c(0,1),
											"admixed"=c(NA,NA))
			}
			admix.props[admix.list$targets,] <- cbind(admix.list$admixture.proportions,1 - admix.list$admixture.proportions)
		}
		pd.coords <- sampling.coords
		coords <- rbind(sampling.coords[[1]],
							sampling.coords[[2]])
	#prune for analysis
		if(drop=="nobody"){
			to.drop <- FALSE
		} else if(drop == "admix.sources"){
			to.drop <- 1:sampled.pops %in% admix.list$sources
		}
		layer.indices <- layer.indices[!to.drop]
		admix.props <- admix.props[!to.drop,]
		allele.counts <- allele.counts[!to.drop,]
		sample.sizes <- sample.sizes[!to.drop,]
		coords <- 	coords[!to.drop,]
	#make parameter list
		par.list <- list("layer.indices" = layer.indices,
					 	 "admix.props" = admix.props,
					 	 "daughter.source" = daughter.source)
		data.list <- list("allele.freqs" = allele.counts/sample.sizes,
						  "sample.sizes" = sample.sizes,
						  "coords"	= coords,
						  "pd.coords" = pd.coords)
		sim.list <- list("par.list" = par.list,
						 "data.list" = data.list)
		return(sim.list)
}

get.daughter.source <- function(sampling.coords){
	daughter.source <- sample(c(1:nrow(sampling.coords[[1]])),1)
	return(daughter.source)
}

# make dataset for use by, e.g., spatialStructure
generate.spChimp.dataset <- function(n.loci,sampling.coords,n.chromo,theta,migration.rate,split.times,admix.list=NULL,drop="nobody"){
	#recover()
	#Allele Counts & Sample sizes
		migration.rate <- check.migration.rate(migration.rate)
		split.times <- check.split.times(split.times)
		daughter.source <- get.daughter.source(sampling.coords)
		data.matrix <- do.call(cbind,
							replicate(n.loci,
								ms(sampling.coords,n.chromo,theta,migration.rate,split.times,admix.list,daughter.source)))
		sampled.pops <- sum(unlist(lapply(sampling.coords,nrow)))
		population.membership <- unlist(lapply(1:sampled.pops,function(i){rep(i,n.chromo)}))
		allele.counts <- matrix(0,nrow=sampled.pops,ncol=n.loci)
		for(i in 1:sampled.pops){
			allele.counts[i,] <- colSums(data.matrix[which(population.membership==i),,drop=FALSE])
		}
		sample.sizes <- matrix(n.chromo,nrow=sampled.pops,ncol=n.loci)
	#Return sim output
		sim.dataset <- make.par.list(sampling.coords,sampled.pops,admix.list,drop,allele.counts,sample.sizes,daughter.source)
		sim.dataset$par.list[["theta"]] <- theta
		sim.dataset$par.list[["migration.rate"]] <- migration.rate
		sim.dataset$par.list[["split.times"]] <- split.times
		
	return(sim.dataset)
}