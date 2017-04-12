get.projection.matrix <- function(mean.sample.sizes){
	k <- length(mean.sample.sizes)
	transformation.matrix <- get.transformation.matrix(mean.sample.sizes)
	qr.transformation.matrix <- qr(t(transformation.matrix))
	projection.matrix <- qr.Q(qr.transformation.matrix)[,1:qr.transformation.matrix$rank]
	stopifnot(qr.transformation.matrix$rank == sum(abs(eigen(transformation.matrix)$values - 1) < 1e-2) )
	return(projection.matrix)
}

get.transformation.matrix <- function(mean.sample.sizes){
	k <- length(mean.sample.sizes)
	transformation.matrix <- diag(k) - matrix(mean.sample.sizes/sum(mean.sample.sizes),nrow=k,ncol=k,byrow=TRUE)
	return(transformation.matrix)
}

project.sample.covariance <- function(sample.covariance,proj.mat){
	proj.samp.cov <- t(proj.mat) %*% sample.covariance %*% proj.mat
	return(proj.samp.cov)
}

get.norm.factor <- function(freqs,n.loci,sample.sizes){
	pseudo.freqs <- rbind(freqs,rep(0.5,n.loci))	
	if(class(sample.sizes)=="matrix"){
		sample.sizes <- rbind(sample.sizes,rep(1,n.loci))
		sum.sample.sizes <- colSums(sample.sizes)
		pseudo.means <- sapply(1:n.loci,
							function(i){
								sum(
									(pseudo.freqs[,i] * sample.sizes[,i]) / 
										sum.sample.sizes[i],na.rm=TRUE)
							}
						)
	} else {
		sample.sizes <- c(sample.sizes,1)
		sum.sample.sizes <- sum(sample.sizes)
		pseudo.means <- apply(pseudo.freqs,2,
							function(x){
								sum(
									(x * c(sample.sizes,1)) / 
										sum.sample.sizes,na.rm=TRUE)
							}
						)
	}
	norm.factor <- sqrt(pseudo.means * (1-pseudo.means))
	norm.factor <- matrix(norm.factor,nrow=nrow(freqs),ncol=n.loci,byrow=TRUE)
	return(norm.factor)
}

get.mean.freqs <- function(freqs,n.loci,sample.sizes){
	if(class(sample.sizes)=="matrix"){
		sum.sample.sizes <- colSums(sample.sizes)
		mean.freqs <- matrix(
						sapply(1:n.loci,
								function(i){
									sum((freqs[,i] * sample.sizes[,i]) / 
										sum.sample.sizes[i],na.rm=TRUE)
								}
						),
					  nrow=nrow(freqs),ncol=n.loci,byrow=TRUE)
	} else {
		sum.sample.sizes <- sum(sample.sizes)		
		mean.freqs <- matrix(
						apply(freqs,2,
								function(x){
									sum((x * sample.sizes) / 
										sum.sample.sizes,na.rm=TRUE)
								}
						),nrow=nrow(freqs),ncol=n.loci,byrow=TRUE)
	}
	return(mean.freqs)
}
	
stdize.freqs <- function(freqs,n.loci,sample.sizes){
	norm.factor <- get.norm.factor(freqs,n.loci,sample.sizes)
	mean.freqs <- get.mean.freqs(freqs,n.loci,sample.sizes)
	mc.freqs <- freqs-mean.freqs
	norm.freqs <- freqs/norm.factor
	std.freqs <- (freqs-mean.freqs)/norm.factor
	stdized.freqs <- list("std.freqs" = std.freqs,
							"norm.freqs" = norm.freqs)
	return(stdized.freqs)
}

stdize.cov <- function(stdized.freqs,proj.mat){
	norm.cov <- cov(t(stdized.freqs$norm.freqs),use="pairwise.complete.obs")
	std.cov <- cov(t(stdized.freqs$std.freqs),use="pairwise.complete.obs")
	proj.std.cov <- project.sample.covariance(std.cov,proj.mat)
	stdized.cov <- list("norm.cov" = norm.cov,
						"std.cov" = std.cov,
						"proj.std.cov" = proj.std.cov)
	return(stdized.cov)
}

process.freq.data <- function(freqs,sample.sizes){
	invars <- apply(freqs,2,function(x){length(unique(x)) == 1})
	freqs <- freqs[,!invars]
	sample.sizes <- sample.sizes[,!invars]
	mean.sample.sizes <- rowMeans(sample.sizes)
	n.loci <- ncol(freqs)
	proj.mat <- get.projection.matrix(mean.sample.sizes)
	stdized.freqs <- stdize.freqs(freqs,n.loci,sample.sizes)
	stdized.cov <- stdize.cov(stdized.freqs,proj.mat)
	obsSigma <- stdized.cov$proj.std.cov
	processed.freq.data <- list("sample.sizes" = mean.sample.sizes,
								"freqs" = freqs,
								"stdized.freqs" = stdized.freqs,
								"stdized.cov" = stdized.cov,
								"obsSigma" = obsSigma,
								"projMat" = proj.mat)
	return(processed.freq.data)
}

viz.output <- function(n.chains,data.block,model.fit,processed.freq.data,sim.dataset){
	output <- lapply(1:n.chains,
						function(i){
							viz.chain.output(i,data.block,model.fit,processed.freq.data,sim.dataset)
						})
	return(output)
}

viz.chain.output <- function(chain.no,data.block,model.fit,processed.freq.data,sim.dataset){
	#recover()
	best <- which.max(get_logposterior(model.fit,inc_warmup=TRUE)[[chain.no]])
	proj.cov <- geoStructure::get.proj.par.cov(model.fit,chain.no,data.block$N)[best,,]
	Z <- extract(model.fit,"sourceCoords",permute=FALSE,inc_warmup=TRUE)[,chain.no,]
	allDist <- fields::rdist(sim.dataset$data.list$coords)
	tmat <- get.transformation.matrix(rowMeans(sim.dataset$data.list$sample.sizes))
	par.cov <- geoStructure::get.par.cov(model.fit,chain.no,data.block$N)[best,,]
	pdf(file=paste0("fit_ch",chain.no,".pdf",sep=""),width=9,height=8)
		par(mar=c(2,2,2,2))
		layout(mat=matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=TRUE))
		plot(data.block$obsSigma, proj.cov,pch=20,
				xlab="observed covariance",
				ylab="parametric fit")
			abline(0,1,col=2)
		plot(allDist,processed.freq.data$stdized.cov$std.cov,col="purple",cex=1.5)
			points(data.block$pDist,
					processed.freq.data$stdized.cov$std.cov[1:data.block$nP,1:data.block$nP],
					col="blue",pch=20,cex=1.5)
			points(data.block$dDist,
					processed.freq.data$stdized.cov$std.cov[(data.block$nP+1):data.block$N,(data.block$nP+1):data.block$N],
					col="red",pch=20,cex=1.5)
			points(allDist,tmat %*% par.cov %*% t(tmat),col=1,pch=20,,cex=0.8)
		posterior <- get_logposterior(model.fit,inc_warmup=TRUE)[[chain.no]]
		n.iter <- length(posterior)
		plot(posterior[(n.iter/2+1):n.iter],type='l')
		plot(data.block$pCoords,pch=19,cex=4,col="blue",
				ylim=range(data.block$pCoords[,2],data.block$dCoords[,2],Z[,2])+c(-0.5,0.5),
				xlim=range(data.block$pCoords[,1],data.block$dCoords[,1],Z[,1])+c(-0.5,0.5))
			text(data.block$pCoords,labels=1:nrow(data.block$pCoords),col="white")
			points(data.block$dCoords+c(0.1,0.1),col="red",cex=4,pch=19)
			text(data.block$dCoords+c(0.1,0.1),labels=1:nrow(data.block$dCoords),col="black")
				points(sim.dataset$data.list$coords[sim.dataset$par.list$daughter.source,,drop=FALSE],col="green",pch=1,cex=5)
			points(Z,col=adjustcolor(1,0.1),pch=20)
	dev.off()
	output <- list("best" = best,
				   "proj.cov" = proj.cov,
				   "Z" = Z,
				   "allDist" = allDist,
				   "tmat" = tmat,
				   "par.cov" = par.cov,
				   "posterior" = posterior)
	return(output)
}

set.sp.param.limits <- function(param){
	lim <- abs(diff(range(param)))/10
	return(lim)
}

make.data.block <- function(pCoords,dCoords,processed.freq.data,earth=FALSE){
	if(earth){
		pDist <- fields::rdist.earth(pCoords)
		dDist <- fields::rdist.earth(dCoords)
	} else {
		pDist <- fields::rdist(pCoords)
		dDist <- fields::rdist(dCoords)
	}
	# pDist <- pDist/sd(c(pDist,dDist))
	# dDist <- dDist/sd(c(pDist,dDist))
	data.block <- list("N" = nrow(pCoords) + nrow(dCoords),
					   "nP" = nrow(pCoords),
					   "nD" = nrow(dCoords),
					   "L" = ncol(processed.freq.data$freqs),
					   "obsSigma" = processed.freq.data$obsSigma,
					   "pCoords" = pCoords,
					   "dCoords" = dCoords,
					   "pDist" = pDist,
					   "dDist" = dDist,
					   "projMat" = processed.freq.data$projMat,
					   "pCentroid" = colMeans(pCoords),
					   "sCoordPrVar" = matrix(c(sd(pDist),0,0,sd(pDist)),nrow=2,ncol=2,byrow=TRUE),
					   "xLow" = min(pCoords[,1]) - set.sp.param.limits(pCoords[,1]),
					   "xHigh" = max(pCoords[,1]) + set.sp.param.limits(pCoords[,1]),
					   "yLow" = min(pCoords[,2]) + set.sp.param.limits(pCoords[,2]),
					   "yHigh" = max(pCoords[,2]) + set.sp.param.limits(pCoords[,2]))
	return(data.block)
}


source("~/Dropbox/spaceCHIMP/code/spaceCHIMP.R")
require(rstan)
load(list.files(pattern="dataset.Robj"))

processed.freq.data <- process.freq.data(sim.dataset$data.list$allele.freqs,
										 sim.dataset$data.list$sample.sizes)
data.block <- make.data.block(pCoords = sim.dataset$data.list$pd.coords[["parental"]],
							  dCoords = sim.dataset$data.list$pd.coords[["daughter"]],
							  processed.freq.data = processed.freq.data,
							  earth=FALSE)

	
model.fit <- stan(model_code=model.block,
				  data=data.block,
				  refresh=1e2,
				  chains=2,
				  iter=5e3,
				  thin=5e3/500)
output <- viz.output(2,data.block,model.fit,processed.freq.data,sim.dataset)
save(data.block,model.fit,processed.freq.data,output,file="output.Robj")

