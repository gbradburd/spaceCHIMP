source("code/ms.spChimp.sim.R")

ms.calls <- generate.ms.command.line.values(diploid.population.size = 1e3,
											locus.size = 1e6,
											per.bp.mu = 1e-7,
											migration.fraction = 0.0005,
											split.times = c(1e4,1e3,5e2))


admix.sources <- c(26,27,29,30)
admix.targets <- c(19,20,24,25)
admix.props <- gtools::rdirichlet(4,c(1.5,1.5))[,1]

sim.dataset <- generate.spChimp.dataset(n.loci=5e3,
										 sampling.coords = list("parental" = as.matrix(expand.grid(1:5,1:5)),
						 										"daughter" = as.matrix(expand.grid(4:6,4:6))),
										 n.chromo = 15,
										 theta = ms.calls$theta,
										 migration.rate = ms.calls$m,
										 split.times = ms.calls$split.times,
										 admix.list = list("sources" = admix.sources,
														   "targets" = admix.targets,
														   "admixture.proportions" = admix.props,
														   "time.point" = 0.0001),
										 drop="nobody")
save(sim.dataset,file="sim.spChimp.dataset.Robj")

samp.cov <- cov(t(sim.dataset$data.list$allele.freqs))										 

pdf(file="test.pdf",width=10,height=5)
par(mfrow=c(1,2))
	plot(sim.dataset$data.list$geoDist,samp.cov)
		points(sim.dataset$data.list$geoDist[1:25,1:25],samp.cov[1:25,1:25],col="blue",pch=20)
		points(sim.dataset$data.list$geoDist[26:34,26:34],samp.cov[26:34,26:34],col="red",pch=20)
	plot(sim.dataset$data.list$coords[1:25,],pch=19,cex=4,col=rainbow(25,start=3/6,end=6/6)[as.numeric(cut(colSums(samp.cov[26:34,1:25]),25))])
		text(sim.dataset$data.list$coords[1:25,],labels=1:25)
			points(sim.dataset$data.list$coords[sim.dataset$par.list$daughter.source,,drop=FALSE],col="green",pch=1,cex=5)
	# plot(c(samp.cov[26:34,1:25]))
		# points(samp.cov[26:34,sim.dataset$par.list$daughter.source],col="red",pch=20)
dev.off()
