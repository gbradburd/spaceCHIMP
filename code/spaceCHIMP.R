model.block <- "
functions {
	vector pairwise_distance(matrix pCoords, matrix sourceCoords, int nP) {
		vector[nP] D;
		for(i in 1:nP) {
            	D[i] = sqrt(dot_self(pCoords[i] - sourceCoords[1]));
        	}
    		return D;
    }
    matrix psCov(int nD, int nP, real alpha0, real alphaD, real alpha2, vector geoDist){
    		matrix[nD,nP] psCov;
    		matrix[1,nP] vecCov;
    		for(i in 1:nP){
    			vecCov[1,i] = alpha0 * exp( -(alphaD* geoDist[i])^alpha2);
    		}
    		for(i in 1:nD){
    			psCov[i] = vecCov[1];
    		}
    		return psCov;
    }
    matrix layerCov(int N, real alpha0, real alphaD, real alpha2, real mu, matrix geoDist){
		matrix[N,N] layerCov;
		for(i in 1:N-1){
			for(j in (i+1):N){
				layerCov[i,j] = alpha0 * ( exp( -(alphaD* geoDist[i, j])^alpha2 ) + mu);
				layerCov[j,i] = layerCov[i,j];
			}
		}
		for(i in 1:N){
			layerCov[i,i] = alpha0 * (1 + mu);
		}
    		return layerCov;
    }
    matrix pdCov(int N, int nP, int nD, vector alpha0, vector alphaD, vector alpha2, real mu, matrix pDist, matrix dDist, vector pSourceDist, vector nugget){
		matrix[N,N] Sigma;
		Sigma[1:nP,1:nP] = layerCov(nP, alpha0[1], alphaD[1], alpha2[1], 0, pDist);
		Sigma[(nP+1):(nP+nD),(nP+1):(nP+nD)] = layerCov(nD, alpha0[2], alphaD[2], alpha2[2], mu, dDist);
		Sigma[(nP+1):(nP+nD),1:nP] = psCov(nD, nP, alpha0[1], alphaD[1], alpha2[1], pSourceDist);
		Sigma[1:nP,(nP+1):(nP+nD)] = (Sigma[(nP+1):(nP+nD),1:nP])';
		for(i in 1:N){
			Sigma[i,i] = Sigma[i,i] + nugget[i];
		}
		return Sigma;
	}
}
data {
	int<lower=2> N; 	  			// total number of samples
	int<lower=2> nP;				// total number of parental population samples
	int<lower=2> nD;				// total number of daugther population samples
	int<lower=N+1> L;	    	// number of loci
	cov_matrix[N-1] obsSigma; 	// observed projected covariance
	matrix[N, N-1] projMat;		// projection matrix
	matrix[nP,2] pCoords;		// observed sampling coordinates of parent population samples
	matrix[nD,2] dCoords;		// observed sampling coordinates of daughter population samples
	matrix[nP,nP] pDist;		// matrix of pairwise distances between parent population samples
	matrix[nD,nD] dDist;		// matrix of pairwise distances between daughter population samples
	vector[2] pCentroid;		// centroid of parental coordinates
	matrix[2,2] sCoordPrVar;	// vcov of multivariate normal prior on source coords
	real xLow;
	real xHigh;
	real yLow;
	real yHigh;
//	matrix[1,2] fx_sCoords;
}
parameters {
	vector<lower=0>[2] alpha0;													// sill of the parametric covariance
	vector<lower=0>[2] alphaD;													// effect of geographic distance in the parametric covariance
	vector<lower=0, upper=2>[2]  alpha2;										// exponential slope parameter in the parametric covariance
	real<lower=0> mu;															// drift shared by all daughter population samples following their founding
  	vector<lower=0>[N] nugget; 													// sample-specific variance (allele sampling error + sample-specific drift)
	real<lower=xLow,upper=xHigh> xSourceCoord;
	real<lower=yLow,upper=yHigh> ySourceCoord;
}
transformed parameters {
	matrix[N-1, N-1] projSigma;			// the mean-centered and projected covariance matrix (dim = N-1 x N-1)
	cov_matrix[N] Sigma;				// this specifies the parametric, admixed covariance matrix
	vector[nP] pSourceDist;				// vector of pairwise distances between parental population samples and estimated source location
	matrix[1,2] sourceCoords;			// coordinates in the parent population layer of the source of the daughter population
//	sourceCoords = fx_sCoords;
	sourceCoords[1,1] = xSourceCoord;
	sourceCoords[1,2] = ySourceCoord;
	pSourceDist = pairwise_distance(pCoords,sourceCoords,nP);
	Sigma = pdCov(N, nP, nD, alpha0, alphaD, alpha2, mu, pDist, dDist, pSourceDist, nugget);
	projSigma = quad_form(Sigma, projMat);
}
model {
	alpha0 ~ normal(0,1);											// prior on alpha0
	alphaD ~ exponential(1);										// prior on alphaD
	alpha2 ~ uniform(0,2);											// prior on alpha2
	mu ~ normal(0,1);												// prior on mu
	nugget ~ normal(0,1);											// prior on nugget
	sourceCoords[1] ~ multi_normal(pCentroid,sCoordPrVar);			// prior on source coords
	L*obsSigma ~ wishart(L,projSigma);								// likelihood function	
}
"