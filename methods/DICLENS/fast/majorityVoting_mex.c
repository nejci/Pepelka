
// compile: mex -largeArrayDims OPTIMFLAGS="/openmp $OPTIMFLAGS" majorityVoting_mex.c
#include "mex.h"
#include "stdlib.h"
#include "time.h"

#define RANDOM_BREAK_TIES 1

/* Random number on range 1..n */
static int rand_int(int n) {
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do {
    rnd = rand();
  } while (rnd >= limit);
  return rnd % n;
}

/* Make random permutation of array 1..n with Fisher-Yates method.*/
void shuffle(double *array, long long n) {
  int i, j, tmp;

  for (i = (int)n - 1; i > 0; i--) {
    j = rand_int(i + 1);
    tmp = (int)array[j];
    array[j] = (int)array[i];
    array[i] = tmp;
  }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	/* INPUTS */
	double *comp;                       /* vector [nClust x 1] */
	double *compSizes;                  /* vector [nComp x 1] */
	double *P;                     /* input matrix [N x nClust] */
	double *nElemInCls;             /* vector [1 x nClust] */


	/* INTERNALS */
	long long N, nClust, numComp, maxCompSize, vote, maxVote, maxVoteInd, numUniq;
	long long i, compId, cI, pI, Prow;  
	mxArray *votes_mat, *compIdx_mat, *compIdxCol_mat, *maxVec_mat, *uniqVec_mat, *posVec_mat;
	double *votes, *compIdx, *compIdxCol, *maxVec, *uniqVec, *posVec;

	// OUTPUTS
	double *finalClusters, *nElemInFinCls;


	mxArray * iterVec_mat;
	double* iterVec;
	mxArray *randIn[1], *randOut[1];

	srand((unsigned int)time(NULL));

	/*===================================================================*/
	/* check for proper number of arguments */
	if(nrhs != 5) {
		mexErrMsgTxt ("5 inputs required.");
	}
	if(nlhs != 2) {
		mexErrMsgTxt ("2 outputs required.");
	}

	// Get components labels
	nClust = (long long) mxGetM(prhs[0]);
	comp = mxGetPr(prhs[0]);

	// Get number and sizes of components in a graph
	numComp = (long long) mxGetM(prhs[1]);
	compSizes = mxGetPr(prhs[1]);

	// Get matrix P
	P = mxGetPr(prhs[2]); // mxGetPr(prhs[0]);

	// Get vector nElemInCls
	nElemInCls = mxGetPr(prhs[3]);

	// Get scalar N
	N = (long long)mxGetScalar(prhs[4]); // number of data points

	/* Check inputs consistency */
	if( nClust != mxGetM(prhs[2])) {
		mexErrMsgTxt("Input sizes inconsistent.");
	}
	if( nClust != mxGetM(prhs[3])) {
		mexErrMsgTxt("Input sizes inconsistent.");
	}
	if( mxGetN(prhs[3]) != 1) {
		mexErrMsgTxt("nElemInCls must be col vector.");
	}

	// MAJORITY VOTE

	// Precompute max component size and create matrix of indices - LUT
	maxCompSize = 0;
	for(i=0;i<numComp;i++){
		if(compSizes[i] > maxCompSize){
			maxCompSize = (long long)compSizes[i];
		}
	}

	compIdx_mat = mxCreateDoubleMatrix(numComp,maxCompSize,mxREAL);
	compIdx = mxGetPr(compIdx_mat);

	compIdxCol_mat = mxCreateDoubleMatrix(numComp,1,mxREAL);
	compIdxCol = mxGetPr(compIdxCol_mat);

	for(i=0;i<nClust;i++){
		compId = (long long)comp[i]-1;
		compIdx[(long long)(compId + numComp*compIdxCol[compId])] = (double)i;
		compIdxCol[compId]++;
	}
	mxDestroyArray(compIdxCol_mat);

	votes_mat = mxCreateDoubleMatrix(numComp,N,mxREAL);
	votes = mxGetPr(votes_mat);

	//#pragma omp parallel for if(numComp>10000) shared(numComp,compSizes,votes,P,nClust,compIdx) private(compId,cI,pI,Prow)
	for(compId=0;compId<numComp;compId++){        
		for(cI=0;cI<compSizes[compId];cI++){
			Prow = (long long)compIdx[compId+numComp*cI];
			for(pI=0;pI<nElemInCls[Prow];pI++){
				votes[(long long)(compId + P[Prow+pI*nClust]*numComp)]++;
			}
		}        
	}
	mxDestroyArray(compIdx_mat);

	// find max value and index of each votes column
	maxVec_mat = mxCreateDoubleMatrix(1,N,mxREAL);
	maxVec = mxGetPr(maxVec_mat);


	iterVec_mat = mxCreateDoubleMatrix(1,numComp,mxREAL);
	iterVec = mxGetPr(iterVec_mat);

	for(cI=0;cI<numComp;cI++){
		iterVec[cI] = (double)cI+1;
	}
	if(RANDOM_BREAK_TIES){
		// if rand() is not sufficient (32-bit), call Matlab
		if(numComp >= RAND_MAX){
			randIn[0] = mxCreateDoubleScalar((double)numComp);
			mexCallMATLAB(1,randOut,1,randIn,"randperm");				
			mxDestroyArray(iterVec_mat);
			iterVec_mat = randOut[0];
			iterVec = mxGetPr(iterVec_mat);
			mxDestroyArray(randIn[0]);
		}
		else{
			shuffle(iterVec, numComp);
		}
	}


	//#pragma omp parallel for shared(numComp,votes,N,finalClusters) private(cI,vI,maxVote,maxVoteInd)
	for(i=0;i<N;i++){
		maxVote = 0;
		maxVoteInd = 0;
		for(cI=0;cI<numComp;cI++){
			vote = (long long)votes[(long long)(iterVec[cI]-1)+i*numComp];
			if(vote > maxVote){
				maxVote = vote;
				maxVoteInd = (long long)(iterVec[cI]-1);
			}
		}
		maxVec[i]=(double)maxVoteInd;
	}
	mxDestroyArray(votes_mat);  

	// Count number of unique elements in maxVec - number of components in finalClusters
	uniqVec_mat = mxCreateDoubleMatrix(1,numComp,mxREAL);
	uniqVec = mxGetPr(uniqVec_mat);
	posVec_mat = mxCreateDoubleMatrix(1,numComp,mxREAL);
	posVec = mxGetPr(posVec_mat);

	for(i=0;i<N;i++){
		uniqVec[(long long)maxVec[i]]++;
	}
	numUniq = 0;
	for(i=0;i<numComp;i++){
		if(uniqVec[i]>0){
			posVec[i] = (double)numUniq;
			numUniq++;
		}
	}
	mxDestroyArray(uniqVec_mat);

	plhs[0] = mxCreateDoubleMatrix(numUniq,N,mxREAL);
	finalClusters = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(numUniq,1,mxREAL);
	nElemInFinCls = mxGetPr(plhs[1]);

	for(i=0;i<N;i++){
		finalClusters[(long long)posVec[(long long)maxVec[i]]+i*numUniq] = 1;
		nElemInFinCls[(long long)(posVec[(long long)maxVec[i]])]++;
	}


	mxDestroyArray(posVec_mat);
	mxDestroyArray(maxVec_mat);
}