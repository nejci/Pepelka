
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
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	/* INPUTS */

	double *P;                     /* input matrix [N x nClust] */
	double *nElemInCls;             /* vector [1 x nClust] */
	double *AM;                      // matrix N x N
	double *Eind_i, *Eind_j, *E_val; // vectors [E_len x 1]

	// OUTPUTS
	double *avgICS, *avgECS, *numClustHist;

	/* INTERNALS */
	long long N, nClust, numComp, maxCompSize, vote, maxVote, maxVoteInd, nFin;
	long long i, compId, cI, pI, Prow, j, cInd1, cInd2,currInd, row, col, cLen1, cLen2, E_len, edgeInd, sparseCurr;  
	double nElem, sumAM;
	double ICS, ECS;

	mxArray *votes_mat, *compIdx_mat, *compIdxCol_mat, *maxVec_mat, *uniqVec_mat, *posVec_mat;
	mxArray *finalClusters_mat, *nElemInFinCls_mat;
	mxArray *sparseIn[6], *sparseOut[1],*compOut[2];
	double *mG_i, *mG_j, *mG_v, *finalClusters, *nElemInFinCls;
	double *votes, *compIdx, *compIdxCol, *maxVec, *uniqVec, *posVec, *comp, *compSizes;

	mxArray * iterVec_mat;
	double* iterVec;
	mxArray *randIn[1], *randOut[1];

	srand((unsigned int)time(NULL));

	/*===================================================================*/
	/* check for proper number of arguments */
	if(nrhs != 6) {
		mexErrMsgTxt ("6 inputs required.");
	}
	if(nlhs != 3) {
		mexErrMsgTxt ("3 outputs required.");
	}

	// Get matrix P
	nClust = (long long) mxGetM(prhs[0]);
	P = mxGetPr(prhs[0]);

	// Get vector nElemInCls
	nElemInCls = mxGetPr(prhs[1]);

	// Get matrix AM
	N = (long long)mxGetM(prhs[2]); // number of data points
	AM = mxGetPr(prhs[2]);

	// Get Eind_i, _j, E_val
	E_len = (long long)mxGetM(prhs[3]);
	Eind_i = mxGetPr(prhs[3]);
	Eind_j = mxGetPr(prhs[4]);
	E_val = mxGetPr(prhs[5]);



	/* Check inputs consistency */
	if( nClust != mxGetM(prhs[1])) {
		mexErrMsgTxt("nElemInCls size inconsistent with P.");
	}
	if( N != mxGetN(prhs[2])) {
		mexErrMsgTxt("AM not square.");
	}
	// Eind_i, _j, val
	if( mxGetM(prhs[3]) != mxGetM(prhs[4]) || mxGetM(prhs[4]) != mxGetM(prhs[5]) ) {
		mexErrMsgTxt("Length of Eind_i, Eind_j and E_val not equal.");
	}
	if( mxGetN(prhs[3]) != 1) {
		mexErrMsgTxt("Eind_i must be column vector.");
	}
	if( mxGetN(prhs[4]) != 1) {
		mexErrMsgTxt("Eind_j must be column vector.");
	}
	if( mxGetN(prhs[5]) != 1) {
		mexErrMsgTxt("E_val must be column vector.");
	}


	// Create output vectors
	plhs[0] = mxCreateDoubleMatrix(1,E_len+1,mxREAL);
	avgICS = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(1,E_len+1,mxREAL);
	avgECS = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(1,E_len+1,mxREAL);
	numClustHist = mxGetPr(plhs[2]);

	
	for(edgeInd=0; edgeInd<E_len; edgeInd++){


		// CREATE GRAPH AND GET CONNECTED COMPONENTS
		// Add an edge from SMST into meta-cluster graph. Make it symmetric.
		sparseIn[0] = mxCreateDoubleMatrix(2*edgeInd, 1, mxREAL); // i
		sparseIn[1] = mxCreateDoubleMatrix(2*edgeInd, 1, mxREAL); // j
		sparseIn[2] = mxCreateDoubleMatrix(2*edgeInd, 1, mxREAL); // s
		sparseIn[3] = mxCreateDoubleScalar((double)nClust); // m
		sparseIn[4] = mxCreateDoubleScalar((double)nClust); // n
		sparseIn[5] = mxCreateDoubleScalar((double)2*edgeInd); // nzmax
		
		mG_i = mxGetPr(sparseIn[0]);
		mG_j = mxGetPr(sparseIn[1]);
		mG_v = mxGetPr(sparseIn[2]);

		for(sparseCurr=0;sparseCurr<edgeInd;sparseCurr++){
			mG_v[sparseCurr] = E_val[sparseCurr];
			mG_v[sparseCurr+edgeInd] = E_val[sparseCurr];
			mG_i[sparseCurr] = Eind_i[sparseCurr];
			mG_i[sparseCurr+edgeInd] = Eind_j[sparseCurr];
			mG_j[sparseCurr] = Eind_j[sparseCurr];
			mG_j[sparseCurr+edgeInd] = Eind_i[sparseCurr];
		}
		// call sparse.m
		mexCallMATLAB(1, sparseOut, 6, sparseIn, "sparse");

		// call components_mex.mexw64
		mexCallMATLAB(2, compOut, 1, sparseOut, "components_mex");

		// Get components labels	
		comp = mxGetPr(compOut[0]);

		// Get number and sizes of components in a graph
		numComp = (long long) mxGetM(compOut[1]);
		compSizes = mxGetPr(compOut[1]);



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
		mxDestroyArray(iterVec_mat);

		// Count number of unique elements in maxVec - number of components in finalClusters
		uniqVec_mat = mxCreateDoubleMatrix(1,numComp,mxREAL);
		uniqVec = mxGetPr(uniqVec_mat);
		posVec_mat = mxCreateDoubleMatrix(1,numComp,mxREAL);
		posVec = mxGetPr(posVec_mat);

		for(i=0;i<N;i++){
			uniqVec[(long long)maxVec[i]]++;
		}
		nFin = 0;
		for(i=0;i<numComp;i++){
			if(uniqVec[i]>0){
				posVec[i] = (double)nFin;
				nFin++;
			}
		}
		mxDestroyArray(uniqVec_mat);

		finalClusters_mat = mxCreateDoubleMatrix(nFin,N,mxREAL);
		finalClusters = mxGetPr(finalClusters_mat);
		nElemInFinCls_mat = mxCreateDoubleMatrix(nFin,1,mxREAL);
		nElemInFinCls = mxGetPr(nElemInFinCls_mat);


		for(i=0;i<N;i++){
			finalClusters[(long long)posVec[(long long)maxVec[i]]+i*nFin] = 1;
			nElemInFinCls[(long long)(posVec[(long long)maxVec[i]])]++;
		}

		mxDestroyArray(posVec_mat);
		mxDestroyArray(maxVec_mat);

		mxDestroyArray(sparseIn[0]);
		mxDestroyArray(sparseIn[1]);
		mxDestroyArray(sparseIn[2]);
		mxDestroyArray(sparseIn[3]);
		mxDestroyArray(sparseIn[4]);
		mxDestroyArray(sparseIn[5]);
		mxDestroyArray(sparseOut[0]);
		mxDestroyArray(compOut[0]);
		mxDestroyArray(compOut[1]);

		// Compute ICS and ECS
		// Precompute indeces location of ones in finalClusters (aka. find)
		//#pragma omp parallel for shared(nFin,N,finalClusters,nElemInFinCls) private(cInd1,currInd,i)
		for(cInd1=0;cInd1<nFin;cInd1++){
			currInd = 0;
			for(i=0;i<N;i++){
				if(finalClusters[cInd1+i*nFin]){
					finalClusters[cInd1+currInd*nFin] = (double)i;
					currInd++;
				}
				if(currInd >= nElemInFinCls[cInd1]){
					break;
				}
			}
		}

		ICS = 0;
		for (cInd1=0;cInd1<nFin;cInd1++){
			cLen1 = (long long)nElemInFinCls[cInd1];
			if (cLen1 > 1){
				nElem = (cLen1 * (cLen1 -1)) / 2.0;
				sumAM = 0;
				for(i=0;i<cLen1-1;i++){
					row = (long long)finalClusters[cInd1+i*nFin];
					for(j=i+1;j<cLen1;j++){					
						col = (long long)finalClusters[cInd1+j*nFin];
						sumAM += AM[row+col*N];
					}
				}
				ICS += sumAM / nElem;
			}
		}
		ICS /= nFin; // avgICS

		/*===================================================================*/
		// Compute ECS
		ECS = 0;
		for (cInd1=0;cInd1<nFin-1;cInd1++){
			for (cInd2=cInd1+1;cInd2<nFin;cInd2++){
				cLen1 = (long long)nElemInFinCls[cInd1];
				cLen2 = (long long)nElemInFinCls[cInd2];
				nElem = (double)(cLen1 * cLen2);
				sumAM = 0;
				if (nElem > 0){
					for(i=0;i<cLen1;i++){
						row = (long long)finalClusters[cInd1+i*nFin];
						for(j=0;j<cLen2;j++){						
							col = (long long)finalClusters[cInd2+j*nFin];
							sumAM += AM[row+col*N];
						}
					}
					ECS += sumAM / nElem;
				}
			}
		}
		ECS /= (nFin * (nFin-1) / 2.0); // avgECS


		// store ICS and ECS
		avgICS[E_len-edgeInd] = ICS; // store from end to start
		avgECS[E_len-edgeInd] = ECS;

		// save info about recent loop iteration: numClust
		numClustHist[E_len-edgeInd] = (double) nFin;

		mxDestroyArray(finalClusters_mat);
		mxDestroyArray(nElemInFinCls_mat);
	}

}