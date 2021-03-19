
#include "mex.h"
#include "omp.h"
#include "math.h"
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Inputs
	double *E;  // relabelled ensemble matrix
	long long M, N;
	long long nClust, maxElem;
	// Outputs
	double *P;
	double *nElemInCls;
	// Internal
	long long i, Erow, label, Ntmp, Ninc, nbytes, destination_offset, bytes_to_copy;

	mxArray *Ptmp_mat, *newSpace_mat;
	double *Ptmp, *newSpace, *newptr;


	/* Check for proper number of input and output arguments */    
	if (nrhs != 2) {
		mexErrMsgTxt("Two inputs argument required.");
	} 
	if(nlhs != 2){
		mexErrMsgTxt("2 outputs required.");
	}

	// E matrix [M X N]
	N = mxGetM(prhs[0]); // number of data points
	M = mxGetN(prhs[0]); // number of ensemble members
	E = mxGetPr(prhs[0]);

	// Get size of ensemble M
	nClust = (long long)mxGetScalar(prhs[1]);    

	// Create temp P matrix, half size of N; if we will need more space, we will allocate another matrix
	Ntmp = (long long)floor(N/2.0);
	Ninc = (long long)ceil(N*0.2); // 20% increase size if out of space
	Ptmp_mat = mxCreateDoubleMatrix(nClust,Ntmp,mxREAL);
	Ptmp = mxGetPr(Ptmp_mat);

	// Create vector for storing number of elements in each cluster
	plhs[1] = mxCreateDoubleMatrix(nClust,1,mxREAL);
	nElemInCls = mxGetPr(plhs[1]);


	//------------------------------------------------
	// Compute P and sum of cols
	maxElem = 0;
	//#pragma omp parallel for shared(E, N, M) private(i,Erow,label)
	for (i=0; i<N*M; i++){
		Erow = i % N;
		if(E[i]>0){
			label = (long long)E[i]-1;

			// Need to add some space to Ptmp?
			if(nElemInCls[label] >= Ntmp){
				newSpace_mat = mxCreateDoubleMatrix(nClust,Ninc,mxREAL);
				newSpace = mxGetPr(newSpace_mat);

				nbytes = (nClust) * (Ntmp + Ninc)* sizeof(double);//size of new array
				destination_offset = nClust * Ntmp; //start offset for copying
				bytes_to_copy = nClust * Ninc * sizeof(double);
				//ptr = mxGetPr(prhs[0]);
				newptr = (double *)mxRealloc(Ptmp, nbytes);//reallocate array
				mxSetPr(Ptmp_mat,newptr);
				//ptr = mxGetPr(prhs[1]);
				memcpy(newptr+destination_offset,newSpace,bytes_to_copy);//actual copy
				mxSetN(Ptmp_mat,Ntmp + Ninc);//fix dimension
				Ptmp = newptr;
				Ntmp += Ninc;
				mxDestroyArray(newSpace_mat);
			}

			Ptmp[label+(long long)nElemInCls[label]*nClust] = (double)Erow;
			nElemInCls[label]++;
			if(nElemInCls[label] > maxElem){
				maxElem = (long long)nElemInCls[label];
			}			
		}
	}  
	

	// Create P matrix - trim Ptmp
	plhs[0] = mxCreateDoubleMatrix(nClust,maxElem,mxREAL);
	P = mxGetPr(plhs[0]);
	memcpy(P,Ptmp,maxElem*nClust*sizeof(double)); // copy contents of Ptmp, only relevant
	mxDestroyArray(Ptmp_mat);

}

