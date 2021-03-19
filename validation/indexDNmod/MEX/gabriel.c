
#include "mex.h"

#define A(i,j) a[(i)+(j)*M]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Inputs
	double *distMat;
	double tol;
	// Outputs
	double *d;
	// Internal
	long long n;
	long long i,j,k;
	
	/* Check for proper number of input and output arguments */
	if (nrhs != 2) {
		mexErrMsgTxt("2 inputs argument required.");
	}
	if(nlhs != 1){
		mexErrMsgTxt("1 outputs required.");
	}

	// distMat matrix [n x n]
	n = mxGetN(prhs[0]);
	if (n != mxGetM(prhs[0])){
		mexErrMsgTxt("distMat matrix is not square!");
	}
	distMat = mxGetPr(prhs[0]);
	
	// tol
	tol = mxGetScalar(prhs[1]);
	
	// d [n x n], initialized to 0
	plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
	d = mxGetPr(plhs[0]);
	
	// copy d, only upper triangle
	for (i=0;i<n;i++){
		for (j=i+1; j<n; j++){
			d[i+j*n] = distMat[i+j*n];
            //mexPrintf("i: %d, j: %d, d: %f\n",i,j,d[i+j*n]);
		} 
	}
	
	// Cycle thru all possible pairs of points
	for (i=0; i<(n-1); i++){
		for (j=(i+1); j<n; j++){
			for (k=0; k<n; k++){
				if (k != i && k != j){
					if (distMat[i+k*n]+distMat[j+k*n]-tol <= distMat[i+j*n]){
						d[i+j*n] = 0.0;
						break;
					}
				}
			}
		}
	}
}