// % function [wcl,pc]= weightCl(E,no_allcl)
// %==========================================================================
// % FUNCTION: wcl = weightCl(E)
// % DESCRIPTION: This function computes weight for each pair of clusters using
// %              their shared members (Jaccard Coefficient)
// %
// % INPUTS:   E = N-by-M matrix of cluster ensemble
// %
// % OUTPUT: wcl = an weighted cluster matrix
// %==========================================================================
// % copyright (c) 2010 Iam-on & Garrett
// % optimization for speed: Nejc Ilc, 2014
// %==========================================================================

// compile: mex -largeArrayDims OPTIMFLAGS="/openmp $OPTIMFLAGS" weightCl.c
#include "mex.h"
#include <math.h>
#include <omp.h>
#include "limits.h"

long long myRound(double x) { return (long long)floor(x+0.5); }

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    long long nCls;               /* input: number of all clusters */
    double *E;                 /* MxN input matrix E */
    long long N,M;                /* size of matrix E */
    double * WCL;              /* output matrix WCL */
    double * PC;               /* output matrix PC */
    //mwSize PC_dims[2];         /* dimensions of PC*/
    long long Ei, WCLi, PCi;
    long long Erow,label, WCLrow, WCLcol;
    long long numPair;
    long long numIntersect, numUnion;
    int PCi_sum;
    double WCLval;
    
    
    /* check for proper number of arguments */
    if(nrhs != 2) {
        mexErrMsgIdAndTxt("LCE:weightCl:nrhs","Two inputs required.");
    }
    if(nlhs < 1) {
        mexErrMsgIdAndTxt("LCE:weightCl:nlhs","At least one output required.");
    }
    if(nlhs > 2) {
        mexErrMsgIdAndTxt("LCE:weightCl:nlhs","Max two output required.");
    }
    
    /* make sure the second input argument is scalar */
    if( !mxIsDouble(prhs[1]) ||
            mxIsComplex(prhs[1]) ||
            mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("LCE:weightCl:notScalar","Second input argument must be a scalar.");
    }
    
    /* get number of rows and columns of first input argument */
    N = mxGetM(prhs[0]); // number of data points
    M = mxGetN(prhs[0]); // number of ensemble members
    
    /* get the value of the scalar input  */
    nCls = (long long)mxGetScalar(prhs[1]);
    
    /* create a pointer to the real data in the input matrix  */
    E = mxGetPr(prhs[0]);
    
    /* create the output matrix WCL */
    plhs[0] = mxCreateDoubleMatrix(nCls,nCls,mxREAL);
    /* get a pointer to the real data in the output matrix */
    WCL = mxGetPr(plhs[0]);
    
    /* create the output matrix WCL */
    //PC_dims[0] = N;
    //PC_dims[1] = nCls;
    //plhs[1] = mxCreateNumericArray(2,PC_dims,mxINT8_CLASS,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N,nCls,mxREAL);
    /* get a pointer to the real data in the output matrix */
    PC = mxGetPr(plhs[1]);
    
    
    // % pc = zeros(N,no_allcl); % matrix indicates if data point belongs to the cluster (1=y, 0=n), row=data, col = cluster
    // % for i=1:N
    // %     pc(i,E(i,:))=1; % pc(i,j) = 1 if data i belongs to cluster j
    // % end
    //     for (i=0; i<N; i++){
    //         for (j=0; j<M; j++){
    //             label = E[i+j*N];
    //             mexPrintf("i:%d, j: %d, E: %d, PCind: %d\n", i,j,label, i+(label-1)*N);
    //             PC[i+(label-1)*N] = 1;
    //         }
    //     }
    
    #pragma omp parallel for shared(E, N, M) private(Ei,Erow,label)
    for (Ei=0; Ei<N*M; Ei++){
        Erow = Ei % N;
        label = (long long) E[Ei];
        //mexPrintf("Ei:%d, Erow: %d, label: %d, PCind: %d\n", Ei,Erow,label, Erow+(label-1)*N);
        PC[Erow+(label-1)*N] = 1;
    }
    
    // % %find number of shared data points for each pair of clusters ==> intersect/union
    // % wcl = zeros(no_allcl,no_allcl);
    // % for i=1:no_allcl-1
    // %     for ii=i+1:no_allcl
    // %         pcSum = pc(:,i)+pc(:,ii);
    // %         tmp = sum(pcSum>0);
    // %         if tmp > 0
    // %             wcl(i,ii) = sum(pcSum==2) / tmp; %intersection/union
    // %         end
    // %     end
    // % end
    // % wcl = wcl + wcl';
    
    //mexPrintf("int: %d, long: %d, double: %d, mwSize: %d\n", sizeof(int),sizeof(long),sizeof(double),sizeof(mwSize));
    //mexPrintf("INT_MAX: %ul\n",UINT_MAX);
    //mexPrintf("myRound(0.5)=%lf\n",myRound(0.49999));
    
    numPair = nCls*(nCls-1)/2;
    
    #pragma omp parallel for
    for (WCLi=0; WCLi < numPair; WCLi++){
        
        //mexPrintf("Num threads %d, thread ID %d.\n", omp_get_num_threads(), omp_get_thread_num());
        // Compute indeces of upper triangular matrix
        WCLcol = myRound(floor(-0.5 + 0.5 * sqrt(1 + 8.0 * WCLi)) + 2);
        WCLrow = myRound(WCLcol * (3.0 - WCLcol) / 2.0 + WCLi)-1;
        WCLcol -= 1;
        
        numIntersect = 0;
        numUnion = 0;
        for (PCi=0; PCi<N; PCi++){
            PCi_sum = (int)(PC[WCLrow*N + PCi] + PC[WCLcol*N + PCi]);
            if (PCi_sum > 0){
                numUnion++;
            }
            if (PCi_sum == 2){
                numIntersect++;
            }
        }
        if(numUnion > 0){
            WCLval = numIntersect / (double)numUnion;
            WCL[WCLcol*nCls + WCLrow] = WCLval;
            WCL[WCLrow*nCls + WCLcol] = WCLval;
        }        
    }
}


