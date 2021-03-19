
// compile: mex -largeArrayDims OPTIMFLAGS="/openmp $OPTIMFLAGS" cts_S_mex.c
#include "mex.h"
#include <math.h>
#include <omp.h>


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    /* INPUTS */
    double *E;           /* input matrix E [MxN] */
    double *WCL;         /* input matrix WCL [nCls x nCls] */
    double *C;           /* input matrix WCL [1 x M] */
    double dc;           /* decay ratio - scalar */
    
    /* OUTPUTS */
    double * S;          /* output matrix S [N x N] */
    
    /* INTERNALS */
    long long N,M;       /* size of matrix E */
    long long nCls, nCls2; /* number of all clusters */
    long long Crow, Ccol;
    mxArray *wCT_mat;
    double *wCT;         /* matrix wCT [nCls x nCls]*/
    double wCT_val, wCT_max, local_max, S_val;
    long long q, i, j, m, WCLi, Ei, Ej;
    
    /*===================================================================*/
    /* check for proper number of arguments */
    if(nrhs != 4) {
        mexErrMsgIdAndTxt("LCE:cts:nrhs","Four inputs required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("LCE:cts:nlhs","One output required.");
    }
    
    /* make sure the fourth input argument is scalar */
    if( !mxIsDouble(prhs[3]) ||
            mxIsComplex(prhs[3]) ||
            mxGetNumberOfElements(prhs[3])!=1 ) {
        mexErrMsgIdAndTxt("LCE:cts:notScalar","dc must be a scalar.");
    }
    
    /* get number of rows and columns of first input argument */
    N = mxGetM(prhs[0]); // number of data points
    M = mxGetN(prhs[0]); // number of ensemble members
    E = mxGetPr(prhs[0]);
    
    /* get number of rows and columns of second input argument */
    nCls = mxGetM(prhs[1]); 
    nCls2 = mxGetN(prhs[1]);
    if( nCls != nCls2 ) {
        mexErrMsgIdAndTxt("LCE:cts:err","Second input (WCL matrix) mot square.");
    }
    WCL = mxGetPr(prhs[1]);

    /* get number of rows and columns of third input argument */
    Crow = mxGetM(prhs[2]); 
    Ccol = mxGetN(prhs[2]);
    if( Crow != 1 ) {
        mexErrMsgIdAndTxt("LCE:cts:err","C has to be row vector.");
    }
    if( Ccol != M ) {
        mexErrMsgIdAndTxt("LCE:cts:err","C has to be row vector with M number of elements.");
    }
    C = mxGetPr(prhs[2]);
    
    /* get the value of the scalar input  */
    dc = mxGetScalar(prhs[3]);

    /* create the output matrix S */
    plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
    /* get a pointer to the real data in the output matrix */
    S = mxGetPr(plhs[0]);
    
    // mexPrintf("M: %d, N:%d, nCls: %d\n",M,N,nCls);
    
    
    //===================================================================
    // COMPUTE wCT
    //===================================================================
    wCT_mat = mxCreateDoubleMatrix(nCls,nCls,mxREAL);
    wCT = mxGetPr(wCT_mat);
    
    //omp_set_num_threads(8);
    
     wCT_max = 0;
     
    #pragma omp parallel for shared(WCL,nCls,wCT,C) private(i,q,j,WCLi,wCT_val)
    //#pragma omp parallel for
    for (i=0;i<nCls-1;i++){ // for each cluster
        //mexPrintf("Num threads %d, thread ID %d.\n", omp_get_num_threads(), omp_get_thread_num());
        // which ensemble member?
        q = 0;
        while ( (i+1) > C[q]){
            q++;
        }
        
        for (j=i+1;j<C[q];j++){ // for each other cluster
            // sum of pairwise minimum values of WCL(i,:) and WCL(j,:)
            wCT_val = 0;
            for (WCLi=0;WCLi<nCls;WCLi++){
                wCT_val += fmin(WCL[WCLi*nCls+i], WCL[WCLi*nCls+j]);
            }
            
//             if (wCT[j*nCls+i]!=0 || wCT[i*nCls+j]!=0){
//                 mexPrintf("!!! i: %d, j: %d\n",i,j);
//             }
            wCT[j*nCls+i] = wCT_val;
            wCT[i*nCls+j] = wCT_val;
        }
    }
    
    for(i=0;i<nCls*nCls;i++){
        if (wCT[i] > wCT_max){
            wCT_max = wCT[i];
        }
    }
    
//     wCT_max = 0;        
//     #pragma omp parallel for
//     for ( i=0;i<nCls*nCls;i++) {
//         if (wCT[i] > wCT_max){
//             #pragma omp critical
//             {
//                 if (wCT[i] > wCT_max){
//                     wCT_max = wCT[i];
//                 }
//             }
//         }
//     }
        
        
    
    //mexPrintf("wCT_max: %lf\n",wCT_max);
    
    if(wCT_max > 0){
       #pragma omp parallel for
       for(i=0;i<nCls*nCls;i++){
           wCT[i] = wCT[i] / wCT_max;
       } 
    }
    
    #pragma omp parallel for
    for(i=0;i<nCls;i++){
        wCT[i*nCls+i] = 1;
    }
    
//     // Display
//     mexPrintf("\nwCT matrix:\n");
//     for(i=0;i<nCls;i++){
//         for(j=0;j<nCls;j++){
//             mexPrintf("%lf ",wCT[j*nCls+i]);
//         }
//         mexPrintf("\n");
//     }
    
    //===================================================================
    // COMPUTE S
    //===================================================================
    //#pragma omp parallel for shared(M,N,E,S,wCT,dc,nCls) private(m,i,j,Ei,Ej,S_val)
    //#pragma omp parallel for
    for (m=0;m<M;m++){
        for (i=0;i<N-1;i++){
            for (j=i+1;j<N;j++){
                Ei = (long long)E[m*N+i]-1;
                Ej = (long long)E[m*N+j]-1;
                
                S_val = S[j*N+i];
                if(Ei == Ej){                    
                    S_val += 1.0/M;
                }
                else {
                    S_val += wCT[Ej*nCls+Ei]*dc/(double)M;
                }
                S[j*N+i] = S_val;
                S[i*N+j] = S_val;

            }
        }    
    }
    
    #pragma omp parallel for shared(N,S) private(i)
    for(i=0;i<N;i++){
        S[i*N+i] = 1;
    }
    
//     // Display
//     mexPrintf("\nS matrix:\n");
//     for(i=0;i<N;i++){
//         for(j=0;j<N;j++){
//             mexPrintf("%lf ",S[j*N+i]);
//         }
//         mexPrintf("\n");
//     }
    
}


