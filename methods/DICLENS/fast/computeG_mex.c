
#include "mex.h"
#include "math.h"
#include "omp.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Inputs
	double *P, *AM, *sep, *nElemInCls;
	long long M, N;
	long long nClust, maxElem; 
	double memoryLimitGB;
	//long long maxElemInCls;
	// Outputs
	double *GSi, *GSj, *GSv;
	// Internal
	long long i,j,p,Cik,Cjl, cLen1, cLen2, P_i, P_j;
	long long numPair, nElem, Gidx, nEdges, edgeIdx;
	double ECS, maxG, simVal;
	mxArray *Gi_mat, *Gj_mat, *Gv_mat;
	double *Gi, *Gj, *Gv;
	double growRatio;
	long long requiredMemBytes, limitMem, sizeG, growSize, newSizeBytes;
	double *newptr;
	int PARALLEL_MODE;

	/* Check for proper number of input and output arguments */
	if (nrhs != 5) {
		mexErrMsgTxt("5 inputs argument required.");
	}
	if(nlhs != 3){
		mexErrMsgTxt("3 outputs required.");
	}

	// P matrix [nClust x maxElem]
	maxElem = mxGetN(prhs[0]); // max. number of data points in cluster
	nClust = mxGetM(prhs[0]); // number of ensemble members
	P = mxGetPr(prhs[0]);

	// AM matrix [N X N]
	N = mxGetM(prhs[1]); // number of data points
	AM = mxGetPr(prhs[1]);

	// Vector - max index in each cluster, first is 0
	M = mxGetN(prhs[2])-1;
	sep = mxGetPr(prhs[2]);

	// Vector - number of elements in each cluster
	nElemInCls = mxGetPr(prhs[3]);

	// Get size of memory in GB - if less memory required, run parallel
	memoryLimitGB = mxGetScalar(prhs[4]); 

	// number of pairs - elements in graph G (triangle)
	numPair = nClust*(nClust-1)/2;

	// Create vectors for graph G

	// Consider memory consumption
	// If numPair results in less than 1 GB, allocate all of it.
	requiredMemBytes = numPair*3*sizeof(double); // Gi, Gj, Gv each takes numPair*sizeof(double)
	limitMem = (long long)(memoryLimitGB * 1024*1024*1024); // Bytes
	growRatio = 0.1; // growth ratio -> sizeG * growRatio 

	PARALLEL_MODE = 1;

	if(requiredMemBytes > limitMem){
		sizeG = (long long)floor(limitMem/(3.0*sizeof(double)));
		PARALLEL_MODE = 0;
	}
	else{
		sizeG = numPair;
	}
	growSize = (long long)ceil((double)sizeG*growRatio);

	Gi_mat = mxCreateDoubleMatrix(sizeG,1,mxREAL);
	Gj_mat = mxCreateDoubleMatrix(sizeG,1,mxREAL);
	Gv_mat = mxCreateDoubleMatrix(sizeG,1,mxREAL);

	Gi = mxGetPr(Gi_mat);
	Gj = mxGetPr(Gj_mat);
	Gv = mxGetPr(Gv_mat);



	//====================================================================
	// PARALLEL FULL MEMORY MODE
	//====================================================================
	if(PARALLEL_MODE){
	#pragma omp parallel for private(p,Cik,Cjl,cLen1,cLen2,nElem,ECS,P_i,P_j,i,j,Gidx) shared(M,sep,nClust,nElemInCls,N,P,AM,Gi,Gj,Gv)
		for (p=0; p < M; p++){
			// Compute indeces of pairs for this partition
			for (Cik=(long long)sep[p]; Cik < sep[p+1]; Cik++ ){ // pth partition
				cLen1 = (long long)nElemInCls[Cik];
				for (Cjl=(long long)sep[p+1]; Cjl < nClust; Cjl++ ){ // other partitions
					cLen2 = (long long)nElemInCls[Cjl];
					nElem = cLen1 * cLen2;                
					if (nElem >0) {
						ECS = 0;
						for(i=0;i<cLen1;i++){
							P_i = (long long)P[Cik+i*nClust];
							for(j=0;j<cLen2;j++){
								P_j = (long long)P[Cjl+j*nClust];
								ECS = ECS + AM[P_j * N + P_i];
							}						
						}
						if(ECS>0){
							Gidx = (long long)(Cjl-1 + Cik*(nClust - Cik/2.0 -1.5));
							Gi[Gidx] = (double)Cik+1;
							Gj[Gidx] = (double)Cjl+1;
							Gv[Gidx] = (double)ECS/(double)nElem;
						}
					}
				}
			}
		}
		// Get max value of G and count non-zero elements, i.e. number of edges
		nEdges = 0;
		maxG = 0;
		for (i=0;i<numPair;i++){
			if (Gv[i]){
				if(Gv[i]>maxG){
					maxG = Gv[i];
				}
				nEdges ++;
			}
		}
		maxG +=1;

		// Create vectors for graph GS (similarity)
		plhs[0] = mxCreateDoubleMatrix(2*nEdges,1,mxREAL);
		plhs[1] = mxCreateDoubleMatrix(2*nEdges,1,mxREAL);
		plhs[2] = mxCreateDoubleMatrix(2*nEdges,1,mxREAL);

		GSi = mxGetPr(plhs[0]);
		GSj = mxGetPr(plhs[1]);
		GSv = mxGetPr(plhs[2]);

		edgeIdx = 0;
		for (i=0;i<numPair;i++){
			if (Gv[i] != 0){
				simVal = maxG - Gv[i];
				GSi[edgeIdx] = Gi[i];
				GSj[edgeIdx] = Gj[i];
				GSv[edgeIdx] = simVal;
				GSi[edgeIdx+1] = Gj[i];
				GSj[edgeIdx+1] = Gi[i];
				GSv[edgeIdx+1] = simVal;
				edgeIdx +=2;
			}
		}
	}
	else{

		//====================================================================
		// SERIAL MEMORY SAVING MODE
		//====================================================================

		nEdges = 0;
		maxG = 0;
		// CANNOT be parallel!
		for (p=0; p < M; p++){
			// Compute indeces of pairs for this partition
			for (Cik=(long long)sep[p]; Cik < sep[p+1]; Cik++ ){ // pth partition
				cLen1 = (long long)nElemInCls[Cik];
				for (Cjl=(long long)sep[p+1]; Cjl < nClust; Cjl++ ){ // other partitions
					cLen2 = (long long)nElemInCls[Cjl];
					nElem = cLen1 * cLen2;                
					if (nElem >0) {
						ECS = 0;
						for(i=0;i<cLen1;i++){
							P_i = (long long)P[Cik+i*nClust];
							for(j=0;j<cLen2;j++){
								P_j = (long long)P[Cjl+j*nClust];
								ECS = ECS + AM[P_j * N + P_i];
							}						
						}
						if(ECS>0){

							// add some space to arrays
							if(nEdges >= sizeG){
								sizeG += growSize;
								newSizeBytes = sizeG*sizeof(double);

								newptr = (double *)mxRealloc(Gi, newSizeBytes); //reallocate array
								mxSetPr(Gi_mat,newptr);
								mxSetM(Gi_mat,sizeG);//fix dimension
								Gi = newptr;

								newptr = (double *)mxRealloc(Gj, newSizeBytes); //reallocate array
								mxSetPr(Gj_mat,newptr);
								mxSetM(Gj_mat,sizeG);//fix dimension
								Gj = newptr;

								newptr = (double *)mxRealloc(Gv, newSizeBytes); //reallocate array
								mxSetPr(Gv_mat,newptr);
								mxSetM(Gv_mat,sizeG);//fix dimension
								Gv = newptr;
							}

							Gi[nEdges] = (double)Cik+1;
							Gj[nEdges] = (double)Cjl+1;
							Gv[nEdges] = (double)ECS/(double)nElem;
							// Get max value of G and count non-zero elements, i.e. number of edges
							if(Gv[nEdges]>maxG){
								maxG = Gv[nEdges];
							}
							nEdges ++;
						}
					}
				}
			}
		}
		maxG +=1;

		// Create vectors for graph GS (similarity)
		plhs[0] = mxCreateDoubleMatrix(2*nEdges,1,mxREAL);
		plhs[1] = mxCreateDoubleMatrix(2*nEdges,1,mxREAL);
		plhs[2] = mxCreateDoubleMatrix(2*nEdges,1,mxREAL);

		GSi = mxGetPr(plhs[0]);
		GSj = mxGetPr(plhs[1]);
		GSv = mxGetPr(plhs[2]);

		edgeIdx = 0;
		for (i=0;i<nEdges;i++){
			simVal = maxG - Gv[i];
			GSi[i] = Gi[i];
			GSj[i] = Gj[i];
			GSv[i] = simVal;
			GSi[i+nEdges] = Gj[i];
			GSj[i+nEdges] = Gi[i];
			GSv[i+nEdges] = simVal;		
		}
	}

	//mexPrintf("PARALLEL_MODE: %d, pairs: %d, edges: %d, ratio: %lf\n",PARALLEL_MODE, numPair, nEdges, nEdges/(double)numPair);
	
	mxDestroyArray(Gi_mat);
	mxDestroyArray(Gj_mat);
	mxDestroyArray(Gv_mat);
}


