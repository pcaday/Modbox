#include "mex.h"
#include <string.h>

// Invert matrices in parallel.
//
//   M = ParallelMatrixInverse(A)
//
//  A is:     n x n x i1 x ... x ik, representing an i1 x ... x ik array of n x n matrices.
//

template <class T>
        void ParallelMatrixInverse2x2(T *M, const T *A, mwSize r)
{
    int l;
    T d;
    
    for (l=0; l<r; l++) {
		d = 1 / (A[0]*A[3]-A[1]*A[2]);
		
        M[0] = A[3]*d;
		M[1] = -A[1]*d;
		M[2] = -A[2]*d;
		M[3] = A[0]*d;
		
        A += 4; M += 4;
    }
}


template <class T>
void ParallelMatrixInverse3x3(T *M, const T *A, mwSize r)
{
    int l;
    T d;
    
	// 0 1 2
	// 3 4 5
	// 6 7 8
    for (l=0; l<r; l++) {
		// Using Cramer's rule...
		
		d = 1 / (A[0]*(A[4]*A[8]-A[5]*A[7])
			    +A[1]*(A[5]*A[6]-A[3]*A[8])
			    +A[2]*(A[3]*A[7]-A[4]*A[6]));
		
		M[0] = +d*(A[4]*A[8]-A[5]*A[7]);
		M[1] = -d*(A[1]*A[8]-A[2]*A[7]);
		M[2] = +d*(A[1]*A[5]-A[2]*A[4]);
		M[3] = -d*(A[3]*A[8]-A[5]*A[6]);
		M[4] = +d*(A[0]*A[8]-A[2]*A[6]);
		M[5] = -d*(A[0]*A[5]-A[2]*A[3]);
		M[6] = +d*(A[3]*A[7]-A[4]*A[6]);
		M[7] = -d*(A[0]*A[7]-A[1]*A[6]);
		M[8] = +d*(A[0]*A[4]-A[1]*A[3]);
		
        A += 9; M += 9;
    }
}




void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
    const mxArray *aMX;
	mxArray *mMX;
    mwSize n;
    const mwSize *aDims;
	mwSize elementsA, aND;
	int i;
	
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    if (nrhs != 1) mexErrMsgTxt("Expected 1 input argument.");
    
    aMX = prhs[0];

	
	/*
	 * Analyze data type of inputs
	 */
	bool is_single, is_complex;
	
    if (mxIsDouble(aMX))
		is_single = false;
    else if (mxIsSingle(aMX))
		is_single = true;
	else mexErrMsgTxt("The inputs must be single or double.");
	
	if (mxIsComplex(aMX))
		//		is_complex = true;
		mexErrMsgTxt("Complex inputs not supported");
    else
		is_complex = false;

	
	/*
	 * Verify the dimensions of the inputs
	 */
	aDims = mxGetDimensions(aMX);	
	aND = mxGetNumberOfDimensions(aMX);
	
    n = aDims[0];
	if (n != aDims[1]) mexErrMsgTxt("The matrices must be square.");

	elementsA = 1;
	for (i = 2; i < aND; i++) elementsA *= aDims[i];

	plhs[0] = mMX = mxCreateNumericArray(aND, aDims, is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
										 is_complex ? mxCOMPLEX : mxREAL);
	
	/*
	 * Actually do the work :)
	 */
	
	switch (n) {
		case 2:
			if (is_single) {
				ParallelMatrixInverse2x2((float *) mxGetData(mMX),
										 (float *) mxGetData(aMX),
										 elementsA);
			} else {
				ParallelMatrixInverse2x2((double *) mxGetData(mMX),
										 (double *) mxGetData(aMX),
										 elementsA);
			}
			break;
		case 3:
			if (is_single) {
				ParallelMatrixInverse3x3((float *) mxGetData(mMX),
										 (float *) mxGetData(aMX),
										 elementsA);
			} else {
				ParallelMatrixInverse3x3((double *) mxGetData(mMX),
										 (double *) mxGetData(aMX),
										 elementsA);
			}
			break;
		default:
			mexErrMsgTxt("Currently only 2x2 and 3x3 matrices supported.");
	}
}