#include "mex.h"
#include <string.h>
#include <algorithm>

using std::max;
using std::min;

// Multiply matrices in parallel.
//
//   M = ParallelMatrixMultiply(A, B)
//
//  A is:     m x n x i1 x ... x ik, representing an i1 x ... x ik array of m x n matrices
//  B is:     n x p x i1 x ... x ik, representing an i1 x ... x ik array of n x p matrices
//

template <class T>
        void ParallelMatrixMultiply(T *M, const T *A, const T *B, mwIndex strideA, mwIndex strideB, mwSize m, mwSize n, mwSize p, mwSize r)
{
    int i, j, k, l;
    T v;
    
    for (l=0; l<r; l++) {
        for (j=0; j<p; j++) {
            for (i=0; i<m; i++) {
                v = 0;
                for (k=0; k<n; k++) {
                    v += A[k*m+i]*B[j*n+k];
                }
                
                *M++ = v;
            }
        }
        
        A += strideA;
        B += strideB;
    }
}



void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
    const mxArray *aMX, *bMX;
	mxArray *mMX;
    mwSize m, n, p, r;
    const mwSize *aDims, *bDims;
	mwSize *mDims;
	mwSize elementsA, elementsB, aND, bND, mND;
	mwIndex strideA, strideB;
	int i;
	
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    if (nrhs != 2) mexErrMsgTxt("Expected 2 input arguments.");
    
    aMX = prhs[0];
    bMX = prhs[1];

	
	/*
	 * Analyze data type of inputs
	 */
	bool is_single, is_complex;
	
    if (mxIsDouble(aMX) && mxIsDouble(bMX))
		is_single = false;
    else if (mxIsSingle(aMX) && mxIsSingle(bMX))
		is_single = true;
	else mexErrMsgTxt("The inputs must be both single or both double.");
	
	if (mxIsComplex(aMX) && mxIsComplex(bMX))
		//		is_complex = true;
		mexErrMsgTxt("Complex inputs not supported");
    else if (~mxIsComplex(aMX) && ~mxIsComplex(bMX))
		is_complex = false;
	else mexErrMsgTxt("The inputs must be both real or both complex.");

	
	/*
	 * Verify the dimensions of the inputs
	 */
	aDims = mxGetDimensions(aMX);
	bDims = mxGetDimensions(bMX);
	
    m = aDims[0];
    n = aDims[1];
	p = bDims[1];

	if (n != bDims[0]) mexErrMsgTxt("The matrices have incompatible sizes.");

	aND = mxGetNumberOfDimensions(aMX);
	bND = mxGetNumberOfDimensions(bMX);

	elementsA = elementsB = 1;
	for (i = 2; i < aND; i++) elementsA *= aDims[i];
	for (i = 2; i < bND; i++) elementsB *= bDims[i];
	
	if (elementsA != elementsB && elementsA != 1 && elementsB != 1)
		mexErrMsgTxt("The inputs do not have the same number of matrices.");
	r = max(elementsA, elementsB);
	
	strideA = (elementsA == 1) ? 0 : m*n;
	strideB = (elementsB == 1) ? 0 : n*p;
	
	// Find which of the arguments has more dimensions, and use its size for the output
	//  matrix.  
	mND = (bND > aND) ? bND : aND;
	mDims = (mwSize *) mxMalloc(sizeof(mwSize) * mND);
	mDims[0] = m; mDims[1] = p;
	memcpy(mDims+2, ((bND > aND) ? bDims : aDims)+2, sizeof(mwSize) * (mND-2));
	
	
	plhs[0] = mMX = mxCreateNumericArray(mND, mDims, is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
										 is_complex ? mxCOMPLEX : mxREAL);
	
	/*
	 * Actually do the work :)
	 */
	
	if (is_single) {
		ParallelMatrixMultiply((float *) mxGetData(mMX),
							   (float *) mxGetData(aMX),
							   (float *) mxGetData(bMX),
							   strideA, strideB,
							   m, n, p, r);
    } else {
		ParallelMatrixMultiply((double *) mxGetData(mMX),
							   (double *) mxGetData(aMX),
							   (double *) mxGetData(bMX),
							   strideA, strideB,
							   m, n, p, r);
	}
	
	mxFree(mDims);
}