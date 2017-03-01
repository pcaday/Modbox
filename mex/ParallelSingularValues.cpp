#include "mex.h"
#include <string.h>
#include <algorithm>
#include <cmath>

using std::max;

// Find singular values in parallel.
//
//   E = ParallelSingularValues(A)
//
//  A is:     n x n x i1 x ... x ik, representing an i1 x ... x ik array of n x n matrices.
//  E is:         n x i1 x ... x ik, representing an i1 x ... x ik array of length-n lists of eigenvalues.
//
//   A must be real; E will be real.


template <class T>
        void ParallelSV2x2(T *Ereal, const T *A, mwSize r)
{
    for (int l=0; l<r; l++) {
        // Find roots of characteristic polynomial of A*A
        //  using the quadratic formula.
        T tr  = A[0]*A[0] + A[1]*A[1] + A[2]*A[2] + A[3]*A[3];
		T det = A[0]*A[3]-A[1]*A[2];
        det *= det;
        
        T disc = sqrt(tr*tr-4*det);
        
        // Output singular values in decreasing order, like MATLAB's
        //  svd function.
        *Ereal++ = sqrt(0.5*(tr+disc));
        *Ereal++ = sqrt(0.5*(tr-disc));

        A += 4;
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

    int oND = max(aND-1,2);
    int* oDims = new int[oND];

    oDims[0] = 1;
    oDims[1] = 1;
    for (i=1;i<aND;i++)
        oDims[i-1] = aDims[i];
    
    plhs[0] = mMX = mxCreateNumericArray(oND, oDims,
                                         is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
										 mxREAL);
	
    delete[] oDims;
    
	
	/*
	 * Actually do the work.
	 */
	
	switch (n) {
		case 2:
			if (is_single) {
				ParallelSV2x2((float *) mxGetData(mMX),
						      (float *) mxGetData(aMX),
										elementsA);
			} else {
				ParallelSV2x2((double *) mxGetData(mMX),
                              (double *) mxGetData(aMX),
										elementsA);
			}
			break;
		default:
			mexErrMsgTxt("Currently only 2x2 matrices supported.");
	}
}