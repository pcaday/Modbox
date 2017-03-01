#include "mex.h"
#include <string.h>
#include <algorithm>

using std::max;

// Take determinants in parallel.
//
//   M = ParallelDet(A)
//
//  A is:     n x n x i1 x ... x ik, representing an i1 x ... x ik array of n x n matrices.
//

template <class T>
        void ParallelDet2x2(T *D, const T *A, mwSize r)
{
    for (int l=0; l<r; l++) {
		*D++ = A[0]*A[3]-A[1]*A[2];
		
        A += 4;
    }
}


template <class T>
void ParallelDet3x3(T *D, const T *A, mwSize r)
{
    // 0 1 2
	// 3 4 5
	// 6 7 8
    for (int l=0; l<r; l++) {
		*D++ = A[0]*(A[4]*A[8]-A[5]*A[7])
			  +A[1]*(A[5]*A[6]-A[3]*A[8])
			  +A[2]*(A[3]*A[7]-A[4]*A[6]);
        A += 9;
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

    int oND = max(aND-2,2);
    int* oDims = new int[oND];

    oDims[0] = 1;
    oDims[1] = 1;
    for (i=0;i<aND-2;i++)
        oDims[i] = aDims[i+2];
    
    plhs[0] = mMX = mxCreateNumericArray(oND, oDims, is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
										 is_complex ? mxCOMPLEX : mxREAL);
	
    delete[] oDims;
    
	
	/*
	 * Actually do the work.
	 */
	
	switch (n) {
		case 2:
			if (is_single) {
				ParallelDet2x2((float *) mxGetData(mMX),
							   (float *) mxGetData(aMX),
										 elementsA);
			} else {
				ParallelDet2x2((double *) mxGetData(mMX),
					   		   (double *) mxGetData(aMX),
										 elementsA);
			}
			break;
		case 3:
			if (is_single) {
				ParallelDet3x3((float *) mxGetData(mMX),
							   (float *) mxGetData(aMX),
										 elementsA);
			} else {
				ParallelDet3x3((double *) mxGetData(mMX),
					   		   (double *) mxGetData(aMX),
										 elementsA);
			}
			break;
		default:
			mexErrMsgTxt("Currently only 2x2 and 3x3 matrices supported.");
	}
}