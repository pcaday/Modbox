#include "mex.h"
#include <string.h>
#include <algorithm>
#include <cmath>

using std::max;

inline bool mxIsScalar(const mxArray *mx)
{
    return (mxGetNumberOfElements(mx) == 1);
}

// Find signatures of symmetric matrices in parallel.
//
//   S = ParallelSignature(A, tol)
//   S = ParallelSignature(A, tol, 1)
//
//    A is: n x n x i1 x ... x ik, representing an i1 x ... x ik array of n x n matrices.
//    S is:         i1 x ... x ik, representing an i1 x ... x ik array of signatures.
//
//   tol is a threshold; eigenvalues with absolute value less
//    than tol will be counted as zeros. Negative tol counts as positive.
//
//   If ParallelSignature is called with a third argument, it will count
//    the number of eigenvalues greater than or equal to tol, instead of
//    the signature.
//
//   A must be real with real eigenvalues. The result is undefined if A has complex eigenvalues.
//   S will be integral (but stored as the same data type as A.)


template <class T, bool countPos>
int SignatureValue(T eigenvalue, T tol)
{
    if (countPos)
        return eigenvalue >= tol;
    else
        return (eigenvalue >= tol) - (eigenvalue <= -tol);
}


template <class T, bool countPos>
void ParallelSignature2x2(T *S, const T *A, mwSize r, T tol)
{
    for (int l=0; l<r; l++) {
        // Find roots of characteristic polynomial with the quadratic
        //  formula.
        T b = -A[0]-A[3];
		T c = A[0]*A[3]-A[1]*A[2];
        T disc = b*b-4*c;
        
        // Force the discriminant to be nonnegative in case numerical error
        //  makes it slightly negative.
        disc = max(disc, T(0.));
        disc = sqrt(disc);

        // Compute the eigenvalues and signature
        T e = 0.5*(-b+disc);
        int s = SignatureValue<T,countPos>(e,tol);
        
        e = 0.5*(-b-disc);
        s += SignatureValue<T,countPos>(e,tol);
        
        A += 4;
        *S++ = s;
    }
}

template <class T>
void ParallelSignature2x2(T *S, const T *A, mwSize r, T tol, bool countPos)
{
    if (countPos)
        ParallelSignature2x2<T,true>(S,A,r,tol);
    else
        ParallelSignature2x2<T,false>(S,A,r,tol);
}





void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
    const mxArray *aMX, *tolMX;
	mxArray *mMX;
    mwSize n;
    const mwSize *aDims;
	mwSize elementsA, aND;
	int i;
    bool countPos;
	
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    if (nrhs < 2 || nrhs > 3) mexErrMsgTxt("Expected 2 or 3 input arguments.");
    
    aMX = prhs[0];
    tolMX = prhs[1];
	countPos = (nrhs == 3);     // If three arguments given, count number
                                //  of eigenvalues >= tol.
    
	/*
	 * Analyze data type of inputs
	 */
	bool is_single, is_complex;
	
    if (mxIsDouble(aMX))
		is_single = false;
    else if (mxIsSingle(aMX))
		is_single = true;
    else mexErrMsgTxt("The first argument must be single or double.");
	
    if (!mxIsNumeric(tolMX))
        mexErrMsgTxt("The tolerance must be a numeric value.");
    
	if (mxIsComplex(aMX) || mxIsComplex(tolMX))
		//		is_complex = true;
		mexErrMsgTxt("Inputs must be real.");
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

    if (!mxIsScalar(tolMX))
        mexErrMsgTxt("The tolerance must be a scalar.");
    
    /*
     * Create output array
     */
    int oND = max(aND-2,2);
    int* oDims = new int[oND];

    oDims[0] = 1;
    oDims[1] = 1;
    for (i=2;i<aND;i++)
        oDims[i-2] = aDims[i];
    
    plhs[0] = mMX = mxCreateNumericArray(oND, oDims, is_single ? mxSINGLE_CLASS : mxDOUBLE_CLASS,
										 mxREAL);
	
    delete[] oDims;
    
	
	/*
	 * Actually do the work.
	 */
	
	switch (n) {
		case 2:
			if (is_single) {
				ParallelSignature2x2((float *) mxGetData(mMX),
                                     (float *) mxGetData(aMX),
								     elementsA,
                                     (float) mxGetScalar(tolMX),
                                     countPos);
			} else {
				ParallelSignature2x2((double *) mxGetData(mMX),
                                     (double *) mxGetData(aMX),
									 elementsA,
                                     (double) mxGetScalar(tolMX),
                                     countPos);
			}
			break;
		default:
			mexErrMsgTxt("Currently only 2x2 matrices supported.");
	}
}