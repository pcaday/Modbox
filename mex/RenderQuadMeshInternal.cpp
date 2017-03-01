#include "mex.h"
#include <cmath>
#include <algorithm>

using namespace std;

#define DEBUG 0

#if DEBUG
    #define mexDebugPrintf mexPrintf
#else
    #define mexDebugPrintf(...)
#endif

#define pi 3.141592653589793238
    
// Controls how smooth the edges are (radius of filter).
static const double Fuzziness = 1;

// For twisted quadrilaterals, this is the smallest we will allow the
//  cross product of the normals to the two twisted edges to be.
// If the
static const double DetThreshold = 0.000001;

// For triangles, if the projection of a vector from one edge to the
//  opposite vertex onto the normal to that edge has length less than this
//  threshold, the three points are considered collinear and the triangle
//  discarded.
static const double CollinearThreshold = 0.000001;

// For polygons, if two vertices are within this distance of each other,
//  they're considered to be in the same location.
static const double CoincidenceThreshold = 0.000001;

// Maximum distance a vertex can be from the pixel center and still be
//  considered at the center
static const double CenterThreshold = 0.00001;

  
// Non-negative-valued modulo, by jthill:
//   http://stackoverflow.com/questions/14997165/fastest-way-to-get-a-positive-modulo-in-c-c?lq=1
//
// Like the remainder (%) operator but returns a nonnegative modulus 
//  even when i is negative.
inline int positive_mod(int i, int n)
{
    const int shift = CHAR_BIT*(sizeof i) - 1;
    int m = i%n;
    return m + (m>>shift & n);
}

// arccos/arcsin/sqrt that automatically select between acos/acosf and asin/asinf
//  depending on data type
template <class T>
inline T arccos(T v)
{
	return acos(v);
}

template <>
inline float arccos<float>(float v)
{
	return acosf(v);
}


template <class T>
inline T arcsin(T v)
{
	return asin(v);
}

template <>
inline float arcsin<float>(float v)
{
	return asinf(v);
}

template <class T>
inline T sqroot(T v)
{
	return sqrt(v);
}

template <>
inline float sqroot<float>(float v)
{
	return sqrtf(v);
}



inline bool mxIsScalar(const mxArray *a)
{
    return mxGetM(a)*mxGetN(a) == 1;
}
  
#define min(a,b) (((a) < (b)) ? (a) : (b));
#define max(a,b) (((a) > (b)) ? (a) : (b));

template<class T>
        inline T EdgeCutoff(T x)
{
    if (x <= -1.) return 1;
    if (x >= 1.) return 0;
	return (arccos(x) - x*sqroot(1-x*x)) * (1/pi);
}

// Approximate the integral of the filter over a wedge shaped region
//  with center at distance d from the filter's center, starting from the
//  line through the filter center and the wedge center and with an 
//  internal angle arccos(b).
template<class T>
        inline T CornerCutoff(T d, T b)
{
    T a, e, f, theta, gamma;
    bool reverse;
	
	reverse = (b < 0);
    if (reverse) {
		d = -d; b = -b;
	}
	
    b = min(b, T(1));
	theta = arccos(b);
	gamma = asin(d*sqrt(1-b*b));
	
	e = cos(theta-gamma);
	f = abs(sin(theta-gamma));
	mexDebugPrintf("CC: %f,%f,%f,%f,%f,%f,%d\n", d,b,e,f,theta,gamma,reverse);
	a = (0.5/pi)*((e-d)*f + arccos(e) - e*sqrt(1-e*e));
    
	return reverse ? (0.5-a) : a;
	
}

// Inline function which adds a weighted contribution to a pixel in an image,
//  used by RenderTriangle
//
// Parameters (all inputs except imageRe/imageIm are input/output)
//
//          (x1,x2): Coordinates of pixel (1-based)
//             mask: Mask from triangle
//             m, n: Dimensions of image
// imageRe, imageIm: Real and imaginary parts of image. imageIm may be NULL.
//     valRe, valIm: Real and imaginary parts of the value of this triangle.
//         wRe, wIm: Real and imaginary parts of the weighting function. Either or both may be NULL.
//         periodic: Flag, true if the image should be treated as periodic.
//
template <class T>
inline void ContributePixel(int x1, int x2, T mask, int m, int n, T *imageRe, T *imageIm,
							T valRe, T valIm, const T *wRe, const T *wIm, bool periodic)
{
	// Switch from 1-based coordinates (x1,x2) to 0-based coordinates (im,jm);
	int i0 = x1 - 1;
	int j0 = x2 - 1;
	
	// For periodic grids, reduce coordinates modulo grid size.
	// Use positive_mod for correct results if i0,j0 happen to be negative.
	if (periodic) {
		i0 = positive_mod(i0, m);
		j0 = positive_mod(j0, n);
	}
	
	// Multiply contribution to this pixel by weight function
	//  (which may be complex) and the quadrilateral's value
	//  (which may also be complex), and add the result to the
	//  image.
	if (!wRe && !wIm) {
		// Fast track, if no weighting function:
		imageRe[j0*m+i0] += mask * valRe;
		if (imageIm)
			imageIm[j0*m+i0] += mask * valIm;
	} else {
		// Slow track: Multiply (wRe + wIm*i) by (valRe + valIm*i)
		T wr = wRe ? wRe[j0*m+i0] : 1.;
		T wi = wIm ? wIm[j0*m+i0] : 0.;
		
		imageRe[j0*m+i0] += mask * (wr * valRe - wi * valIm);
		if (imageIm)
			imageIm[j0*m+i0] += mask * (wr * valIm + wi * valRe);
	}
}

template <class T, int N>
void RenderConvexPolygon(T *imageRe, T *imageIm, T x1[N], T x2[N], 
						 T valRe, T valIm, const T *wRe, const T *wIm, int m, int n, bool periodic)
{
	T n1[N], n2[N], c[N], edgeposv1[N], edgeposv2[N], edgeposc[N];
	T triarea;
	
	// For triangles:
	// Take cross product between 0-1 edge and 0-2 edge to get the area
	//  (used later).
	if (N == 3)
		triarea = 0.5 * abs((x1[1]-x1[0])*(x2[2]-x2[0]) - (x2[1]-x2[0])*(x1[2]-x1[0]));

	// Initialize the bounding box
	int min_x1 = INT_MAX, min_x2 = INT_MAX;
	int max_x1 = INT_MIN, max_x2 = INT_MIN;

	// Report on vertices
	for (int l = 0; l < N; l++)
	{
		mexDebugPrintf("Vertex %d: (%f,%f)\n", l, x1[l], x2[l]);
	}
	
	// Round up the normals.
	for (int l = 0; l < N; l++)
	{
		T d1 = x1[(l+1)%N] - x1[l];
        T d2 = x2[(l+1)%N] - x2[l];
        T sqedgelen = d1*d1 + d2*d2;
        
		if (sqedgelen <= CoincidenceThreshold*CoincidenceThreshold) {
            // Two vertices in the same location
            mexDebugPrintf("Coincident vertices: %d, %d\n", l, (l+1)%N);
			if (N > 3) {
				// Eliminate the duplicate vertex.
				for (int k = l+1; k < N-1; k++) {
					x1[k] = x1[k+1];
					x2[k] = x2[k+1];
				}
				// Render the resulting (N-1)-gon
				RenderConvexPolygon<T,N-1>(imageRe, imageIm, x1, x2, valRe, valIm, wRe, wIm, m, n, periodic);
			}
            return;
        }
		
		// Get the affine functional representing projected normalized distance along this edge.
		//  With these functionals, one vertex is at 0 and the other is at 1. 
		T isqedgelen = 1 / sqedgelen;
		edgeposv1[l] = d1 * isqedgelen;
		edgeposv2[l] = d2 * isqedgelen;
		edgeposc[l]  = x1[l]*edgeposv1[l] + x2[l]*edgeposv2[l];
		
		// Get the affine functional representing distance from the edge, divided by the
		//  filter radius (Fuzziness)
		T iedgelen = sqrt(isqedgelen) * (1 / Fuzziness);
		n1[l] = -d2 * iedgelen;
		n2[l] = +d1 * iedgelen;
		c[l] = n1[l]*x1[l] + n2[l]*x2[l];
		
		// Adjust the bounding box.
		min_x1 = min(min_x1, (int)  ceil(x1[l]-Fuzziness));
        min_x2 = min(min_x2, (int)  ceil(x2[l]-Fuzziness));
        max_x1 = max(max_x1, (int) floor(x1[l]+Fuzziness));
        max_x2 = max(max_x2, (int) floor(x2[l]+Fuzziness));
		
		mexDebugPrintf("Edge halfplane  %d: (%f,%f).x - %f\n", l, n1[l], n2[l], c[l]);
		mexDebugPrintf("Edge positioner %d: (%f,%f).x - %f\n", l, edgeposv1[l], edgeposv2[l], edgeposc[l]);
	}
	
	// Check if the normals are pointing inwards by seeing if vertex 2 is on the
	//  inside of the halfplane for the 0-1 edge. If not, flip all normals.
	T fval = n1[0]*x1[2] + n2[0]*x2[2] - c[0];
	if (fval < 0) {
		mexDebugPrintf("Flipping normals.\n");
		// Wrong way!
		for (int l = 0; l < N; l++) {
			n1[l] = -n1[l];
			n2[l] = -n2[l];
			c[l] = -c[l];
		}
	}
	// Record the orientation too for corner work later.
	T orientation = (fval < 0) ? -1 : 1;
	
	// Clip the bounding box to the image domain if not periodic
	if (!periodic) {
        min_x1 = max(min_x1,1);   min_x2 = max(min_x2,1);
        max_x1 = min(max_x1,m);   max_x2 = min(max_x2,n);
    }
	
	// OK, now loop through the pixels!
	for (int i = min_x1; i <= max_x1; i++) {
		T edgedist[N];
		T edgepos[N];
		
		// Restart the edgedist, edgepos variables to prevent too much accumulated
		//  numerical error.
		for (int l = 0; l < N; l++) {
			edgedist[l] = n1[l]*i + n2[l]*min_x2 - c[l];
			edgepos[l] = edgeposv1[l]*i + edgeposv2[l]*min_x2 - edgeposc[l];
		}
		
        for (int j = min_x2; j <= max_x2; j++) {
			T mask = 1;
			bool has_corner[N];
			
			mexDebugPrintf("Pixel (%d,%d): edgedists %f, %f, %f, %f, edgepos %f, %f, %f, %f\n", i, j, edgedist[0], edgedist[1], edgedist[2], edgedist[3], edgepos[0], edgepos[1], edgepos[2], edgepos[3]);
			// First check: does the filter overlap the triangle at all?
			//  or is the filter completely inside the triangle?
			bool all_in = true;
			for (int l = 0; l < N; l++) {
				if (edgedist[l] <= -1.) goto nextpixel;		// If we're completely out of the triangle, skip to the next pixel.
															// Yes, that's a goto statement.
				all_in &= (edgedist[l] >= 1.);
			}
			
			if (!all_in) {
				// Second check: corners?
				for (int l = 0; l < N; l++) {
					T a1 = x1[l] - i;
					T a2 = x2[l] - j;
					
					T corner_dist = a1*a1 + a2*a2;
					
					has_corner[l] = (corner_dist < Fuzziness*Fuzziness);
					
					if (has_corner[l]) {
						// Compute corner contribution.
						T b1 = x1[(l+1)%N]   - x1[l];
						T b2 = x2[(l+1)%N]   - x2[l];
						T c1 = x1[(l+N-1)%N] - x1[l];
						T c2 = x2[(l+N-1)%N] - x2[l];
						
						T bnorm_sq = b1*b1 + b2*b2;
						T cnorm_sq = c1*c1 + c2*c2;
						T b_dp = a1*b1 + a2*b2;
						T c_dp = a1*c1 + a2*c2;
						
						T tb = (b_dp + sqrt(b_dp*b_dp + bnorm_sq*(Fuzziness*Fuzziness - corner_dist))) / bnorm_sq;
						T tc = (c_dp + sqrt(c_dp*c_dp + cnorm_sq*(Fuzziness*Fuzziness - corner_dist))) / cnorm_sq;
						
						// Compute area of triangle between vertices a, b, c.
						//  Already done for triangles.
						if (N != 3)
							triarea = 0.5 * abs(b1*c2 - b2*c1);
						
						// Triangular contribution to corner area
						T corner_tri = tb * tc * triarea * (1. / (pi*Fuzziness*Fuzziness));
						
						// Now compute arc contribution. Find boundary points:
						T bb1 = a1 - tb*b1;
						T bb2 = a2 - tb*b2;
						T bc1 = a1 - tc*c1;
						T bc2 = a2 - tc*c2;
						
						// Take dot product: cosine of angle between them
						// Note b is the edge going to vertex l+1;
						//      c is the edge from vertex l.
						T dp = (bb1*bc1 + bb2*bc2) * (1. / (Fuzziness*Fuzziness));
						
						// Also take cross product, for its sign.
						T cp = (bb1*bc2 - bb2*bc1) * orientation;
						
						// Get cosine of half the angle between them, using half-angle formula:
						T dpm = sqrt((1 + dp) * 0.5);
						
						// Adjust correctly in case the original angle was >= 180 degrees.
						if (cp < 0)
							dpm = -dpm;
						
						// Then we get the arc contribution.
						T corner_arc = EdgeCutoff(dpm);
						
						mexDebugPrintf(" Corner %d: triangle (%f) + arc (%f)\n", l, corner_tri, corner_arc);
						mexDebugPrintf("           bpoints (%f,%f), (%f,%f)\n", bb1, bb2, bc1, bc2);
						mexDebugPrintf("           disp    (%f,%f), (%f,%f)\n", b1, b2, c1, c2);
						mexDebugPrintf("           dist    %f, %f\n", tb, tc);
						mexDebugPrintf("           a       (%f,%f)\n", a1, a2);
						mexDebugPrintf("           cp(%c)   %f\n", orientation == 1 ? '+' : '-', cp);
						mexDebugPrintf("           dp      %f -> %f\n", dp, dpm);
						
						
						mask += corner_arc + corner_tri;
					}
				}
				
				// Edge contributions.
				for (int l = 0; l < N; l++) {
					if (edgedist[l] < 1) {
						if ((edgepos[l] >= 0 && edgepos[l] <= 1)		// Remove this edge contribution if it comes
							|| has_corner[l] || has_corner[(l+1)%N]) {  // from an actual edge of the poly or if one
																		// of its vertices is in sight.
							T econtrib = EdgeCutoff(edgedist[l]);
							mask -= econtrib;
							
							mexDebugPrintf(" Edge %d: removing %f\n", l, econtrib);
						}
					}
				}

				// If the mask is untouched, only false edges intersect the filter
				//  so the filter does not intersect the triangle. (We already know the
				//  filter's not completely inside the triangle).
				if (mask == 1) goto nextpixel;
			}
			
			// Add in this pixel.
			mexDebugPrintf(" Mask: %f\n", mask);
			if (mask > 0)
				ContributePixel(i, j, mask, m, n, imageRe, imageIm, valRe, valIm, wRe, wIm, periodic);
			
			
		nextpixel:
			// Update the edgedist, edgepos variables for the next round
			for (int l = 0; l < N; l++) {
				edgedist[l] += n2[l];
				edgepos[l] += edgeposv2[l];				
			}
		}
	}
}

// No partial specialization...

template <>
void RenderConvexPolygon<double,2>(double *imageRe, double *imageIm, double x1[2], double x2[2], 
								   double valRe, double valIm, const double *wRe, const double *wIm, int m, int n, bool periodic)
{
	// Two-gons are invisible :)
}

template <>
void RenderConvexPolygon<float,2>(float *imageRe, float *imageIm, float x1[2], float x2[2], 
								  float valRe, float valIm, const float *wRe, const float *wIm, int m, int n, bool periodic)
{
	// Two-gons are invisible :)
}

// Render a triangle with vertices chosen from the lists ax1[], ax2[].
//  i0, i1, i2 are the indices of the triangle's vertices in ax1[], ax2[].
//  valRe, valIm is the (complex) value of the triangle.
//  wRe, wIm are the (complex) weight functions for the image
//  m, n are the dimensions of the image and weight functions
//  periodic is true if the triangle should wrap around the image edges.
template<class T>
        void RenderTri(T *imageRe, T *imageIm, const T ax1[], const T ax2[], int i0, int i1, int i2, T valRe, T valIm, const T *wRe, const T *wIm, int m, int n, bool periodic)
{
    T x1[3], x2[3];
    
    // Round up the coordinates of the vertices.
    x1[0] = ax1[i0];    x2[0] = ax2[i0];
    x1[1] = ax1[i1];    x2[1] = ax2[i1];
    x1[2] = ax1[i2];    x2[2] = ax2[i2];
    
	// Pass them off to RenderConvexPolygon
	RenderConvexPolygon<T,3>(imageRe, imageIm, x1, x2, valRe, valIm, wRe, wIm, m, n, periodic);
}


template<class T>
        void RenderQuad(T *imageRe, T *imageIm, const T ix1[4], const T ix2[4], T valRe, T valIm, const T *wRe, const T *wIm, int m, int n, bool periodic, bool unwrap)
{
    T n11, n12, n21, n22, dp, dp1, dp2, d, det, n1, n2, p1, p2;
    int i, diff;
    bool side[2];
    int straddles = 0;
    T x1[4] = {ix1[0], ix1[1], ix1[2], ix1[3]};
    T x2[4] = {ix2[0], ix2[1], ix2[2], ix2[3]};

	// Look for any NaN vertices
	int nanvertex = -1;
	bool hasnan = false;

	for (i = 0; i <= 3; i++) {
        if (std::isnan(x1[i]) || std::isnan(x2[i])) {
			if (hasnan) {
				mexDebugPrintf("Two NaN vertices; returning\n");
                return;	
			}
			nanvertex = i;
			hasnan = true;
		}
    }
	
    // Unwrap the quad's coordinates if asked for.
    if (unwrap) {
        // Choose a non-NaN vertex.
        int i0 = (nanvertex == 0) ? 1 : 0;
    
        // Unwrap the rest of the vertices based on the first non-NaN
        //  vertex.
        for (i = i0 + 1; i <= 3; i++) {
            if (i == nanvertex) continue;
            
            diff = lround((x1[i] - x1[i0]) / m) * m;
            
            mexDebugPrintf("Unwrapping: adjusting %f by %d\n", x1[i], diff);
            x1[i] -= diff;
            
            diff = lround((x2[i] - x2[i0]) / n) * n;
            
            mexDebugPrintf("Unwrapping: adjusting %f by %d\n", x2[i], diff);
            x2[i] -= diff;
        }
    }

    // Draw a triangle if we have a NaN vertex.
    if (hasnan) {
		mexDebugPrintf("NaN vertex (%d): triangle\n", nanvertex);
		RenderTri(imageRe, imageIm, x1, x2, (nanvertex+1)%4, (nanvertex+2)%4, (nanvertex+3)%4,
				  valRe, valIm, wRe, wIm, m, n, periodic);
		return;
    }
    
    // Look through vertices to classify quadrilateral as
    //  convex/nonconvex/twisted (self-intersecting).
    for (i = 0; i <= 3; i++) {
        n1 = -(x2[(i+1)%4]-x2[i]);
        n2 = x1[(i+1)%4]-x1[i];
        if ((n1*n1+n2*n2) < CoincidenceThreshold*CoincidenceThreshold) {
            // Two vertices in the same location: it's a triangle.
            mexDebugPrintf("Duplicate vertices (%d,%d): triangle\n", i, (i+1)%4);
            RenderTri(imageRe, imageIm, x1, x2, (i+1)%4, (i+2)%4, (i+3)%4,
					  valRe, valIm, wRe, wIm, m, n, periodic);
            return;
        }
        
        // Take the dot product of vertex i with the normal for the
        //  edge between it and vertex i + 1.
        dp = n1*x1[i] + n2*x2[i];
        
        // Check which side of the edge the two vertices not on this edge
        //  are on.
        side[0] = (x1[(i+2)%4]*n1 + x2[(i+2)%4]*n2) > dp;
        side[1] = (x1[(i+3)%4]*n1 + x2[(i+3)%4]*n2) > dp;
        
        // Vertices on opposite sides: nonconvex or twisted quadrilateral.
        // If there are edges like this, there should be exactly two of
        //  them.
        // Record these edges in a bitmask.
        straddles |= (side[0] != side[1]) << i;
    }
    
    mexDebugPrintf("Straddle mask: %d\n", straddles);
    
    // Now split into triangles.
    switch (straddles) {
        case 0:
            // Convex quadrilateral: render it!
			RenderConvexPolygon<T,4>(imageRe, imageIm, x1, x2, valRe, valIm, wRe, wIm, m, n, periodic);
			break;
        case 3:
        case 12:
            // Nonconvex: split on the v1-v3 edge
            RenderTri(imageRe, imageIm, x1, x2, 0, 1, 3, valRe, valIm, wRe, wIm, m, n, periodic);
            RenderTri(imageRe, imageIm, x1, x2, 2, 1, 3, valRe, valIm, wRe, wIm, m, n, periodic);
            break;
        case 6:
        case 9:
            // Nonconvex: split on the v0-v2 edge
            RenderTri(imageRe, imageIm, x1, x2, 1, 0, 2, valRe, valIm, wRe, wIm, m, n, periodic);
            RenderTri(imageRe, imageIm, x1, x2, 3, 0, 2, valRe, valIm, wRe, wIm, m, n, periodic);
            break;
        case 5:
        case 10:
            int a0, a1, b0, b1;
            // Twisted: find intersection point. This case shouldn't happen
            //  very often so it's not super time-critical.
            // Here we find where the two twisted edges intersect and
            //  break up the quadrilateral into two triangles sharing
            //  that intersection point.
            
            // Get the vertex numbers a0,a1,b0,b1 of the two twisted edges,
            //  a and b.
            if (straddles == 5) {
                a0 = 0; a1 = 1; b0 = 2; b1 = 3;
            } else {
                a0 = 1; a1 = 2; b0 = 3; b1 = 0;
            }

            // Get normals to the two twisted edges and assemble into
            //  a matrix.
            n11 = -(x2[a1] - x2[a0]);  n12 = x1[a1] - x1[a0];
            n21 = -(x2[b1] - x2[b0]);  n22 = x1[b1] - x1[b0];
            
            // v1 is always in the first twisted edge; v3 is always in
            //  the second; use v1 and v3 to take dot products with the
            //  the normals.
            dp1 = x1[1]*n11 + x2[1]*n12;
            dp2 = x1[3]*n21 + x2[3]*n22;

            mexDebugPrintf("Matrix:\n%f %f\n%f %f\n", n11, n12, n21, n22);

            // Get and check the determinant.
            det = n11*n22 - n12*n21;
            if (det > -DetThreshold && det < DetThreshold) {
                // The two twisted edges are nearly parallel; don't worry
                //  about drawing this quad.
                mexDebugPrintf("Nearly collinear twisted edges; not drawing.\n");
                return;
            }
            
            // Solve matrix equation n*p = dp to find intersection point
            d = 1. / det;
            p1 = (+n22*dp1 - n12*dp2) * d;
            p2 = (-n21*dp1 + n11*dp2) * d;
            
            mexDebugPrintf("Twisted quad: intersection point (%f,%f)\n", p1, p2);
            
            // Create new vertex list with intersection point p added.
            T nx1[5] = {x1[0], x1[1], x1[2], x1[3], p1};
            T nx2[5] = {x2[0], x2[1], x2[2], x2[3], p2};
            
            // Render triangles.
            RenderTri(imageRe, imageIm, nx1, nx2, a1, b0, 4, valRe, valIm, wRe, wIm, m, n, periodic);
            RenderTri(imageRe, imageIm, nx1, nx2, a0, b1, 4, valRe, valIm, wRe, wIm, m, n, periodic);

            break;
    }
}

template<class T>
        void RenderQuadMeshInternal(T *imageRe, T *imageIm, const T *x1, const T *x2, const T *vRe, const T *vIm, const T *wRe, const T *wIm, int k, int l, int m, int n, bool periodic, bool unwrap)
{
    int i, j;
    
    for (i=0; i<k; i++) {
        for (j=0; j<l; j++) {
            // Extract the vertices of this quad
            T rx1[4] = {x1[j*(k+1)+i], x1[j*(k+1)+i+1], x1[(j+1)*(k+1)+(i+1)], x1[(j+1)*(k+1)+i]};
            T rx2[4] = {x2[j*(k+1)+i], x2[j*(k+1)+i+1], x2[(j+1)*(k+1)+(i+1)], x2[(j+1)*(k+1)+i]};
            // And render it!
            T vr = vRe[j*k+i];
            T vi = vIm ? vIm[j*k+i] : 0;
            RenderQuad(imageRe, imageIm, rx1, rx2, vr, vi, wRe, wIm, m, n, periodic, unwrap);
        }
    }
}


/*
 * Render an (anti-aliased) quadrilateral mesh and add it to the image
 *
 *  O = RenderQuadMeshInternal(O, x1, x2, v, p)
 *  O = RenderQuadMeshInternal(O, x1, x2, v, p, w)
 *
 *      O: image (m x n array)
 * x1, x2: (k+1) by (l+1) matrices with the coordinates of the mesh points.
 *      v: k by l matrix with the values to assign to each quadrilateral
 *      p: Flag controlling wraparound:
 *          p = 0: no wraparound (nonperiodic grid)
 *          p = 1: wraparound (periodic grid)
 *          p = 2: wraparound + unwrap coordinates of each quad
 *      w: optional m x n weight function.
 *          If w is given, the rendered mesh is multiplied by w before
 *          being added to the image O.
 *
 * O, v, and w may be complex.
 *
 * NOTE: the coordinates x1, x2 are 1-indexed (MATLAB-style) --
 *         so x1 = 1 for pixels at the left edge of the image, and
 *            x1 = m for pixels at the right edge of the image.
 *
 */

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    const mxArray *oMX, *x1MX, *x2MX, *vMX, *perMX, *wMX;
    mxArray *outMX;
    unsigned int k, l, m, n;
    bool periodic, unwrap;
    int per;
    void *wRe, *wIm;
    
    if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
    if (nrhs < 5 || nrhs > 6) mexErrMsgTxt("Wrong number of input arguments (expecting 5).");
    
    oMX = prhs[0];
    x1MX = prhs[1];
    x2MX = prhs[2];
    vMX = prhs[3];
    perMX = prhs[4];
    
    if (mxIsComplex(x1MX) || mxIsComplex(x2MX))
        mexErrMsgTxt("x1 and x2 must be real.");
    
    m = mxGetM(oMX);
    n = mxGetN(oMX);
    
    k = mxGetM(x1MX);
    l = mxGetN(x1MX);
    
    if (mxGetM(x2MX) != k || mxGetN(x2MX) != l)
        mexErrMsgTxt("x1 and x2 have different sizes.");
    k--;
    l--;
    if (mxGetM(vMX) != k || mxGetN(vMX) != l)
        mexErrMsgTxt("v must have size one less in each dimension than x1 and x2.");
    
    if (!mxIsScalar(perMX) || mxIsComplex(perMX) || (!mxIsNumeric(perMX) && !mxIsLogical(perMX)))
        mexErrMsgTxt("The periodic flag must be a real scalar.");
    
    per = mxGetScalar(perMX);
    periodic = (per > 0);
	unwrap = (per > 1);
    
    if (nrhs > 5) {
        wMX = prhs[5];
        if (mxGetM(wMX) != m || mxGetN(wMX) != n)
            mexErrMsgTxt("The weight function must have the same dimensions as the image.");
        wRe = mxGetData(wMX);
        wIm = mxGetImagData(wMX);
    } else {
        wRe = NULL;
        wIm = NULL;
        wMX = NULL;
    }
    
    plhs[0] = outMX = mxDuplicateArray(oMX);
    if (mxIsDouble(oMX) && mxIsDouble(x1MX) && mxIsDouble(x2MX) && mxIsDouble(vMX) && (!wMX || mxIsDouble(wMX))) {
        RenderQuadMeshInternal((double *) mxGetData(outMX),
                (double *) mxGetImagData(outMX),
                (double *) mxGetData(x1MX),
                (double *) mxGetData(x2MX),
                (double *) mxGetData(vMX),
                (double *) mxGetImagData(vMX),
                (double *) wRe,
                (double *) wIm,
                k, l, m, n, periodic, unwrap);
    } else if (mxIsSingle(oMX) && mxIsSingle(x1MX) && mxIsSingle(x2MX) && mxIsSingle(vMX) && (!wMX || mxIsSingle(wMX))) {
        RenderQuadMeshInternal((float *) mxGetData(outMX),
                (float *) mxGetImagData(outMX),
                (float *) mxGetData(x1MX),
                (float *) mxGetData(x2MX),
                (float *) mxGetData(vMX),
                (float *) mxGetImagData(vMX),
                (float *) wRe,
                (float *) wIm,
                k, l, m, n, periodic, unwrap);
    } else
        mexErrMsgTxt("Inputs must be all single or all double.");
}