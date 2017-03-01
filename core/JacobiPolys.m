function p = JacobiPolys(N, alpha, beta)
% Return the Jacobi polynomials P_n^(alpha,beta) for 0 <= n <= N.
%
% The polynomials are returned as columns in an (N+1) x (N+1) upper
%  triangular matrix, where the n'th polynomial is stored in the (n + 1)'st
%  column.
% 
% Each column lists the coefficients for that polynomial,
%  arranged from lowest order to highest,
%  e.g. the column [2; -1; 1] stands for x^2 - x + 2.
%
% An example of the output is
%   JacobiPolys(2,4,0) = [ 1     2     1
%                          0     3     7
%                          0     0     7 ]
%
%  representing the polynomials 1, 3x+2, and 7x^2+7x+1.
%   
%
% Reference:
%        Section 2.2 of Y. Shkolinsky, "Prolate spheroidal wave
%                 functions on a disc-- Integration and approximation of
%                 two-dimensional bandlimited functions."
%
%        Applied and Computational Harmonic Analysis 22 (2007), 235-256.
%	     http://www.sciencedirect.com/science/article/pii/S1063520306000960
%

N = floor(N);
assert(N >= 0, 'JacobiPoly:Npos', 'N must be nonnegative');

p = zeros(N+1);

% 0th polynomial
p(1,1) = 1;

if N < 1, return; end

% 1st polynomial
p(1,2) = (alpha - beta) / 2;
p(2,2) = (alpha + beta) / 2 + 1;

for n = 1:(N-1)
	% Compute (n+1)'st polynomial.
	a1n = 2 * (n+1) * (n+alpha+beta+1) * (2*n+alpha+beta);
	a2n = -(2*n + alpha + beta + 1) * (alpha^2 - beta^2);
	a3n = (2*n + alpha + beta) * (2*n + alpha + beta + 1) * (2*n + alpha + beta + 2);
	a4n = 2 * (n + alpha) * (n + beta) * (2*n + alpha + beta + 2);

	p(:,n+2) = (a3n*[0; p(1:end-1,n+1)] - a2n*p(:,n+1) - a4n*p(:,n)) / a1n;
	%                  ^-- multiplication by x
end


end