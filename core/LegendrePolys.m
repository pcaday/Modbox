% Generate Legendre polynomials.
%
% P = LegendrePoly(N,d) computes the Legendre polynomials of
%  degrees 0 <= n <= N in dimension d.
%
% It returns the polynomials as an (N+1)x(N+1) matrix. The j'th column of
% the matrix holds the coefficients of the Legendre polynomial of degree
% (j-1), arranged from lowest order term to highest. That is, P(i,j) is the
% coefficient of x^(i-1) in the Legendre polynomial P_{j-1,d} of degree j-1
% and dimension d.
%
% Reference: K. Atkinson and W. Han, "Spherical Harmonics and
%              Approximations on the Unit Sphere: An Introduction." Lecture
%              Notes in Mathematics 2044. Berlin-Heidelberg:
%              Springer-Verlag 2012. DOI 10.1007/978-3-642-25983-8_2

function P = LegendrePolys(N, d)

% LegendreCache is a cell array of matrices. LegendreCache{d} represents
% the output of LegendrePoly for dimension d, and whatever size the matrix
% LegendreCache{d} is.
persistent LegendreCache

assert(d >= 2 & (floor(d) == d), 'LegendrePoly:badD', 'd must be an integer greater than 1');
assert(N >= 0 & (floor(N) == N), 'LegendrePoly:badN', 'N must be a nonnegative integer');

% Check if we have cached the Legendre polynomials for this dimension, up
% to this value of N.
if length(LegendreCache) >= d && size(LegendreCache{d}, 1) >= N+1
    % If so, retrieve the answer from the cache
    P = LegendreCache{d}(1:N+1,1:N+1);
    fprintf('LegendrePoly: (%d,%d) cached\n', N, d);
else
    % Otherwise, compute it.
    
    P = zeros(N+1);
    
    % Order zero polynomial
    P(1,1) = 1;
    
    if N < 1, return; end
    
    % Order 1 polynomial
    P(2,2) = 1;
    
    for n = 2:N
        P(1:n+1,n+1) = (2*n+d-4)/(n+d-3) * [0; P(1:n,n)]...
            - (n-1)/(n+d-3) * P(1:n+1,n-1);
    end
    
    % Save for later.
    LegendreCache{d} = P;
end

end