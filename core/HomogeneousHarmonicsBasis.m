% Generate an orthonormal basis for the space of homogeneous harmonics of
% degree n in d dimensions.
%
%   [H,I] = HomogeneousHarmonicsBasis(n, d)
%
% Inputs
%       n = degree of spherical harmonics (>= 0)
%       d = dimension                     (>= 1)
%
% Outputs
%       H = matrix with basis for spherical harmonics of degree n,
%            dimension d. Each harmonic is represented by a column of H,
%            The i'th entry in each column are coefficients for the
%            monomial given in the corresponding row of I.
%       I = companion matrix, giving the monomials used in H.
%            Each row of I contains the degrees in x_1,...,x_d.
%
%   E.g. if hypothetically
%
%      H = [ 0    -1             I = [ 0  0  2
%            1     3 ]                 1  0  1 ]
%
%   Then the rows of I represent monomials:
%      I = [ 0     0     2        -->  z^2
%            1     0     1 ]      -->  xz
%
%   and the columns of H represent polynomials
%      column 1 of H = [0  1]' --> xz
%      column 2 of H = [-1 3]' --> -z^2 + 3xz
%
% Credits
%     This function uses formula (2.164) from
%           K. Atkinson and W. Han, "Spherical Harmonics and
%              Approximations on the Unit Sphere: An Introduction." Lecture
%              Notes in Mathematics 2044. Berlin-Heidelberg:
%              Springer-Verlag 2012. DOI 10.1007/978-3-642-25983-8_2
%
%    Combining Atkinson and Han's (2.164) with (2.159) and proposition 2.42,
%      we get that, given bases y_{m,d-1,k} (k=1,2,...) for Y_{m,d-1}, a
%      basis for spherical harmonics of degree n in dimension d is
%
%        y_{n,m,d,k} = c(n,m,d) * P        (x_d) * Y        (x ,...,x   )
%                                  n-m,d+2m         m,d-1,k   1      d-1
%
%    for m = 0,1,...,n, and all values of k.
%    Here c(n,m,d) is a constant depending on n, m, d, and P_{n-m,d+2m} is
%    the Legendre polynomial of degree n-m, and dimension d+2m (also known
%    as an ultraspherical polynomial).
%
%    Note that while Y_{m,d-1,k} may be homogeneous, the Legendre
%    polynomial P is not homogeneous in general. However, P only contains
%    terms of orders n-m, n-m-2, n-m-4,... We can effectively make P
%    homogeneous of degree n-m by multiplying each term of order n-m-2k by
%    the factor (x_1^2 + ... + x_d^2)^k. This won't affect the value of
%    y_{n,m,d,k} on the unit sphere.
%
%    To multiply Y_{m,d-1,k} by P and do the homogenizing process we just
%    talked about, we use a basis consisting of the monomials in variables
%    x1,...,x_{d-1} where each variable is raised to a power between 0 and
%    n. (This is an (n+1)^(d-1) dimensional space.) We don't have to keep
%    track of the degree of the last variable xd, because we will choose
%    its power to make the polynomial homogeneous of degree n.

function [H, I] = HomogeneousHarmonicsBasis(n, d)

assert(d >= 1 & (floor(d) == d), 'HomogeneousHarmonicsBasis:badD', 'd must be a positive integer.');
assert(n >= 0 & (floor(n) == n), 'HomogeneousHarmonicsBasis:badN', 'n must be a nonnegative integer.');

persistent HomHarmCache

if n == 0,
    % For n = 0, harmonics are just constants. The norm of a constant c is
    % |c|^2 * (surface area of (n-1) sphere), so choose
    %                c = (surface area)^(-1/2)
    H = [SurfaceAreaOfNSphere(d-1)^(-1/2)];                     %#ok<NBRAK>
    I = zeros(1, d);
elseif all(size(HomHarmCache) >= [n d]) && ~isempty(HomHarmCache{n,d}),
    cache = HomHarmCache{n,d};
    fprintf('HomogeneousHarmonicBasis: (%d,%d) cached\n', n, d);
    H = cache.H;
    I = cache.I;
elseif d == 1,
    if n == 1,
        H = [sqrt(1/2)];                                        %#ok<NBRAK>
        I = [1];                                                %#ok<NBRAK>
    else
        H = [];                         % If d = 1 and n >= 2, there are no
        I = [];                         % harmonic polynomials...
    end
else
    sz = (n+1)^(d-1);               % dimension of the basis B of monomials
    inds = (n+1).^(0:(d-2));        %
    % Prepare matrix M representing multiplication by (x_1^2+...x_d^2).
    % Multiplication by x_d^2 does nothing (because the power of x_d is
    %   chosen to be whatever it needs to be to make the result homogeneous
    %   of degree n).
    % Multiplication by x_i^2, i < d, is equal to incrementing the power of
    %   x_i by 2. In the linear order of the monomials, that's the same as
    %   adding 2*(n+1)^(i-1) = 2*inds(i) to the index.
    M = spdiags(ones(sz, d), [0 -2*inds], sz, sz);
    
    % Create a cell array to store the parts of the basis.
    newBasisParts = cell(n+1,1);
    
    for m = 0:n
        % Compute the constant c_{n,m,d}
        c = 2^(-m-(d-2)/2) / gamma(m + (d-1)/2) * ...
            sqrt((2*n+d-2)*factorial(n+m+d-3)/factorial(n-m));
        % Compute the Legendre polynomial of degree n-m and dimension
        %  d+2*m.
        P = LegendrePolys(n-m, d+2*m);
        P = P(:,end);
        % Prepare matrix A representing multiplication by P,
        %  plus homogenizing.
        A = sparse(sz, sz);          % initialize to all zeros
        for k=0:(length(P)-1)/2
            A = A + P(end-2*k) * M^k;
        end
        [Hp,Ip] = HomogeneousHarmonicsBasis(m, d-1);
        % "Stretch" Hp. Hps represents the same polynomials as Hp, except
        % that the coefficients are relative to the full monomial basis.
        Hps = zeros(sz,size(Hp,2));
        Hps(Ip*inds'+1,:) = Hp;
        
        % Multiply Hps by c and P and homogenize, then finish off by adding
        %  the newly generated polynomials to the list.
        newBasisParts{m+1} = A*Hps*c;
    end
    
    % Combine the parts of the bases from m = 0...n into one big basis.
    H = [newBasisParts{:}];
    
    % Degrees in each variable
    degs = rem(floor(bsxfun(@rdivide, (0:(sz-1))', inds)), n+1);
    % Remaining degree for the variable x_d to make a monomial of degree d.
    degd = n - sum(degs, 2);
    
    % Set up I
    I = [degs degd];
    
    % Remove unused monomials from H and I.
    empties = ~any(H,2);
    H(empties,:) = [];
    I(empties,:) = [];
    
    % Cache H and I for later use.
    cache.H = H;
    cache.I = I;
    HomHarmCache{n,d} = cache;
end

end