function [gammas, phis, phis_bessel, chis] = PSWFRadial(N, c, K, renorm, dtol, tol, samples)
% Compute radial components of generalized prolate spheroidal wave functions,
%   using Shkolinsky's algorithm (expansion in Jacobi polynomials).
%
%   References: Y. Shkolinsky, "Prolate spheroidal wave functions on a disc--
%                 Integration and approximation of two-dimensional bandlimited
%                 functions."
%               Applied and Computational Harmonic Analysis 22 (2007), 235-256.
%				http://www.sciencedirect.com/science/article/pii/S1063520306000960
%
%               D. Slepian, "Prolate spheroidal wave functions, Fourier
%                 analysis and uncertainty IV: Extensions to higher
%                 dimensions; generalized prolate spheroidal wave
%                 functions."
%               Bell System Technical J. Nov 1964
%
%          Note: the polynomial basis functions T_{N,n} have the scaling
%          from Shkolinsky's paper. The T_{N,n} in Slepian's paper differ
%          by a factor of
%
%                h_{N,n} = sqrt(2*(2n+N+1)) * ((n+N) choose n).
%
%          Shkolinsky's scaling makes the T_{N,n} not only orthogonal but
%           orthonormal, while Slepian's scaling has the advantage that 0
%           <= T_{N,n}(x) <= 1 for x in (0,1].
%
%
%  [gammas, phis] = PSWFRadial(N, c, K, renorm, dtol, tol, samples);
%  [gammas, phis, phis_bessel] = PSWFRadial(...)
%  [gammas, phis, phis_bessel, chis] = PSWFRadial(...)
%
% Inputs
%         N: degree
%         c: bandlimit
%         K: maximum index of radial functions (total K+1 radial functions)
%              *NOTE* the Jacobi polynomials lose precision for K > 40
%              approximately. For larger values of K, phis_bessel (see below)
%              may be a better choice than phis.
%              Default: 10.
%
%    renorm: Boolean flag. If true, phis and phis_bessel will be divided by
%              x^(N+1/2) (useful for the GeneralizedPSWFs routine).
%              Default: false.
%
%      dtol: tolerance for the coefficients d_j^{N,n}. Any d_j^{N,n} below
%              dtol will be considered zero.
%              Default: 0.
%
%       tol: error tolerance. Default: 0.
%
%   samples: If nonzero, the number of sample points (minus 1) at which
%             to sample the output functions. If zero, the output functions
%             are not sampled.
%            Default: 0.
%
% Outputs
%      gammas: vector of eigenvalues (sorted high to low) - the gamma_{n,N}
%                in Slepian and Shkolinsky's papers.
%
%        phis: cell array of corresponding eigenfunctions (radial parts of
%                generalized PSWFs). These take an r value (radius)
%                and return the value of the function at that radius.
%              The phis are computed using Jacobi polynomials, and can be
%                unstable for larger K (> 40 approximately) and are 
%                restricted to the domain [0, 1].
%              If samples = 0 (see above), these are function handles,
%                otherwise they are griddedInterpolant objects
%                (see MATLAB documentation). Both can be called the same
%                way (phis{i}(r_values)), and can handle vector inputs.
%
%              These are the phi_{n,N} in Slepian and Shkolinsky's papers.
%
% phis_bessel: same as phis, but computed in a different way, using Bessel
%                functions. They appear more stable and are not restricted
%                to the domain [0, 1] like the Jacobi polynomials, but are
%                slower to compute.
%
%        chis: vector of eigenvalues of the phis w.r.t. the matrix used
%                internally to solve for the coefficients used in
%                constructing the phis. 

% Substitute default arguments if needed
if nargin < 3, K = 10; end
if nargin < 4, renorm = false; end
if nargin < 5, dtol = 0; end
if nargin < 6, tol = 0; end
if nargin < 7, samples = 0; end

CheckTypes({N,c,K,renorm,dtol,tol,samples}, {'scalar','scalar','scalar','logical','scalar','scalar','scalar'});

%validateattributes(N, 'numeric', 'scalar', 1)
%validateattributes(c, 'numeric', 'scalar', 2)
%validateattributes(K, 'numeric', 'scalar', 3)
%validateattributes(renorm, 'numeric', 'scalar', 4)
%validateattributes(dtol, 'numeric', 'scalar', 5)
%validateattributes(tol, 'numeric', 'scalar', 6)


doBessel = false;
if nargout > 2, doBessel = true; end

% Prepare the matrix omtxs for substituting 1 - 2x^2.
omtxs = zeros(2*K+1, K+1);
omtxs(1,1) = 1;
for n = 1:K
    omtxs(:,n+1) = omtxs(:,n) - 2*[0; 0; omtxs(1:end-2,n)];
end

% Prepare the coefficients for scaling P_n (h_{N,n} / (n+N choose n))
ns = (0:K)';                            % make it a column vector

facs = (ns + N) ./ ns;
facs(1) = 1;
npncn = cumprod(facs);                  % npncn = (n+N choose n)

pscale = diag(sqrt(2*(2*ns+N+1)));      % Skholinsky's normalization --
                                        % makes the T_{N,n} orthonormal
%pscale = diag(1./npncn);
                                        % Slepian's normalization --
                                        % ensures 0 <= T_{N,n} <= 1
Jp = JacobiPolys(K,N,0);
Tp = omtxs * Jp * pscale;


% Prepare the matrix for the polynomial coefficients

kappa = (N + 2*ns + 1/2) .* (N + 2*ns + 3/2);
h = sqrt(2 * (2*ns+N+1)) .* npncn;
gammap1 = -((ns+N+1).^2) ./ ((2*ns+N+1) .* (2.*ns+N+2)) .* (h ./ [h(2:end); 1]);
gamma0 = (2*ns.*(ns+1) + N * (2*ns+N+1)) ./ ((2*ns+N) .* (2*ns+N+2));
if (N == 0)
    gamma0(1) = 0.5;            % The denominator in gamma0 is 0 at n = 0 for N = 0
end                             %  ... fix up the value.

% It's symmetric, so copy gammap1 to gammam1.
gammam1 = [1; gammap1(1:end-1)];

B = full(spdiags([-c^2*gammap1 -kappa-c^2*gamma0 -c^2*gammam1],...
    [-1 0 1], K+1, K+1));
[D, chis] = eig(B);
chis = diag(chis);


D = fliplr(D);                  % These are the coefficients d_j^{N,n}
                                %  using Shkolinsky's scaling. (d_j^{N,n} =
                                %  D(j+1,n+1). )
D(abs(D) < dtol) = 0;
D_Sl = bsxfun(@times, D, h(:)); % These are the coefficients d_j^{N,n}
                                %  using Slepian's scaling.

p = Tp*D;                       % These are the coefficients of the
                                %  polynomials.

chis = flipud(chis);

% Calculate eigenvalues (gammas) of the corresponding eigenfunctions from
% the coefficients d_j^{N,n}. (Equation (46) in Slepian)
gammas = c^(N+1/2) / (2^(N+1) * gamma(N+2)) * D_Sl(1,:) ./ sum(D_Sl,1);
gammas = gammas';

% Output the function handle array for the eigenfunctions (phis) found.
phis = cell(K+1,1);
deg = size(p,1)-1;

if renorm
    for n = 1:(K+1)
        phis{n} = @(x) JacobiPolyEvalRenorm(x, p(:,n), N, deg);
    end
else
    for n = 1:(K+1)
        phis{n} = @(x) JacobiPolyEval(x, p(:,n), N, deg);
    end
end


if doBessel
    % Calculate function handles which compute an alternate expression for
    % phi using a Bessel function representation. These are valid on
    %  [0,inf), not just [0,1]. We can also get these from the coefficients
    %  d_j^{N,n}, which is cool.
    %
    % This is equation (44) in Slepian's paper.
    
    phis_bessel = cell(K+1,1);
    
    js = 0:K;
    facs = (N + js) ./ js;
    facs(1) = 1;
    npjcj = cumprod(facs);               % npjcj = (n+j) choose j
    
    if renorm,
        for n = 1:(K+1)
            zeroVal = phis{n}(0);
            phis_bessel{n} = @(x) BesselPSWFEvalRenorm(x,...
                D_Sl(:,n).' ./ npjcj, N, K, c, gammas(n), zeroVal);
        end
    else
        for n = 1:(K+1)
            zeroVal = phis{n}(0);
            phis_bessel{n} = @(x) BesselPSWFEval(x,...
                D_Sl(:,n).' ./ npjcj, N, K, c, gammas(n), zeroVal);
        end
    end
end


% Remove radial PSWFs with eigenvalues less than tol.
phis(abs(gammas)<tol) = [];
if doBessel
    phis_bessel(abs(gammas)<tol) = [];
end
chis(abs(gammas)<tol) = [];
gammas(abs(gammas)<tol) = [];


% If requested, sample the phis.
% In a short test, 'spline' interpolation was actually faster than 'linear'
%  after the first call to the interpolant.
if samples > 0,
    samplePoints = linspace(0,1,samples+1);
    for j = 1:length(phis)
        phis{j} = griddedInterpolant(samplePoints, phis{j}(samplePoints), 'spline');
        phis_bessel{j} = griddedInterpolant(samplePoints, phis_bessel{j}(samplePoints), 'spline');
    end
end

end





% These four functions implement the radial PSWF functions (the phi's).
% They are referenced in the anonymous function handles returned by
%  PSWFRadial.

% Evaluate a phi using Jacobi polynomial coefficients (renormalized)
%
% We can do this in one line with
%   vals = bsxfun(@power, xs, (0:deg)) * coeffs;
% but this is faster (compare to vander.m built in to MATLAB).
%
function vals = JacobiPolyEvalRenorm(xs, coeffs, ~, deg)
xpows = ones(size(xs));
vals = zeros(size(xs));
for i = 0:deg
    vals = vals + xpows * coeffs(i+1);
    xpows = xpows .* xs;
end
end

% Evaluate a phi using Jacobi polynomial coefficients (not renormalized)
%
% We can do this in one line with
%   vals = bsxfun(@power, xs, ((0:deg) + N + 1/2)) * coeffs;
% but this is faster.
%
function vals = JacobiPolyEval(xs, coeffs, N, deg)
xpows = xs .^ (N+1/2);
vals = zeros(size(xs));
for i = 0:deg
    vals = vals + xpows * coeffs(i+1);
    xpows = xpows .* xs;
end 
end

% Evaluate a phi using Bessel functions (renormalized)
%
% The old code was: 
%    phis_bessel{n} = @(x) sum(D_Sl(:,n)' .* besselj(N+2*js+1, c*x) ./ ...
%                (sqrt(c) * x^(N+1) * npjcj)) / gammas(n);
function vals = BesselPSWFEvalRenorm(xs, coeffs, N, K, c, gamma, zeroVal)
vals = zeros(size(xs));
for j = 0:K
    vals = vals + coeffs(j+1) * besselj(N + 2*j + 1, c*xs);
end
vals = vals ./ ((gamma * sqrt(c)) * xs.^(N+1));

% The formula above is undefined at x = 0, so we substitute zeroVal
%  (which we get using the Jacobi polynomials) instead.
vals(xs == 0) = zeroVal;
end

% Evaluate a phi using Bessel functions (not renormalized)
%
% The old code was:
%            phis_bessel{n} = @(x) sum(D_Sl(:,n)' .* besselj(N+2*js+1, c*x) ./ ...
%                (sqrt(c*x) * npjcj)) / gammas(n);
function vals = BesselPSWFEval(xs, coeffs, N, K, c, gamma, zeroVal)
vals = zeros(size(xs));
for j = 0:K
    vals = vals + coeffs(j+1) * besselj(N + 2*j + 1, c*xs);
end
vals = vals ./ (gamma * sqrt(c*xs));

% The formula above is undefined at x = 0, so we substitute zeroVal
%  (which we get using the Jacobi polynomials) instead.
vals(xs == 0) = zeroVal;
end