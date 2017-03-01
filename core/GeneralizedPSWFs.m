% Compute generalized prolate spheroidal wave functions (PSWFs)
%

function [psis, lambdas] = GeneralizedPSWFs(c, d, tol, K, useBessel, samples)

if nargin < 2, d = 2; end
if nargin < 3, tol = 1e-8; end
if nargin < 4, K = 50; end
if nargin < 5, useBessel = false; end
if nargin < 6, samples = 0; end

if d == 1
    Nmax = 1;           % Only do N = 0 and N = 1 for 1D.
else
    Nmax = 100;
end

p = d - 2;

psis = [];
lambdas = [];

for N = 0:Nmax
    % Renormalize (divide phi by r^(N+(p+1)/2)) except in 1D.
    doRenorm = d > 1;
    
    % Compute the radial parts of the PSWFs
    [gammas, renorm_phis, renorm_phis_bessel] = PSWFRadial(N + p/2, c, K, doRenorm, 0, tol, samples);
    k = length(renorm_phis);

    if useBessel, renorm_phis = renorm_phis_bessel; end
    
    if isempty(gammas), break; end
    
    if d > 1
        betas = gammas*c^(-(p+1)/2);
        alphas = (-1i)^N*(2*pi)^(1+p/2)*betas;
    
        % Divide phi by r^(N+(p+1)/2).
        % The r^N factor cancels out the r^N growth in the homogeneous
        %  harmonic polynomials.
        %
        % Note this perfectly cancels the x^(N+1/2) factor in the radial PSWF
        %  leaving only the Jacobi polynomials.
        [H, I] = HomogeneousHarmonicsBasis(N, d);
        
        % Generate all products of renormalized radial parts (renorm_phis)
        %  and spherical harmonics of order N (columns of H)
        h = size(H,2);
        newGPSWFs = cell(h,k);
        for i = 1:h
            for j = 1:k
                newGPSWFs{i,j} = ...
                    @(x) renorm_phis{j}(sqrt(sum(x.*x,2))) ...
                    .* EvalCompressedNdPoly(H(:,i), I, x);
            end
        end

        fprintf('GeneralizedPSWFs: N = %d, %d radial parts, %d harmonics\n', N, k, h);
        
        repalphas = repmat(alphas', h, 1);
        
        psis = [psis; newGPSWFs(:)];       %#ok<AGROW>
        lambdas = [lambdas; repalphas(:)]; %#ok<AGROW>
    else
        alphas = (1i)^N*(pi/2)^(1/2)*gammas;
        
        newGPSWFs = cell(size(renorm_phis));
        % Extend to the whole real line. If we're using Bessel
        %  functions, the radial PSWFs aren't defined all the way to -1,
        %  so we extend them to even (N = 0) or odd (N = 1) functions
        %  on the real line.
        if N == 1,
            for j = 1:k
                newGPSWFs{j} = @(x) sign(x) .* renorm_phis{j}(abs(x));
            end
        else
            for j = 1:k
                newGPSWFs{j} = @(x) renorm_phis{j}(abs(x));
            end            
        end
        psis = [psis; newGPSWFs(:)]; %#ok<AGROW>
        lambdas = [lambdas; alphas]; %#ok<AGROW>
    end
    
    if N == Nmax && d > 1
        warning('GeneralizedPSWFs:ReachedNMax',...
            'Maximum N reached; some PSWFs will be missing. Increase Nmax in GeneralizedPSWFs.m to find them.');
    end
end

end