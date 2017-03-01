function [g,gmax] = ComputeUnscaledG(xi,wedgeCreator,b,valueClass)

global GlobalGmaxMinXiRadius

% Get the nu-adapted basis.
nuBasis = wedgeCreator.NuBasis(b);
nuBasis = cast(nuBasis, valueClass);

% Change coordinates to nu basis
xiInNuBasis = ParallelMatrixMultiply(nuBasis.', xi);

n = size(nuBasis,1);

% Allocate g.
sz = size(xi);
sz(1) = [];
smSizeG = [prod(sz) n*(n-1)/2];
sizeG = [sz n*(n-1)/2];
g = zeros(smSizeG,valueClass);

% Invert first (xi') coordinate.
xiInNuBasis(1,:) = 1 ./ xiInNuBasis(1,:);

% Compute and store g_ij = xi''_i * xi''_j / xi'
idx = 1;
for j = 2:n
    for i = 2:j
        g(:,idx) = xiInNuBasis(i,:) .* xiInNuBasis(j,:) .* xiInNuBasis(1,:);
        idx = idx + 1;
    end
end

% Compute maximum ||g||_2 over the box's support if requested.
if nargout > 1,
    boxSupp = logical(wedgeCreator.Wedge(b));
    
    % Do a cutoff to remove small |xi| from the support
    % (to avoid getting spuriously high gmax's)
    grid = wedgeCreator.grid;
    xiCF = grid.coordFuns();
    xiCF{1} = xiCF{1} / max(xiCF{1}(:));
    xiCF{2} = xiCF{2} / max(xiCF{2}(:));
    xiR = max(abs(xiCF{1}), abs(xiCF{2}));
    boxSupp(xiR < GlobalGmaxMinXiRadius) = false;
    
    gmax = sqrt(max(sum(g.*g, 2) .* boxSupp(:)));  % Maximum 2-norm 
end

g = reshape(g, sizeG);

end