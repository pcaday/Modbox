% Apply an FIO to a plane wave with a given frequency.
%
%   Af = ApplyFIOPatchSingleFrequency(patch, xi~, outGrid, inGridT [, diffeo])
%   Af = ApplyFIOPatchSingleFrequency(patch, xi~, outGrid, inGridT)
%
% Inputs
%      patch: FIOPatch object
%        xi~: Frequency of plane wave
%    outGrid: y grid to generate the output on
%    inGridT: x~ grid for the input
%     diffeo: Change of coordinates to use (defaults to identity).
% Outputs
%         Af: The FIO applied to the plane wave.
%
%
%   Af = ApplyFIOPatchSingleFrequency(ppatch, xi~, outGrid)
%
% Same as above, but ppatch is a ProcessedFIOPatch, and inGridT and diffeo
%  are taken from ppatch.

function Afxi = ApplyFIOPatchSingleFrequency(patch, xiT, outGrid, inGridT, diffeo)
    global GlobalProcessNewMaslov
    global GlobalProcessMaslovSign
    global GlobalSignatureTol
    
    if nargin == 3,
        % We're given a ProcessedFIOPatch, not an FIOPatch.
        inGridT = patch.inGrid;
        diffeo = patch.diffeo;
        patch = patch.patch;
    elseif nargin == 4,
        diffeo = [];
    end
    
    n = patch.dims;
    
    edgeGridT = inGridT.PixelEdgeGrid;
    
    xT = inGridT(:);
    xTEdge = edgeGridT(:);
    
    if ~isempty(diffeo),
        [x,xi] = diffeo.ForwardTransform(xT,xiT);
        [xEdge,xiEdge] = diffeo.ForwardTransform(xTEdge,xiT);
    else
        x = xT;
        xi = xiT;
        xEdge = xTEdge;
        xiEdge = xiT;
    end
    
    warp = patch.CanonicalTransformationY(xEdge,xiEdge);
    a = patch.PrincipalSymbol(x,xi);
    
    [dxy, dxiy] = patch.DCanonicalTransformationY(x,xi);
    if isempty(dxy), error('DCanonicalTransformationY must be implemented.'); end
        
    if ~isempty(diffeo),
        % Get dy/dx~ and dy/dxi~ with the chain rule.
        [dxdxt,dxidxt,dxidxit] = diffeo.DInverseCotangentT(xT, xiT);
        dxty =  ParallelMatrixMultiply(dxy, dxdxt)...
              + ParallelMatrixMultiply(dxiy, dxidxt);
        dxity = ParallelMatrixMultiply(dxiy, dxidxit);
    else
        dxty = dxy;
        dxity = dxiy;
    end
    
    % Divide principal symbol by half-density factor.
    a = a ./ sqrt(abs(ParallelDet(dxty)));
    
    % Compute and multiply by Maslov factor.
    dxitxt = -ParallelMatrixMultiply(ParallelMatrixInverse(dxty), dxity);
    if GlobalProcessNewMaslov
        sgn = ParallelSignature(dxitxt, GlobalSignatureTol, 1);
        a = a .* 1i.^(-GlobalProcessMaslovSign*sgn);
    else
        sgn = ParallelSignature(dxitxt, GlobalSignatureTol);
        a = a.* exp(-GlobalProcessMaslovSign*pi*0.25i*sgn);
    end
    
    a(~isfinite(a)) = 0.;
    
    Afxi_x = exp(1i * xT.' * (xiT(:))) .* a;
    Afxi_x = reshape(Afxi_x, size(inGridT));
    Afxi_x = Function.WithValues(inGridT, Afxi_x);
    
    Afxi = Function.Zeros(outGrid, class(inGridT.mins));
    Afxi.f = complex(Afxi.f);
%    if any(outGrid.periodic),
%        wrap = 2;
%    else
%        wrap = 0;
%    end
    warp = reshape(warp, [2 size(inGridT)+1]);
    Pushforward(Afxi, Afxi_x, warp);  % , wrap);
    
    if nargout == 0,
        Afxi.plot; colorbar;
        plottitle = sprintf('FIO applied to plane wave \\xi~ = (%s)  ', ...
            strjoin(num2str(xiT), ','));
        normtitle(plottitle, Afxi)
    end
end