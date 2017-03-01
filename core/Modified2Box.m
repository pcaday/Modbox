% Apply an FIO to a function, approximately.
%
%         Au = Modified2Box(pA, u, outGrid, valueClass, zeroPadding)
%   [Au,Aus] = Modified2Box(...)
%
% Inputs
%          pA: a ProcessedFIO object (created by the ProcessedFIO.Process
%               factory method) representing the FIO.
%           u: the input function, as a Function object
%     outGrid: the Grid object representing the codomain of the FIO.
%  valueClass: (optional) class of the output to use: 'single' or 'double'.
%                Default: 'single'.
% zeroPadding: (optional) amount of zero padding for FFTs.
%                e.g. a factor of 3.0 with a 50x50 input grid would
%                 zero-pad all inputs to 150x150 before taking FFTs.
%                Default: 1.0 (no padding).
%
% Outputs
%          Au: the computed result of applying the FIO pA to the function u
%               as a Function object.
%         Aus: a cell array which holds the contributions to Au due to each
%               patch.
%
% Grid Notes
%    - The grid of the input function u must match the input grid used when
%       creating pA. To apply the same FIO to another size of input,
%       the original FIO must be re-processed with the corresponding grid.
%
%    - The grid of the output, outGrid, does not need to match the output
%       grid used by pA. However, if pA's output grid is periodic in one or
%       more dimensions, outGrid must also be periodic in those dimensions
%       (only) and the periods should match.
%

% This version of the algorithm approximates the pseudodifferential
%  amplitude b(x,x',xi) by b(x',x',xi).
% This approximation is accurate to leading order.

function [Au,Aus] = Modified2Box(processedFIO, u, outGrid, valueClass, zeroPadding)
global GlobalPSWFCache
global GlobalAlphaThreshold
global GlobalM2BVerbosity
global GlobalFreqNoiseFloor
global GlobalDFTRoutine
global GlobalIDFTRoutine
global GlobalUsePseudopolarFFT

assert(u.grid == processedFIO.inGrid, 'The input function does not match the processed FIO''s grid. Reprocess the FIO with the matching grid.');

if nargin < 4, valueClass = 'single'; end
if nargin < 5, zeroPadding = 1; end

% Make a copy of u, then cast its values to the desired class.
%  (The copy is made so the caller's u is not modified; it's passed by
%   reference.)
u = u.copy();
u.f = cast(u.f, valueClass);

patchList = processedFIO.patches(:);

n = processedFIO.fio.dims;
d = n*(n-1)/2;
Au = Function.Zeros(outGrid, valueClass);
if nargout > 1,
    Aus = cell(size(patchList));
end

% Make sure A is complex. Otherwise, RenderQuadMesh won't compute the
%  imaginary part of pushforwards.
Au.f = complex(Au.f);

% msgChars will hold the number of characters from the last informational
%  message.
msgChars = 0;

% Loop over patches.
patchNum = 0;
for pPatch = patchList.'
    patchNum = patchNum + 1;
    
    % Pull u back by the diffeomorphism.
    xt = pPatch.inGrid(:);
    xt = cast(xt, valueClass);
    x = pPatch.localPatch.diffeo.InverseTransform(xt);
    ut = u.Pullback(pPatch.inGrid, x);
    clear x;
    
    % Compute and store its DFT.
    utHat = GlobalDFTRoutine(ut, zeroPadding);
    
    % Apply noise floor to the DFT
    utHat.f(abs(utHat.f) < GlobalFreqNoiseFloor) = 0;
    
    % Prepare some variables for later:
    %  ukHat: the contents of the k'th box of u, in the frequency domain.
    % ukmHat: temporary variable for tensor product terms.
    ukHat = utHat.copy();
    ukmHat = ukHat.copy();
    
    boxes = pPatch.wedgeCreator.Count();
    
    deg = pPatch.localPatch.globalPatch.degree;
    
    % Prepare the |xi|^m multiplier for the degree of the FIO.
    utHatGrid = utHat.grid;
    xit = utHatGrid(:);
    xit = cast(xit, valueClass);
    xitR = sum(xit.*xit, 1) .^ (deg / 2);
    xitR = reshape(xitR, size(utHatGrid));
    
    % If the degree is negative, set the (0,0) entry(ies) of the multiplier
    %  to 0, rather than leaving them infinite.
    if deg < 0,
        if GlobalUsePseudopolarFFT,
            xitR(1,:) = 0;
        else
            subs = DFTCenterSubs(utHatGrid);
            xitR(subs{:}) = 0;
        end
    end
    
    contrib = Function.Zeros(ut.grid, valueClass);
    clear utHatGrid
    
    % TEST
    contribTotal = contrib.copy;
    
    % For each box in the patch:
    for b = 1:boxes
        % Get the amplitude
        amp = pPatch.amplitude{b};
        
        % Multiply by box cutoff in the frequency domain.
        ukHat.f = utHat.f .* pPatch.wedgeCreator.Wedge(b);

        c = pPatch.expCs(b);
        if ~isfinite(c),
            % Skip NaN and Inf c's
            continue;
        elseif (c == 0),
            % c = 0: no exponential term.
            % Just multiply by amplitude and wedge cutoff function.
            
            % Map NaN's to zero.
            amp(~isfinite(amp)) = 0;
            amp = reshape(amp, size(ut.f));
                
            % Multiply by the wedge function in frequency domain
            if deg ~= 0, ukHat.f = ukHat.f .* xitR; end
        
            % Return to spatial domain and multiply by amplitude.
            up = GlobalIDFTRoutine(ukHat, ut.grid);
            up.f = up.f .* amp;
                
            contrib.f = up.f;
        else
            % Get PSWFs to approximate the second-order exponential term.
            [~,psis,lambdas] = GlobalPSWFCache.GetPSWFs(d,c);
            
            % Get the g~ function and scale it to get g.
            [g,gmax] = ComputeUnscaledG(xit, pPatch.wedgeCreator, b, valueClass);
            g = g / (gmax * pPatch.cExtras(b));
            g(~isfinite(g) | ukHat.f(:) == 0) = 0;
        
            % For each tensor product term:
            contrib.f(:) = 0;
            nSignificant = 0;
            
            for j = 1:length(lambdas)
                % Prepare the alpha function (x' multiplier).
                alpha = psis{j}(pPatch.expFs(:,:,b)) .* lambdas(j);
                alpha = alpha .* amp;
                alpha(~isfinite(alpha)) = 0.;
                
                % Skip if the alpha function is negligible.
                if max(abs(alpha(:))) < GlobalAlphaThreshold, continue; end
                alpha = reshape(alpha, size(ut.f));
                
                % Tracing output.
                if GlobalM2BVerbosity > 2 && n == 2,
                    subplot(2,3,6,'replace');
                    uh = abs(ukHat);
                    uh.f = max(log(uh.f), -10);
                    uh.plot; colorbar;
                    title('Original log \hat u_k');
                end
                
                % Multiply by the theta function (xi multiplier)
                %  in frequency domain.
                ukmHat.f = ukHat.f .* reshape(psis{j}(g), size(ukHat.f));
                
                % Multiply by |xi|^d if non-zero degree
                if deg ~= 0, ukmHat.f = ukmHat.f .* xitR; end
        
                % Return to spatial domain and multiply by alpha.
                up = GlobalIDFTRoutine(ukmHat, ut.grid);
                up.f = up.f .* alpha;
                contrib.f = contrib.f + up.f;
                
                % Count this as a significant term in the output
                nSignificant = nSignificant + 1;
                
                % Debugging output
                if GlobalM2BVerbosity > 2 && n == 2,
                    subplot(2,3,1,'replace');
                    imagesc(abs(alpha).'); axis xy; colormap cool; colorbar;
                    title \alpha
                    subplot(2,3,2,'replace');
                    imagesc(abs(reshape(psis{j}(g), size(ukmHat.f))).'); axis xy; colorbar;
                    title \hat\theta
                    subplot(2,3,3,'replace');
                    uh = abs(ukmHat);
                    uh.f = max(log(uh.f), -10);
                    uh.plot; colorbar;
                    title('Multiplied log \hat u_k');
                    suptitle(sprintf('Patch %d, box %d, term %d', patchNum, b, j));
                            
                    drawnow
                end
                clear up
                clear alpha
            end
            
            % Debugging output
            if GlobalM2BVerbosity > 0,
                if GlobalM2BVerbosity == 1,
                    ss = repmat('\b', 1, msgChars);
                else
                    ss = '\n';
                end
                ss = sprintf(ss);
                msgChars = fprintf('%sPatch %d, box %d: %d significant terms', ... 
                    ss, patchNum, b, nSignificant) - length(ss);
            end
            
            clear psis
            clear lambdas
        end
        
        clear g
        clear amp
        
        % Debugging output
        if GlobalM2BVerbosity >= 1 && any(isnan(contrib.f(:)))
            warning('Modified2Box: NaN encountered in output.');
        end
            
        % TEST
        contribTotal.f = contribTotal.f + contrib.f;
        
        % Pushforward (with Jacobian factor)
        Pushforward(Au, contrib, pPatch.warps{b}); %, processedFIO.fio.wrapFlag);
        
        % Tracing output
        if GlobalM2BVerbosity > 1 && n == 2,
            subplot(2,3,4,'replace');
            Au.plot; colorbar;
            title('Real part')
            subplot(2,3,5,'replace');
            imagAf = imag(Au);
            imagAf.plot; colorbar;
            title('Imaginary part')
            hold on
            DrawIrregularGrid(pPatch.warps{b});
            suptitle(sprintf('Patch %d, box %d', patchNum, b));

            drawnow
        end
    end
    
    clear xit
    clear xitR
    clear contrib
    
    % After each patch is done, save the current cumulative sum.
    if nargout > 1, Aus{patchNum} = Au.copy(); end
    
    if GlobalM2BVerbosity > 1, suptitle(''); end
end

% Currently Afs is a cumulative sum of the patch contributions.
% Undo the cumulative sum to retrieve the individual contributions.
if nargout > 1,
    for patchNum = length(patchList):-1:2
        Aus{patchNum}.f = Aus{patchNum}.f - Aus{patchNum - 1}.f;
    end
end

end
