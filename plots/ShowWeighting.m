% Shows the weighting of the different SRD contributions.
%
%   SHOWWEIGHTING(pFIO, wcr)
%
% pFIO: a ProcessedFIO object.
%  wcr: WedgeCreator for the original (x,xi) grid.
%
% For best results, set amplitude to identically 1 in the FIO's patches
%  and set
%     GlobalProcessNoDyDxFactor = true
%     GlobalProcessMaslovSign   = 0
%     GlobalFunctionCxPlot      = false
%
% Creates one figure for each of the patches, numbered 51...50+k,
%  where k is the number of SRDs. Each figure shows the amplitude
%  pulled back to original (x,xi) coordinates. The scroll bar controls
%  the xi value (box number).
%
% Finally, a similar figure, figure 50, is created, with the sum of all the
%  amplitudes.
%
% If the amplitude is identically 1, then this will show the cutoffs for
%  each SRD, and the sum of the cutoffs (which should be 1).
function ShowWeighting(pFIO, wcr)
    xGrid = pFIO.inGrid;
    xi = wcr.Centers();
    
    x = xGrid.AllPoints;
    x = repmat(x, [ones(1,ndims(x)) wcr.Count]);
    xi = reshape(xi, [size(xi,1) ones(1,ndims(x)-2), wcr.Count]);
    xi = repmattosize(xi, size(x));
    
    npatches = length(pFIO.patches);
    a_pull = cell(1,npatches);
    a_sum = 0;
    
    for k = 1:npatches
        xTGrid = pFIO.patches(k).inGrid;
        wcrT = pFIO.patches(k).wedgeCreator;
        
        a = reshape([pFIO.patches(k).amplitude{:}], ...
            [size(xTGrid) wcrT.Count]);
        
        % For testing:
%         %a = mod(1:wcrT.Count(), 2);   % xi stripes
%         a = (1:wcrT.Count()) / 10;
%         a = reshape(a, 1, 1, []);
%         a = repmat(a, [size(xTGrid) 1]);
%         xTCF = xTGrid.coordFuns();
%         ax = xTCF{1} .* xTCF{2};
%         ax = repmattosize(ax, size(a));
%         a = a + ax;
        
        I = CotangentInterpolator(a, xTGrid, wcrT);
                
        [xt,xit] = pFIO.patches(k).localPatch.diffeo.ForwardTransform(x,xi);
        
        a_pull{k} = I(xt,xit);
        a_pull{k}(isnan(a_pull{k})) = 0;
        a_sum = a_sum + a_pull{k};
    end
    
    figure(50)
    clf
    ShowCotangentFun(a_sum, wcr, xGrid, [0 max(a_sum(:))], false);
    
    for k = 1:npatches
        figure(50+k)
        clf
        %range = {};
        range = {[0 max(a_pull{k}(:))]};
        ShowCotangentFun(a_pull{k}, wcr, xGrid, range{:}, false);
    end
end