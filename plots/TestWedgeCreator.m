function TestWedgeCreator(wcr, whichTest)

if nargin < 2, whichTest = 'slider'; end



switch whichTest
    case 'sum'
        f = zeros(size(wcr.grid));
        for b = 1:wcr.Count
            f = f + wcr.Wedge(b);
        end
        claimedSum = wcr.Sum;
        
        subplot 221
        imagesc(f.'); axis xy; colorbar; colormap winter
        title('Sum of wedge functions');
        subplot 222
        imagesc(claimedSum.'); axis xy; colorbar;
        title('Sum from the Sum method');
        subplot 223
        imagesc((f-claimedSum).'); axis xy; colorbar;
        title('Discrepancy');
    case 'slider'
        slider(@(box) DisplayWedge(box, @identityFilter), [0.51 wcr.Count()+0.49])
        title('Wedges');
    case 'supportslider'
        slider(@(box) DisplayWedge(box, @supportFilter), [0.51 wcr.Count()+0.49])
        title('Wedges');
    case 'logslider'
        slider(@(box) DisplayWedge(box, @logFilter), [0.51 wcr.Count()+0.49])
        title('Wedges');
end




function DisplayWedge(box, filter)
    box = round(box);
    fun = Function.WithValues(wcr.grid, filter(wcr.Wedge(box)));
    fun.plot; colormap gray; axis image
    nuBasis = wcr.NuBasis(box);
    
    hold on
    
    limx = xlim;
    limy = ylim;
    center = [mean(limx) mean(limy)];
    r = min(diff(limx), diff(limy)) / 2;
    plot([0 r*nuBasis(1,1)]+center(1), [0 r*nuBasis(2,1)]+center(2), 'ro-');
    plot([0 r*nuBasis(1,2)]+center(1), [0 r*nuBasis(2,2)]+center(2), 'go-');
    
    AddImageXYTracker
    colorbar
    
    hold off
end


end


function y = identityFilter(x)
    y = x;
end

function y = supportFilter(x)
    y = double(x ~= 0);
end

function y = logFilter(x)
    y = log10(abs(x));
end
