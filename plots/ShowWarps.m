function ShowWarps(pPatch)
    % Find min and max x and y values.
    
    nwarps = length(pPatch.warps);
    mins = [inf; inf];
    maxes = [-inf; -inf];
    
    for i = 1:nwarps
        warp = pPatch.warps{i};
        mins = min(mins, min(warp(:,:), [], 2));
        maxes = max(maxes, max(warp(:,:), [], 2));
    end
    
    arrow_center = (mins+maxes) / 2;
    arrow_len = min(maxes-mins) * 0.25;

    xlim([mins(1) maxes(1)])
    ylim([mins(2) maxes(2)])
    
    slider(@ShowOneWarp, [0.51 nwarps+0.49])
    
    function ShowOneWarp(sliderVal)
        bt = round(sliderVal);
        cla
        DrawIrregularGrid(pPatch.warps{bt});
        
        nu = pPatch.wedgeCreator.Center(bt);
        pts = [arrow_center arrow_center + nu*arrow_len];
        plot(pts(1,:), pts(2,:), 'g-');
    end
end