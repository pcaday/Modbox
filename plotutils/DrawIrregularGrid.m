function DrawIrregularGrid(gridPts, N, doLegend)
dims = ndims(gridPts) - 1;
if nargin < 3 || isempty(doLegend), doLegend = false; end

switch dims
    case 2
        if nargin < 2 || isempty(N), N = 15; end

        Nx = size(gridPts, 2);
        Ny = size(gridPts, 3);
        % Choose evenly-spaced indices to draw grid lines at.
        Xi = round(linspace(1, Nx, min(N,Nx)));
        Yi = round(linspace(1, Ny, min(N,Ny)));
        
        % And plot them!
        hold on
        for i = Yi
            hx = plot(squeeze(gridPts(1,:,i)), squeeze(gridPts(2,:,i)),'r-');
        end
        for i = Xi
            hy = plot(squeeze(gridPts(1,i,:)), squeeze(gridPts(2,i,:)),'b-');
        end
        
        if doLegend,
            legend([hx hy], {'x_1' 'x_2'});
        end
    case 3
        if nargin < 2 || isempty(N), N = 6; end

        Nx = size(gridPts, 2);
        Ny = size(gridPts, 3);
        Nz = size(gridPts, 4);
        % Choose evenly-spaced indices to draw grid lines at.
        Xi = round(linspace(1, Nx, min(N,Nx)));
        Yi = round(linspace(1, Ny, min(N,Ny)));
        Zi = round(linspace(1, Nz, min(N,Nz)));
        
        % And plot them!
        hold on
        for i = Yi
            for j = Zi
                hx = plot3(squeeze(gridPts(1,:,i,j)), squeeze(gridPts(2,:,i,j)), squeeze(gridPts(3,:,i,j)), 'r-');
            end
        end
        for i = Zi
            for j = Xi
                hy = plot3(squeeze(gridPts(1,j,:,i)), squeeze(gridPts(2,j,:,i)), squeeze(gridPts(3,j,:,i)), 'g-');
            end
        end
        for i = Xi
            for j = Yi
                hz = plot3(squeeze(gridPts(1,i,j,:)), squeeze(gridPts(2,i,j,:)), squeeze(gridPts(3,i,j,:)), 'b-');
            end
        end
        
        if doLegend,
            legend([hx hy hz], {'x_1' 'x_2', 'x_3'});
        end
    otherwise
        error('DrawIrregularGrid:dims', 'Only 2D and 3D grids supported.');
end

end