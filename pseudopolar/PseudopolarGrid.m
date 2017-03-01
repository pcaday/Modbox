% Pseudopolar grid support.
%
% The grid is organized with angles starting at 0 and increasing
% counterclockwise to 2*pi. (i.e., the ordinary way)
classdef PseudopolarGrid < IrregularGrid
    properties
        nr
        nAngles
        dr
        rGV
        anglesGV
    end
    
    
    methods
        function ppGrid = PseudopolarGrid(grid, oversampling)
            if nargin < 2, oversampling = 1; end
            
            % Check inputs
            assert(grid.dims == 2, 'Only 2D grids supported.');
            [a,b] = size(grid);
            assert(a == b, 'Grid must be square for pseudopolar DFT.');

            % Get sizes
            N = double(ceil(a*oversampling));    % double is here because circshift can't handle single.
            M = N;
            
            % Construct the list of frequencies in the PP grid
            % #radii:  N+1
            % #angles: 4*M
            L = grid.maxes - grid.mins + grid.ds;
            ns = grid.ns;
            maxfreq = ns * pi ./ L;
            points = PseudopolarGrid.GeneratePseudopolarGrid(maxfreq, M, N);
            
            % Superclass constructor to create the IrregularGrid
            ppGrid = ppGrid@IrregularGrid(points);
            
            % Store information about the pseudopolar grid in the
            % properties.
            maxfreq = max(maxfreq);
            ppGrid.dr = maxfreq/N;
            ppGrid.rGV = (0:N)*maxfreq/N;
            xis = squeeze(points(:,end,:));
            ppGrid.anglesGV = mod(atan2(xis(2,:), xis(1,:)), 2*pi).';
            
            ppGrid.nr = N + 1;
            ppGrid.nAngles = 4*M;
        end
    end    
    
    methods (Static)
        function [freqs,r,x] = GeneratePseudopolarGrid(maxfreq, M, N)
            x = 1 - 2*(0:(M-1))/M ;  % linspace(1, -1, M);
            r = linspace(0, 1, N+1).';
            tx = bsxfun(@times, x, r);
            tr = repmat(r, size(x));
            
            freqs = fcat([tx -tr -tx  tr] * maxfreq(1), ...
                [tr  tx -tr -tx] * maxfreq(2));
            freqs = circshift(freqs, [0 0 double(ceil((M-1)/2))]);
        end
    end
end