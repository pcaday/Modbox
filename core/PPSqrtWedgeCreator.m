% Square-root wedge creator adapted for pseudopolar DFT.
classdef PPSqrtWedgeCreator < WedgeCreator
    properties
        ratio
        N
        empty
        mask
        angles
        dthetas
        centers
        smoothness
        aspect
        rMin
        rMinInside
        rMaxInside
        rMax
    end
    
    methods
        % swc = PPSqrtWedgeCreator(grid)
        % swc = PPSqrtWedgeCreator(grid, ratio)
        % swc = PPSqrtWedgeCreator(grid, ratio, rmin, rmax)
        % swc = PPSqrtWedgeCreator(grid, ratio, rmin, rmax, smoothness)
        % swc = PPSqrtWedgeCreator(grid, ratio, rmin, rmax, smoothness, class)
        %
        %  Create a WedgeCreator that creates wedges, the total number
        %   of which is proportional to sqrt(n), where n is the geometric
        %   average of the grid's dimensions.
        %
        %  Inputs
        %      grid: the frequency grid to use
        %     ratio: [optional] ratio between sqrt(n) and number of wedges.
        %                 The number of wedges is approximately
        %                 ratio*sqrt(n).
        % rmin,rmax: [optional] inside and outside radii for the
        %                 pseudo-radial cutoff
        %smoothness: [optional] smoothness of the angular cutoffs,
        %                 ranging from 0 (hard cutoff) to 0.5 (smoothest).
        %     class: [optional] string representing class of the wedge
        %                 functions, 'single' or 'double'. Default: single.
        function swc = PPSqrtWedgeCreator(grid, ratio, rMin, rMax, smoothness, valueClass)
            if nargin < 2 || isempty(ratio), ratio = 6; end
            if nargin < 3 || isempty(rMin), rMin = -inf; end
            if nargin < 4 || isempty(rMax), rMax = inf; end
            if nargin < 5 || isempty(smoothness), smoothness = 0.25; end
            if nargin < 6 || isempty(valueClass), valueClass = 'single'; end
            
            swc.ratio = ratio;
            swc.grid = grid;
            swc.smoothness = smoothness;
            
            assert(grid.dims == 2, 'Only 2D grids currently supported.');
            
            n = size(grid, 1) - 1;
            
            % Calculate total number of wedges;
            %  then, round it to a multiple of 8 for symmetry.
            swc.N = ceil(sqrt(n) * ratio);
            swc.N = 8 * ceil(swc.N / 8);
            
            swc.empty = zeros(grid.realSize, valueClass);
            
            % Create the (pseudo-)radial mask.
            xiR = repmattosize(linspace(0,1,n+1).', size(grid));
            rMinInside_ = 1.5*rMin;
            rMaxInside_ = 0.75*rMax;
            swc.mask = FullWindow(xiR, rMin, rMinInside_, rMaxInside_, rMax, 'smooth');
            
            swc.rMin = rMin;
            swc.rMax = rMax;
            swc.rMinInside = rMinInside_;
            swc.rMaxInside = rMaxInside_;
            
            % Make sure we always mask the zero frequency components.
            % swc.mask(1,:) = 0;
            
            % Get the angles of all the rays in the pseudopolar DFT.
            xis = grid(end,:);
            swc.angles = atan2(xis(2,:), xis(1,:));

            % Get the aspect ratio of the original grid.
            M = size(grid, 2) / 4;
            swc.aspect = grid(end, floor(M/2)+1);

            % Get the box centers
            dth = 2*pi/swc.N;
            sz = swc.aspect;
            ths = (0:swc.N-1) * dth;
            nus = [sz(1)*cos(ths); sz(2)*sin(ths)];
            nuNorms = sqrt(nus(1,:).^2 + nus(2,:).^2);
            swc.centers = bsxfun(@rdivide, nus, nuNorms);
            
            % Get the spacing between boxes
            thetas = atan2(swc.centers(2,:), swc.centers(1,:));
            swc.dthetas = mod(thetas - circshift(thetas, [0 1]), 2*pi);
        end
        
        % N = wedgeCreator.Count()
        %
        %   Return the total number of wedges, N.
        %
        function N = Count(wcr)
            N = wcr.N;
        end
        
        % f = wedgeCreator.Wedge(i)
        %
        %   Get the values of the frequency cutoff function f for wedge i.
        %   f will be an array of the cutoff function values on the
        %   grid used to construct the SqrtWedgeCreator.
        function f = Wedge(wcr, i)
            sm = wcr.smoothness;
            
            a = wcr.Center(i);
            a = atan2(a(2), a(1));
            rot_angle = mod(wcr.angles - a + pi, 2*pi) - pi;
            
            dthp = wcr.dthetas(mod(i,wcr.N)+1);
            dthm = wcr.dthetas(i);
            
            cutoff = FullWindow(rot_angle,  -dthm * (0.5 + sm), ...
                                            -dthm * (0.5 - sm), ...
                                             dthp * (0.5 - sm), ...
                                             dthp * (0.5 + sm), 'smooth');
            
            f = repmattosize(cutoff, size(wcr.grid));
            f = f .* wcr.mask;
        end
        
        % nu = wedgeCreator.Center(i)
        %
        %   Get a unit vector, nu, pointing toward the center of wedge i.
        %
        function nu = Center(wcr, i)
            nu = wcr.centers(:,i);
        end
        
        % nus = wedgeCreator.Centers
        %
        %   Get an n x k matrix of unit vectors
        %   pointing toward the center of the wedges.
        %
        %   The i'th column is the vector pointing to wedge i.
        %
        function nus = Centers(wcr)
            nus = wcr.centers;
        end
        
        % mask = wedgeCreator.Sum
        %
        %   Returns the sum of all the wedges (a wedge sum? :p)
        %
        function mask = Sum(wcr)
            mask = wcr.mask;
        end
        
        % boxes = wedgeCreator.NearestBoxes(xi)
        %
        %   Find the nearest box numbers to the given xis
        %
        %    xi is an n x i1 x ... x ik array
        % boxes is a  m x i1 x ... x ik array of box numbers
        function boxes = NearestBoxes(wcr, xi)
            N_ = wcr.N;
            ns = wcr.grid.ns(:);
            th = atan2(xi(2,:) / ns(2), xi(1,:) / ns(1));
            dth = 2*pi/N_;
            bx = th / dth;
            N_ = uint16(N_);
            bx_above = mod(uint16( ceil(bx)), N_) + 1;
            bx_below = mod(uint16(floor(bx)), N_) + 1;
            boxes = [bx_above; bx_below];
            boxes = reshape(boxes, size(xi));
        end
        
        function lpCeil = LowPassCeiling(wcr)
        % xiLP = wedgeCreator.LowPassCeiling()
        %   
        %   Return the largest frequencies not completely handled
        %  by this WedgeCreator. In other words, if the sum of the wedges
        %  does not add up to 1 at xi, then
        %
        %                          |xi_i| <= |xiLP_i|
        %
        %  for all i.
            maxFreq = max(wcr.grid.rGV);
            lpCeil = maxFreq * wcr.rMinInside;
            lpCeil = max(lpCeil, 0);    % We always cut off the zero
                                        % frequency regardless of rMin.
        end
    end
end