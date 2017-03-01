classdef SqrtWedgeCreator2D < WedgeCreator
    properties
        ratio
        N
        empty
        mask
        smoothness = -1         % Smoothness value: 0 to 0.5,
                                %  or -1 for rendering wedges as antialiased
                                %  triangles.
        rInsideMin
        rInsideMax
        rMin
        rMax
    end
    
    properties (SetAccess = protected)
        xiAngles
    end
    
    methods
        % swc = SqrtWedgeCreator2D(grid)
        % swc = SqrtWedgeCreator2D(grid, ratio)
        % swc = SqrtWedgeCreator2D(grid, ratio, rmin, rmax)
        % swc = SqrtWedgeCreator2D(grid, ratio, rmin, rmax, smoothness)        
        % swc = SqrtWedgeCreator2D(grid, ratio, rmin, rmax, smoothness, class)
        %
        %  Create a WedgeCreator that creates wedges, the total number
        %   of which is proportional to sqrt(n), where n is the geometric
        %   average of the grid's dimensions.
        %
        %  Inputs
        %       grid: the frequency grid to use
        %      ratio: [optional] ratio between sqrt(n) and number of wedges.
        %                  The number of wedges is approximately
        %                  ratio*sqrt(n).
        %  rmin,rmax: [optional] inside and outside radii for the
        %                  pseudo-radial cutoff
        % smoothness: [optional] smoothness of the angular cutoffs,
        %                  ranging from 0 (hard cutoff) to 0.5 (smoothest).
        %                  Can also be -1, where the wedge cutoff functions
        %                  are created as antialiased triangles (default).
        %      class: [optional] string representing class of the wedge
        %                  functions, 'single' or 'double'. Default: single.
        function swc = SqrtWedgeCreator2D(grid, ratio, rMin, rMax, smoothness, valueClass)
            if nargin < 2 || isempty(ratio), ratio = 4; end
            if nargin < 3 || isempty(rMin), rMin = -inf; end
            if nargin < 4 || isempty(rMax), rMax = inf; end
            if nargin < 5 || isempty(smoothness), smoothness = -1; end
            if nargin < 6 || isempty(valueClass), valueClass = class(grid.ns); end
            
            swc.ratio = ratio;
            swc.grid = grid;
            
            assert(grid.dims == 2, 'Expected a 2D grid.');
            
            n = prod(grid.realSize) .^ (1 / grid.dims);
            
            % Calculate total number of wedges;
            %  then, round it to a multiple of 8 for symmetry.
            swc.N = ceil(sqrt(n) * ratio);
            swc.N = 8 * ceil(swc.N / 8);
            
            swc.empty = zeros(grid.realSize, valueClass);
            
            maxFreq = abs(grid.mins);
            
            xiCF = grid.coordFuns();
            xiCF{1} = xiCF{1} / maxFreq(1);
            xiCF{2} = xiCF{2} / maxFreq(2);
            xiR = max(abs(xiCF{1}), abs(xiCF{2}));
            xiR = cast(xiR, valueClass);
            
            rInsideMin_ = 1.5*rMin;
            rInsideMax_ = 0.75*rMax;
            swc.mask = FullWindow(xiR, rMin, rInsideMin_, rInsideMax_, rMax, 'smooth');
            
            swc.rMin = rMin;
            swc.rMax = rMax;
            swc.rInsideMin = rInsideMin_;
            swc.rInsideMax = rInsideMax_;
            
            swc.smoothness = smoothness;
            if smoothness >= 0,
                % Precompute angles for later.
                swc.xiAngles = mod(atan2(xiCF{2},xiCF{1}), 2*pi);
            end
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
            global GlobalWedgeThreshold
            
            smoothness_ = wcr.smoothness;
            
            dth = 2*pi/wcr.N;
            ns = wcr.grid.ns;
            center = ceil((ns+1)/2);
            if smoothness_ == -1,
                x1 = center(1) + ns(1)*[NaN 0 cos((i-1.5)*dth) cos((i-0.5)*dth)];
                x2 = center(2) + ns(2)*[NaN 0 sin((i-1.5)*dth) sin((i-0.5)*dth)];
                f = OldRenderQuadInternal(wcr.empty, x1, x2);
            else
                smoothness_ = max(smoothness_, 0.0000001);
                    
                if i > 1,
                    f = FullWindow(wcr.xiAngles, (i-1.5-smoothness_)*dth, (i-1.5+smoothness_)*dth, (i-0.5-smoothness_)*dth, (i-0.5+smoothness_)*dth, 'smooth');
                else
                    f = HalfWindow(wcr.xiAngles, (0.5+smoothness_)*dth, (0.5-smoothness_)*dth, 'smooth')...
                      + HalfWindow(wcr.xiAngles, 2*pi-(0.5+smoothness_)*dth, 2*pi-(0.5-smoothness_)*dth, 'smooth');
                end
            end
            
            f = f .* wcr.mask;
            f(f < GlobalWedgeThreshold) = 0;
        end
        
        % 
        
        % nu = wedgeCreator.Center(i)
        %
        %   Get a unit vector, nu, pointing toward the center of wedge i.
        %
        function nu = Center(wcr, i)
            dth = 2*pi/wcr.N;
            sz = -wcr.grid.mins;   % Note: assuming we have a real DFT grid here!
            th = (i-1) * dth;
            nu = [sz(1)*cos(th); sz(2)*sin(th)];
            nu = nu / norm(nu);
        end
        
        % nus = wedgeCreator.Centers
        %
        %   Get an n x k matrix of unit vectors
        %   pointing toward the center of the wedges.
        %
        %   The i'th column is the vector pointing to wedge i.
        %
        function nus = Centers(wcr)
            dth = 2*pi/wcr.N;
            sz = -wcr.grid.mins;   % Note: assuming we have a real DFT grid here!
            ths = (0:wcr.N-1) * dth;
            nus = [sz(1)*cos(ths); sz(2)*sin(ths)];
            nuNorms = sqrt(nus(1,:).^2 + nus(2,:).^2);
            nus = bsxfun(@rdivide, nus, nuNorms);
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
        %   |xi_i| <= |xiLP_i| for all i.
            maxFreq = abs(wcr.grid.mins);
            lpCeil = maxFreq * wcr.rInsideMin;
        end

    end
end