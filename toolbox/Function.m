% Function class, representing a function sampled on a grid.
%
% The functions can be scalar-, vector-, or matrix-valued.
%
classdef Function < handle & matlab.mixin.Copyable
    properties
        grid = []        % Grid for this function (class Grid)
        f = []           % Function values (n-dimensional array, with values
                         %   corresponding to the points of the grid).
        m = [1]                                       %#ok<NBRAK> 
                         % Dimension(s) of the codomain (a vector of any size)
                         %   Default: [1] (scalar function)
        description = '' % Optional human-readable description of this
                         %  function.
        defaultInterpolationMethod = 'cubic';
        defaultExtrapolationMethod = 'zero';
        % Default interpolation and extrapolation method used in the
        %  sampleAt method - can be any of the methods
        %  griddedInterpolant knows about:
        %    interpolation: 'nearest', 'linear', 'pchip' (1D), 'cubic',
        %                     'spline'
        %    extrapolation: 'nearest', 'linear', 'pchip' (1D), 'cubic',
        %                     'spline', 'none'
        %
        % Also, the extrapolation method can be 'zero', indicating the
        %  function should be assumed to be zero outside its domain.
        %
        % For R2012b and below, griddedInterpolant doesn't do extrapolation
        %  with the 'nearest', 'linear', or 'cubic' methods, and does with
        %  the 'spline' and 'pchip' methods.
        % In this case, defaultExtrapolationMethod does not affect the
        %  extrapolation method, except that 'zero' extrapolation does work
        %  if the interpolation method is 'nearest', 'linear', or 'cubic'.
        
        periodic = false % True to treat this function as a periodic
                         %  function by default for sampling (see SampleAt
                         %  method)
    end
    
    properties (Access = protected)
        interpolant = [] % Instance of the griddedInterpolant class which
                         %  performs interpolation of the function.
    end
    
    % Constructors
    methods (Static)
        function fun = Zeros(newGrid, valueClass, newM)
        % Create a zero function.
        %  With no arguments, this is the same as the default constructor.
        %
        % valueClass is an optional string giving the class of the new
        %   function [default: 'double'].
        % m is the dimension(s) of the function's codomain
        %   [default: 1 (scalar function)].
        
            % Substitute default class if not given
            if nargin < 2 || isempty(valueClass)
                valueClass = 'double';
            end
            % Substitute default newM if not given or empty
            if nargin < 3 || isempty(newM)
                newM = [1];                                 %#ok<NBRAK>
            end
            
            CheckTypes({valueClass,newM}, {'char','integer'}, {2 3})
            assert(isa(newGrid, 'Grid') || isa(newGrid, 'IrregularGrid'),...
                'Argument 1 is not a grid.');
            
            fun = Function();
            
            fun.grid = newGrid;
            fun.m = newM;
            fun.f = zeros([realSize(newGrid) newM], valueClass);
        end
        
        function fun = Constant(newGrid, value, newM)
        % Create a constant function.
        %
        %   f = Function.Constant(grid, c);
        %   f = Function.Constant(grid, c, m);
        %
        % grid: Grid object for the function
        %    c: Constant value
        %    m: (optional) Dimension(s) of the function's codomain.
        %         Defaults to 1 (scalar function).
        
            % Substitute default newM if not given or empty
            if nargin < 3 || isempty(newM)
                newM = [1];                                 %#ok<NBRAK>
            end
            
            CheckTypes({newM}, {'integer'}, {3})
            assert(isa(newGrid, 'Grid') || isa(newGrid, 'IrregularGrid'),...
                'Argument 1 is not a grid.');
            
            fun = Function();
            
            fun.grid = newGrid;
            fun.m = newM;
            fun.f = repmat(value,[realSize(newGrid) newM]);
        end
        
        function fun = WithValues(newGrid, newValues)
        % Create a function, initialized with specific values.
        % The dimension(s) of the function's codomain are automatically
        %  determined.
            CheckTypes({newValues}, {'numeric'}, 2)
            assert(isa(newGrid, 'Grid') || isa(newGrid, 'IrregularGrid'),...
                'Argument 1 is not a grid.');
            
            n = newGrid.dims;
            assert(ndims(newValues) >= n,...
                'Function:WithValues:NotEnoughDimensions',...
                'Not enough dimensions in the initialization values');
            sz = size(newValues);
            if (n == 1 && sz(2) == 1), sz(2) = []; end   % Allow initialization via column vector
            assert(all(sz(end-n+1:end) == realSize(newGrid)),...
                'Function:WithValues:WrongSize',...
                'Initialization values are the wrong size for the grid.');
            newM = sz(1:end-n);
            if isempty(newM)
                newM = [1];                              %#ok<NBRAK>
            end
            
            fun = Function();
            fun.grid = newGrid;
            fun.m = newM;
            fun.f = newValues;
        end
        
        function fun = WithHandle(newGrid, handle)
        % Create a function, initialized using a function handle.
        % The dimension(s) of the function's codomain are automatically
        %  determined.
        %
        %     f = Function.WithHandle(grid, function_handle)
        %
        % function_handle will be passed the complete set of grid points
        %  as an n x (size of grid) array. 'coord' can be used to extract
        %  given coordinates.
        %
        
        % In 1D, if function_handle returns a row vector, the result is
        %  transposed to a column vector, to allow for simpler
        %  initialization, i.e.
        %    
        %     function_handle = @exp
        %
        %  instead of
        %
        %     function_handle = @(x) exp(coord(x,1));
        %
        %
            CheckTypes({newGrid,handle}, {'Grid','function_handle'})
            
            vals = handle(newGrid.AllPoints());
        %    if newGrid.dims == 1, vals = vals(:); end   
            fun = Function.WithValues(newGrid, vals);
        end
        
        function fun = RandomScalar(newGrid)
            CheckTypes({newGrid}, {'Grid'});
            
            fun = Function();
            fun.grid = newGrid;
            fun.m = 1;
            fun.f = rand(size(newGrid));
        end
    end
    
    
    %
    %
    methods
        function ft = transpose(f)
            ft = f.copy();
            ft.f = ft.f.';
        end
        
        function ft = ctranspose(f)
            ft = f.copy();
            ft.f = ft.f';
        end
        
        function n = norm(f, p)
            if nargin < 2, p = 2; end
            
            n = norm(f.f(:), p);
            if isfinite(p),
                dx = prod(f.grid.ds);
                n = n * (dx ^ (1 / p));
            end
        end
        
        function r = real(f)
        % Return the real part of the Function, as a new Function
            r = f.copy();
            r.f = real(r.f);
        end
        
        function r = imag(f)
        % Return the imaginary part of the Function, as a new Function
            r = f.copy();
            r.f = imag(r.f);
        end
        
        function b = isreal(f)
        % Check if this Function has an imaginary component.
            b = isreal(f.f);
        end
        
        function g = abs(f)
        % Return the absolute value of this Function, as a new Function
            g = f.copy();
            g.f = abs(g.f);
        end
        
        function l = log(f)
        % Return the log of this Function, as a new Function.
            l = f.copy();
            l.f = log(l.f);
        end
        
        function l = plog(f)
        % Return a 'pseudo-polar log', (f/|f|) * log(|f|/M), where
        %  M = min {|f(x)| : f(x) nonzero}.
            l = f.copy();
            lf = log(abs(l.f));
            lf = lf - min(lf(isfinite(lf)));
            l.f = sign(l.f) .* lf;
        end
        
        function m = modify(f, modifier)
        % Modify this Function using a function handle.
        %  (convenience function).
        %
        %    m = f.modify(function_handle);
        %
        % The function handle is passed the raw data of f (i.e. f.f)
        %  and returns the raw data of m.
            m = f.copy();
            m.f = modifier(f.f);
        end
        
        function fc = Component(f, inds)
            k = length(inds);
            assert(k == length(f.m), 'Component has the wrong number of indices.');
            
            subs = [repmat({':'}, [1 f.grid.dims]) num2cell(k(:).')];
            fc = f.copy();
            fc.f = f.f(subs{:});
            fc.m = 1;
        end
        
        function fs = Slice(f, inds)
        % Return a slice of this Function along one or more variables.
        %
        %   slice = f.Slice(inds)
        %
        % inds is a length-n cell array (where n is the # of dimensions in the
        %  domain).
        % inds{i} may be either a scalar representing a value for variable i
        %  or [] indicating that variable i is free.
        % At least one variable must be free.
        %
        % For example, if f is a 3D Function on the domain [-1,1]^3, then
        %    f.Slice({0.5, [], -0.5})
        % returns a 1D function on the y grid [-1,1], containing the values
        %  of f at x = 0.5, z = -0.5.
        %
        % If inds{i} is not on a whole grid point, it will be rounded
        %  to the nearest grid point.
        % If it is outside the grid entirely, it will be set to the nearest
        %  grid edge, and a warning issued.
            
            sliceVars = cellfun(@isempty, inds);
            assert(any(sliceVars), 'Function:Slice:noSliceVars', 'at least one variable must be free.');
            
            inds(sliceVars) = {0};
            try
                inds = [inds{:}];
            catch
                assert('Function:Slice:wrongInput', 'input must be a cell array containing scalars and empty matrices.');
            end
            
            assert(length(inds) == f.grid.dims, 'Function:Slice:wrongSize', 'input does not have the correct number of entries.');
            inds = round(f.grid.toGridSpace(inds(:))).';
            
            if any((inds <= 0 | inds > f.grid.ns) & ~sliceVars),
                warning('Function:Slice:outOfGrid', 'requested slice is outside grid, choosing nearest slice.');
                inds(inds <= 0) = 1;
                inds(inds > f.grid.ns) = f.grid.ns(inds > f.grid.ns);
            end
            
            inds = num2cell(inds);
            inds(sliceVars) = {':'};
            gs = f.grid.Slice(find(sliceVars));           %#ok<FNDSB>  <-- Code Analyzer is confused.
            fsVals = f.f(inds{:});
            sz = [size(fsVals) ones(1,f.grid.dims)];
            sz(~sliceVars) = [];
            fsVals = reshape(fsVals, sz);
            fs = Function.WithValues(gs, fsVals);
        end
            
        function spy(f, varargin)
        % View the sparsity pattern (support) of this Function
            assert(f.grid.dims == 2, 'spy only works for 2D Functions.');

            spy(f.f.', varargin{:})
            axis normal   % Undo the axis scaling constraints that 'spy' does
            axis xy
        end
        
        
        function values = SampleAt(fun, points, varargin)
        % Sample the function (using interpolation) at given points.
        %  Uses MATLAB's built-in griddedInterpolant class to do the
        %  hard work.
        %
        %   f.SampleAt(points, ['InterpolationMethod', inMethod,]
        %                      ['ExtrapolationMethod', exMethod,]
        %                      ['Periodic', periodic]);
        %
        % Inputs
        %   points: an n x i1 x i2 x ... x ik array of points
        %            (where n is the dimension of the domain)
        % inMethod: (optional) the interpolation method to use, overriding
        %            the function's defaultInterpolationMethod.
        % exMethod: (optional) the extrapolation method to use, overriding
        %            the function's defaultExtrapolationMethod. Only
        %            applies to non-periodic functions!
        % periodic: (optional) whether to treat the function as periodic
        %            for sampling, overriding the value of the grid's
        %            periodic property.
        %
        % Outputs
        %   values: an m1 x ... x mK x i1 x i2 x ... x ik array of function
        %            values at the input points. Here m1 x ... mK are the
        %            dimensions of the codomain of f.
        %
        % Currently only applies to scalar functions -- needs to be fixed!
            persistent InterpolantsHaveExtrapMethod
            
            assert(isa(fun.grid, 'Grid'), 'SampleAt requires a regular grid.');
            
            p = inputParser;
            addParamValue(p, 'InterpolationMethod', fun.defaultInterpolationMethod);
            addParamValue(p, 'ExtrapolationMethod', fun.defaultExtrapolationMethod);
            addParamValue(p, 'Periodic', fun.grid.periodic, @islogical);
            
            parse(p, varargin{:});
            isPeriodic = p.Results.Periodic;
            method = p.Results.InterpolationMethod;
            exMethod = p.Results.ExtrapolationMethod;
            
            zeroExtrap = strcmp(exMethod, 'zero');
            if zeroExtrap, exMethod = 'none'; end
                
            % Check the dimensions of the array of input points.
            n = fun.grid.dims;
            originalShape = size(points);
            assert(n == 1 || n == originalShape(1), 'Function:SampleAt:DimMismatch',...
                'The input points must have the same dimension as the function''s domain.');
            newShape = originalShape;
            if (n == newShape(1))
                % Slice off the first dimension (unless it's not 1 and our
                %   domain is 1D)
                newShape(1) = [];
            end
            if length(newShape) < 2
                % Add a trailing 1 if necessary to make the shape vector's
                % size at least 2.
                newShape = [newShape 1];
            end
            
            % Do the interpolants already exist?
            if ~isempty(fun.interpolant)
                % Yep! Update the interpolation method
                %  and function values
                fun.interpolant.Method = method;
                fun.interpolant.Values = fun.f;
            else
                % Nope! Create new ones.
                %sz1 = fun.grid.NPoints();
                %for i = 1:length(fun.m(:))
                fun.interpolant = griddedInterpolant(...
                    fun.grid.gridVectors(), fun.f, method);
                %end
            end
            
            % Set the extrapolation method if supported by MATLAB.
            if isempty(InterpolantsHaveExtrapMethod),
                % Cache whether we have the ExtrapolationMethod property
                %  because isprop is slow!
                InterpolantsHaveExtrapMethod = isprop(fun.interpolant, 'ExtrapolationMethod');
            end
            if InterpolantsHaveExtrapMethod,
                fun.interpolant.ExtrapolationMethod = exMethod;
            end
            
            % Reshape the points into a matrix
            points = reshape(points, n, []);
            
            % If the periodic option was specified, wrap the points
            if numel(isPeriodic) == 1,
                isPeriodic = repmat(isPeriodic, [1 n]);
            end
            for i = 1:n
                if isPeriodic(i),
                    xiMin = fun.grid.mins(i);
                    xiMax = fun.grid.maxes(i);
                    points(i,:) = xiMin + mod(points(i,:)-xiMin, xiMax-xiMin);
                end
            end
            
            % Sample the function at the points
            values = fun.interpolant(points.');
            % Reshape the values to match the inputs
            values = reshape(values, newShape);
            
            % If zero extrapolating, the points landing outside the
            %  function's grid currently have NaN values. Convert those to
            %  zero.
            if zeroExtrap,
                values(isnan(values)) = 0;
            end
        end
        
        function nf = Pullback(fun, yGrid, yGridImages)
        % Compute pullback (coordinate changes)
        %
        % Inputs
        %          yGrid: new grid (class Grid)
        %    yGridImages: array of images of the yGrid's points in f's own
        %                  grids. Can be a matrix or an (n+1)-dimensional
        %                  array.
        %
        %                 Matrix Form:
        %                 ------------
        %                 yGridImages(:,i) is the image of the i'th point
        %                  in yGrid, using linear indexing.
        %
        %                 (n+1)-D Array Form:
        %                 -------------------
        %                 yGridImages(:,i1,i2,...,in) is the image of the
        %                  point represented by (i1,...,in).
        %
        % Outputs
        %             nf: fun with coordinate changes
        %                   (Function, on the grid specified by yGrid)
        
            CheckTypes({yGrid,yGridImages}, {'Grid','numeric'});
            
            n = fun.grid.dims;
            if ismatrix(yGridImages),
                assert(size(yGridImages, 2) == yGrid.NPoints,...
                    'Function:CoordChange:wrongRowCount',...
                    'yGridImages doesn''t have the right number of columns.');
                yGridImages = reshape(yGridImages, [n realSize(yGrid)]);
            else
                assert(all([n realSize(yGrid)] == size(yGridImages)), ...
                    'Function:CoordChange:dimNotConsistent', ...
                    'yGrid and yGridImages don''t have matching dimensions.');
            end
            
            % Function.SampleAt should create the interpolant if not already created,
            %  then ask it to interpolate.
            interpolatedVals = fun.SampleAt(yGridImages);
            nf = Function.WithValues(yGrid, interpolatedVals);
        end
        
        function nf = Resample(fun, newGrid, varargin)
        % Sample the function on a new grid.
        %
        %    nf = f.Resample(newGrid, ...)
        %
        % Inputs:
        %     newGrid: Grid object
        %
        %     Additional inputs are passed to SampleAt.
        %
        % Outputs:
        %          nf: Function object on the new grid.
        %
            CheckTypes({newGrid}, {'Grid'});
            
            interpolatedVals = fun.SampleAt(newGrid.AllPoints(), varargin{:});
            nf = Function.WithValues(newGrid, interpolatedVals);
        end
        
        
        function fHat = DFT(f, zeroPadding)
        % Compute the DFT of the function.
        %
        %  fHat = f.DFT([zeroPadding])
        %
        %    zeroPadding is an optional argument. If specified, the
        %     dimensions of f are multiplied by zeroPadding and rounded up
        %     before applying the FFT.
        %
        %  fHat is returned as a new Function on a new grid.
        %     If zeroPadding > 1, its grid is bigger than f's
        %     by a factor of zeroPadding (up to rounding).
        
            assert(isa(f.grid, 'Grid'), 'DFT requires a regular grid.');
            
            % Default zero-padding ratio is 1 (no padding).
            if nargin < 2, zeroPadding = 1; end
            assert(zeroPadding >= 1, 'Function:FFT:zpSmall',...
                'zeroPadding must be at least 1.0');
            
            % Check types on input
            CheckTypes({f, zeroPadding}, {'Function', 'scalar'});
            
            % Create a grid for the FFT function
            fHatGrid = f.grid.DFTGrid(zeroPadding);
            paddedNs = size(fHatGrid);  % not realSize -- fftn will complain
            
            % Compute the (shifted) FFT.
            fftvals = fftshift(fftn(f.f, paddedNs));
            
            % Create the new function.
            fHat = Function.WithValues(fHatGrid, fftvals);
        end
        
        
        
        function f = IDFT(fHat, fGrid)
        % Compute the inverse DFT of the function.
        %
        %  f = fHat.IDFT(fGrid)
        %
        %  fGrid is the Grid object to use for the inverse FFT, f.
        %
        %  f is returned as a new Function on the grid fGrid.
        %
        %     If fHat's grid is bigger than fGrid (for instance, from
        %  zero-padding in the DFT method), f will be trimmed down to the
        %  size of fGrid.
        
            assert(isa(fHat.grid, 'Grid'), 'IDFT requires a regular grid.');
            
            % Check data types
            CheckTypes({fHat,fGrid}, {'Function','Grid'});
            
            % Perform inverse FFT
            ifftvals = ifftn(ifftshift(fHat.f));
            
            % If the original function was zero-padded before applying the FFT,
            %  trim down the inverse FFT to the function's original size.
            if any(realSize(fGrid) ~= realSize(fHat.grid))
                % The dimension here is arbitrary, so construct a cell array
                %  with the subscripts we want, then do the subscripting by
                %  converting the cell array to a comma-separated list.
                % This last idea is from Nikola Sprljan,
                %  http://www.mathworks.com/matlabcentral/newsreader/view_thread/111351
                % I believe it's also in the source of a MATLAB built-in
                % function.
                sz = realSize(fGrid);
                subs = arrayfun(@(i) 1:sz(i), 1:length(sz), 'UniformOutput', false);
                ifftvals = ifftvals(subs{:});
            end
            
            % Create the function itself
            f = Function.WithValues(fGrid, ifftvals);
        end
		
        
        function h = plot(fun, vrange, cx, swapXY)
        % Plot this function:
        %
        %     f.PLOT
        %     f.PLOT([min max])
        %     f.PLOT([], cxflag[, swapXY])
        %     f.PLOT([min max], cxflag[, swapXY])
        %
        % In the last two forms, if cxflag is true, the function is plotted
        %  as a complex function (using compleximagesc), while if cxflag
        %  is false, just the real part is displayed.
        %
        % If cxflag is omitted, the value of the boolean global variable
        %  GlobalFunctionCxPlot is used instead (if unset, it defaults
        %  to false).
        %
        % swapXY is an optional flag: if true, the x- and y-axes are
        %  swapped; if false or omitted, they remain the same.
            global GlobalFunctionCxPlot
            
            if nargin < 3, cx = GlobalFunctionCxPlot; end
            if isempty(cx), cx = false; end
            
            if nargin < 4, swapXY = false; end
            
            irregular = ~isa(fun.grid, 'Grid');
            
            assert(prod(fun.m) == 1, 'Function:plot:vectorValued',...
                'Plotting vector-valued functions not implemented.');
            switch fun.grid.dims
                case 1
                    if cx,
                        warning('Function:plot:complex1d',...
                            'Complex 1D plotting not implemented.');
                        cx = false;
                    end
                    
                    xlimits = [fun.grid.mins, fun.grid.maxes];
                    X = fun.grid.coordFuns{1};
                    Y = real(fun.f);
                    
                    if swapXY,
                        hh = plot(Y,X);
                    
                        ylim(xlimits);
                        if nargin > 1 && ~isempty(vrange),
                            xlim(vrange);
                        end
                    else
                        hh = plot(X,Y);
                    
                        xlim(xlimits);
                        if nargin > 1 && ~isempty(vrange),
                            ylim(vrange);
                        end
                    end
                case 2
                    ranges = {};
                    if ~irregular,
                        xrange = [fun.grid.mins(1) fun.grid.maxes(1)];
                        yrange = [fun.grid.mins(2) fun.grid.maxes(2)];
                        ranges = {xrange, yrange};
                    end
                    % Transpose the function to switch x and y coordinates.
                    %   x is the first index (rows) which MATLAB shows as y
                    %   and y is the second index (column index) which
                    %   MATLAB shows as x.
                    if cx,
                        imfunc = @compleximagesc;
                        im = fun.f;
                    else
                        imfunc = @imagesc;
                        im = real(fun.f);
                    end
                    
                    if swapXY,
                        im = im.';
                        ranges = fliplr(ranges);
                    end
                    
                    if nargin > 1 && ~isempty(vrange),
                        hh = imfunc(ranges{:}, im.', vrange);
                    else
                        hh = imfunc(ranges{:}, im.');
                    end
                    axis xy;
                    
                otherwise
                    assert(false, 'Function:plot:tooManyD', ...
                        'Plotting only implemented for 1D and 2D functions.');
            end
            
            % Add a disclaimer title if the function is not real,
            %  stating that we're only showing the real part.
            if ~cx && any(imag(fun.f(:)))
                title('Real part of function.')
            end
            
            % Only return an output argument if one is requested
            %  (to avoid 'ans' output)
            % This technique is from someone else (don't remember who!)
            if nargout > 0,
                h = hh;
            end
        end
        
        function h = semilogy(fun, vrange, cx)
            % Log plot of a function.
            %     f.SEMILOGY
            %     f.SEMILOGY([min max])
            %     f.SEMILOGY([], cxflag[, swapXY])
            %     f.SEMILOGY([min max], cxflag[, swapXY])
            %
            %  1D log plots are drawn with
            %  semilogy; 2D log plots are drawn by displaying the log
            %  of the function as a (complex) image.
            %
            % For 2D log plots, [min max] can be replaced by a single value
            %  s, in which case max will be automatically
            %  determined, and min = max - s.
            %
            assert(prod(fun.m) == 1, 'Function:semilogy:vectorValued',...
                'Plotting vector-valued functions not implemented.');

            global GlobalFunctionCxPlot
            
            if nargin < 2, vrange = []; end
            if nargin < 3, cx = GlobalFunctionCxPlot; end
            if isempty(cx), cx = false; end
            
            if (fun.grid.dims == 1),
                xlimits = [fun.grid.mins, fun.grid.maxes];
                X = fun.grid.coordFuns{1};
                Y = real(fun.f);
                
                hh = semilogy(X,Y);
                
                xlim(xlimits);
                if ~isempty(vrange),
                    ylim(vrange);
                end
            else
                if cx,
                    fmod = fun.plog;
                    absmax = max(abs(fmod.f(isfinite(fmod.f))));
                else
                    fmod = fun.abs.log;
                    absmax = max(fmod.f(isfinite(fmod.f)));
                end
                
                if isempty(vrange),
                    vrange = [absmax - 20 absmax];
                elseif length(vrange) == 1,
                    vrange = [absmax-vrange absmax];
                else
                    vrange = log(abs(vrange));
                end
                
                fmod.plot(vrange, cx);
            end
            
            if nargout > 0, h = hh; end
        end
        
        function ap = animate(fun,tvar,scale,varargin)
            % Produce an animated plot of a time-dependent Function.
            %
            %    f.animate
            %    f.animate(tvar)
            %    f.animate(tvar,scale)
            %    f.animate(tvar,scale,...)
            %    f.animate([],scale,...)
            %    f.animate(tvar,[],...)
            %    f.animate([],[],...)
            %    ap = f.animate(...)
            %
            % tvar, if present, indicates the index of the time variable.
            %  If omitted or [], the last variable is assumed to be time.
            %
            % scale indicates the units of the time axis, in seconds.
            %  If omitted or [], it defaults to 1 second.
            %
            % The function is displayed as with the PLOT method. Any extra
            %  arguments are passed to PLOT.
            %
            % animate optionally returns a structure ap that can be
            %  passed to the animatestop method. Otherwise, the animation
            %  continues looping until the axis is closed or the plot is
            %  overwritten.
            %
            ndim = fun.grid.dims;
            assert(ndim > 1, 'Function:TimePlot:only1D',...
                'Time plot only supported for >1D Functions.');
            
            if nargin < 2 || isempty(tvar), tvar = ndim; end;
            if nargin < 3 || isempty(scale), scale = 1.; end;
            
            
            dt = fun.grid.ds(tvar);
            tFrame = dt * scale;
            tFrame = max(1,round(tFrame*1000)) * 0.001;   % Round to nearest millisecond to avoid warnings from MATLAB
            
            % Stuff the data for the timer task to find.
            tData.tvar = tvar;
            tData.fun = fun;
            tData.i = 1;
            tData.hAxes = gca;
            tData.plotArgs = varargin;
            tData.ih = [];
            
            % Create a timer to produce the animation,
            t = timer('BusyMode', 'queue', 'ExecutionMode', 'fixedSpacing', ...
                'Period', tFrame, 'TimerFcn', @Function.AnimateTask,...
                'UserData', tData);          
            
            start(t);
            
            % Return a structure with the timer if output requested.
            if nargout > 0,
                ap.timer = t;
            end
        end
        
        
    end
    
    methods (Static)
        function animatestop(ap)
            % Stop a time animation created by the animate method.
            % This is a static method:
            %
            %  Function.animatestop(ap)
            %
            % ap is the return value of the animate method.
            try
                stop(ap.timer);
                delete(ap.timer);
            catch
                error('Function:animatestop:badAP', 'Invalid input or error stopping animation.');
            end
        end
        
        function AnimateTask(t, ~)
            % Helper function (timer callback) for the Function.animate method.
            
            % t is a timer object.
            % Get data from t:
            tData = get(t, 'UserData');
            
            tvar = tData.tvar;
            fun = tData.fun;
            i = tData.i;
            hAxes = tData.hAxes;
            ih = tData.ih;
            
            if (~isempty(ih) && ~ishandle(ih)),
                % Our image was deleted, probably because someone
                % wanted to overwrite the plot. Abort!
                stop(t);
                delete(t);
                return;
            end
            
            nt = fun.grid.ns(tvar);
            ndim = fun.grid.dims;
            
            % Prepare the slice:
            slice = cell(1,ndim);
            slice{tvar} = fun.grid.fromGridSpace(i, tvar);
            
            % Plot the slice on the already prepared axes.
            if ishandle(hAxes),
                oldAx = gca;
                axes(hAxes);
                if ndim == 4,
                    hold off
                    %fun.Slice(slice).vol3d(tData.plotArgs{:});
                    %view(3);
                    fs = fun.Slice(slice).abs;
                    if isempty(tData.plotArgs),
                        tData.plotArgs = {max(fs.f(:)) * 0.5};
                    end
                    if ~isempty(tData.ih), delete(tData.ih); end
                    tData.ih = fs.isosurface(tData.plotArgs{:});
                else
                    tData.ih = fun.Slice(slice).plot(tData.plotArgs{:});            
                end
                title(sprintf('t = %f', slice{tvar}));
                axes(oldAx);
                
                tData.i = mod(i, nt) + 1;
                set(t, 'UserData', tData);
            else
                % Oops! Our axes got deleted. Stop the animation.
                stop(t);
                delete(t);
            end
        end
        
    end
    
    methods
        function h = isosurface(fun, isovalue)
            assert(prod(fun.m) == 1, 'Function:isosurface:vectorValued',...
                'This function isn''t a scalar function.');
            assert(isa(fun.grid, 'Grid'), 'Function:isosurface:irregular',...
                'isosurface requires a regular grid.');
            assert(fun.grid.dims == 3, 'Function:isosurface:not3D',...
                'isosurface only works for 3D Functions');
            
            % Following MATLAB's doc page for isosurface...
            
            cfs = fun.grid.coordFuns();
            X = cfs{1};
            Y = cfs{2};
            Z = cfs{3};
            % We have to permute the first and second indices because
            %  isosurface orders the indices (Y,X,Z) [meshgrid-style].
            V = permute(fun.f, [2 1 3]);
            X = permute(X, [2 1 3]);
            Y = permute(Y, [2 1 3]);
            Z = permute(Z, [2 1 3]);
            % Set the axis limits (may not want to keep this here)
            xlim([fun.grid.mins(1), fun.grid.maxes(1)]);
            ylim([fun.grid.mins(2), fun.grid.maxes(2)]);
            zlim([fun.grid.mins(3), fun.grid.maxes(3)]);
            
            h = patch(isosurface(X,Y,Z,V,isovalue));
            isonormals(X,Y,Z,V,h)
            set(h, 'FaceColor', 'blue', 'EdgeColor', 'none')
            camlight
            lighting gouraud
        end
        
        function h = vol3d(fun, varargin)
            % 3D volumetric plot, using `vol3d` from the MATLAB File
            % Exchange.
            %
            %   f.VOL3D
            %   f.VOL3D(...)
            %   h = f.VOL3D(...)
            %
            % Any arguments are passed to vol3d.
            hh = vol3d('CData', fun.f, varargin{:});
            
            if nargout > 0, h = hh; end;
        end
        
        function h = sliceplot(fun,varargin)
            % Create a slice plot of a 3D Function using `slice`.
            %
            %   f.SLICEPLOT(Sx,Sy,Sz)
            %   f.SLICEPLOT(XI,YI,ZI)
            %   f.SLICEPLOT(...)
            %   h = f.SLICEPLOT(...)
            %
            % SLICEPLOT provides the arguments X,Y,Z,V for slice;
            %  all others can be provided by the user (except an axes to
            %  to plot into).
            
            assert(prod(fun.m) == 1, 'Function:slice:vectorValued',...
                'This function isn''t a scalar function.');
            assert(isa(fun.grid, 'Grid'), 'Function:slice:irregular',...
                'slice requires a regular grid.');
            assert(fun.grid.dims == 3, 'Function:slice:not3D',...
                'slice only works for 3D Functions');
            
            gv = fun.grid.gridVectors;
            [X,Y,Z] = meshgrid(gv{:});
            hh = slice(X,Y,Z, permute(fun.f, [2 1 3]), varargin{:});
            set(hh, 'EdgeColor', 'none');
            
            if nargout > 0, h = hh; end;
        end
        
        function h = tsliceplot(fun,varargin)
            % Create a slice plot of a 3D Function using `slice`.
            %
            %   f.TSLICEPLOT(Sx,Sy,Sz)
            %   f.TSLICEPLOT(XI,YI,ZI)
            %   f.TSLICEPLOT(...)
            %   h = f.TSLICEPLOT(...)
            %
            % Calls Function.sliceplot; additionally, the alpha channel is
            %  set so that the slices' transparency are proportional to
            %  the absolute value of the real part of f.
            %
            % In particular, the slices are transparent where f = 0.
            hh = sliceplot(fun,varargin{:});
            for i = 1:numel(hh)
                set(hh(i), 'FaceAlpha', 'flat', 'AlphaData', abs(get(hh(i), 'CData')));
            end
            if nargout > 0, h = hh; end;
        end
    end
end