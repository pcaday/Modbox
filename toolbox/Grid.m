classdef Grid < handle
    % Read-write properties.
    properties (SetObservable)
        mins        % Minimum values in each coordinate (n-vector)
        maxes       % Maximum values in each coordinate (n-vector)
        ns          % Number of grid points in each coordinate (n-vector)
        periodic    % Whether the grid is periodic in each coordinate (boolean n-vector)
    end
    
    % Read-only properties
    properties (SetAccess = protected)
        dims        % # of dimensions of grid.
        ds          % Distance between grid points in each coordinate (n-vector)
        gridVecs    % 1 x dims cell array of grid vectors for each coordinate
    end
    
    
    methods
        function grid = Grid(mins, maxes, ns, periodic)
            % Create a Grid object.
            %
            %  grid = Grid.WithSpacing(mins, maxes, ns);
            %  grid = Grid.WithSpacing(mins, maxes, ns, periodic);
            %
            % Inputs
            %      mins: vector of minimum values in each dimension
            %     maxes: vector of maximum values in each dimension
            %        ns: nonnegative integer vector: number of grid points in each dimension
            %  periodic: [optional] logical vector specifying for each dimension
            %             whether the grid is periodic in that dimension.
            %            May be a scalar to specify all dimensions at once.
            %            Defaults to nonperiodic (false).
            
            % Check data types of input.
            CheckTypes({mins, maxes, ns}, {'numeric','numeric','numeric'});
            
            % Check that mins, maxes, and ns all have the same length.
            grid.dims = length(mins);
            
            if nargin < 4 || isempty(periodic),
                periodic = false(1, grid.dims);
            elseif length(periodic) == 1,
                periodic = repmat(periodic, 1, grid.dims);
            end
            
            assert(all(grid.dims == [length(maxes) length(ns) length(periodic)]),...
                'Grid:Grid:dimMismatch',...
                'The inputs must have the same lengths.');
            
            % Assign mins, maxes, and periodic flags
            grid.mins = mins;
            grid.maxes = maxes;
            grid.periodic = periodic;
            
            % Add listener for ns to adjust its value and ds after setting
            grid.addlistener('ns', 'PostSet', @Grid.NSListener);
            % Add listeners for mins, maxes so ds will be updated
            %  when they change.
            grid.addlistener('mins', 'PostSet', @Grid.DSListener);
            grid.addlistener('maxes', 'PostSet', @Grid.DSListener);
            
            % Assign ns, and trigger the update to ns, ds, and the grid
            % vectors.
            grid.ns = ns;
        end
        
        function newGrid = copy(grid)
        % Make a copy of this grid.
        %
        %   newGrid = oldGrid.copy();
        %
        
        %  Use the Grid constructor to make
        %  sure the listeners are copied too (matlab.mixin.Copyable
        %  doesn't copy the property set listeners)
            newGrid = Grid(grid.mins, grid.maxes, grid.ns);
        end
    end
    
    methods (Static)
        function grid = FromGridVectors(gvs)
            % Create a Grid object from a cell array of its grid vectors.
            %
            %  grid = Grid.FromGridVectors({gv1, gv2, ... gvN});
            %
            assert(iscell(gvs) && ~isempty(gvs), 'Grid:FromGridVectors:badGVs',...
                'Expected a non-empty cell array of grid vectors.');
            gridMins = cellfun(@min, gvs);
            gridMaxes = cellfun(@max, gvs);
            gridNs = cellfun(@length, gvs);
            grid = Grid(gridMins, gridMaxes, gridNs);
        end
        
        function grid = WithSpacing(mins, maxes, ds, periodic)
            % Create a Grid object with given grid spacing in each dimension.
            %
            %  grid = Grid.WithSpacing(mins, maxes, ds);
            %  grid = Grid.WithSpacing(mins, maxes, ds, periodic);
            %
            % Inputs
            %      mins: vector of minimum values in each dimension
            %     maxes: vector of maximum values in each dimension
            %        ds: vector of grid spacing in each dimension
            %  periodic: [optional] logical vector specifying for each dimension
            %             whether the grid is periodic in that dimension.
            %            May be a scalar to specify all dimensions at once.
            %            Defaults to nonperiodic (false).
            %
            % If maxes - mins isn't an even multiple of the grid spacing
            %  ds, then ds is rounded to a nearby value.
            if nargin < 4, periodic = false(size(mins)); end
            
            grid = Grid(mins, maxes, 2*ones(size(mins)), periodic);
            grid.setDS(ds);
        end
    end
    
    methods (Static)
        function DSListener(source, event) %#ok<INUSL>
            % Listener for changes in ds
            event.AffectedObject.updateDS();
            event.AffectedObject.updateVectors();
        end
        
        function NSListener(source, event) %#ok<INUSL>
            % Listener for changes in ns
            event.AffectedObject.updateNS();
            event.AffectedObject.updateVectors();
        end
    end
    
    methods (Access = protected)
        function updateNS(grid)
            % Make sure ns is at least two in each dimension, and round
            %  up to the nearest integer.
            grid.ns = max(ceil(grid.ns), 2);
            grid.updateDS();
        end
        
        function updateDS(grid)
            grid.ds = (grid.maxes - grid.mins) ./ (grid.ns - 1);
            grid.updateVectors();
        end
        
        function updateVectors(grid)
           % Create cell array of grid vectors for each dimension.
            grid.gridVecs = cell(1, grid.dims);
            for i = 1:grid.dims
                grid.gridVecs{i} = grid.mins(i) + grid.ds(i) .* (0:(grid.ns(i)-1));
                %grid.gridVecs{i} = linspace(grid.mins(i),grid.maxes(i),grid.ns(i));
            end 
        end
    end
    
    methods
        function [varargout] = subsref(grid, ss)
            % Override of built-in subsref (see SUBSREF).
            %
            % Allows access to the coordinates of grid points using
            %  parenthesis indexing. For example, if G is a 2D Grid, then
            %  G(3,5) returns the coordinates of the grid point at (3,5),
            %  as a length-2 column vector.
            %
            % Possible usages:
            %  G(i1,...,in)
            %      Coordinates of one or more grid points.
            %       i1,...,in may be integers, vectors, ':', and can
            %       involve END.
            %      If i1,...,in specify an array of grid points of size
            %        sz1 x ... x szn, then G(i1,...,in) is an array of size
            %        n x sz1 x ... x szn.
            %  G(:)
            %      This returns the coordinates of all grid points as
            %       an n x numel(G) array. Use G.AllPoints() to get all
            %       grid points as an array which preserves G's structure.
            switch ss(1).type
                case '()'
                    % Return coords of the i'th grid point,
                    %  so grid(i1,...,in) will return the real coordinates
                    %  corresponding to the grid point with indices i1...in.
                    n = grid.dims;
                    % Only allow subscripts with n coordinates, where
                    %  n is the dimension of the grid, and (:) subscripts
                    %  And don't allow sequences of subscripts
                    %     [like grid(1,2,3)(4,5,6)]
                    %   -- actually MATLAB doesn't support this so this
                    %      case should never happen.
                    assert(length(ss) == 1, 'Grid:subsref:multiSub',...
                        'Multiple subscripts not supported.'); % not supported by MATLAB either, so this should NEVER happen!
                    assert(nargout <= 1, 'Grid:subsref', ...
                        'Too many output arguments.');
                    nsubs = length(ss(1).subs);
                    if nsubs == 1 && ischar(ss(1).subs{1})...
                            && ss(1).subs{1} == ':',
                        % Handle G(:) case
                        [coords{1:n}] = ndgrid(grid.gridVecs{:});
                        p = fcat(coords{:});
                        p = reshape(p, n, []);
                    elseif nsubs == n,
                        % Loop through the subscripts
                        subs = [ss(1).subs];
                        for i = 1:n
                            % Replace any ':' inputs by the complete list of
                            % indices for that coordinate.
                            if ischar(subs{i}) && strcmp(subs{i}, ':');
                                subs{i} = 1:grid.ns(i);
                            end
                            % Convert grid indices to the actual coordinates
                            subs{i} = grid.mins(i) + grid.ds(i) * (subs{i}-1);
                            
                            % Alternatively,
                            % subs{i} = grid.gridVecs{i}(subs{i});
                        end
                        % Generate the grid
                        [coords{1:n}] = ndgrid(subs{:});
                        % Concatenate them to get coordinate lists.
                        p = fcat(coords{:});
                        
                        % Remove trailing singleton dimensions from the output
                        %  (like MATLAB's standard subscripting seems to do)
                        % Make sure to leave at least 1 dimension since
                        %  MATLAB's arrays have to have at least 2 dimensions
                        %  total.
                        % Don't need to do this any more!
                        %lasttokeep = max([1 find(cellfun('length', subs)~=1)]);
                        %i = size(p);
                        %i((lasttokeep+1):n) = [];
                        %p = reshape(p, i);
                    else
                        error(...
                            'Grid:subsref:reshapeSub',...
                            'Only subscripts with the same number of entries as the grid dimension are supported, or a single colon (grid(:)).');
                    end
                    varargout = {p};
                otherwise
                    %try
                    % For { } and .'s, use MATLAB's built-in functionality.
                    if nargout > 0,
                        varargout = cell(1, nargout);
                        [varargout{:}] = builtin('subsref', grid, ss);
                    else
                        % No-argument case is tricky!
                        % If an output argument is returned anyway, catch
                        %  it in ans and return it. We fill
                        %  ans with a bogus value, and if it is changed
                        %  by the built-in subsref, we return it.
                        ans = {'VOID!'};
                        builtin('subsref', grid, ss);
                        if ~(iscell(ans) && length(ans) == 1 && ...
                                ischar(ans{1}) && strcmp(ans{1}, 'VOID!')) %#ok<*NOANS>
                            varargout = {ans};
                        end
                    end
                    %catch ME
                    %    throw(ME)
                    %end
            end
        end
        
        % Overload 'end' to allow expressions like grid(2:end,:)
        %  'end' tells MATLAB how "big" G is in each dimension.
        function e = end(grid,k,m)
            if k == m,
                % Return the product of the lengths of the remaining
                %  dimensions, k, k+1, ..., dims. This is what MATLAB does
                %  too.
                e = prod(grid.ns(k:end));
            else
                % Return the length of dimension k
                e = grid.ns(k);
            end
        end
        
        function n = NPoints(grid)
        % Return total number of points in grid.
            n = prod(grid.ns);
        end
        
        function p = AllPoints(grid)
        % Return (n+1)-D array of all the points in the grid.
        % This is the same as grid(:,:,...:) but easier to use if
        %  the number of dimensions of the grid is unknown.
            n = grid.dims;
            [coords{1:n}] = ndgrid(grid.gridVecs{:});
            p = fcat(coords{:});
        end
        
        function d = ndims(grid)
        % Return number of dimensions of grid, MATLAB-style
        %  (returns 2 dimensions for 1D grids). For the real number of
        %  dimensions, use the `dims` property.
            d = max(grid.dims, 2);
        end
        
        function varargout = size(grid, dim)
        % Return size of grid, MATLAB-style
        %   (for 1D grids of length n, returns n x 1)
            if nargin < 2,
                if nargout <= 1,
                    % Size vector
                    s = grid.ns;
                    if length(s) == 1, s = [s 1]; end
                    varargout{1} = s;
                else
                    m = min(nargout, grid.dims);
                    varargout(1:m) = num2cell(grid.ns(1:m));
                    if nargout > grid.dims
                        varargout(m+1:nargout) = {1};
                    else
                        varargout{nargout} = prod(grid.ns(nargout:end));
                    end
                end
            else
                % Size in the given dimension
                if dim > grid.dims
                    s = 1;
                else
                    s = grid.ns(dim);
                end
                varargout{1} = s;
            end
        end
        
        function s = realSize(grid)
        % Returns the (real) size of the grid.
        %  Same as size, but returns [n] for a 1D grid of length n,
        %  instead of [n 1].
            s = grid.ns;
        end
        
        function vectors = gridVectors(grid)
        % Return a cell array of grid vectors for this grid (as used by
        %  ndgrid or meshgrid, for instance).
            vectors = grid.gridVecs;
        end
        
        function funs = coordFuns(grid)
        % Return a 1-by-n cell array of coordinate functions.
        %  These aren't Function objects, just plain arrays.
            funs = cell(1, grid.dims);
            [funs{:}] = ndgrid(grid.gridVectors{:});
        end
        
        function funs = pixelEdgeCoordFuns(grid)
        % Return a 1-by-n cell array of "pixel edge" coordinate functions,
        %  where each grid point is treated as a "pixel" with dimensions
        %  grid.ds. The edges are halfway between adjacent grid points.
        
        %  These aren't Function class objects, just plain arrays.
            edgeGridVecs = cell(1, grid.dims);
            
            for i = 1:grid.dims
                edgeGridVecs{i} = grid.mins(i) + grid.ds(i) .* ((0:grid.ns(i)) - 0.5);
            end 

            funs = cell(1, grid.dims);
            [funs{:}] = ndgrid(edgeGridVecs{:});
        end
        
        function f = eq(grid1, grid2)
        % Compare two Grids for equality.
        %
        % Only works for scalar inputs, but
        %  Grids aren't meant to work in arrays and won't work well
        %  in arrays because () subscripting is overloaded.
            if ~isa(grid1, 'Grid') || ~isa(grid2, 'Grid')
                f = false;
            else
                f = isequal(grid1.mins, grid2.mins) ...
                    && isequal(grid1.maxes, grid2.maxes) ...
                    && isequal(grid1.ns, grid2.ns) ...
                    && isequal(grid1.periodic, grid2.periodic);
            end
        end
        
        function singleGrid = single(grid)
        % Convert a Grid object to single-precision.
            singleGrid = Grid(single(grid.mins),...
                            single(grid.maxes),...
                            single(grid.ns),...
                            logical(grid.periodic));
        end

        function doubleGrid = double(grid)
        % Convert a Grid object to double-precision.
            doubleGrid = Grid(double(grid.mins),...
                            double(grid.maxes),...
                            double(grid.ns),...
                            logical(grid.periodic));
        end
        
        function newGrid = convertClass(grid, newClass)
            % Convert the grid's data to a new class
            %  (currently single and double supported).
            switch (newClass)
                case 'single'
                    newGrid = single(grid);
                case 'double'
                    newGrid = double(grid);
                otherwise
                    error('Unknown class for conversion.');
            end
        end

        function xl = toGridSpace(grid, xg, coord)
        % Convert from global coordinates to grid entries.
        %
        %   xl = grid.toGridSpace(xg)
        %
        %  xg and xl are n x anything arrays, where n = # of dimensions,
        %   representing points in n-D space.
        %
        %
        %   xl = grid.toGridSpace(xg, i)
        %
        %  With this syntax, xl is an array of values for the
        %   i'th variable, rather than containing all coordinates.
        %
            mins_ = grid.mins;
            ds_ = grid.ds;
            if nargin < 3,
                for d = grid.dims:-1:1
                    xl(d,:) = 1 + (xg(d,:) - mins_(d)) / ds_(d);
                end
            else
                xl = 1 + (xg - mins_(coord)) / ds_(coord);
            end
            xl = reshape(xl,size(xg));
        end
        
        function xg = fromGridSpace(grid, xl, coord)
        % Convert from grid entries to global coordinates.
        %
        %   xg = grid.fromGridSpace(xl)
        %
        %  xg and xl are n x anything arrays, where n = # of dimensions,
        %   representing points in n-D space.
        %
        %
        %   xg = grid.fromGridSpace(xl, i)
        %
        %  With this syntax, xl is an array of values for the
        %   i'th variable, rather than containing all coordinates.
        %
            mins_ = grid.mins;
            ds_ = grid.ds;
            if nargin < 3,
                for d = grid.dims:-1:1
                    xg(d,:) = mins_(d) + (xl(d,:) - 1) * ds_(d);
                end
                xg = reshape(xg,size(xl));
            else
                xg = mins_(coord) + (xl - 1) * ds_(coord);
            end
        end

        function in = IsInside(grid, x)
        % Tests whether point(s) are inside the grid.
        %
        %   in = grid.isInside(x)
        %
        % x is an n x anything array, where n is the number of dimensions
        %  of the grid.
        %
        % The output is a logical array of the same size as x, minus the
        %  leading dimension.
            sz = size(x);
            in = true(1, prod(sz(2:end)));
            for i = 1:grid.dims
                in = in & x(i,:) >= grid.mins(i) ...
                        & x(i,:) <= grid.maxes(i);
            end
            
            in = reshape(in, [sz(2:end) 1]);
        end
        
        function xWrap = WrapToGrid(grid, x)
        % Wrap any periodic coordinates of the given point(s) so they lie
        %  inside the grid.
        %
        %   xWrap = grid.WrapToGrid(x);
        %
        % x is an n x anything array, where n is the number of dimensions
        %  of the Grid.
        % The output is an array of the same size as x containing the
        %  wrapped points.
            xWrap = x;
            for i = 1:grid.dims
                if grid.periodic(i),
                    mn = grid.mins(i);
                    mx = grid.maxes(i);
                    xWrap(i,:) = mod(xWrap(i,:) - mn, mx - mn) + mn;
                end
            end
        end
        
        function eGrid = PixelEdgeGrid(grid)
        % Return the "pixel edge" grid of this grid.
        %
        %   eGrid = grid.PixelEdgeGrid();
        %
        % The pixel edge grid points are halfway between the original
        %  grid points, and surrounding the original grid. The pixel
        %  edge grid size is one more than the original in each coordinate.
            halfDS = grid.ds * 0.5;
            eGrid = grid.copy();
            eGrid.ns = eGrid.ns + 1;
            eGrid.mins = eGrid.mins - halfDS;
            eGrid.maxes = eGrid.maxes + halfDS;
        end
        
        function fftg = DFTGrid(grid, zeroPadding)
        % Return the corresponding frequency grid for a spatial grid.
        %
        % Assumes the value of a function on a single grid cell
        %  is represented by its lower-left (in 2D) corner, or in general,
        %  by the grid point at the corner of the cell with the most
        %  negative values in each coordinate.
            if nargin < 2, zeroPadding = 1; end
            
            Ns = realSize(grid);
            paddedNs = ceil(Ns * zeroPadding);
            
            Ls = grid.maxes - grid.mins + grid.ds;
            paddedLs = Ls .* paddedNs ./ Ns;
            
            omegaMins = - floor(paddedNs/2) * (2*pi) ./ paddedLs;
            omegaMaxes = floor((paddedNs-1)/2) * (2*pi) ./ paddedLs;
            
            fftg = Grid(omegaMins, omegaMaxes, paddedNs);
        end
        
        function sg = Slice(grid, dims)
        % Get a slice of this Grid. A "slice" is a grid consisting of one
        %  or more of the coordinates of the original grid.
        %
        % The dimensions to use are specified by the vector dims.
        %
        % Slice can also be used to permute the grid coordinates:
        %  e.g. grid.Slice([2 1]) for a 2D grid swaps the variables.
            sg = Grid(grid.mins(dims), grid.maxes(dims), grid.ns(dims),...
                grid.periodic(dims));
        end
        
        function setDS(grid, ds)
        % Set the grid spacing (ds).
        %
        % Any ds entry does not go evenly into the width of the grid in
        %  its dimension will be rounded.
            assert(length(ds) == length(grid.ds), 'Grid:setDS:badLength',...
                'ds is the wrong length');
            assert(all(grid.maxes > grid.mins), 'Grid:setDS:zeroDim',...
                'The grid has zero size in one or more dimensions.');
            grid.ns = round(1 + (grid.maxes - grid.mins) ./ ds);
        end
        
        function setDSExactly(grid, ds)
        % Set the grid spacing (ds).
        %
        % The grid will be expanded if necessary to accomodate a whole
        %  number of grid points in each direction. If this occurs, the
        %  expansion will occur symmetrically around the center of the
        %  grid.
        %
        % This method logically should not be used for periodic grids.
            assert(length(ds) == length(grid.ds), 'Grid:setDS:badLength',...
                'ds is the wrong length');
            assert(all(grid.maxes > grid.mins), 'Grid:setDS:zeroDim',...
                'The grid has zero size in one or more dimensions.');

            np = ceil((grid.maxes - grid.mins) ./ ds);
            center = (grid.maxes + grid.mins) * 0.5;
            hw = np .* ds * 0.5;
            
            grid.mins  = center - hw;
            grid.maxes = center + hw;
            grid.ns    = 1 + np;
        end
    end
end