classdef IrregularGrid < matlab.mixin.Copyable
    properties
        dims
        points        % n x i1 x ... x in array, containing the grid points
        ns
    end
    
    methods
        function grid = IrregularGrid(points)
            sz = size(points);
            grid.dims = sz(1);
            grid.ns = sz(2:end);
            grid.points = points;
        end
        
        function apoints = AllPoints(grid)
            apoints = grid.points;
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
        % Returns the size of grid.
        %  Same as size, but returns [n] for a 1D grid of length n,
        %  instead of [n 1].
            s = grid.ns;
        end
        
        function funs = coordFuns(grid)
        % Return a 1-by-n cell array of coordinate functions.
        %  These aren't Function class objects, just plain arrays.
            funs = num2cell(grid.points, 2:grid.dims+1);
        end

        function [varargout] = subsref(grid, ss)
        % Overloaded subsref to retrieve grid points using subscripting
        %  e.g. grid(3,5) or grid(:,2) or grid(:).
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
                    assert(length(ss) == 1, 'IrregularGrid:subsref:multiSub',...
                        'Multiple subscripts not supported.'); % not supported by MATLAB either, so this should NEVER happen!
                    assert(nargout <= 1, 'IrregularGrid:subsref', ...
                        'Too many output arguments.');
                    
                    subs = ss(1).subs;
                    assert(~isempty(subs), 'IrregularGrid:subsref:nosubs',...
                        'No subscripts provided.');
                    
                    % Add ':' to the beginning and subsref the points.
                    subs = [':' subs];
                    
                    varargout = {grid.points(subs{:})};
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
        
        function d = ndims(grid)
        % Return number of dimensions of grid, MATLAB-style
        %  (returns 2 dimensions for 1D grids). For the real number of
        %  dimensions, access grid.dims
            d = max(grid.dims, 2);
        end
        
        function f = eq(grid1, grid2)
        % Compare two IrregularGrids (only works for scalar inputs, but
        %  IrregularGrids aren't meant to work in arrays and won't work well
        %  in arrays because () subscripting is overloaded.)
            if ~isa(grid1, 'IrregularGrid') || ~isa(grid2, 'IrregularGrid')
                f = false;
            else
                f = isequal(grid1.points, grid2.points);
            end
        end
        
        function singleGrid = single(grid)
        % Convert to single-precision floating point
            singleGrid = IrregularGrid(single(grid.points));
        end

        function doubleGrid = double(grid)
        % Convert to double-precision floating point
            doubleGrid = IrregularGrid(double(grid.points));
        end
        
        function Draw(grid, varargin)
        % Draw the grid.
        %
        %  grid.Draw
        %  grid.Draw(N)
        %  grid.Draw(N, doLegend)
        %
        % N is the number of grid lines in each coordinate to draw.
        % doLegend is a Boolean flag indicating whether to draw a legend or
        %  not.
            DrawIrregularGrid(grid.points, varargin{:});
        end
    end
end