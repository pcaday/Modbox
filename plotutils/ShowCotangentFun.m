% Display a function on the cotangent bundle (currently only for 2D)
%
%    ShowCotangentFun(f, wcr, xSize, ...)
%
%      f: I1 x ... x In x M or I x M array, where the I indices
%          are for the x coordinates, and the M index is the box
%          number in the frequency domain.
%         If f is I x M, the xSize argument must be provided.
%         f may also be a cell array of such arrays, each representing
%          a function, in which case each function will be placed
%          in a subplot.
%         
%    wcr: the WedgeCreator object that holds that xi values we have
%          data from f on.
%  xGrid: (optional) x grid for the function.
%
%    Any additional arguments are passed to the plot method of the Function
%    class.
%
function ShowCotangentFun(f, wcr, xGrid, varargin)
    global GlobalFunctionCxPlot
    
    % Wrap f in a cell array if needed.
    if ~iscell(f), f = {f}; end
    
    % Deduce x grid if not provided.
    if nargin < 3 || isempty(xGrid)
        xSize = size(f{1});
        xSize(end) = [];
        xGrid = Grid(ones(size(xSize)), xSize, xSize);
    end
    
    % If no range provided, choose one automatically.
    if nargin < 4 || isempty(varargin{1}),
        % Check whether this is a complex plot.
        cx = GlobalFunctionCxPlot;
        if nargin >= 5,
            cx = varargin{2};
        end
        % Choose an appropriate range.
        minv = +inf;
        maxv = -inf;
        for k = 1:length(f)
            v = f{k}(isfinite(f{k}));
            if cx,
                v = abs(v);
            else
                v = real(v);
            end
            minv = min(min(v), minv);
            maxv = max(max(v), maxv);
        end
        if cx, minv = 0; end
        % Store the range in the arguments to Function.plot
        varargin{1} = [minv maxv];
    end
    
    switch wcr.grid.dims
        case 2
            % Flatten the x coordinates if necessary.
            for k = 1:length(f)
                sz = size(f{k});
                f{k} = reshape(f{k}, [], sz(end));
            end
            
            slider(@SCF2DCallback, [0.51 wcr.Count()+0.49])
        otherwise
            error('Only 2D currently supported.');
    end
    
    colormap gray; 
    
    function SCF2DCallback(box)
        box = round(box);
        nuBasis = wcr.NuBasis(box);
        nfuns = length(f);
        
        for kk = 1:nfuns
            subplot(1,nfuns,kk);
            
            fun = Function.WithValues(xGrid, reshape(f{kk}(:,box), size(xGrid)));
            fun.plot(varargin{:}); colorbar;

            hold on

            % Old code -- only works for cartesian DFT grids
%             r = min(abs(wcr.grid.mins));
%             plot([0 r*nuBasis(1,1)], [0 r*nuBasis(2,1)], 'ro-');
%             plot([0 r*nuBasis(1,2)], [0 r*nuBasis(2,2)], 'go-');

            % New code -- works for cartesian and polar DFT grids
            limx = xlim;
            limy = ylim;
            center = [mean(limx) mean(limy)];
            r = min(diff(limx), diff(limy)) / 2;
            plot([0 r*nuBasis(1,1)]+center(1), [0 r*nuBasis(2,1)]+center(2), 'ro-');
            plot([0 r*nuBasis(1,2)]+center(1), [0 r*nuBasis(2,2)]+center(2), 'go-');

            AddImageXYTracker
            
            hold off
        end
    end
end