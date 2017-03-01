% Sound speed "well" function from de Hoop, Uhlmann, Vasy, and Wendt's
%  paper:
%
%    c = c0 + kappa * exp(-|x-x0|^2/sigma^2)
%
% Return value is an Lstruct, ready to pass to EuclideanConformalMetric.
%
function Lstruct = SoundSpeedWell(c0,x0,kappa,sigma,dims,grid)
    if nargin < 1 || isempty(c0), c0 = 2; end
    if nargin < 2 || isempty(x0), x0 = [0; 0]; end
    if nargin < 3 || isempty(kappa), kappa = -0.4; end
    if nargin < 4 || isempty(sigma), sigma = 0.5; end
    if nargin < 5 || isempty(dims), dims = 2; end
    if nargin < 6, grid = Grid([0 0], [1 1], [2 2], [false false]); end
    
    sig2 = sigma * sigma;
    
    Lstruct.dims = dims;
    Lstruct.lambda = @WellLambda;
    Lstruct.dlambda = cell(1,dims);
    for i = 1:dims
        Lstruct.dlambda{i} = @(x) WellDLambda(x,i);
        for j = 1:dims
            Lstruct.d2lambda{i,j} = @(x) WellD2Lambda(x,i,j);
        end
    end
    
    
    
    %%%
    %
    % Helper functions
    %
    %
    function v = WellLambda(x)
        x = grid.WrapToGrid(x);
        ep = expPart(x);
        v = -log(C(x,ep));
    end

    function v = WellDLambda(x, i)
        x = grid.WrapToGrid(x);
        ep = expPart(x);
        v = -DC(x,i,ep) ./ C(x,ep);
    end

    function v = WellD2Lambda(x, i, j)
        x = grid.WrapToGrid(x);
        ep = expPart(x);
        c = C(x,ep);
        v = -D2C(x,i,j,ep) ./ c ...
            + DC(x,i,ep) .* DC(x,j,ep) ./ (c.*c);
    end

    function c = C(~,ep)
        c = c0 + kappa * ep;
    end

    function dc = DC(x, i, ep)
        x_i = (x(i,:) - x0(i));
        % Reshape to correct size
        sz = size(x);
        x_i = reshape(x_i, [sz(2:end) 1]);
        
        dc = -2*kappa/sig2 * ep .* x_i;
    end
       
    function d2c = D2C(x, i, j, ep)
        x_ij = (x(i,:) - x0(i)) .* (x(j,:) - x0(j));
        % Reshape to correct size
        sz = size(x);
        x_ij = reshape(x_ij, [sz(2:end) 1]);
        
        d2c = 4*kappa/(sig2*sig2) * ep .* x_ij;
        if i == j,
            d2c = d2c - 2*kappa/sig2 * ep;
        end
    end

    function ep = expPart(x)
        % Compute |x-x0|^2
        if dims == 2,
            % Do 2D as a special case for speed.
            d = (x(1,:) - x0(1)).^2 + (x(2,:) - x0(2)).^2;
            sz = size(x);
            d = reshape(d, [sz(2:end) 1]);
        else
            d = bsxfun(@minus, x, x0);
            d = fsqueeze(sum(d.*d, 1));
        end
        
        % Compute exponential
        ep = exp(-d/sig2);
    end
end