% Constant sound speed.
%
% Return value is an Lstruct, ready to pass to EuclideanConformalMetric.
%
function Lstruct = ConstantSoundSpeed(c0,dims)
    if nargin < 1 || isempty(c0), c0 = 1; end
    if nargin < 2 || isempty(dims), dims = 2; end
    
    Lstruct.dims = dims;
    Lstruct.lambda = @ConstLambda;
    Lstruct.dlambda = cell(1,dims);
    for i = 1:dims
        Lstruct.dlambda{i} = @ZeroFunc;
        for j = 1:dims
            Lstruct.d2lambda{i,j} = @ZeroFunc;
        end
    end
    
    nlog_c0 = -log(c0);
    
    %%%
    %
    % Helper functions
    %
    %
    function v = ConstLambda(x)
        sz = size(x);
        v = repmat(nlog_c0, [sz(2:end) 1]);
    end

    function v = ZeroFunc(x)
        sz = size(x);
        v = zeros([sz(2:end) 1]);
    end
end