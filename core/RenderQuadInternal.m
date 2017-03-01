%
% A wrapper for RenderQuadMeshInternal that makes it behave like
%  the (old) RenderQuadInternal MEX function.
%
function out = RenderQuadInternal(in, x1, x2, v)
    valueClass = class(in);
    if nargin < 4, v = ones(1, valueClass); end
    
    rx1 = reshape(x1([1 2 4 3]),2,2);
    rx2 = reshape(x2([1 2 4 3]),2,2);
    out = RenderQuadMeshInternal(in, rx1, rx2, v, zeros(1,valueClass));
end
