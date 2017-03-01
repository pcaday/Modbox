function RenderQuad(O, x1, x2, v)
    if nargin < 4, v = 1; end
    
    rx1 = reshape(x1([1 2 4 3]),2,2);
    rx2 = reshape(x2([1 2 4 3]),2,2);
    RenderQuadMesh(O, rx1, rx2, v, 0);
end
