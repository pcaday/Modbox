% Approximate dy/dx(x,nu).
%
%   dxy = ApproximateDXY(y, xGrid, yGrid)
%
% Inputs
%        y: values of y(x,nu) for a fixed nu and all x in the grid.
%    xGrid: the domain (x) grid
%    yGrid: the codomain (y) grid
%
% Outputs
%      dxy: values of dy/dx(x,nu).
%           dxy(i,j,:,...,:) will contain the values of
%            dy_i/dx_j at each point in the x grid.
%
function dxy = ApproximateDXY(y, xGrid, yGrid)
    dxy_edge = Gradient(y,xGrid,yGrid);
    dxy = InterpolateCornersToCenters(dxy_edge, size(y,1), 2);
end